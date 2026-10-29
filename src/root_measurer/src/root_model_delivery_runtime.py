"""Existing pipeline codeobjects with fixed span backend and exact receipt reuse."""
import copy,hashlib,json,os,sys
from pathlib import Path
from types import FunctionType
import numpy as np
import root_automatic_model_measurement as automatic
import root_automatic_source_profile as source_profile
from root_span_model_pipeline import measure_span_tube
from root_exact_query_sections import ExactQuerySectionMeasurer,_exact_query
from root_graph_source_length import pack

P=Path(__file__).resolve().parents[1]
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def bound_function(fn,**namespace):
    result=FunctionType(fn.__code__,dict(fn.__globals__,**namespace),fn.__name__,fn.__defaults__);result.__kwdefaults__=fn.__kwdefaults__
    assert result.__code__ is fn.__code__
    return result
def key(point,normal):
    p,n=_exact_query(point,normal)
    return json.dumps([list(map(pack,p)),list(map(pack,n))],separators=(',',':'))

def measure_integrated(argv,*,stage_cache=None,previous_measurement=None,output_guard=None,input_loader=None):
    """argv remains the existing mesh-only CLI schema. Cache inputs are optional."""
    sys.path.insert(0,str(P/'scripts'))
    import measure_automatic_model as entry
    source=Path(argv[argv.index('--mesh')+1]);source_sha=sha(source);sessions=[];cached_sections={};cache_proof=None
    def read_bound(ref):
        if not isinstance(ref,dict) or set(ref)!={'path','sha256','bytes'}:raise ValueError('FROZEN_CACHE_FILE_REFERENCE_REQUIRED')
        path=Path(ref['path'])
        if path.stat().st_size!=ref['bytes'] or sha(path)!=ref['sha256']:raise ValueError('FROZEN_CACHE_SHA_BYTES_CHANGED')
        return path,json.loads(path.read_text(encoding='utf-8-sig'))
    if previous_measurement is not None:
        previous_path,old=read_bound(previous_measurement)
        if old['input_file_binding']['sha256']!=source_sha:raise ValueError('SECTION_CACHE_REQUIRES_IDENTICAL_ORIGINAL_SOURCE')
        ref=dict(path=str(previous_path),sha256=sha(previous_path),bytes=previous_path.stat().st_size)
        candidates=[('/automatic_generation_evidence/original_triangle_sparse_radius_checks/%d/section'%i,r['section']) for i,r in enumerate(old.get('automatic_generation_evidence',{}).get('original_triangle_sparse_radius_checks',[]))]
        candidates += [('/radius_curve/rows/%d/original_mesh_section'%i,r['original_mesh_section']) for i,r in enumerate((old.get('radius_curve') or {}).get('rows',[]))]
        for pointer,r in candidates:
            if r.get('source_sha256')==source_sha and r.get('exact_query_point') and r.get('index_evidence',{}).get('all_original_faces_classified') is True and r.get('query_was_rounded_to_binary64') is False:
                cached_sections.setdefault(key(r['exact_query_point'],r['exact_query_normal']),(copy.deepcopy(r),dict(reference=ref,pointer=pointer)))
    cache_data=None
    if stage_cache is not None:
        (a,first),(b,second)=map(read_bound,stage_cache)
        # The failed-stage receipt binds the exact input through its frozen contract.
        first_cp=P.parents[1]/first['contract']['workspace_relative_path'];first_contract=json.loads(first_cp.read_text())
        source_ref=first_contract['source'];assert source_ref['sha256']==source_sha and source_ref['bytes']==source.stat().st_size and first['voxelization']['occupied_cells']==len(first['source_cells'])
        assert first['contract']['sha256']==sha(first_cp) and second['new_kernel']['integration_binding']['only_existing_span_thinning_backend_bound'] is True
        second_cp=P.parents[1]/second['contract']['workspace_relative_path'];second_contract=json.loads(second_cp.read_text());assert sha(second_cp)==second['contract']['sha256'] and second_contract['previous_failed_diagnostic']['sha256']==sha(a)
        for r in second_contract['sources']:
            path=P.parents[1]/r['workspace_relative_path'];assert sha(path)==r['sha256'],'STAGE_CACHE_PRODUCER_CHANGED'
        cells=np.asarray(first['source_cells'],np.int64);cells.flags.writeable=False
        cache_data=dict(cells=cells,voxelization=first['voxelization'],kernel=second['new_kernel'],pitch=first_contract['pitch'])
        cache_proof=dict(cells_receipt_sha256=sha(a),kernel_receipt_sha256=sha(b),source_sha256=source_sha,existing_cells_and_kernel_reused=True)
    def factory(phase):
        class CachedSections:
            def __init__(self,mesh,*,execution_domain='local',max_queries=5,max_candidate_faces=100000,max_work=10000000):
                nv,nf=len(mesh.vertices),len(mesh.faces)
                if not 0<nv<=1200000 or not 0<nf<=2000000:raise ValueError('UNCHANGED_COMPLETE_SOURCE_GUARDS')
                if execution_domain=='local':
                    if nv>3000 or nf>1000:raise ValueError('UNCHANGED_LOCAL_SOURCE_GUARDS')
                elif execution_domain!='registered_slurm_complete_source' or not os.environ.get('SLURM_JOB_ID','').isdigit():raise ValueError('REGISTERED_SOURCE_INDEX_REQUIRED')
                if type(max_queries)is not int or not 1<=max_queries<=5 or type(max_candidate_faces)is not int or not 1<=max_candidate_faces<=100000 or type(max_work)is not int or not 0<=max_work<=10000000:raise ValueError('UNCHANGED_QUERY_CANDIDATE_WORK_GUARDS')
                self.mesh=mesh;self.kw=dict(execution_domain=execution_domain,max_queries=max_queries,max_candidate_faces=max_candidate_faces,max_work=max_work);self.engine=None;self.queries=0
                self.array_sha=hashlib.sha256(mesh.vertices.tobytes()+mesh.faces.tobytes()).hexdigest();self.stats=dict(phase=phase,requested=0,cached=0,fresh=0,new_index_built=False,new_logical_work=0);sessions.append(self.stats)
            def section(self,point,normal):
                if self.queries>=self.kw['max_queries']:raise ValueError('UNCHANGED_QUERY_INDEX_CAP')
                k=key(point,normal);self.queries+=1;self.stats['requested']+=1
                hit=cached_sections.get(k)
                if hit and self.kw['max_candidate_faces']==100000 and self.kw['max_work']==10000000 and hit[0]['source_arrays_sha256']==self.array_sha:
                    self.stats['cached']+=1;r=copy.deepcopy(hit[0]);r['reused_exact_Source_section_receipt']=hit[1];r['geometry_reexecuted_for_this_call']=False;r['fresh_complete_original_Source_plane_query']=False;return r
                if self.engine is None:self.engine=ExactQuerySectionMeasurer(self.mesh,**self.kw);self.stats['new_index_built']=True
                self.stats['fresh']+=1;r=self.engine.section(point,normal);self.stats['new_logical_work']+=r.get('index_evidence',{}).get('counted_index_work',0)+r.get('contour_pair_comparison_upper_bound',0);r['geometry_reexecuted_for_this_call']=True;return r
        return CachedSections
    def voxel(v,f,pitch,*args,**kwargs):
        if cache_data is None:return automatic.voxelize_closed_chains(v,f,pitch,*args,**kwargs)
        assert pitch==cache_data['pitch'] and hashlib.sha256(v.tobytes()+f.tobytes()).hexdigest()==cache_data['voxelization']['source_arrays_sha256']
        return cache_data['cells'],copy.deepcopy(cache_data['voxelization'])
    def skeleton(solid,**kwargs):
        if cache_data is None:return measure_span_tube(solid,**kwargs)
        assert solid.pitch==cache_data['pitch'] and np.array_equal(solid.cells,cache_data['cells'])
        return copy.deepcopy(cache_data['kernel'])
    auto_fn=bound_function(automatic.measure_model,measure_tube_solid=skeleton,voxelize_closed_chains=voxel,ExactQuerySectionMeasurer=factory('automatic_node_probe'))
    def model(*args,**kwargs):
        r=auto_fn(*args,**kwargs);r['integration_backend']=dict(name='fixed_existing_Source_span_current_neighborhood_Lee',original_automatic_codeobject_exact=True,original_module_globals_modified=False,stage_cache=cache_proof);return r
    edge_fn=bound_function(source_profile.original_mesh_edge_profile,ExactQuerySectionMeasurer=factory('automatic_edge_profile'))
    def edges(*args,**kwargs):
        start=len(sessions);r=edge_fn(*args,**kwargs);stats=sessions[start:];r['requested_section_stations']=r['actual_fresh_whole_Source_queries'];r['actual_fresh_whole_Source_queries']=sum(s['fresh'] for s in stats);r['reused_exact_Source_section_receipts']=sum(s['cached'] for s in stats);r['Source_index_execution_stats']=copy.deepcopy(stats)
        j=0
        for b in r['index_batches']:
            if b['new_complete_Source_index_instantiated']:
                s=stats[j];j+=1;b['actual_fresh_whole_Source_queries']=s['fresh'];b['reused_exact_Source_receipts']=s['cached'];b['new_complete_Source_index_instantiated']=s['new_index_built']
        return r
    namespace=dict(measure_model=model,original_mesh_edge_profile=edges)
    if output_guard is not None:namespace['_output']=output_guard
    if input_loader is not None:namespace['load_model_input']=input_loader
    result=bound_function(entry.main,**namespace)(argv)
    # Runtime counters are returned separately; the produced measurement is immutable.
    return result,dict(original_entry_codeobject_exact=True,fixed_span_backend=True,stage_cache=cache_proof,Source_queries=sum(s['fresh'] for s in sessions),exact_Source_receipts_reused=sum(s['cached'] for s in sessions),index_execution_stats=sessions,original_module_globals_modified=False)
