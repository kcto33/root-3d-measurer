"""All automatically generated edge stations, explicit fresh whole-index batches."""
import hashlib,math,os
from fractions import Fraction as F
import numpy as np
from root_exact_query_sections import ExactQuerySectionMeasurer
from root_graph_source_length import pack
from root_labelled_path_cover import labelled_path_cover

def original_mesh_edge_profile(mesh,points,adjacency,*,execution_domain='local',max_work=10000000):
    if type(max_work)is not int or not 0<=max_work<=10000000:raise ValueError('AGGREGATE_WORK_GUARD_ONLY_TIGHTENED')
    if len(points)!=len(adjacency) or not 0<len(points)<=1200000:raise ValueError('COMPLETE_AUTOMATIC_GRAPH_REQUIRED')
    if execution_domain=='local' and len(points)>256:raise ValueError('UNCHANGED_LOCAL_256_NODE_GUARD')
    if execution_domain!='local' and (execution_domain!='registered_slurm_complete_source' or not os.environ.get('SLURM_JOB_ID','').isdigit()):raise ValueError('REGISTERED_COMPLETE_SOURCE_REQUIRED')
    pairs=sorted((i,j) for i,row in enumerate(adjacency) for j in row if i<j)
    if len(pairs)>1000:raise ValueError('EXPLICIT_1000_STATION_GUARD_NO_PARTIAL_PROFILE')
    cover=labelled_path_cover(adjacency,expected_nodes=len(points),expected_edges=len(pairs),execution_domain='local' if execution_domain=='local' else 'registered_slurm_complete_graph')
    before=hashlib.sha256(mesh.vertices.tobytes()+mesh.faces.tobytes()).hexdigest()
    xyz=[tuple(F(float(x)) for x in p) for p in points]
    if any(len(p)!=3 for p in xyz):raise ValueError('THREE_AUTOMATIC_COORDINATES_REQUIRED')
    # Freeze every station and its batch before measuring or inspecting success.
    stations=[]
    for eid,(a,b) in enumerate(pairs):
        p=tuple((x+y)/2 for x,y in zip(xyz[a],xyz[b]));n=tuple(y-x for x,y in zip(xyz[a],xyz[b]))
        stations.append(dict(generated_edge_id=eid,generated_node_ids=[a,b],batch_id=eid//5,exact_point=list(map(pack,p)),exact_normal=list(map(pack,n)),position=list(map(float,p)),query_representable=any(n) and all(max(abs(v.numerator).bit_length(),v.denominator.bit_length())<=4096 for v in p+n)))
    rows=[];batches=[];work=0;exhausted=False;nv,nf=len(mesh.vertices),len(mesh.faces)
    for start in range(0,len(stations),5):
        selected=stations[start:start+5];engine=None;batch_work=0;actual=0
        for station in selected:
            if not station['query_representable']:
                result=dict(status='MISSING_ZERO_OR_UNREPRESENTABLE_AUTOMATIC_EDGE_NORMAL',radius=None,area=None)
            elif exhausted or work+3*nv+nf>max_work:
                exhausted=True;result=dict(status='MISSING_DECLARED_TOTAL_SOURCE_PROFILE_WORK_GUARD',radius=None,area=None)
            else:
                if engine is None:
                    engine=ExactQuerySectionMeasurer(mesh,execution_domain=execution_domain,max_queries=5)
                actual+=1
                try:result=engine.section(station['exact_point'],station['exact_normal'])
                except (ValueError,RuntimeError,AssertionError,OverflowError) as exc:result=dict(status='MISSING_ORIGINAL_TRIANGLE_SECTION_ENTRY',radius=None,area=None,error_type=type(exc).__name__,reason=str(exc))
                used=result.get('index_evidence',{}).get('counted_index_work',3*nv+nf)+result.get('contour_pair_comparison_upper_bound',0)
                work+=used;batch_work+=used
                if work>max_work:
                    exhausted=True;result=dict(status='MISSING_DECLARED_TOTAL_SOURCE_PROFILE_WORK_GUARD',radius=None,area=None,retained_overlimit_outcome=result)
            rows.append(dict(station,available=result.get('radius') is not None,area_equivalent_radius=result.get('radius'),section_area=result.get('area'),status=result['status'],original_mesh_section=result,physical_radius=None))
        batches.append(dict(batch_id=start//5,fixed_station_ids=[s['generated_edge_id'] for s in selected],query_denominator=len(selected),actual_fresh_whole_Source_queries=actual,new_complete_Source_index_instantiated=engine is not None,registered_index_max_queries=5,logical_query_work=batch_work))
    assert len(rows)==len(stations) and all(b['actual_fresh_whole_Source_queries']<=5 for b in batches)
    paths=[]
    for path in cover['paths']:
        arc=0.;samples=[]
        for a,b,eid in zip(path['original_node_ids'],path['original_node_ids'][1:],path['original_edge_ids']):
            edge_length=math.dist(points[a],points[b]);row=rows[eid]
            samples.append(dict(generated_edge_id=eid,arc_position=arc+edge_length/2,position=row['position'],radius=row['area_equivalent_radius'],status=row['status']));arc+=edge_length
        paths.append(dict(path_id=path['path_id'],component_id=path['component_id'],closed=path['closed'],sampled_length=arc,original_generated_node_ids=path['original_node_ids'],samples=samples))
    assert hashlib.sha256(mesh.vertices.tobytes()+mesh.faces.tobytes()).hexdigest()==before
    return dict(status='COMPLETE_DECLARED_AUTOMATIC_EDGE_STATION_OUTCOMES',station_denominator=len(rows),available_radii=sum(r['available'] for r in rows),rows=rows,paths=paths,index_batches=batches,
        actual_fresh_whole_Source_queries=sum(b['actual_fresh_whole_Source_queries'] for b in batches),complete_original_mesh_classification_per_query=True,total_logical_query_work=work,max_declared_work=max_work,
        all_stations_fixed_before_Source_outcomes=True,station_choice_uses_GT_or_success=False,all_generated_edges_owned_once=cover['all_original_edges_owned_exactly_once'],
        source_sha256=mesh.source_sha256,original_Source_array_SHA_before_after=before,geometry_or_Source_modified=False,radius_definition='sqrt(area/pi) of the original-triangle simple disk at each exact generated-edge midpoint and chord-normal plane',
        continuous_radius_between_stations=None,failed_stations_not_interpolated=True,physical_radius=None,actual_GT_accuracy=None)
