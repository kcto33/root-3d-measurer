"""Supervised mesh compatibility path, preserving original triangle geometry."""
from collections import Counter
from fractions import Fraction as F
from types import FunctionType
import hashlib
import json
import math
import time
from pathlib import Path
import numpy as np

from root_triangle_measurement_input import read_triangle_ply_bytes
import root_closed_triangle_voxelization as voxel
from root_reconstruction_measurement import ReconstructedSolid
from root_graph_measurement_ledger import graph_measurement_ledger
from root_mesh_resolution import resolution_schedule,qualify_resolution
from root_model_delivery_exports import file_reference,write_json


def closed_contact_chains(vertices,faces,groups):
    """Oriented geometric chain closure permits zero-width edge contacts.

    This is not manifold repair or nonintersection/physical-solid admission.
    Original IDs and triangles remain untouched, including contact surfaces.
    """
    try:return voxel._validate(vertices,faces,groups)
    except ValueError as exc:
        if str(exc)!='GEOMETRIC_ORIENTED_EDGE_PAIR_CLOSURE_REQUIRED':raise
    f=np.asarray(faces);groups=[(0,len(f))] if groups is None else groups
    start=0;proof=[]
    for begin,end in groups:
        if begin!=start or not begin<end<=len(f):raise ValueError('EXACT_FACE_PARTITION_REQUIRED')
        count=Counter();net=Counter()
        for row in f[begin:end]:
            a,b,c=[vertices[int(i)] for i in row]
            u=tuple(x-y for x,y in zip(b,a));v=tuple(x-y for x,y in zip(c,a))
            if not any((u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0])):raise ValueError('DEGENERATE_FACE')
            for x,y in [(a,b),(b,c),(c,a)]:
                key=(x,y) if x<y else (y,x);count[key]+=1;net[key]+=1 if x<y else -1
        if any(net.values()) or any(n%2 for n in count.values()):raise ValueError('NONZERO_GEOMETRIC_CHAIN_BOUNDARY')
        proof.append(dict(face_range=[begin,end],geometric_oriented_boundary=0,
                          contact_geometric_edges=sum(n>2 for n in count.values()),
                          source_IDs_merged=False,source_manifold_or_nonoverlap_certified=False))
        start=end
    if start!=len(f):raise ValueError('EXACT_FACE_PARTITION_REQUIRED')
    return f,list(groups),proof


def sample_mesh(mesh,pitch,origin,*,max_output_cells=5000000):
    shift=np.asarray(origin,float)
    # Translation is solely the sampling lattice frame. All final measurement
    # and display use original source arrays/coordinates.
    fn=FunctionType(voxel.voxelize_closed_chains.__code__,dict(voxel.voxelize_closed_chains.__globals__,_validate=closed_contact_chains),argdefs=voxel.voxelize_closed_chains.__defaults__)
    fn.__kwdefaults__=voxel.voxelize_closed_chains.__kwdefaults__
    cells,receipt=fn(mesh.vertices,mesh.faces,pitch,max_output_cells=max_output_cells,max_abs_winding=1,grid_origin=shift)
    receipt.update(grid_origin=shift.tolist(),original_source_sha256=mesh.source_sha256,
                   source_coordinate_translation_for_sampling_only=True)
    return cells,receipt


class LazyRationalVertices:
    def __init__(self,vertices):self.vertices=vertices;self.cache={}
    def __getitem__(self,index):
        index=int(index)
        if index not in self.cache:self.cache[index]=tuple(F(float(x)) for x in self.vertices[index])
        return self.cache[index]


class PersistentTriangleSections:
    """One conservative complete-plane index; original exact disk predicates."""
    def __init__(self,mesh,*,max_queries,max_candidate_faces=100000):
        from root_exact_query_sections import _exact_indexed_kind
        if not 0<max_queries<=250000:raise ValueError('FIXED_MESH_QUERY_BUDGET_REQUIRED')
        self.mesh=mesh;self.max_queries=max_queries;self.max_faces=max_candidate_faces;self.queries=0
        kind=_exact_indexed_kind();self.kernel=kind.__new__(kind)
        self.kernel.mesh=mesh;self.kernel.vertices=LazyRationalVertices(mesh.vertices)
        self.kernel.faces=mesh.faces;self.kernel.max_segments=10000
        self.before=hashlib.sha256(mesh.vertices.tobytes()+mesh.faces.tobytes()).hexdigest()

    def section(self,point,normal):
        from root_exact_query_sections import _exact_query
        from root_local_triangle_section import _selection
        p=np.asarray(point,float);n=np.asarray(normal,float);ep,en=_exact_query(point,normal)
        if self.queries>=self.max_queries:raise ValueError('MESH_QUERY_BUDGET_EXHAUSTED')
        self.queries+=1;vertices=self.mesh.vertices
        with np.errstate(over='ignore',invalid='ignore'):
            ds=sum((vertices[:,i]-p[i])*n[i] for i in range(3))
            pad=64*np.finfo(float).eps*sum((np.abs(vertices[:,i])+abs(p[i]))*abs(n[i]) for i in range(3))
        positive=ds>pad;negative=ds<-pad;f=self.mesh.faces
        ids=np.flatnonzero(~(positive[f].all(axis=1)|negative[f].all(axis=1)))
        evidence=dict(source_faces=len(f),candidate_faces=len(ids),all_original_faces_classified=True,
                      uncertain_coplanar_tangent_faces_retained=True,infinite_plane=True,source_geometry_crop=False,
                      persistent_index=True,query_index=self.queries,source_sha256=self.mesh.source_sha256)
        if len(ids)>self.max_faces:return dict(radius=None,status='MESH_SECTION_CANDIDATE_BUDGET',index_evidence=evidence)
        distances=[None]*len(vertices)
        for i in np.unique(f[ids]):distances[int(i)]=sum(x*(y-z) for x,y,z in zip(en,self.kernel.vertices[i],ep))
        k=self.kernel;k._exact_indexed_face_ids=ids.tolist();k._exact_indexed_distances=distances;k.faces=f[ids]
        self.selected=None
        def select(instance,loops,qp,qn):
            work=sum(map(len,loops))
            if work>2500 or any(len(loop)>256 for loop in loops):return dict(path=None,status='MESH_SECTION_CONTOUR_WORK_BUDGET')
            result=_selection(instance,loops,qp,qn);self.selected=result['path'];return result
        from types import MethodType
        k._select_local_disk=MethodType(select,k)
        out=k.section(point,normal);out['index_evidence']=evidence
        out['original_triangle_radius_qualified']=out.get('radius') is not None
        if out.get('radius') is not None and self.selected:
            # Exact same-plane polygon first moments, no moved query plane.
            drop=max(range(3),key=lambda i:abs(en[i]));keep=[i for i in range(3) if i!=drop]
            xy=[tuple(q[i] for i in keep) for q in self.selected];area=F(0);moment=[F(0),F(0)]
            for a,b in zip(xy,xy[1:]+xy[:1]):
                w=a[0]*b[1]-a[1]*b[0];area+=w
                for j in range(2):moment[j]+=(a[j]+b[j])*w
            center=[F(0)]*3
            for i,v in zip(keep,[x/(3*area) for x in moment]):center[i]=v
            center[drop]=ep[drop]-sum(en[i]*(center[i]-ep[i]) for i in keep)/en[drop]
            out['centroid']=list(map(float,center))
        return out

    def unchanged(self):
        return hashlib.sha256(self.mesh.vertices.tobytes()+self.mesh.faces.tobytes()).hexdigest()==self.before


def segment_surface_hits(mesh,start,end):
    """Exact rational intersections; conservative finite segment AABB index."""
    p=np.asarray(start,float);d=np.asarray(end,float)-p
    triangles=mesh.vertices[mesh.faces]
    low=np.minimum(p,p+d);high=np.maximum(p,p+d)
    # Outward rounding keeps every possibly intersecting original triangle.
    ids=np.flatnonzero(np.all(triangles.max(axis=1)>=np.nextafter(low,-np.inf),axis=1)&
                       np.all(triangles.min(axis=1)<=np.nextafter(high,np.inf),axis=1))
    if len(ids)>100000:raise ValueError('MESH_SEGMENT_CANDIDATE_BUDGET')
    p=tuple(F(float(x)) for x in start);last=tuple(F(float(x)) for x in end);d=tuple(b-a for a,b in zip(p,last))
    sub=lambda a,b:tuple(x-y for x,y in zip(a,b))
    cross=lambda a,b:(a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
    dot=lambda a,b:sum(x*y for x,y in zip(a,b))
    hits={}
    for i in ids:
        a,b,c=[tuple(F(float(x)) for x in row) for row in triangles[i]]
        e1=sub(b,a);e2=sub(c,a);pv=cross(d,e2);det=dot(e1,pv);tv=sub(p,a)
        if det==0:
            if dot(cross(e1,e2),tv)==0:raise ValueError('MESH_SEGMENT_COPLANAR_CONTACT')
            continue
        u=dot(tv,pv)/det;qv=cross(tv,e1);v=dot(d,qv)/det
        if not (u>=0 and v>=0 and u+v<=1):continue
        t=dot(e2,qv)/det
        if 0<=t<=1:hits.setdefault(t,[]).append(int(i))
    return [dict(t=float(t),exact_t=[str(t.numerator),str(t.denominator)],position=[float(x+t*y) for x,y in zip(p,d)],source_face_ids=ids)
            for t,ids in sorted(hits.items())]


def original_mesh_endpoints(mesh,graph,pitch,engine):
    from root_branch_section_measurement import BranchOwnership
    from root_native_entity_measurement import _tangent
    from root_local_triangle_section import _contains
    from root_fragment_triangle_section import _boundary
    points=np.asarray(graph['points'],float);ownership=BranchOwnership(points,graph['adjacency']);records=[];changes={};routes={}
    for branch in ownership.ledger['branches']:
        if branch['closed']:continue
        ids=branch['source_node_ids']
        for oriented in [ids,ids[::-1]]:
            node=oriented[0]
            if len(graph['adjacency'][node])!=1:continue
            path=points[oriented];arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(path,axis=0),axis=1))]
            inward=_tangent(path,arc,0,pitch,False)
            r=dict(node_id=node,status='endpoint_unresolved',original_position=points[node].tolist(),physical_tip_qualified=False)
            try:
                if inward is None:raise ValueError('MESH_NO_UNIQUE_TERMINAL_TANGENT')
                n=-inward;cut=engine.section(points[node],n);contour=engine.selected
                if cut.get('radius') is None or not contour or cut.get('centroid') is None:raise ValueError('MESH_NO_UNIQUE_TERMINAL_DISK')
                center=np.asarray(cut['centroid']);drop=int(np.argmax(np.abs(n)));keep=[i for i in range(3) if i!=drop]
                xy=[tuple(v[i] for i in keep) for v in contour];q=tuple(F(float(center[i])) for i in keep)
                if not _contains(xy,q) or _boundary(xy,q):raise ValueError('MESH_TERMINAL_CENTROID_NOT_INTERIOR')
                witnesses=[center]+[np.asarray(list(map(float,v))) for v in contour]
                if any(ownership.nearest_owners(w)!={branch['branch_id']} for w in witnesses):raise ValueError('MESH_MIXED_TERMINAL_DISK')
                route=[center.copy()];initial_radius=cut['radius'];travel=0.;trace_closed=False
                for step_number in range(256):
                    step=min(pitch,max(pitch/16,cut['radius']/2));proposal=center+step*n
                    value=engine.section(proposal,n)
                    if value.get('radius') is None:
                        if value['status'] not in ['NO_POSITIVE_TRIANGLE_SECTION','MISSING_UNIQUE_QUERY_CONTAINING_DISK']:
                            raise ValueError('MESH_TERMINAL_TRACE_'+value['status'])
                        trace_closed=True;break
                    proposed=np.asarray(value['centroid']);contour=engine.selected
                    if not contour or any(ownership.nearest_owners(np.asarray(list(map(float,w))))!={branch['branch_id']} for w in contour):
                        raise ValueError('MESH_TERMINAL_TRACE_MIXED_BRANCH')
                    drop=int(np.argmax(np.abs(n)));keep=[i for i in range(3) if i!=drop]
                    xy=[tuple(v[i] for i in keep) for v in contour];qp=tuple(F(float(proposed[i])) for i in keep)
                    if not _contains(xy,qp) or _boundary(xy,qp):raise ValueError('MESH_TERMINAL_TRACE_CENTROID_NOT_INTERIOR')
                    delta=proposed-center;distance=math.hypot(*delta)
                    if delta@n<step/4 or distance>max(2*pitch,cut['radius']):raise ValueError('MESH_TERMINAL_TRACE_NONLOCAL_OR_BACKWARD')
                    if segment_surface_hits(mesh,center,proposed):raise ValueError('MESH_TERMINAL_TRACE_OUTSIDE_SURFACE')
                    route.append(proposed);center=proposed;travel+=distance;cut=value
                    if travel>max(50*pitch,20*initial_radius):raise ValueError('MESH_TERMINAL_TRACE_LENGTH_GUARD')
                    if len(route)>1:
                        new=route[-1]-route[max(0,len(route)-5)];norm=math.hypot(*new)
                        if norm:
                            new/=norm
                            if new@n<.5:raise ValueError('MESH_TERMINAL_TRACE_REVERSAL')
                            n=new
                if not trace_closed:raise ValueError('MESH_TERMINAL_TRACE_STEP_GUARD')
                u=np.cross(n,np.eye(3)[int(np.argmin(np.abs(n)))]);u/=math.hypot(*u);v=np.cross(n,u)
                directions=[n]+[math.cos(math.pi/18)*n+s*math.sin(math.pi/18)*a for a in [u,v] for s in [-1,1]]
                exits=[]
                for direction in directions:
                    hits=segment_surface_hits(mesh,center,center+4*pitch*direction)
                    if not hits or hits[0]['t']<=0:raise ValueError('MESH_NO_LOCAL_ORIGINAL_BOUNDARY_EXIT')
                    first=hits[0]
                    if ownership.nearest_owners(first['position'])!={branch['branch_id']}:raise ValueError('MESH_FOREIGN_TERMINAL_EXIT')
                    exits.append(first)
                spread=max(math.dist(a['position'],b['position']) for a in exits for b in exits)
                if spread>pitch:raise ValueError('MESH_TERMINAL_DIRECTION_UNSTABLE')
                neighbor=oriented[1];end=exits[0]['position']
                neighbor_cut=engine.section(points[neighbor],n)
                if neighbor_cut.get('radius') is None:raise ValueError('MESH_ENDPOINT_NEIGHBOR_NOT_INTERIOR')
                changed_route=[points[neighbor]]+route+[np.asarray(end)]
                for at,(a,b) in enumerate(zip(changed_route,changed_route[1:])):
                    crossings=segment_surface_hits(mesh,a,b)
                    if any(hit['t']<1-64*np.finfo(float).eps or at+2<len(changed_route) for hit in crossings):raise ValueError('MESH_CHANGED_EDGE_CROSSES_ORIGINAL_BOUNDARY')
                changes[node]=end;routes[node]=route+[np.asarray(end)];r.update(status='geometric_surface_estimate',position=end,endpoint_method='original_triangle_owned_disk_trace_and_exact_local_ray',
                    source_sha256=mesh.source_sha256,section_centroid=route[0].tolist(),direction=n.tolist(),directional_exits=exits,boundary_spread=spread,
                    supported_trace_length=travel,trace_points=len(route),
                    original_surface_unchanged=True,precision_certificate=False)
            except (ValueError,RuntimeError,AssertionError) as exc:r['reason']=str(exc)
            records.append(r)
    changed=points.copy()
    for node,q in changes.items():changed[node]=q
    # If both ends share one edge, qualify their joint changed edge as well.
    for a,b in ownership.ledger['source_edge_ids']:
        if a in changes and b in changes:
            try:
                if any(64*np.finfo(float).eps<h['t']<1-64*np.finfo(float).eps for h in segment_surface_hits(mesh,changed[a],changed[b])):
                    raise ValueError('MESH_BOTH_ENDS_CHANGED_EDGE_CROSSES_BOUNDARY')
            except ValueError as exc:
                for node in [a,b]:
                    changed[node]=points[node]
                    next(r for r in records if r['node_id']==node).update(status='endpoint_unresolved',reason=str(exc))
                    routes.pop(node,None)
    output=changed.tolist();adj=[[] for unused in output]
    for a,b in ownership.ledger['source_edge_ids']:
        interior=[]
        if a in routes:interior.extend(reversed(routes[a][:-1]))
        if b in routes:interior.extend(routes[b][:-1])
        chain=[a]
        for q in interior:
            if math.dist(output[chain[-1]],q)<=1e-12:continue
            chain.append(len(output));output.append(list(map(float,q)));adj.append([])
        chain.append(b)
        for u,v in zip(chain,chain[1:]):adj[u].append(v);adj[v].append(u)
    return dict(points=output,adjacency=adj),records


def proxy_axis(cells,pitch,origin,progress=None):
    from root_tube_measurement import measure_tube_solid
    from root_span_model_pipeline import measure_span_tube
    from root_cubical_graph_quotient import quotient_graph
    from root_mesh_body_index import indexed_mesh_body as indexed_body_measurement
    solid=ReconstructedSolid(cells,pitch,{})
    extent=math.prod(map(int,cells.max(axis=0)-cells.min(axis=0)+5))
    if extent>350000000:
        from root_mesh_thinning_storage import guarded_lee_sparse_audit,byte_plan
        low=cells.min(axis=0)-1;shape=cells.max(axis=0)-low+2
        plan=byte_plan(shape,len(cells),8192*1024**2)
        if not plan['admitted']:raise ValueError('MESH_LEE_BYTE_PLAN_EXCEEDS_BUDGET')
        def storage_adapter(volume,**kwargs):
            return guarded_lee_sparse_audit(volume,source_cells=cells-low,**kwargs)
        fn=FunctionType(measure_tube_solid.__code__,dict(measure_tube_solid.__globals__,guarded_skeletonize=storage_adapter),argdefs=measure_tube_solid.__defaults__)
        fn.__kwdefaults__=measure_tube_solid.__kwdefaults__
        kernel=fn(solid,max_dense_cells=math.prod(map(int,shape+2)))
    else:
        kernel=(measure_span_tube if extent<=1000000 else measure_tube_solid)(solid,max_dense_cells=350000000)
    if kernel['representation']['status']!='PASS_SUPPORTED_GEOMETRIC_GRAPH_TOPOLOGY':
        p,adj,receipt=quotient_graph(np.asarray(kernel['source_skeleton_cells'],np.int64),pitch)
        if receipt['skeleton_voxel_H2']:raise ValueError('MESH_PROXY_CAVITIES_NOT_A_PURE_GRAPH')
        kernel.update(points=p.tolist(),adjacency=adj)
    body=indexed_body_measurement(solid,kernel['points'],kernel['adjacency'],progress=progress)
    graph=dict(points=(np.asarray(body['points'])+np.asarray(origin)).tolist(),adjacency=body['adjacency'])
    if not graph['points']:raise ValueError('EMPTY_MESH_AXIS')
    return graph


def measure_mesh(mesh_path,output,*,target_pitch,target_reason,coordinate_units=None,scale_mm=None,progress=None,resume_target_from=None):
    from root_native_entity_measurement import subdivide_measurement_graph,_tangent
    start=time.monotonic();source=file_reference(mesh_path)
    if source['bytes']>268435456:raise ValueError('MESH_SOURCE_BYTE_BUDGET')
    mesh=read_triangle_ply_bytes(Path(mesh_path).read_bytes(),source['sha256'],max_vertices=1200000,max_faces=2000000)
    def emit(phase,**data):
        if progress:progress(phase,elapsed_seconds=time.monotonic()-start,**data)
    graph=None;rows=[];failures=[];runs=[];qualification=None;effective=target_pitch;ends=[]
    try:
        cache={};reuse=None
        if resume_target_from is not None:
            previous=json.loads(Path(resume_target_from).read_text(encoding='utf-8'))
            if previous['schema_version']!='root-mesh-delivery-measurement-v001' or previous['input_file_binding']['sha256']!=mesh.source_sha256:
                raise ValueError('MESH_TARGET_CACHE_SOURCE_MISMATCH')
            if previous['stages']['declared_target_pitch']!=target_pitch or previous['radius_curve']['rows']:
                raise ValueError('MESH_TARGET_CACHE_MUST_BE_PRE_ENDPOINT_PRE_RADIUS')
            previous_prefix=Path(resume_target_from).with_name(Path(resume_target_from).name[:-len('_measurement.json')]).resolve()
            experiment=next(p for p in previous_prefix.parents if (p/'EXPERIMENT_LOG.md').is_file())
            delivery=json.loads(previous_prefix.with_name(previous_prefix.name+'_delivery.json').read_text(encoding='utf-8'))
            if delivery['measurement']['sha256']!=file_reference(resume_target_from)['sha256']:raise ValueError('MESH_TARGET_CACHE_MEASUREMENT_CHANGED')
            execution_key=hashlib.sha256(str(previous_prefix).encode('utf-8')).hexdigest()[:12]
            execution=json.loads((experiment/'snapshot/native_model_delivery_run001/executions'/(previous_prefix.name+'_'+execution_key+'.json')).read_text(encoding='utf-8'))
            producer=next(r for r in execution['source_files'] if r['path'].endswith('src/root_mesh_delivery_measurement.py'))
            archived=experiment/'snapshot/native_model_delivery_run001/implementation_run009/src/root_mesh_delivery_measurement.py'
            if file_reference(archived)['sha256']!=producer['sha256']:raise ValueError('MESH_TARGET_CACHE_UNQUALIFIED_PRODUCER')
            if any(r['status']!='endpoint_unresolved' or r['original_position']!=previous['candidate_centerline']['points'][r['node_id']] for r in previous['geometric_endpoints']):
                raise ValueError('MESH_TARGET_CACHE_ENDPOINTS_ALREADY_CHANGED')
            old=next(r for r in previous['stages']['runs'] if r['name']=='target' and r['pitch']==target_pitch)
            if old['origin']!=[0.,0.,0.] or previous['provenance'].get('input_graph_or_GT_used') is not False:
                raise ValueError('MESH_TARGET_CACHE_FRAME_OR_INFERENCE_MISMATCH')
            cache[(target_pitch,(0.,0.,0.))]=(previous['candidate_centerline'],old['voxelization'])
            reuse=file_reference(resume_target_from)
        for attempt in range(3):
            graphs={}
            for level in resolution_schedule(effective,target_reason):
                emit('mesh_resolution',attempt=attempt+1,level=level['name'],pitch=level['pitch'])
                key=(level['pitch'],tuple(level['origin']))
                if key not in cache:
                    cells,receipt=sample_mesh(mesh,level['pitch'],level['origin'])
                    sampled_axis=proxy_axis(cells,level['pitch'],level['origin'],progress=lambda value:emit('mesh_axis_body',**value)) if len(cells) else dict(points=[],adjacency=[])
                    cache[key]=(sampled_axis,receipt)
                graphs[level['name']]=cache[key][0]
                if level['name']=='target':graph=graphs['target']
                runs.append(dict(attempt=attempt+1,**level,voxelization=cache[key][1]))
            graph=graphs['target'];qualification=qualify_resolution(graphs,effective)
            if qualification['passed']:break
            effective/=2
        if not qualification['passed']:raise ValueError('RESOLUTION_CHALLENGE_FAILED')
        terminal_count=sum(len(a)==1 for a in graph['adjacency'])
        endpoint_engine=PersistentTriangleSections(mesh,max_queries=max(1,terminal_count*258))
        graph,ends=original_mesh_endpoints(mesh,graph,effective,endpoint_engine)
        graph=subdivide_measurement_graph(graph,effective)
        # All original-surface stations/normals fixed before radius outcomes.
        ledger=graph_measurement_ledger(graph['points'],graph['adjacency'],coordinate_units='model_units')
        points=np.asarray(graph['points'],float);stations=[]
        for branch in ledger['branches']:
            ids=branch['source_node_ids'];path=points[ids];arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(path,axis=0),axis=1))]
            for at,(a,b) in enumerate(zip(ids,ids[1:])):
                n=_tangent(path,arc,(arc[at]+arc[at+1])/2,effective,branch['closed'])
                stations.append(dict(node_ids=[int(a),int(b)],branch_id=branch['branch_id'],position=((points[a]+points[b])/2).tolist(),normal=n.tolist() if n is not None else None))
        from root_branch_section_measurement import BranchOwnership
        ownership=BranchOwnership(graph['points'],graph['adjacency'])
        engine=PersistentTriangleSections(mesh,max_queries=max(1,len(stations)))
        for i,s in enumerate(stations):
            value=engine.section(s['position'],s['normal']) if s['normal'] is not None else dict(radius=None,status='MISSING_UNIQUE_BRANCH_NORMAL')
            if value.get('radius') is not None:
                witnesses=[np.asarray(list(map(float,q))) for q in engine.selected]
                if any(ownership.nearest_owners(q)!={s['branch_id']} for q in witnesses):
                    value=dict(value,radius=None,status='MISSING_MIXED_BRANCH_SECTION')
            rows.append(dict(s,radius=value.get('radius'),area=value.get('area'),status=value['status'],available=value.get('radius') is not None,
                             section_centroid=value.get('centroid'),original_surface_sha256=mesh.source_sha256))
            if i%50==0 or i+1==len(stations):emit('radius',processed=i+1,total=len(stations))
        assert engine.unchanged()
    except (ValueError,RuntimeError,AssertionError,MemoryError,OverflowError) as exc:
        failures.append(dict(code='MESH_PIPELINE_BLOCKED',reason=str(exc),position=None));emit('pipeline_blocked',reason=str(exc))
    if graph is not None and not graph['points']:graph=None
    ledger=graph_measurement_ledger(graph['points'],graph['adjacency'],coordinate_units='model_units') if graph else None
    if ledger:
        # Original finite-end surface qualification is explicit, never inherited
        # from the sampled proxy's end coordinates.
        resolved={r['node_id'] for r in ends if r['status']=='geometric_surface_estimate'}
        for node in ledger['endpoint_node_ids']:
            if not any(r['node_id']==node for r in ends):ends.append(dict(node_id=node,status='endpoint_unresolved',original_position=graph['points'][node],reason='ORIGINAL_MESH_TERMINAL_AXIS_NOT_QUALIFIED',physical_tip_qualified=False))
            if node not in resolved:failures.append(dict(code='GEOMETRIC_ENDPOINT_UNRESOLVED',node_id=node,position=graph['points'][node]))
    failures.extend(dict(code=r['status'],position=r['position'],edge_nodes=r['node_ids']) for r in rows if not r['available'])
    assert file_reference(mesh_path)==source
    result=dict(schema_version='root-mesh-delivery-measurement-v001',status='MESH_MEASUREMENT_BLOCKED' if any(f['code']=='MESH_PIPELINE_BLOCKED' for f in failures) else 'MESH_OUTPUT_WITH_EXPLICIT_GAPS',
        case_label=Path(mesh_path).parent.name+' / '+Path(mesh_path).stem,input_file_binding=dict(source,format='triangle_ply'),
        coordinate_units=coordinate_units or 'model_units',scale_mm_per_model_unit=scale_mm,candidate_centerline=graph,
        radius_curve=dict(sampling='edge_midpoint',rows=rows,definition='sqrt(area/pi) on unchanged original triangle surface',stations_fixed_before_radius_outcomes=True),
        length=None if not ledger else dict(value=ledger['measurement_graph_length'],definition='sum of unique final graph edges',endpoint_status=('geometric_surface_estimates' if all(r['status']=='geometric_surface_estimate' for r in ends) else 'endpoint_unresolved') if ends else 'not_applicable_closed_graph',complete_endpoint_length_qualified=False),
        geometric_endpoints=ends,failures=failures,model_bounds=dict(minimum=mesh.vertices.min(axis=0).tolist(),maximum=mesh.vertices.max(axis=0).tolist()),
        provenance=dict(mode='reused_mesh_target_new_challenges' if resume_target_from else 'fresh',input_kind='triangle_ply',mesh_sampling_is_proxy=True,input_graph_or_GT_used=False,
                        reused_target=file_reference(resume_target_from) if resume_target_from else None),
        stages=dict(resolution=qualification,target_reason=target_reason,declared_target_pitch=target_pitch,runs=runs),
        functional_acceptance=False,precision_acceptance=False,worker_seconds=time.monotonic()-start,original_entity_unchanged=True)
    write_json(output,result)
    return result,mesh


def write_mesh_delivery(report_path,prefix,experiment,mesh):
    from root_model_delivery_exports import export_tables
    prefix=Path(prefix);experiment=Path(experiment);report=json.loads(Path(report_path).read_text(encoding='utf-8'))
    exports,paths,coverage=export_tables(report,prefix)
    folder=experiment/'results/model_delivery_local/geometry'/('mesh_'+mesh.source_sha256[:24])
    manifest=folder/'surface.json'
    if manifest.is_file():
        surface=json.loads(manifest.read_text(encoding='utf-8'))
        if surface['source_mesh_sha256']!=mesh.source_sha256:raise ValueError('MESH_DISPLAY_BINDING_CHANGED')
        for ref in surface['buffers'].values():
            actual=file_reference(folder/ref['path'],folder)
            if actual!=ref:raise ValueError('MESH_DISPLAY_BUFFER_CHANGED')
    else:
        folder.mkdir(parents=True,exist_ok=True);buffers={}
        for name,array in [('positions',mesh.vertices.astype('<f4')),('indices',mesh.faces.astype('<u4').reshape(-1))]:
            path=folder/(name+'.bin')
            with path.open('xb') as stream:array.tofile(stream)
            buffers[name]=file_reference(path,folder)
        surface=dict(schema_version='original-triangle-display-v001',source_mesh_sha256=mesh.source_sha256,
                     vertices=len(mesh.vertices),triangles=len(mesh.faces),all_original_faces_displayed=True,
                     display_float32_only=True,measurement_geometry_modified=False,buffers=buffers)
        write_json(manifest,surface)
    surface=dict(surface,buffers={k:file_reference(folder/ref['path'],prefix.parent) for k,ref in surface['buffers'].items()})
    from root_model_parameters import write_parameters
    parameters=write_parameters(report,prefix,mesh=mesh)
    result=dict(schema_version='root-native-delivery-v001',input_kind='triangle_ply',measurement=file_reference(report_path,prefix.parent),
        case_label=report['case_label'],summary=dict(status=report['status'],candidate_length=report.get('length'),paths=len(paths),coverage=coverage,
        automatic_nodes=len((report.get('candidate_centerline') or {}).get('points',[]))),exports=exports,surface=surface,parameters=parameters,
        provenance=report['provenance'],functional_acceptance='NOT_YET_PASSED',precision_acceptance='NOT_YET_PASSED')
    delivery=prefix.with_name(prefix.name+'_delivery.json');write_json(delivery,result)
    viewer=prefix.with_name(prefix.name+'_viewer.html')
    with viewer.open('x',encoding='utf-8') as stream:stream.write(Path(__file__).with_name('model_delivery_viewer.html').read_text(encoding='utf-8'))
    return result


def run_mesh_worker(args):
    from root_local_network_policy import disable_measurement_network,network_policy_state
    disable_measurement_network()
    prefix=Path(args.output_prefix).resolve();experiment=next(p for p in prefix.parents if (p/'EXPERIMENT_LOG.md').is_file())
    progress_path=prefix.with_name(prefix.name+'_progress.jsonl')
    def progress(phase,**data):
        event=dict(phase=phase,**data)
        with progress_path.open('a',encoding='utf-8') as stream:stream.write(json.dumps(event)+'\n')
        print(json.dumps(event),flush=True)
    report=prefix.with_name(prefix.name+'_measurement.json')
    result,mesh=measure_mesh(args.mesh,report,target_pitch=args.mesh_target_pitch,target_reason=args.mesh_target_reason,
                           coordinate_units=args.coordinate_units,scale_mm=args.scale_mm,progress=progress,resume_target_from=args.mesh_resume_target_from)
    delivery=write_mesh_delivery(report,prefix,experiment,mesh)
    progress('complete',status=result['status'],network_policy=network_policy_state())
    return delivery
