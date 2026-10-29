"""Portable native-entity input and local measurement, with no PLY round trip or GT."""
import hashlib
import json
import math
import os
import sys
import time
from fractions import Fraction
from itertools import product
from types import FunctionType
from pathlib import Path, PurePosixPath

import numpy as np

from root_graph_measurement_ledger import graph_measurement_ledger
from root_model_delivery_exports import file_reference, write_json
from root_reconstruction_measurement import load_reconstruction,ReconstructedSolid
from root_local_network_policy import disable_measurement_network,network_policy_state

PROJECT = Path(__file__).resolve().parents[1]
WORKSPACE = PROJECT.parents[1]
SCHEMA = 'root-native-entity-measurement-v001'


def digest_cells(solid):
    return hashlib.sha256(solid.cells.tobytes()).hexdigest()


def load_entity(path, *, stage='regularized_negative_indices'):
    metadata = Path(path).resolve()
    if metadata.is_dir():
        metadata = metadata/'entity_prediction.json'
    ref = file_reference(metadata)
    data = json.loads(metadata.read_text(encoding='utf-8-sig'))
    caches, crops = {}, []
    for field in data.get('fields', []):
        item = field['cache']
        name = PurePosixPath(item['path'].replace('\\', '/')).name
        local = metadata.parent/name
        if not local.is_file():
            raise ValueError('ENTITY_PACKAGE_MISSING_BOUND_CROP: '+name)
        actual = file_reference(local)
        if actual['sha256'] != item['sha256'] or actual['bytes'] != item['bytes']:
            raise ValueError('ENTITY_CROP_SHA_BYTES_MISMATCH: '+name)
        caches[item['path']] = (local, item['sha256'])
        crops.append(actual)
    solid = load_reconstruction(metadata, ref['sha256'], caches, stage=stage)
    bound = dict(metadata=ref, crops=crops, stage=stage, pitch=solid.pitch,
                 occupied_cells=len(solid.cells), cells_sha256=digest_cells(solid),
                 source_cell_count=data['all_source_negative_cells_accounted'],
                 declared_boundary=data.get('mesh'),
                 observation_weights_available=False, known_empty_exterior_available=False,
                 coordinates='canonical cell center=(integer_index+0.5)*pitch',
                 model_modified=False, ply_revoxelization_performed=False)
    return solid, bound, data


def check_entity_unchanged(binding):
    for ref in [binding['metadata']]+binding['crops']:
        current = file_reference(ref['path'])
        if current['sha256'] != ref['sha256'] or current['bytes'] != ref['bytes']:
            raise ValueError('ENTITY_CHANGED_DURING_MEASUREMENT')


def _bound_workspace_file(ref):
    rel = Path(ref['workspace_relative_path'])
    if rel.is_absolute() or '..' in rel.parts:
        raise ValueError('WORKSPACE_RELATIVE_REFERENCE_REQUIRED')
    path = WORKSPACE/rel
    current = file_reference(path)
    if current['sha256'] != ref['sha256'] or current['bytes'] != ref['bytes']:
        raise ValueError('INHERITED_SOURCE_SHA_BYTES_MISMATCH')
    return path


def inherit_bound_measurement(witness, binding):
    """M1 only: old section results remain at their old positions and normals."""
    witness = Path(witness).resolve()
    case = json.loads(witness.read_text(encoding='utf-8-sig'))['case']
    source = case['source_frame']
    if binding['stage'] != 'regularized_negative_indices':
        raise ValueError('INHERITANCE_REQUIRES_MATCHING_REGULARIZED_ENTITY_STAGE')
    if source['entity_prediction']['sha256'] != binding['metadata']['sha256']:
        raise ValueError('INHERITED_ENTITY_METADATA_MISMATCH')
    if source['pitch'] != binding['pitch']:
        raise ValueError('INHERITED_COORDINATE_FRAME_MISMATCH')
    body_path = _bound_workspace_file(source['body'])
    radius_path = _bound_workspace_file(source['radius_profile'])
    body = json.loads(body_path.read_text(encoding='utf-8-sig'))
    profile = json.loads(radius_path.read_text(encoding='utf-8-sig'))
    if body['points'] != profile['points']:
        raise ValueError('INHERITED_RADIUS_NOT_ON_FINAL_INHERITED_AXIS')
    if body.get('source_solid_modified') or body.get('GT_axis_or_radius_used'):
        raise ValueError('INHERITED_RESULT_NOT_MODEL_ONLY')
    graph = dict(points=body['points'], adjacency=body['adjacency'])
    rows = [dict(node_id=r['node_id'], position=r['position'], normal=r['normal'],
                 radius=r['new_radius'], available=r['new_radius'] is not None,
                 status=r['new_status']) for r in profile['rows']]
    radius = dict(sampling='node_halfedge', rows=rows, definition=profile['radius_definition'],
                  inherited_positions_and_planes_unchanged=True,
                  reference_coverage_not_claimed=True)
    provenance = dict(mode='inherited',
                      reason='M1 reuses SHA-bound existing full-entity body and radius results',
                      sources=[file_reference(witness),file_reference(body_path),file_reference(radius_path)],
                      fresh_automatic_measurement=False, new_geometry_queries=0,
                      stage=binding['stage'])
    return graph, radius, provenance


def _tangent(path, arc, position, pitch, closed):
    extent = arc[-1]
    if extent <= 0:
        return None
    if closed:
        width = min(10*pitch, extent)
        locations = np.mod(position+np.linspace(-width/2, width/2, 21), extent)
    else:
        locations = np.linspace(max(0.,position-5*pitch),min(extent,position+5*pitch),21)
    q = np.column_stack([np.interp(locations,arc,path[:,a]) for a in range(3)])
    q -= q[10]
    derivative = np.array([math.fsum(2*(i-10)/77*float(q[i,a]) for i in range(21))
                           for a in range(3)])
    noise=64*np.finfo(float).eps*max(float(np.abs(path).max()),pitch,1e-30)
    derivative[np.abs(derivative)<=noise]=0.
    norm = math.hypot(*derivative)
    return derivative/norm if norm > 0 and math.isfinite(norm) else None


def axis_aligned_section(solid, point, normal, slab_cache=None):
    """Same-plane square union, including grid-face cuts; no duplicate shared faces."""
    normal=np.asarray(normal,float)
    active=np.flatnonzero(normal)
    if len(active)!=1:
        return None
    axis=int(active[0]);other=[a for a in range(3) if a!=axis]
    h=Fraction(solid.pitch)
    coordinates=[Fraction(float(x))/h for x in point]
    level=coordinates[axis]
    lower=level.numerator//level.denominator
    levels=[lower-1,lower] if level.denominator==1 else [lower]
    if slab_cache is None:slab_cache={}
    if axis not in slab_cache:
        slabs={}
        for cell in solid.cells:
            slabs.setdefault(int(cell[axis]),set()).add(tuple(int(cell[a]) for a in other))
        slab_cache[axis]=slabs
    squares=set().union(*(slab_cache[axis].get(k,set()) for k in levels))
    q=[coordinates[a] for a in other]
    def incident(v):
        lo=v.numerator//v.denominator
        return [lo-1,lo] if v.denominator==1 else [lo]
    expected={(a,b) for a in incident(q[0]) for b in incident(q[1])}
    hits=expected&squares
    base=dict(radius=None,area=None,exact_axis_plane=True,plane_was_moved=False,
              overlapping_face_squares_deduplicated=True)
    if not hits:return dict(base,status='QUERY_OUTSIDE_SECTION')
    component=set();stack=[next(iter(hits))]
    while stack:
        at=stack.pop()
        if at in component:continue
        component.add(at)
        x,y=at
        stack.extend(v for v in [(x-1,y),(x+1,y),(x,y-1),(x,y+1)]
                     if v in squares and v not in component)
    if not hits<=component:
        return dict(base,status='QUERY_ON_AMBIGUOUS_SECTION_COMPONENT_CONTACT')
    if not expected<=component:
        return dict(base,status='QUERY_ON_SECTION_OUTER_BOUNDARY')
    boundary={}
    for x,y in component:
        corners=[(x,y),(x+1,y),(x+1,y+1),(x,y+1)]
        for a,b in zip(corners,corners[1:]+corners[:1]):
            edge=tuple(sorted((a,b)))
            boundary[edge]=boundary.get(edge,0)+1
    adjacency={}
    for (a,b),count in boundary.items():
        if count==1:
            adjacency.setdefault(a,set()).add(b);adjacency.setdefault(b,set()).add(a)
    area=float(len(component)*h*h)
    result=dict(base,area=area,selected_polygons=len(component),
                exact_area_square_count=len(component))
    touching=any((x+dx,y+dy) in squares and (x+dx,y+dy) not in component
                 for x,y in component for dx,dy in [(-1,-1),(-1,1),(1,-1),(1,1)])
    if touching or any(len(v)!=2 for v in adjacency.values()):
        return dict(result,status='NONMANIFOLD_OR_CONTACT_SECTION')
    seen=set();loops=0
    for start in adjacency:
        if start in seen:continue
        loops+=1;todo=[start]
        while todo:
            at=todo.pop()
            if at in seen:continue
            seen.add(at);todo.extend(adjacency[at]-seen)
    if loops!=1:return dict(result,status='NONDISK_SECTION_HOLES_NOT_FILLED')
    center=list(map(float,point));witnesses=[]
    for i,a in enumerate(other):
        center[a]=float((Fraction(sum(v[i] for v in component),len(component))+Fraction(1,2))*h)
    for square in sorted(component):
        witness=list(map(float,point))
        for i,a in enumerate(other):witness[a]=float((square[i]+Fraction(1,2))*h)
        witnesses.append(witness)
    return dict(result,status='CLOSED_DISK_AXIS_ALIGNED_CUBICAL_SECTION',radius=math.sqrt(area/math.pi),
                section_centroid_evidence=dict(centroid=center,_interior_polygon_witnesses=witnesses,
                    definition='area centroid of the unchanged axis-aligned cube-union plane section'))


class NativeSectionIndex:
    """Enumerate a complete query-connected plane slab; retain original cut math."""
    def __init__(self,solid,kernel_kind):
        self.solid=solid;self.kernel_kind=kernel_kind
        self.occupied={tuple(map(int,c)):i for i,c in enumerate(solid.cells)}
        self.steps=list(product((-1,0,1),repeat=3))
        self.maximum=float(np.abs(solid.cells).max()) if len(solid.cells) else 0.
        self.extreme=int(np.argmax(np.max(np.abs(solid.cells),axis=1))) if len(solid.cells) else None
        self.last_stats=None

    def section(self,point,normal):
        p=np.asarray(point,float)/self.solid.pitch;n=np.asarray(normal,float)
        size=math.hypot(*n)
        if not size>0:raise ValueError('NONZERO_SECTION_NORMAL_REQUIRED')
        n=n/size
        tol=256*np.finfo(float).eps*max(1.,float(abs(p).max()),self.maximum)
        width=.5*float(abs(n).sum())+tol
        origin=np.floor(p).astype(np.int64)
        seeds=[]
        for step in self.steps:
            key=tuple(int(origin[a])+step[a] for a in range(3))
            if key in self.occupied and all(key[a]-tol<=p[a]<=key[a]+1+tol for a in range(3)):
                seeds.append(key)
        if not seeds:
            return dict(radius=None,area=None,status='QUERY_OUTSIDE_SECTION',
                        index_evidence=dict(source_cells=len(self.solid.cells),selected_cells=0))
        pending=list(seeds);tested=set();selected=set()
        while pending:
            key=pending.pop()
            if key in tested:continue
            tested.add(key)
            if abs(sum((key[a]+.5-p[a])*n[a] for a in range(3)))>width:continue
            selected.add(key)
            for step in self.steps:
                nxt=tuple(key[a]+step[a] for a in range(3))
                if nxt in self.occupied and nxt not in tested:pending.append(nxt)
        # Kernel diagonal face-path checks may need occupied zero-area neighbors.
        halo=set(selected)
        for key in selected:
            for axis in range(3):
                for sign in (-1,1):
                    nxt=list(key);nxt[axis]+=sign;nxt=tuple(nxt)
                    if nxt in self.occupied:halo.add(nxt)
        ids={self.occupied[key] for key in halo}
        if self.extreme is not None:ids.add(self.extreme)
        ordered=sorted(ids)
        subset=ReconstructedSolid(self.solid.cells[ordered],self.solid.pitch,
                                  dict(query_index_subset=True,source_material_modified=False))
        result=self.kernel_kind(subset).section(point,normal)
        self.last_stats=dict(source_cells=len(self.solid.cells),selected_cells=len(ordered),
                             plane_component_cells=len(selected),original_tolerance_extremum_retained=True,
                             full_source_global_components_enumerated=False,
                             query_connected_plane_cells_complete=True,source_modified=False)
        result['index_evidence']=self.last_stats
        return result


def indexed_body_measurement(solid,points,adjacency,progress=None):
    import root_branch_section_measurement as original
    kind=original.owned_section_class()
    class Engine:
        def __init__(self,source):self.index=NativeSectionIndex(source,kind)
        def section(self,point,normal):return self.index.section(point,normal)
    fn=original.branch_section_centroid_measurement
    bound=FunctionType(fn.__code__,dict(fn.__globals__,owned_section_class=lambda:Engine),
                       fn.__name__,fn.__defaults__)
    bound.__kwdefaults__=fn.__kwdefaults__
    return bound(solid,points,adjacency,progress=progress)


def geometric_endpoints(solid,graph,*,_ownership=None,_only_endpoints=None,_retry=True):
    """Trace supported section centroids to finite model ends, never living tips."""
    from root_reconstruction_measurement import cubical_ray_exit
    from root_section_centroid import _covered
    source=np.asarray(graph['points'],float);source_adj=graph['adjacency']
    ledger=graph_measurement_ledger(source,source_adj,coordinate_units='model_units')
    from root_branch_section_measurement import BranchOwnership
    from root_cubical_disk_section import interior_boundary_section_class
    ownership=BranchOwnership(source,source_adj) if _ownership is None else _ownership
    section_index=NativeSectionIndex(solid,interior_boundary_section_class())
    slab_cache={}
    records=[];routes={};occupied=set(map(tuple,solid.cells));h=solid.pitch
    def supported(a,b):
        low=np.floor(np.minimum(a,b)/h).astype(int)-1
        high=np.floor(np.maximum(a,b)/h).astype(int)
        if math.prod(map(int,high-low+1))<=4096:
            boxes=np.asarray([p for p in product(*(range(int(x),int(y)+1) for x,y in zip(low,high)))
                              if p in occupied],np.int64).reshape(-1,3)
        else:boxes=solid.cells
        return bool(len(boxes)) and _covered(boxes,h,np.asarray(a),np.asarray(b))
    def section(point,normal):
        value=axis_aligned_section(solid,point,normal,slab_cache)
        return section_index.section(point,normal) if value is None else value
    def bounded_exit(point,direction,maximum,branch_id):
        try:ray=cubical_ray_exit(solid,point,direction)
        except ValueError as exc:
            if str(exc)!='RAY_START_NOT_IN_OCCUPIED_CELL':raise
            inward=np.nextafter(point,point-direction*h)
            ray=cubical_ray_exit(solid,inward,direction)
            if math.dist(ray['boundary_point'],point)>1e-8*h:raise
            ray=dict(ray,boundary_point=point.tolist(),boundary_distance=0.)
        end=np.asarray(ray['boundary_point'],float)
        if ray['boundary_distance']>maximum or not supported(point,end):return None
        if ownership.nearest_owners((point+end)/2)!={branch_id}:return None
        return end
    for branch in ledger['branches']:
        ids=branch['source_node_ids'];path=source[ids]
        identities=[ownership.node_branches[i] for i in ids if i<len(ownership.node_branches)]
        owners=set.intersection(*identities) if identities else set()
        owner_id=next(iter(owners)) if len(owners)==1 else branch['branch_id']
        arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(path,axis=0),axis=1))]
        for at in [0,len(ids)-1]:
            node=ids[at]
            if len(source_adj[node])!=1:continue
            if _only_endpoints is not None and node not in _only_endpoints:continue
            record=dict(node_id=node,status='endpoint_unresolved',original_position=source[node].tolist(),
                        physical_tip_qualified=False)
            interior=path[1:-1] if len(path)>=4 else path
            interior_arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(interior,axis=0),axis=1))]
            tangent=_tangent(interior,interior_arc,0. if at==0 else interior_arc[-1],solid.pitch,False)
            if tangent is None:
                record['reason']='NO_SUPPORTED_OUTWARD_TANGENT';records.append(record);continue
            direction=-tangent if at==0 else tangent
            try:
                initial=section(source[node],tangent)
                evidence=initial.get('section_centroid_evidence')
                record['initial_section_evidence']=dict(position=source[node].tolist(),
                    normal=tangent.tolist(),status=initial['status'],radius=initial.get('radius'),
                    unique_disk=initial.get('radius') is not None and evidence is not None,
                    witness_count=0 if evidence is None else len(evidence['_interior_polygon_witnesses']))
                if initial.get('radius') is None or evidence is None:
                    record['reason']='ENDPOINT_SOURCE_SECTION_NOT_A_UNIQUE_DISK';records.append(record);continue
                if any(ownership.nearest_owners(w)!={owner_id}
                       for w in evidence['_interior_polygon_witnesses']):
                    record['reason']='ENDPOINT_SECTION_HAS_MIXED_BRANCH_SUPPORT';records.append(record);continue
                origin=np.asarray(evidence['centroid'],float)
                neighbor=source_adj[node][0]
                if not supported(source[neighbor],origin):
                    record['reason']='ENDPOINT_SOURCE_CENTRE_NOT_CONNECTED';records.append(record);continue
                route=[origin];radius=initial['radius'];travel=0.;finished=False
                protected=path[arc>10*h] if at==0 else path[arc<arc[-1]-10*h]
                for step_number in range(256):
                    step=min(h,max(h/16,radius/2))
                    proposal=route[-1]+step*direction
                    cut=section(proposal,direction)
                    if cut.get('radius') is None:
                        if cut['status'] not in ('QUERY_OUTSIDE_SECTION','NO_POSITIVE_AREA_SECTION'):
                            if cut['status'] in ('NONDISK_SECTION_HOLES_NOT_FILLED','NONMANIFOLD_OR_CONTACT_SECTION'):
                                end=bounded_exit(route[-1],direction,4*h,owner_id)
                                if end is not None:
                                    if math.dist(route[-1],end)>1e-12:route.append(end)
                                    record.update(endpoint_method='short_boundary_ray_after_section_refusal',terminal_section_status=cut['status'])
                                    finished=True;break
                            record['reason']='ENDPOINT_TRACE_'+cut['status'];break
                        end=bounded_exit(route[-1],direction,2*h,owner_id)
                        if end is None:
                            record['reason']='ENDPOINT_EXIT_NOT_LOCALLY_SUPPORTED';break
                        if math.dist(route[-1],end)>1e-12:route.append(end)
                        finished=True;break
                    ev=cut.get('section_centroid_evidence')
                    if ev is None or any(ownership.nearest_owners(w)!={owner_id}
                                         for w in ev['_interior_polygon_witnesses']):
                        record['reason']='ENDPOINT_TRACE_MIXED_OR_UNDEFINED_SECTION';break
                    new=np.asarray(ev['centroid'],float)
                    delta=new-route[-1];distance=math.hypot(*delta)
                    if float(delta@direction)<step/4 or distance>max(2*h,radius):
                        record['reason']='ENDPOINT_TRACE_NONLOCAL_OR_BACKWARD';break
                    if not supported(route[-1],new):
                        end=bounded_exit(route[-1],direction,4*h,owner_id)
                        if end is not None:
                            if math.dist(route[-1],end)>1e-12:route.append(end)
                            record.update(endpoint_method='short_boundary_ray_after_disconnected_centroid',terminal_section_status=cut['status'])
                            finished=True;break
                        record['reason']='ENDPOINT_TRACE_OUTSIDE_ENTITY';break
                    if len(protected) and float(np.linalg.norm(protected-new,axis=1).min())<1.5*h:
                        record['reason']='ENDPOINT_TRACE_REJOINS_EXISTING_AXIS';break
                    route.append(new);travel+=distance;radius=cut['radius']
                    if travel>max(50*h,20*initial['radius']):
                        record['reason']='ENDPOINT_TRACE_LENGTH_GUARD';break
                    recent=np.asarray(route[-5:]);new_direction=recent[-1]-recent[0]
                    norm=math.hypot(*new_direction)
                    if norm>0:
                        new_direction/=norm
                        if float(new_direction@direction)<.5:
                            record['reason']='ENDPOINT_TRACE_DIRECTION_REVERSAL';break
                        direction=new_direction
                if finished:
                    routes[node]=route
                    record.update(status='geometric_surface_estimate',position=route[-1].tolist(),
                                  extension=math.fsum(math.dist(a,b) for a,b in zip(route,route[1:])),
                                  trace_points=len(route),direction=direction.tolist(),
                                  centered_source_section=origin.tolist(),
                                  region='finite delivered cube union',
                                  ray_status='SOURCE_SUPPORTED_SECTION_TRACE_AND_BOUNDARY_EXIT',
                                  original_entity_unchanged=True)
                elif 'reason' not in record:record['reason']='ENDPOINT_TRACE_STEP_GUARD'
                if not finished:
                    record.update(last_supported_position=route[-1].tolist(),last_direction=direction.tolist(),
                                  supported_trace_length=travel,trace_points=len(route))
                    try:record['boundary_ray_diagnostic']=cubical_ray_exit(solid,route[-1],direction)
                    except ValueError as exc:record['boundary_ray_diagnostic']=dict(reason=str(exc))
            except (ValueError,RuntimeError) as exc:record['reason']=str(exc)
            records.append(record)
    for a,b in ledger['source_edge_ids']:
        if a in routes and b in routes and not supported(routes[a][0],routes[b][0]):
            for node in [a,b]:
                routes.pop(node,None)
                for r in records:
                    if r['node_id']==node:r.update(status='endpoint_unresolved',reason='ENDPOINT_CENTRES_NOT_SOURCE_CONNECTED')
    points=source.tolist();adj=[[] for _ in points]
    for node,route in routes.items():points[node]=route[-1].tolist()
    for a,b in ledger['source_edge_ids']:
        interior_points=[]
        if a in routes:interior_points.extend(reversed(routes[a][:-1]))
        if b in routes:interior_points.extend(routes[b][:-1])
        chain=[a]
        for value in interior_points:
            if math.dist(points[chain[-1]],value)<=1e-12:continue
            node=len(points);points.append(list(map(float,value)));adj.append([]);chain.append(node)
        chain.append(b)
        for u,v in zip(chain,chain[1:]):adj[u].append(v);adj[v].append(u)
    result=dict(points=points,adjacency=adj)
    if _retry and any(r['status']!='geometric_surface_estimate' for r in records):
        from root_native_endpoint_retry import retry_from_interior
        result,records=retry_from_interior(solid,result,records,ownership)
    return result,records


def subdivide_measurement_graph(graph,max_length):
    points=[list(p) for p in graph['points']];adj=[[] for _ in points]
    before=graph_measurement_ledger(points,graph['adjacency'],coordinate_units='model_units')
    for a,b in before['source_edge_ids']:
        first=np.asarray(points[a]);last=np.asarray(points[b]);length=math.dist(first,last)
        count=max(1,math.ceil(length/max_length-1e-12));previous=a
        for i in range(1,count):
            node=len(points);points.append((first+(last-first)*(i/count)).tolist());adj.append([])
            adj[previous].append(node);adj[node].append(previous);previous=node
        adj[previous].append(b);adj[b].append(previous)
    after=graph_measurement_ledger(points,adj,coordinate_units='model_units')
    assert before['source_components']==after['source_components']
    assert before['source_cycle_rank']==after['source_cycle_rank']
    assert math.isclose(before['measurement_graph_length'],after['measurement_graph_length'],rel_tol=1e-12,abs_tol=1e-12)
    return dict(points=points,adjacency=adj,subdivided_without_geometry_change=True)


def edge_radius_profile(solid, points, adjacency, progress=None):
    from root_branch_section_measurement import BranchOwnership
    from root_cubical_disk_section import interior_boundary_section_class
    points = np.asarray(points,float)
    ledger = graph_measurement_ledger(points,adjacency,coordinate_units='model_units')
    owner = BranchOwnership(points,adjacency)
    engine = NativeSectionIndex(solid,interior_boundary_section_class())
    stations = []
    for branch in ledger['branches']:
        ids = branch['source_node_ids']
        path = points[ids]
        lengths = np.linalg.norm(np.diff(path,axis=0),axis=1)
        arc = np.r_[0,np.cumsum(lengths)]
        for at,(a,b) in enumerate(zip(ids,ids[1:])):
            normal = _tangent(path,arc,(arc[at]+arc[at+1])/2,solid.pitch,branch['closed'])
            stations.append(dict(node_ids=[int(a),int(b)],branch_id=branch['branch_id'],
                                 position=((points[a]+points[b])/2).tolist(),
                                 normal=None if normal is None else normal.tolist()))
    rows = [];slab_cache={}
    for at,station in enumerate(stations):
        if station['normal'] is None:
            result = dict(radius=None,status='MISSING_UNIQUE_BRANCH_NORMAL')
        else:
            result = axis_aligned_section(solid,station['position'],station['normal'],slab_cache)
            if result is None:result = engine.section(station['position'],station['normal'])
            if result.get('radius') is not None:
                witnesses = result['section_centroid_evidence']['_interior_polygon_witnesses']
                if any(owner.nearest_owners(p) != {station['branch_id']} for p in witnesses):
                    result = dict(result,radius=None,status='MISSING_MIXED_BRANCH_SECTION')
        center = result.get('section_centroid_evidence',{}).get('centroid')
        rows.append(dict(station,radius=result.get('radius'),available=result.get('radius') is not None,
                         status=result['status'],area=result.get('area'),section_centroid=center,
                         centroid_offset=None if center is None else math.dist(center,station['position'])))
        if progress is not None and (at % 50 == 0 or at+1 == len(stations)):
            progress('radius', processed=at+1,total=len(stations),
                     available=sum(r['available'] for r in rows))
    return dict(sampling='edge_midpoint',rows=rows,
                definition='sqrt(area/pi) of one branch-associated disk of the unchanged cube union',
                tangent_definition='21 sample local line derivative over at most 10 native pitches',
                stations_fixed_before_radius_outcomes=True)


def measure_native(entity, output, *, stage='regularized_negative_indices',
                   coordinate_units=None, scale_mm=None, inherit_from=None,
                   resume_axis_from=None,backend='auto', max_dense_cells=350000000, progress=None,
                   _axis_order=None,_stop_after_endpoints=False,_allow_axis_candidate=True):
    start = time.monotonic()
    def emit(phase, **data):
        if progress is not None:
            progress(phase,elapsed_seconds=time.monotonic()-start,**data)
    emit('load_entity')
    solid,binding,metadata = load_entity(entity,stage=stage)
    emit('entity_loaded', occupied_cells=len(solid.cells),pitch=solid.pitch,
         occupied_shape=(solid.cells.max(axis=0)-solid.cells.min(axis=0)+1).tolist())
    if scale_mm is not None and (not math.isfinite(scale_mm) or scale_mm <= 0):
        raise ValueError('POSITIVE_FINITE_MM_SCALE_REQUIRED')
    requested_backend=backend
    if backend=='auto':
        # Domain size is not a geometric reason to use a different thinning
        # estimator. The small-span path creates spurious leaves on curved
        # solids; use the same topology-guarded compiled kernel at every size.
        backend='existing_guarded_lee'
    if _axis_order is None:
        _axis_order=(tuple(json.loads((PROJECT/'configs/native_model_delivery_v001.json').read_text(encoding='utf-8'))['compiled_lee_axis_order'])
                     if backend=='existing_guarded_lee' else (0,1,2))
    graph = radius = None
    failures = []
    provenance = dict(mode='fresh',backend=backend,requested_backend=requested_backend,
                      execution_network_policy=network_policy_state(),
                      graph_policy='native_cube_h1_quotient_junction_ports_v002',
                      endpoint_policy='owned_trace_interior_and_short_boundary_challenge_v003',
                      implementation_sources=[file_reference(Path(__file__).with_name(name),WORKSPACE) for name in
                          ['root_native_entity_measurement.py','root_cubical_graph_quotient.py',
                           'root_native_junction_fit.py','root_native_endpoint_retry.py','root_native_axis_order.py',
                           'root_local_network_policy.py','root_cubical_sections.py','root_section_centroid.py']],
                      fresh_automatic_measurement=True,
                      input_graph_or_GT_used=False,ply_revoxelization_performed=False)
    provenance['algorithm_axis_order']=list(_axis_order)
    stages = {};endpoint_records=[]
    try:
        if inherit_from is not None:
            graph,radius,provenance = inherit_bound_measurement(inherit_from,binding)
            emit('inherited_results_bound',nodes=len(graph['points']))
        elif resume_axis_from is not None:
            previous_path=Path(resume_axis_from).resolve()
            previous=json.loads(previous_path.read_text(encoding='utf-8'))
            if previous.get('axis_stage_reference'):
                ref=previous['axis_stage_reference'];stage_path=Path(ref['path']).resolve()
                if not stage_path.is_relative_to(PROJECT/'experiments'):
                    raise ValueError('AXIS_STAGE_REFERENCE_OUTSIDE_PROJECT_EXPERIMENTS')
                current=file_reference(stage_path)
                if current['sha256']!=ref['sha256'] or current['bytes']!=ref['bytes']:
                    raise ValueError('AXIS_STAGE_CACHE_CHANGED')
                previous=json.loads(stage_path.read_text(encoding='utf-8'))
            old=previous['entity_binding']
            if previous.get('schema_version')!=SCHEMA or previous['provenance']['mode']=='inherited':
                raise ValueError('RESUME_REQUIRES_PRIOR_AUTOMATIC_NATIVE_AXIS')
            if any(old[k]!=binding[k] for k in ['stage','pitch','cells_sha256','occupied_cells']):
                raise ValueError('RESUME_ENTITY_STAGE_OR_COORDINATES_MISMATCH')
            if old['metadata']['sha256']!=binding['metadata']['sha256']:
                raise ValueError('RESUME_ENTITY_METADATA_MISMATCH')
            if previous['coordinate_units']!=(coordinate_units or solid.provenance['units']):
                raise ValueError('RESUME_COORDINATE_UNIT_MISMATCH')
            if previous['provenance'].get('backend')!=backend:
                raise ValueError('RESUME_BACKEND_MISMATCH')
            if previous['provenance'].get('graph_policy')!=provenance['graph_policy'] and previous.get('stages',{}).get('skeleton',{}).get('representation',{}).get('status')!='PASS_SUPPORTED_GEOMETRIC_GRAPH_TOPOLOGY':
                raise ValueError('RESUME_GRAPH_POLICY_MISMATCH_REQUIRES_FRESH_AXIS')
            graph=previous['candidate_centerline']
            if not graph or not graph['points']:raise ValueError('RESUME_HAS_NO_AXIS')
            if graph.get('subdivided_without_geometry_change') and not previous.get('axis_stage_only'):
                raise ValueError('RESUME_REQUIRES_UNMODIFIED_BODY_AXIS_CHECKPOINT')
            # A previous geometry-endpoint update must not be silently applied twice.
            for endpoint in previous.get('geometric_endpoints',[]):
                graph['points'][endpoint['node_id']]=endpoint['original_position']
            stages=previous.get('stages',{})
            provenance.update(mode='reused_axis_new_sections',fresh_automatic_measurement=False,
                              algorithm_axis_order=previous['provenance'].get('algorithm_axis_order',[0,1,2]),
                              fresh_section_measurement=True,reused_axis=file_reference(previous_path))
            emit('axis_stage_reused',nodes=len(graph['points']))
        else:
            emit('skeleton',backend=backend)
            if backend == 'existing_span':
                from root_span_model_pipeline import measure_span_tube
                kernel = measure_span_tube(solid,max_dense_cells=max_dense_cells)
            elif backend == 'existing_guarded_lee':
                from root_native_axis_order import original_frame_kernel
                kernel=original_frame_kernel(solid,_axis_order,max_dense_cells)
            else:
                raise ValueError('UNKNOWN_NATIVE_BACKEND')
            if kernel['representation']['status']!='PASS_SUPPORTED_GEOMETRIC_GRAPH_TOPOLOGY':
                from root_cubical_graph_quotient import quotient_graph
                cells=np.asarray(kernel['source_skeleton_cells'],np.int64)
                positions,adjacency,receipt=quotient_graph(cells,solid.pitch)
                receipt['legacy_representation']=kernel['representation']
                receipt['status']='PASS_SUPPORTED_GEOMETRIC_GRAPH_TOPOLOGY' if receipt['skeleton_voxel_H2']==0 else 'UNRESOLVED_CAVITIES_NOT_A_PURE_GRAPH'
                ids=receipt['source_voxel_ids']
                kernel.update(points=positions.tolist(),adjacency=adjacency,
                              source_skeleton_cells=cells[ids].tolist(),representation=receipt)
                kernel['ledger']=graph_measurement_ledger(positions,adjacency,coordinate_units='caller_solid_coordinate_units')
                kernel['cubical_boundary_radii']=[kernel['cubical_boundary_radii'][i] for i in ids]
                emit('graph_quotient',cycle_rank=receipt['graph_cycle_rank'],cavities=receipt['skeleton_voxel_H2'],
                     filled_relations_removed=receipt['pivot_chords_removed'])
            stages['skeleton'] = {k:kernel[k] for k in ['status','thinning','representation']}
            graph = dict(points=kernel['points'],adjacency=kernel['adjacency'])
            emit('skeleton_complete',nodes=len(graph['points']),
                 representation=kernel['representation']['status'])
            if kernel['representation']['status'] == 'PASS_SUPPORTED_GEOMETRIC_GRAPH_TOPOLOGY':
                from root_body_measurement_api import measure_centerline_body
                emit('body_correction')
                body = indexed_body_measurement(
                    solid,kernel['points'],kernel['adjacency'],
                    progress=lambda value:emit('body_correction',**value))
                checked = measure_centerline_body(solid,coordinate_units=coordinate_units or
                                                 solid.provenance['units'],kernel_graph=kernel,
                                                 frozen_body=body,include_radius_profile=False)
                assert checked['body_geometry'] is body
                graph = dict(points=body['points'],adjacency=body['adjacency'])
                stages['body_correction'] = dict(status=body['status'],moved_nodes=body['moved_nodes'])
            else:
                failures.append(dict(code='UNSUPPORTED_SOURCE_GRAPH_REPRESENTATION',position=None))
        if inherit_from is None and graph and graph['points']:
            from root_native_junction_fit import correct_junctions,METHOD
            if stages.get('skeleton',{}).get('representation',{}).get('status')=='PASS_SUPPORTED_GEOMETRIC_GRAPH_TOPOLOGY' and stages.get('junction_fit',{}).get('method')!=METHOD:
                emit('junction_correction')
                graph,junction_receipt=correct_junctions(solid,graph,progress=emit)
                stages['junction_fit']=junction_receipt
            axis_path=Path(output).with_name(Path(output).stem+'_axis.json')
            write_json(axis_path,dict(schema_version=SCHEMA,axis_stage_only=True,
                entity_binding=binding,coordinate_units=coordinate_units or solid.provenance['units'],
                candidate_centerline=graph,provenance=provenance,stages=stages))
            provenance['axis_stage_reference']=file_reference(axis_path)
            emit('geometric_endpoints')
            graph,endpoint_records=geometric_endpoints(solid,graph)
            if _allow_axis_candidate and not _stop_after_endpoints and backend=='existing_guarded_lee' and resume_axis_from is None:
                stages['axis_selection']=dict(policy='retain qualified primary012; try only fixed120 on unresolved endpoints',
                    selected_order=list(_axis_order),candidate_attempted=False,GT_used=False)
            if (_allow_axis_candidate and not _stop_after_endpoints and inherit_from is None
                and resume_axis_from is None and backend=='existing_guarded_lee'
                and any(r['status']!='geometric_surface_estimate' for r in endpoint_records)):
                from root_native_axis_order import endpoint_candidate_admission
                order=tuple(json.loads((PROJECT/'configs/native_model_delivery_v001.json').read_text(encoding='utf-8'))['compiled_lee_endpoint_candidate_axis_order'])
                emit('axis_candidate',order=list(order))
                trial_path=Path(output).with_name(Path(output).stem+'_axis120_attempt.json')
                primary=dict(candidate_centerline=graph,geometric_endpoints=endpoint_records,entity_binding=binding,stages=stages)
                primary_path=Path(output).with_name(Path(output).stem+'_primary_endpoint_attempt.json')
                write_json(primary_path,dict(primary,schema_version=SCHEMA,status='PRIMARY_ENDPOINT_ATTEMPT_ONLY',
                    provenance=dict(provenance),radius_curve=None,coordinate_units=coordinate_units or solid.provenance['units']))
                trial,_=measure_native(entity,trial_path,stage=stage,coordinate_units=coordinate_units,
                    scale_mm=scale_mm,backend='existing_guarded_lee',max_dense_cells=max_dense_cells,
                    progress=lambda phase,**data:emit('axis_candidate_progress',order=list(order),candidate_phase=phase,**{k:v for k,v in data.items() if k!='elapsed_seconds'}),
                    _axis_order=order,_stop_after_endpoints=True,_allow_axis_candidate=False)
                admission=endpoint_candidate_admission(primary,trial)
                selection=dict(policy='retain qualified primary012; try only fixed120 on unresolved endpoints',
                    primary_axis_reference=provenance['axis_stage_reference'],primary_attempt=file_reference(primary_path),candidate=file_reference(trial_path),
                    admission=admission,selected_order=list(_axis_order),candidate_attempted=True,GT_used=False)
                if admission['qualified']:
                    graph=trial['candidate_centerline'];endpoint_records=trial['geometric_endpoints']
                    stages=trial['stages'];selection['selected_order']=list(order)
                    provenance.update(algorithm_axis_order=list(order),primary_axis_reference=provenance['axis_stage_reference'],
                                      axis_stage_reference=trial['axis_stage_reference'])
                    emit('axis_candidate_selected',order=list(order))
                stages['axis_selection']=selection
            if not _stop_after_endpoints:
                graph=subdivide_measurement_graph(graph,solid.pitch)
                emit('radius',processed=0,total=sum(map(len,graph['adjacency']))//2)
                radius = edge_radius_profile(solid,graph['points'],graph['adjacency'],progress=emit)
    except (ValueError,RuntimeError,AssertionError,MemoryError,OverflowError) as exc:
        failures.append(dict(code='NATIVE_PIPELINE_BLOCKED',position=None,
                             error_type=type(exc).__name__,reason=str(exc)))
        emit('pipeline_blocked',error_type=type(exc).__name__,reason=str(exc))
    if radius:
        failures.extend(dict(code=r['status'],position=r['position'],node_id=r.get('node_id'),
                             edge_nodes=r.get('node_ids')) for r in radius['rows'] if not r['available'])
    units = coordinate_units or solid.provenance['units']
    ledger = graph_measurement_ledger(graph['points'],graph['adjacency'],coordinate_units=units) \
        if graph and graph['points'] else None
    if ledger:
        resolved={r['node_id'] for r in endpoint_records if r['status']=='geometric_surface_estimate'}
        failures.extend(dict(code='GEOMETRIC_ENDPOINT_UNRESOLVED',node_id=n,
                             position=graph['points'][n]) for n in ledger['endpoint_node_ids'] if n not in resolved)
        failures.extend(dict(code='ISOLATED_COMPONENT_NO_UNIQUE_AXIS',node_id=n,
                             position=graph['points'][n]) for n in ledger['isolated_node_ids'])
    check_entity_unchanged(binding)
    assert digest_cells(solid) == binding['cells_sha256']
    result = dict(schema_version=SCHEMA,
                  status='NATIVE_AXIS_ENDPOINT_ATTEMPT_ONLY' if _stop_after_endpoints else
                  'INHERITED_OUTPUT_WITH_EXPLICIT_GAPS' if inherit_from else
                  'NATIVE_OUTPUT_WITH_EXPLICIT_GAPS' if ledger and radius and not any(
                      f['code'] in ['UNSUPPORTED_SOURCE_GRAPH_REPRESENTATION','NATIVE_PIPELINE_BLOCKED'] for f in failures) else 'NATIVE_MEASUREMENT_BLOCKED',
                  case_label=metadata.get('root',Path(entity).stem),
                  input_file_binding=dict(sha256=binding['metadata']['sha256'],
                                          bytes=binding['metadata']['bytes'],format='native_entity_package'),
                  entity_binding=binding,coordinate_units=units,scale_mm_per_model_unit=scale_mm,
                  candidate_centerline=graph,radius_curve=radius,
                  length=None if ledger is None else dict(value=ledger['measurement_graph_length'],
                         definition='sum of unique final graph edges',
                         endpoint_status='not_applicable_closed_graph' if not ledger['endpoint_node_ids'] else
                         'geometric_surface_estimates' if endpoint_records and all(
                             r['status']=='geometric_surface_estimate' for r in endpoint_records) else 'endpoint_unresolved',
                         complete_endpoint_length_qualified=False),
                  model_bounds=dict(minimum=(solid.cells.min(axis=0)*solid.pitch).tolist(),
                                    maximum=((solid.cells.max(axis=0)+1)*solid.pitch).tolist()),
                  provenance=provenance,stages=stages,failures=failures,geometric_endpoints=endpoint_records,
                  axis_stage_reference=provenance.get('axis_stage_reference'),
                  functional_acceptance=False,precision_acceptance=False,
                  worker_seconds=time.monotonic()-start,original_entity_unchanged=True)
    write_json(output,result)
    emit('measurement_written',path=Path(output).name)
    return result,solid


def run_native_worker(args):
    disable_measurement_network()
    from root_model_delivery_exports import write_native_delivery
    prefix = Path(args.output_prefix).resolve()
    experiment = next((p for p in prefix.parents if (p/'EXPERIMENT_LOG.md').is_file()),None)
    if experiment is None or not prefix.is_relative_to(experiment/'results'):
        raise ValueError('REGISTERED_EXPERIMENT_RESULTS_REQUIRED')
    measurement = prefix.with_name(prefix.name+'_measurement.json')
    log = prefix.with_name(prefix.name+'_progress.jsonl')
    if any(p.exists() for p in [measurement,log,prefix.with_name(prefix.name+'_delivery.json')]):
        raise ValueError('OUTPUT_PREFIX_ALREADY_EXISTS')
    prefix.parent.mkdir(parents=True,exist_ok=True)
    def progress(phase,**data):
        event = dict(phase=phase,**data)
        with log.open('a',encoding='utf-8') as stream:
            stream.write(json.dumps(event,ensure_ascii=False,allow_nan=False)+'\n')
        print(json.dumps(event,ensure_ascii=False,allow_nan=False),flush=True)
    result,solid = measure_native(args.entity,measurement,stage=args.entity_stage,
                                  coordinate_units=args.coordinate_units,
                                  scale_mm=args.scale_mm,inherit_from=args.inherit_from,
                                  resume_axis_from=args.resume_axis_from,
                                  backend=args.native_backend,progress=progress)
    progress('exports_and_original_boundary')
    data,out,html = write_native_delivery(measurement,prefix,experiment,solid,
        morphology_options=dict(vertical_axis=getattr(args,'vertical_axis','z'),down_sign=getattr(args,'down_sign',-1)))
    check_entity_unchanged(result['entity_binding'])
    progress('complete',delivery=str(out),viewer=str(html),summary=data['summary'])
    return data
