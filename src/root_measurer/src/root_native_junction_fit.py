"""Fit finite junction positions from source-section centroids on incident ports."""
import itertools
import math
import numpy as np
from root_branch_section_measurement import BranchOwnership
from root_cubical_disk_section import interior_boundary_section_class
from root_graph_measurement_ledger import graph_measurement_ledger
from root_section_centroid import _covered

METHOD='source_owned_port_line_intersection_v001'


def correct_junctions(solid,graph,progress=None):
    from root_native_entity_measurement import NativeSectionIndex,axis_aligned_section,_tangent
    points=np.asarray(graph['points'],float);adj=graph['adjacency'];h=solid.pitch
    ledger=graph_measurement_ledger(points,adj,coordinate_units='model_units')
    ownership=BranchOwnership(points,adj);engine=NativeSectionIndex(solid,interior_boundary_section_class());slabs={}
    occupied=set(map(tuple,solid.cells));ports={i:[] for i,a in enumerate(adj) if len(a)>=3}
    paths={}
    for branch in ledger['branches']:
        ids=branch['source_node_ids'];path=points[ids];arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(path,axis=0),axis=1))]
        paths[branch['branch_id']]=(ids,path,arc)
        if branch['closed']:continue
        for end,reverse in [(ids[0],False),(ids[-1],True)]:
            if end not in ports:continue
            q=path[::-1] if reverse else path;s=arc[-1]-arc[::-1] if reverse else arc
            maximum=min(20*h,.4*float(s[-1]));centers=[];radii=[];locations=[]
            if maximum>3*h:
                for at in np.linspace(3*h,maximum,12):
                    pos=np.array([np.interp(at,s,q[:,a]) for a in range(3)])
                    normal=_tangent(q,s,float(at),h,False)
                    if normal is None:continue
                    section=axis_aligned_section(solid,pos,normal,slabs)
                    if section is None:section=engine.section(pos,normal)
                    ev=section.get('section_centroid_evidence')
                    if section.get('radius') is None or ev is None:continue
                    if any(ownership.nearest_owners(w)!={branch['branch_id']} for w in ev['_interior_polygon_witnesses']):continue
                    centers.append(ev['centroid']);radii.append(section['radius']);locations.append(float(at))
            record=dict(branch_id=branch['branch_id'],at_end=reverse,source_length=float(s[-1]),
                        fixed_query_count=12,valid_centroids=len(centers),qualified=False)
            if len(centers)>=3 and max(locations)-min(locations)>=2*h:
                centers=np.asarray(centers);center=centers.mean(axis=0);_,_,vectors=np.linalg.svd(centers-center,full_matrices=False)
                direction=vectors[0];normal_projector=np.eye(3)-np.outer(direction,direction)
                residual=float(np.sqrt(np.mean(np.sum(((centers-center)@normal_projector)**2,axis=1))))
                record.update(qualified=residual<=.5*h,center=center.tolist(),direction=direction.tolist(),
                              residual=residual,anchor_arc=max(locations),radius=float(np.median(radii)))
            ports[end].append(record)
    candidates={};receipts=[]
    for node,rows in ports.items():
        receipt=dict(node_id=node,ports=rows,status='UNCHANGED_UNQUALIFIED_PORTS')
        if len(rows)==len(adj[node]) and all(r['qualified'] for r in rows):
            projectors=[np.eye(3)-np.outer(r['direction'],r['direction']) for r in rows]
            matrix=sum(projectors);values=np.linalg.eigvalsh(matrix)
            if values[0]>1e-6:
                candidate=np.linalg.solve(matrix,sum(a@np.array(r['center']) for a,r in zip(projectors,rows)))
                residual=max(float(np.linalg.norm(a@(candidate-np.array(r['center'])))) for a,r in zip(projectors,rows))
                displacement=math.dist(candidate,points[node]);bound=min(4*h,.25*min(r['source_length'] for r in rows))
                if displacement<=bound and residual<=h and tuple(np.floor(candidate/h).astype(int)) in occupied:
                    candidates[node]=candidate;receipt.update(status='CANDIDATE',position=candidate.tolist(),
                                                              displacement=displacement,residual=residual,bound=bound)
                else:receipt.update(status='UNCHANGED_POSITION_OR_RESIDUAL_GUARD',displacement=displacement,residual=residual)
        receipts.append(receipt)
    def construct(active):
        result=points.copy();influences={i:set() for i in range(len(points))}
        for node in active:result[node]=candidates[node];influences[node].add(node)
        for branch in ledger['branches']:
            if branch['closed']:continue
            ids,path,arc=paths[branch['branch_id']];total=arc[-1]
            ends=[i for i in [ids[0],ids[-1]] if i in active]
            anchors={n:next(r['anchor_arc'] for r in ports[n] if r['branch_id']==branch['branch_id']) for n in ends}
            if len(ends)==2 and sum(anchors.values())>=total:
                for i,node in enumerate(ids):
                    result[node]=(1-arc[i]/total)*candidates[ids[0]]+(arc[i]/total)*candidates[ids[-1]]
                    influences[node].update(ends)
                continue
            for end in ends:
                reverse=end==ids[-1];route=ids[::-1] if reverse else ids;s=total-arc[::-1] if reverse else arc
                at=anchors[end];anchor=np.array([np.interp(at,s,points[route,a]) for a in range(3)])
                for position,node in zip(s,route):
                    if position>at:break
                    result[node]=(1-position/at)*candidates[end]+position/at*anchor;influences[node].add(end)
        return result,influences
    active=set(candidates);rolled=set()
    while active:
        result,influences=construct(active);bad=set()
        for a,b in ledger['source_edge_ids']:
            if not influences[a] and not influences[b]:continue
            lo=np.floor(np.minimum(result[a],result[b])/h).astype(int)-1;hi=np.floor(np.maximum(result[a],result[b])/h).astype(int)
            if math.prod(map(int,hi-lo+1))<=4096:
                boxes=np.array([idx for idx in itertools.product(*(range(int(x),int(y)+1) for x,y in zip(lo,hi))) if idx in occupied],np.int64).reshape(-1,3)
            else:boxes=solid.cells
            if math.dist(result[a],result[b])<=1e-10*h or not len(boxes) or not _covered(boxes,h,result[a],result[b]):
                bad.update(influences[a]|influences[b])
        if not bad:break
        active-=bad;rolled.update(bad)
    result,_=construct(active)
    for receipt in receipts:
        if receipt['node_id'] in active:receipt['status']='MOVED_WITH_SOURCE_SUPPORTED_PORTS'
        elif receipt['node_id'] in rolled:receipt['status']='ROLLED_BACK_SOURCE_SEGMENT_GUARD'
    if progress:progress('junction_correction_complete',junctions=len(ports),moved=len(active),rolled_back=len(rolled))
    return dict(points=result.tolist(),adjacency=adj),dict(method=METHOD,junctions=receipts,moved_junctions=len(active),
        moved_nodes=int(np.count_nonzero(np.any(result!=points,axis=1))),source_modified=False,
        all_node_and_edge_ids_preserved=True,short_links_never_merged=True,GT_used=False)
