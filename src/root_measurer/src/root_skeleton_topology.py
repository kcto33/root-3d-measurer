"""Recover serialized native skeleton adjacency; preserve physical edges once."""
from collections import deque
import numpy as np
from scipy.spatial import cKDTree

def branches_from_adjacency(points,adjacency,radii):
    components={};label=0
    for start in range(len(points)):
        if start in components:continue
        queue=[start];components[start]=label
        while queue:
            u=queue.pop()
            for v in sorted(adjacency[u]):
                if v not in components:components[v]=label;queue.append(v)
        label+=1
    degree=np.array([len(x) for x in adjacency]);used=set();paths=[]
    def walk(start,first):
        path=[start,first];used.add(tuple(sorted((start,first))));previous,current=start,first
        # Stop at the first junction itself, before traversing its next edge.
        while degree[current]==2 and current!=start:
            choices=sorted(v for v in adjacency[current] if v!=previous)
            if not choices:break
            nxt=choices[0];edge=tuple(sorted((current,nxt)))
            if edge in used:break
            used.add(edge);path.append(nxt);previous,current=current,nxt
        return path
    for start in range(len(points)):
        if degree[start]==2:continue
        for first in sorted(adjacency[start]):
            if tuple(sorted((start,first))) not in used:paths.append(walk(start,first))
    for start in range(len(points)):
        for first in sorted(adjacency[start]):
            if tuple(sorted((start,first))) not in used:paths.append(walk(start,first))
    expected={tuple(sorted((u,v))) for u,ns in enumerate(adjacency) for v in ns}
    assert used==expected
    assert sum(len(path)-1 for path in paths)==len(expected)
    branches=[]
    for i,path in enumerate(paths):
        p=points[path]
        branches.append({'branch_id':f'native_branch{i:04d}','nodes':p.tolist(),
                         'radii':[radii[j] for j in path],'arc_length':float(np.linalg.norm(np.diff(p,axis=0),axis=1).sum()),
                         'component_id':components[path[0]]})
    return branches,{'components':label,'physical_edges':len(expected),'all_edges_once':True,
                     'degree1':int(sum(degree==1)),'degree3plus':int(sum(degree>=3))}

def recover_native(decoded):
    values={};allpoints=[]
    for b in decoded['branches']:
        for p,r in zip(b['nodes'],b['radii']):
            key=tuple(np.round(p,12));allpoints.append(key)
            if r is not None and np.isfinite(r):
                if key in values:assert abs(values[key]-r)<1e-9
                values[key]=float(r)
    allpoints.extend(tuple(np.round(p,12)) for p in decoded['branchpoints']+decoded['endpoints'])
    points=np.array(sorted(set(allpoints)),float);tree=cKDTree(points);distance,index=tree.query(points,k=2)
    pitch=float(distance[:,1].min());vectors=(points[index[:,1]]-points)[np.isclose(distance[:,1],pitch,rtol=0,atol=1e-9)]
    x=vectors[0]/np.linalg.norm(vectors[0])
    y=next(v/np.linalg.norm(v) for v in vectors if abs(np.dot(v/np.linalg.norm(v),x))<1e-7)
    basis=np.array([x,y,np.cross(x,y)]).T;coordinates=(points-points[0])@basis/pitch;grid=np.rint(coordinates).astype(int)
    residual=float(abs(coordinates-grid).max());assert residual<1e-6,'NOT_NATIVE_LATTICE'
    lookup={tuple(c):i for i,c in enumerate(grid)};assert len(lookup)==len(points)
    adjacency=[set() for p in points]
    for i,c in enumerate(grid):
        for offset in np.ndindex(3,3,3):
            delta=np.array(offset)-1
            if not delta.any():continue
            j=lookup.get(tuple(c+delta))
            if j is not None:adjacency[i].add(j)
    radii=[values.get(tuple(p)) for p in points]
    branches,audit=branches_from_adjacency(points,adjacency,radii)
    diag=decoded['diagnostics']
    assert len(points)==diag['skeleton_voxels'],'SKELETON_NODE_IDENTITY_FAILED'
    assert audit['degree1']==diag['endpoint_count'] and audit['degree3plus']==diag['branchpoint_count'],'NATIVE_DEGREE_IDENTITY_FAILED'
    assert audit['components']==decoded['component_count'],'NATIVE_COMPONENT_IDENTITY_FAILED'
    audit.update(nodes=len(points),pitch=pitch,lattice_residual=residual,missing_radius_nodes=sum(v is None for v in radii),
                 identification='verified voxel lattice and original26-neighbor adjacency, not arbitrary endpoint distance',
                 existing_radius_unchanged=True,no_added_coordinates=True)
    return {'status':'OK','branches':branches,'branchpoints':points[[i for i,x in enumerate(adjacency) if len(x)>=3]].tolist(),
            'endpoints':points[[i for i,x in enumerate(adjacency) if len(x)==1]].tolist(),'component_count':audit['components'],
            'total_length':sum(b['arc_length'] for b in branches),'native_adjacency_audit':audit}
