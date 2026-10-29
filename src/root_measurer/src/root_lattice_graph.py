"""Remove digital cell shortcuts without cutting face edges or genuine loops.

A diagonal edge is redundant only if strictly shorter neighbor edges provide
a path inside that same unit lattice cell. This is local graph interpretation,
not a claim that all remaining voxel edges represent biological root length.
"""
import numpy as np

def reduce_lattice_graph(points,adjacency,pitch):
    points=np.asarray(points,float)
    if points.ndim!=2 or points.shape[1]!=3 or pitch<=0:raise ValueError('LATTICE_INPUT')
    grid_float=(points-points.min(axis=0))/pitch;grid=np.rint(grid_float).astype(np.int64)
    if abs(grid_float-grid).max()>1e-6:raise ValueError('NON_LATTICE_COORDINATES')
    if len(set(map(tuple,grid)))!=len(points):raise ValueError('DUPLICATE_LATTICE_POINTS')
    out=[set(ns) for ns in adjacency];edges={tuple(sorted((i,j))) for i,ns in enumerate(out) for j in ns}
    rank={}
    for i,j in edges:
        if i==j or i not in out[j]:raise ValueError('ASYMMETRIC_OR_SELF_EDGE')
        delta=abs(grid[i]-grid[j]);r=int(delta@delta)
        if delta.max()>1 or r<1 or r>3:raise ValueError('NON_NEIGHBOR_EDGE')
        rank[(i,j)]=r
    removed=[]
    ordered=sorted(edges,key=lambda e:(rank[e],tuple(grid[e[0]]),tuple(grid[e[1]])))
    for i,j in ordered:
        r=rank[(i,j)]
        if r==1:continue
        lo=np.minimum(grid[i],grid[j]);hi=np.maximum(grid[i],grid[j]);seen={i};todo=[i]
        while todo and j not in seen:
            u=todo.pop()
            for v in sorted(out[u]):
                e=tuple(sorted((u,v)))
                if v in seen or rank[e]>=r or not np.all((lo<=grid[v])&(grid[v]<=hi)):continue
                seen.add(v);todo.append(v)
        if j in seen:
            out[i].remove(j);out[j].remove(i);removed.append((i,j))
    def components(a):
        seen=set();count=0
        for start in range(len(a)):
            if start in seen:continue
            count+=1;todo=[start];seen.add(start)
            while todo:
                u=todo.pop()
                for v in a[u]-seen:seen.add(v);todo.append(v)
        return count
    before=components(adjacency);after=components(out)
    if before!=after:raise AssertionError('LOCAL_REDUCTION_COMPONENT_LOSS')
    remaining=len(edges)-len(removed)
    audit={'method':'strictly-shorter-edge path inside same unit lattice cell','pitch':float(pitch),'nodes':len(points),
        'edges_before':len(edges),'edges_after':remaining,'diagonal_shortcuts_removed':len(removed),'face_edges_removed':0,
        'components_before':before,'components_after':after,'cycle_rank_before':len(edges)-len(points)+before,
        'cycle_rank_after':remaining-len(points)+after,'coordinates_unchanged':True,'components_joined':0,'global_tree_enforcement':False,
        'biological_structure_qualified':False}
    return out,audit
