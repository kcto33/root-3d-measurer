"""Exact foreground26/background6 cubical topology using empty-row intervals.

Every empty grid cell in the padded finite domain belongs to one interval.
Adjacent rows connect iff their integer intervals overlap. This changes only
storage, not occupancy, adjacency, cavity size or the 350m-cell guard.
"""
import math
import numpy as np
from root_reconstruction_measurement import _indices
from root_cubical_surface import component_labels
from root_volume_topology import cubical_euler


def sparse_solid_topology(cells, max_domain_cells=350000000, max_intervals=4000000):
    cells=_indices(cells)
    counts=cubical_euler(cells)
    if not len(cells):
        return dict(counts,foreground_components=0,enclosed_empty_components=0,tunnel_dimension=0,
                    enclosed_empty_sizes=[],component_occupied_counts=[],allocated_dense_cells=0,
                    interval_nodes=0,physical_topology_qualified=False,
                    unknown_as_empty_for_derived_topology=True)
    lo=cells.min(axis=0)-1;shape=cells.max(axis=0)-lo+2
    domain=math.prod(map(int,shape))
    if domain>max_domain_cells:raise ValueError('SPARSE_TOPOLOGY_DOMAIN_GUARD')
    axis=int(np.argmax(shape));order=[axis]+[i for i in range(3) if i!=axis]
    xyz=(cells-lo)[:,order];sx,sy,sz=map(int,shape[order])
    capacity=len(cells)+sy*sz+1
    if capacity>max_intervals:raise ValueError('SPARSE_TOPOLOGY_INTERVAL_GUARD')
    rowkeys=xyz[:,1]*sz+xyz[:,2]
    permutation=np.lexsort((xyz[:,0],rowkeys));keys=rowkeys[permutation];xs=xyz[permutation,0]
    starts=np.r_[0,np.flatnonzero(np.diff(keys))+1];ends=np.r_[starts[1:],len(keys)]
    occupied_rows={int(keys[a]):xs[a:b] for a,b in zip(starts,ends)}
    parent=np.arange(capacity,dtype=np.int64);sizes=np.zeros(capacity,np.int64)
    next_id=1
    def find(a):
        while parent[a]!=a:
            parent[a]=parent[parent[a]];a=int(parent[a])
        return a
    def join(a,b):
        a,b=find(a),find(b)
        if a==b:return
        if sizes[a]<sizes[b]:a,b=b,a
        parent[b]=a;sizes[a]+=sizes[b];sizes[b]=0
    def overlap(a,b):
        i=j=0
        while i<len(a) and j<len(b):
            al,ah,ai=a[i];bl,bh,bi=b[j]
            if al<=bh and bl<=ah:join(ai,bi)
            if ah<bh:i+=1
            else:j+=1
    previous_slab=[None]*sz;total_background=0
    for y in range(sy):
        last=None
        for z in range(sz):
            row=[];left=0
            occupied=occupied_rows.get(y*sz+z,())
            for x in occupied:
                x=int(x)
                if left<x:
                    row.append((left,x-1,next_id));sizes[next_id]=x-left
                    total_background+=x-left;next_id+=1
                left=x+1
            if left<sx:
                row.append((left,sx-1,next_id));sizes[next_id]=sx-left
                total_background+=sx-left;next_id+=1
            for a,b,i in row:
                if a==0 or b==sx-1 or y==0 or y==sy-1 or z==0 or z==sz-1:
                    join(0,i)
            above=previous_slab[z]
            if above is not None:overlap(row,above)
            if last is not None:overlap(row,last)
            previous_slab[z]=row;last=row
    if total_background!=domain-len(cells):raise AssertionError('BACKGROUND_CELL_ACCOUNTING')
    exterior=find(0)
    cavity_roots=[i for i in range(1,next_id) if parent[i]==i and i!=exterior]
    cavities=sorted((int(sizes[i]) for i in cavity_roots),reverse=True)
    foreground,labels=component_labels(cells)
    h1=foreground+len(cavities)-counts['euler']
    if h1<0:raise AssertionError('NEGATIVE_CUBICAL_H1')
    return dict(counts,foreground_components=foreground,enclosed_empty_components=len(cavities),
                tunnel_dimension=h1,enclosed_empty_sizes=cavities,
                component_occupied_counts=sorted(np.bincount(labels).tolist(),reverse=True),
                allocated_dense_cells=0,interval_nodes=next_id-1,storage_axis=axis,
                finite_padded_shape=shape.tolist(),domain_cells=domain,
                background_cells_accounted=total_background,foreground_connectivity=26,
                background_connectivity=6,unknown_as_empty_for_derived_topology=True,
                physical_topology_qualified=False,no_geometry_modified=True)
