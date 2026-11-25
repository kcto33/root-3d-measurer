"""Lee0.26 rules/order on occupied-neighbor storage, without a dense empty box.

Algorithm reference (BSD-3-Clause): scikit-image v0.26.0,
src/skimage/morphology/_skeletonize_lee_cy.pyx.in, Ignacio Arganda-Carreras
and the scikit-image contributors. Existing pinned Euler LUT/octants and
whole-vanished-component sequential recovery remain authoritative.

This is a storage implementation, not a new skeleton or endpoint estimator.
The dense measuring path always includes a blank border, so all six 3D
directions apply even when occupied coordinates span only one plane.
"""
from functools import lru_cache
import hashlib,time
import numpy as np
from root_guarded_thinning import LUT,OCTANTS,deletable_current

OFFSETS=np.asarray([(p,r,c) for p in (-1,0,1) for c in (-1,0,1) for r in (-1,0,1)],np.int64)
FACE_BITS=(12,14,16,10,22,4)
CENTER=1<<13
REVERSE=(np.uint32(1)<<np.arange(26,-1,-1,dtype=np.uint32))
ADJ=[]
for q in OFFSETS:
    mask=0
    for j,p in enumerate(OFFSETS):
        if j!=13 and np.max(np.abs(q-p))<=1:mask|=1<<j
    ADJ.append(mask)

@lru_cache(maxsize=200000)
def simple_point(code):
    remaining=int(code)&~CENTER
    if not remaining:return True # upstream octree returns True for zero labels
    pending=remaining&-remaining;remaining^=pending
    while pending:
        bit=pending&-pending;pending^=bit;j=bit.bit_length()-1
        found=remaining&ADJ[j];remaining^=found;pending|=found
    return not remaining

def _table(cells,max_dense_cells):
    p=np.asarray(cells)
    if p.ndim!=2 or p.shape[1]!=3 or not np.issubdtype(p.dtype,np.integer):raise ValueError('SPARSE_LEE_INTEGER_CELLS_REQUIRED')
    n=len(p)
    if n>=np.iinfo(np.int32).max:raise ValueError('SPARSE_LEE_NEIGHBOR_INDEX_GUARD')
    # Twelve bytes per formerly admitted dense cell conservatively covers the
    # original bool/padded-result/int-label neighborhood. Keep the same caller
    # resource parameter; account for arrays actually allocated instead.
    estimated=n*(27*4+24+4+1+32)
    if estimated>int(max_dense_cells)*12:raise ValueError('SPARSE_LEE_STORAGE_GUARD')
    order=np.lexsort((p[:,2],p[:,1],p[:,0]));p=p[order].astype(np.int64,copy=True)
    if not n:return p,np.empty((0,27),np.int32),dict(estimated_bytes=estimated)
    low=p.min(0)-2;shape=p.max(0)-low+3
    if np.prod(shape.astype(object))>np.iinfo(np.int64).max:raise ValueError('SPARSE_LEE_INTEGER_DOMAIN_GUARD')
    stride=np.array([int(shape[1])*int(shape[2]),int(shape[2]),1],np.int64);keys=(p-low)@stride
    if np.any(np.diff(keys)<=0):raise ValueError('SPARSE_LEE_UNIQUE_CELLS_REQUIRED')
    neighbors=np.full((n,27),n,np.int32)
    for j,delta in enumerate(OFFSETS):
        target=keys+int(delta@stride);where=np.searchsorted(keys,target);hit=where<n;ii=np.flatnonzero(hit);hit[ii]=keys[where[ii]]==target[ii]
        neighbors[hit,j]=where[hit]
    return p,neighbors,dict(estimated_bytes=estimated,neighbor_table_bytes=neighbors.nbytes,occupied_nodes=n,
        empty_box_cells=int(np.prod(shape.astype(object))),allocated_dense_domain_cells=0,caller_dense_limit_unchanged=int(max_dense_cells),
        allocation_byte_limit=int(max_dense_cells)*12)

def _euler_ok(codes):
    total=np.zeros(len(codes),np.int32)
    for octant in OCTANTS:
        index=np.ones(len(codes),np.uint16)
        for j,neighbor in enumerate(octant):index|=(((codes>>neighbor)&1)<<(7-j)).astype(np.uint16)
        total+=LUT[index]
    return total==0

def _run(cells,max_dense_cells,sequential=False,progress=None):
    start=time.monotonic();p,neighbors,plan=_table(cells,max_dense_cells);n=len(p)
    if not n:return p,dict(passes=0,removed=0,storage=plan)
    alive=np.ones(n+1,bool);alive[n]=False;codes=np.zeros(n,np.uint32)
    for j in range(27):codes|=(neighbors[:,j]!=n).astype(np.uint32)<<j
    passes=0;removed=0
    while True:
        changed=0;passes+=1
        for face in FACE_BITS:
            ids=np.flatnonzero(alive[:n]&((codes&(1<<face))==0));cc=codes[ids]
            if not len(ids):continue
            if sequential:
                unique,inverse=np.unique(cc,return_inverse=True)
                keep=np.fromiter((deletable_current(int(c)) for c in unique),bool,len(unique))[inverse]
                ids=ids[keep];predicate=deletable_current
            else:
                ids=ids[(np.bitwise_count(cc)!=2)&_euler_ok(cc)];cc=codes[ids]
                unique,inverse=np.unique(cc,return_inverse=True)
                keep=np.fromiter((simple_point(int(c)) for c in unique),bool,len(unique))[inverse]
                ids=ids[keep];predicate=simple_point
            # ids are in the same (plane,row,column) collection/recheck order as
            # upstream Cython. Recheck only simple-point in compiled Lee mode.
            for i in ids:
                if predicate(int(codes[i])):
                    alive[i]=False;js=neighbors[i];valid=alive[js]
                    codes[js[valid]]&=~REVERSE[valid];changed+=1
        removed+=changed
        if progress is not None:progress(dict(passes=passes,remaining=int(alive[:n].sum()),removed=removed,seconds=time.monotonic()-start))
        if not changed:break
        if passes>n+1:raise AssertionError('SPARSE_LEE_CONVERGENCE_GUARD')
    return p[alive[:n]],dict(passes=passes,removed=removed,seconds=time.monotonic()-start,storage=plan,
        estimator='same_Lee026_global_order' if not sequential else 'same_sequential_current_neighbor_guard')

def sparse_raw_lee(cells,*,max_dense_cells=350000000,progress=None):
    return _run(cells,max_dense_cells,progress=progress)

def sparse_guarded_lee(cells,*,max_dense_cells=350000000,receipt=None,_verified_connected=False,progress=None):
    from root_cubical_surface import component_labels
    source=np.asarray(cells,np.int64);before=hashlib.sha256(source.tobytes()).hexdigest()
    if not len(source):
        if receipt is not None:receipt.update(status='PASS_SPARSE_LEE_EMPTY',source_components=0,derived_components=0)
        return source.copy()
    if _verified_connected:count=1;labels=np.zeros(len(source),np.int32)
    else:count,labels=component_labels(source)
    nodes,raw=_run(source,max_dense_cells,progress=progress)
    # Map the sorted result back to source component IDs without a dense label.
    lo=source.min(0)-1;shape=source.max(0)-lo+2;stride=np.array([int(shape[1])*int(shape[2]),int(shape[2]),1],np.int64)
    keys=(source-lo)@stride;order=np.argsort(keys);nk=(nodes-lo)@stride;at=np.searchsorted(keys[order],nk)
    if len(at) and (np.any(at>=len(source)) or np.any(keys[order[at]]!=nk)):raise AssertionError('SPARSE_LEE_ADDED_MATERIAL')
    present=set(map(int,labels[order[at]]));missing=sorted(set(range(count))-present);repairs=[];parts=[nodes]
    for component in missing:
        part=source[labels==component];recovered,audit=_run(part,max_dense_cells,sequential=True,progress=progress)
        parts.append(recovered);repairs.append(dict(source_component=component,source_voxels=len(part),recovered_derived_voxels=len(recovered),audit=audit))
    nodes=np.concatenate(parts);nodes=nodes[np.lexsort((nodes[:,2],nodes[:,1],nodes[:,0]))]
    after=component_labels(nodes)[0]
    if after!=count:raise AssertionError('SPARSE_LEE_COMPONENT_PRESERVATION_FAILED')
    if hashlib.sha256(source.tobytes()).hexdigest()!=before:raise AssertionError('SPARSE_LEE_SOURCE_CHANGED')
    audit=dict(status='PASS_NONEMPTY_COMPONENT_MEASUREMENT_GUARD',storage='occupied27_neighbor_global_order_Lee026',source_components=count,derived_components=after,
        missing_compiled_components=missing,fallback_components=repairs,source_material_unchanged=True,occupied_material_added=0,components_joined=0,
        genuine_cycles_forced_to_tree=False,physical_centerline_qualified=False,new_deletion_conditions=False,raw_passes=raw,
        dense_empty_box_allocated=False,dense_compiled_execution_replaced_only_for_storage=True)
    if receipt is not None:receipt.update(audit)
    return nodes
