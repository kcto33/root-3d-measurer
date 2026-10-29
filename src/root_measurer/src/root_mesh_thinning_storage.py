"""Same compiled Lee, sparse component audit; no new deletion conditions."""
import hashlib
import math
import numpy as np
from skimage.morphology import skeletonize
from root_cubical_surface import component_labels
from root_guarded_thinning import sequential_lee_thinning


def byte_plan(shape,occupied,budget_bytes):
    shape=np.asarray(shape,dtype=np.int64)
    if shape.shape!=(3,) or np.any(shape<=0):raise ValueError('POSITIVE_3D_STORAGE_SHAPE_REQUIRED')
    # Source+padded result coexist. Audit hashes/containment use bounded slabs,
    # not full-volume copies, and sparse component adjacency is released before
    # the compiled call. Reserve covers source mesh/cache/Python allocations.
    dense=math.prod(map(int,shape));padded=math.prod(map(int,shape+2))
    slab=min(dense,32*int(shape[1])*int(shape[2]))
    sparse_upper=1024*int(occupied)
    reserve=1536*1024**2
    estimate=max(dense+padded+3*slab,dense+sparse_upper)+reserve
    return dict(domain_cells=dense,source_bool_bytes=dense,padded_bool_bytes=padded,
                sparse_audit_upper_bytes=sparse_upper,reserve_bytes=reserve,
                estimated_peak_bytes=estimate,budget_bytes=budget_bytes,
                admitted=estimate<=budget_bytes,allocation_model='two dense bool arrays; sparse26 component audit;32-plane slab checks')


def digest(volume):
    h=hashlib.sha256()
    for begin in range(0,volume.shape[0],32):h.update(np.ascontiguousarray(volume[begin:begin+32]).data)
    return h.hexdigest()


def guarded_lee_sparse_audit(volume,method='lee',*,max_dense_cells=None,receipt=None,
                             budget_bytes=8192*1024**2,source_cells=None):
    if method!='lee' or volume.dtype!=np.bool_ or volume.ndim!=3:raise ValueError('SAME_BOOLEAN_3D_LEE_REQUIRED')
    if not volume.flags.c_contiguous:raise ValueError('C_CONTIGUOUS_VOLUME_REQUIRED_FOR_BYTE_PLAN')
    cells=np.argwhere(volume) if source_cells is None else np.asarray(source_cells,np.int64)
    if cells.ndim!=2 or cells.shape[1]!=3 or len(np.unique(cells,axis=0))!=len(cells):raise ValueError('UNIQUE_SPARSE_SOURCE_CELLS_REQUIRED')
    if len(cells) and not np.all(volume[tuple(cells.T)]):raise ValueError('SPARSE_AUDIT_SOURCE_CELL_MISMATCH')
    if np.count_nonzero(volume)!=len(cells):raise ValueError('SPARSE_AUDIT_ALL_SOURCE_CELLS_REQUIRED')
    plan=byte_plan(volume.shape,len(cells),budget_bytes)
    if not plan['admitted']:raise ValueError('MESH_LEE_BYTE_PLAN_EXCEEDS_BUDGET')
    if not len(cells):
        result=skeletonize(volume,method='lee')
        if receipt is not None:receipt.update(status='PASS_SAME_COMPILED_LEE_EMPTY_SOURCE',source_components=0,derived_components=0,
            source_material_unchanged=True,existing_compiled_Lee_reused=True,new_deletion_rules=False,byte_plan=plan)
        return result
    before=digest(volume);count,labels=component_labels(cells)
    result=skeletonize(volume,method='lee')
    nodes=np.argwhere(result)
    shape=volume.shape;keys=np.ravel_multi_index(cells.T,shape);order=np.argsort(keys);sorted_keys=keys[order]
    node_keys=np.ravel_multi_index(nodes.T,shape);where=np.searchsorted(sorted_keys,node_keys)
    if np.any(where>=len(cells)) or np.any(sorted_keys[where]!=node_keys):raise AssertionError('LEE_ADDED_MATERIAL')
    present=set(map(int,labels[order[where]]));missing=sorted(set(range(count))-present);repairs=[]
    for component in missing:
        part=cells[labels==component];low=part.min(axis=0);high=part.max(axis=0)+1;part_shape=high-low
        # Same sequential recovery as old guard; do not run an oversized repair.
        if math.prod(map(int,part_shape+2))*12>budget_bytes-plan['source_bool_bytes']-plan['padded_bool_bytes']-plan['reserve_bytes']:
            raise ValueError('MESH_MISSING_COMPONENT_RECOVERY_BYTE_BUDGET')
        local=np.zeros(tuple(part_shape),bool);local[tuple((part-low).T)]=True
        recovered,audit=sequential_lee_thinning(local,max_dense_cells=math.prod(map(int,part_shape+2)))
        box=tuple(slice(int(a),int(b)) for a,b in zip(low,high));result[box]|=recovered
        repairs.append(dict(source_component=component,source_voxels=len(part),recovered_derived_voxels=int(recovered.sum()),audit=audit))
    for begin in range(0,volume.shape[0],32):
        if np.any(result[begin:begin+32]&~volume[begin:begin+32]):raise AssertionError('LEE_ADDED_MATERIAL')
    after_count,_=component_labels(np.argwhere(result))
    if after_count!=count or digest(volume)!=before:raise AssertionError('LEE_COMPONENT_OR_SOURCE_CHANGED')
    audit=dict(status='PASS_SAME_COMPILED_LEE_SPARSE_COMPONENT_AUDIT',source_components=count,derived_components=after_count,
        missing_compiled_components=missing,fallback_components=repairs,original_compiled_output_unchanged=not missing,
        source_material_unchanged=True,new_deletion_rules=False,existing_compiled_Lee_reused=True,
        existing_sparse_component_algorithm_reused=True,byte_plan=plan,
        universal_topological_or_biological_qualification=False)
    if receipt is not None:receipt.update(audit)
    return result
