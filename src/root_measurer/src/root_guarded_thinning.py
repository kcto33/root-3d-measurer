"""Measurement-only guard against Lee thinning deleting whole components.

The Euler LUT/octant convention follows scikit-image 0.26.0 (BSD-3-Clause),
_skeletonize_lee_cy.pyx.in. The installed library is never modified. Its original
six face directions are retained; all deletion conditions are rechecked against
the current neighborhood, including endpoints and isolated points. This is not
a certification of anatomical centerlines or biological length.
"""
from functools import lru_cache
import hashlib
import numpy as np
from scipy import ndimage
from skimage.morphology import skeletonize as compiled_lee
from skimage.morphology._skeletonize_lee_cy import fill_Euler_LUT

LUT=np.asarray(fill_Euler_LUT(),dtype=np.int32)
if hashlib.sha256(LUT.tobytes()).hexdigest()!='02a9bf53102d0d16aadf12586421741d67b820140dc9b926197043e73496bd10':
    raise RuntimeError('UNQUALIFIED_EULER_LUT_RUNTIME')
LUT.flags.writeable=False
OCTANTS=((2,1,11,10,5,4,14),(0,9,3,12,1,10,4),
    (8,7,17,16,5,4,14),(6,15,7,16,3,12,4),
    (20,23,19,22,11,14,10),(18,21,9,12,19,22,10),
    (26,23,17,14,25,22,16),(24,25,15,16,21,22,12))
FACE_ORDER=((1,-1),(1,1),(2,1),(2,-1),(0,1),(0,-1))
POWERS=(1<<np.arange(27,dtype=np.int64))
CONNECT26=np.ones((3,3,3),bool)


def neighborhood_code(patch):
    # Upstream convention orders (plane, column, row), with center index 13.
    return int(np.dot(np.asarray(patch,dtype=np.int64).transpose(0,2,1).reshape(-1),POWERS))


def euler_deletion_invariant(code):
    total=0
    for octant in OCTANTS:
        index=1
        for j,neighbor in enumerate(octant):
            if (code>>neighbor)&1:index|=1<<(7-j)
        total+=int(LUT[index])
    return total==0


@lru_cache(maxsize=200000)
def deletable_current(code):
    if not ((code>>13)&1):return False
    # A current endpoint or singleton must survive, even if it was not one at
    # initial parallel-candidate collection. No source component may vanish.
    if code.bit_count()-1<=1:return False
    if not euler_deletion_invariant(code):return False
    flat=np.array([(code>>i)&1 for i in range(27)],bool);flat[13]=False
    patch=flat.reshape(3,3,3).transpose(0,2,1)
    return ndimage.label(patch,structure=CONNECT26)[1]==1


def _input(volume,max_dense_cells):
    source=np.asarray(volume)
    if source.ndim!=3 or source.dtype!=np.bool_:raise ValueError('BOOLEAN_3D_MEASUREMENT_VOLUME_REQUIRED')
    if any(n==0 for n in source.shape):raise ValueError('NONZERO_DOMAIN_REQUIRED')
    if int(np.prod(np.array(source.shape,dtype=object)+2))>max_dense_cells:
        raise ValueError('UNCHANGED_MEASUREMENT_DENSE_GUARD')
    return source


def sequential_lee_thinning(volume,max_dense_cells=350000000):
    source=_input(volume,max_dense_cells);before=hashlib.sha256(source.tobytes()).digest()
    image=np.pad(source,1);removed=0;passes=0
    while True:
        changed=False;passes+=1
        for axis,step in FACE_ORDER:
            core=image[1:-1,1:-1,1:-1]
            neighbor=[slice(1,-1)]*3;neighbor[axis]=slice(1+step,image.shape[axis]-1+step)
            boundary=core&~image[tuple(neighbor)]
            candidates=[]
            for p in np.argwhere(boundary)+1:
                a,b,c=map(int,p);code=neighborhood_code(image[a-1:a+2,b-1:b+2,c-1:c+2])
                if deletable_current(code):candidates.append((a,b,c))
            for a,b,c in candidates:
                code=neighborhood_code(image[a-1:a+2,b-1:b+2,c-1:c+2])
                if deletable_current(code):image[a,b,c]=False;removed+=1;changed=True
        if not changed:break
        if passes>int(np.count_nonzero(source))+1:raise AssertionError('THINNING_CONVERGENCE_BOUND')
    result=image[1:-1,1:-1,1:-1].copy()
    assert not np.any(result&~source) and np.count_nonzero(source)-np.count_nonzero(result)==removed
    assert hashlib.sha256(source.tobytes()).digest()==before
    old_count=ndimage.label(source,structure=CONNECT26)[1];new_count=ndimage.label(result,structure=CONNECT26)[1]
    assert old_count==new_count,'DERIVED_THINNING_COMPONENTS_CHANGED'
    return result,dict(status='COMPLETE_SEQUENTIAL_CURRENT_NEIGHBORHOOD_THINNING',
        source_components=int(old_count),derived_components=int(new_count),passes=passes,
        derived_voxels_removed=removed,source_material_unchanged=True,occupied_material_added=0,
        current_endpoint_Euler_and_simple_recheck=True,full_runtime_Betti_check=False,
        physical_centerline_qualified=False,biological_length_qualified=False)


def guarded_skeletonize(volume,method='lee',*,max_dense_cells=350000000,receipt=None):
    """Keep compiled Lee exact except for entirely vanished source components."""
    if method!='lee':raise ValueError('ONLY_REGISTERED_LEE_METHOD_SUPPORTED')
    source=_input(volume,max_dense_cells);source_hash=hashlib.sha256(source.tobytes()).hexdigest()
    result=compiled_lee(source,method='lee')
    labels,count=ndimage.label(source,structure=CONNECT26)
    present=set(map(int,np.unique(labels[result])))-{0}
    missing=sorted(set(range(1,count+1))-present);repairs=[]
    if missing:
        boxes=ndimage.find_objects(labels)
        result=result.copy()
        for component in missing:
            box=boxes[component-1];part=labels[box]==component
            repaired,audit=sequential_lee_thinning(part,max_dense_cells=max_dense_cells)
            result[box]|=repaired
            repairs.append(dict(source_component=component,source_voxels=int(np.count_nonzero(part)),
                recovered_derived_voxels=int(np.count_nonzero(repaired)),audit=audit))
    assert not np.any(result&~source) and hashlib.sha256(source.tobytes()).hexdigest()==source_hash
    after_count=ndimage.label(result,structure=CONNECT26)[1];assert after_count==count
    audit=dict(status='PASS_NONEMPTY_COMPONENT_MEASUREMENT_GUARD',source_components=int(count),
        derived_components=int(after_count),missing_compiled_components=missing,fallback_components=repairs,
        original_compiled_output_unchanged=not missing,source_material_unchanged=True,
        occupied_material_added=0,components_joined=0,genuine_cycles_forced_to_tree=False,
        universal_topological_or_biological_qualification=False)
    if receipt is not None:receipt.update(audit)
    return result
