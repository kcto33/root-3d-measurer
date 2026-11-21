"""Exact estimator equality, including ordering, fallback and resource bounds."""
import itertools
import numpy as np
import pytest
from skimage.morphology import skeletonize
from root_guarded_thinning import guarded_skeletonize,sequential_lee_thinning
from root_sparse_lee import sparse_raw_lee,sparse_guarded_lee,_run
from root_reconstruction_measurement import ReconstructedSolid
from root_native_axis_order import original_frame_kernel
from root_sparse_solid_topology import sparse_solid_topology

def dense(cells,guard=False):
    low=cells.min(0)-1;shape=cells.max(0)-low+2;volume=np.zeros(shape,bool);volume[tuple((cells-low).T)]=True
    before=volume.copy();audit={}
    result=guarded_skeletonize(volume,receipt=audit) if guard else skeletonize(volume,method='lee')
    assert np.array_equal(before,volume)
    return np.argwhere(result)+low,audit

@pytest.mark.parametrize('seed',range(24))
def test_sparse_raw_matches_compiled_every_cell(seed):
    rng=np.random.default_rng(seed);volume=rng.random((7,8,9))<(.15+.03*seed);cells=np.argwhere(volume)
    cells=cells[rng.permutation(len(cells))]-[20,30,40];copy=cells.copy()
    expected,_=dense(cells);actual,_=sparse_raw_lee(cells)
    assert np.array_equal(expected,actual)
    assert np.array_equal(copy,cells)

@pytest.mark.parametrize('shape',[(2,2,10),(4,4,12),(3,4,11),(5,5,5)])
def test_guarded_recovery_matches_old_sequential_rule(shape):
    cells=np.indices(shape).reshape(3,-1).T;expected,audit=dense(cells,True);proof={}
    actual=sparse_guarded_lee(cells,receipt=proof)
    assert np.array_equal(expected,actual)
    assert bool(audit['fallback_components'])==bool(proof['fallback_components'])
    volume=np.ones(shape,bool);out,unused=sequential_lee_thinning(volume)
    sequential,unused=_run(cells,350000000,sequential=True)
    assert np.array_equal(np.argwhere(out),sequential)

def connected_path(n=30):
    return np.unique(np.concatenate([np.c_[np.arange(n),np.zeros(n,int),np.zeros(n,int)],
        np.c_[np.full(n,n-1),np.arange(n),np.zeros(n,int)],
        np.c_[np.full(n,n-1),np.full(n,n-1),np.arange(n)]]),axis=0)

@pytest.mark.parametrize('order',list(itertools.permutations(range(3))))
def test_large_connected_storage_preserves_full_kernel(order):
    cells=connected_path();solid=ReconstructedSolid(cells,.01,{})
    reference=original_frame_kernel(solid,order,2000000)
    actual=original_frame_kernel(solid,order,5000)
    for key in ['source_skeleton_cells','points','adjacency','cubical_boundary_radii','ledger','representation']:
        assert reference[key]==actual[key]
    assert actual['thinning']['global_dense_allocation_avoided']
    assert actual['thinning']['peak_component_dense_cells']==0
    assert np.array_equal(solid.cells,cells)

def test_billion_empty_domain_no_dense_allocation():
    cells=np.tile(np.arange(1001)[:,None],(1,3));proof={}
    actual=sparse_guarded_lee(cells,receipt=proof)
    assert np.array_equal(actual,cells)
    assert proof['raw_passes']['storage']['empty_box_cells']>1000000000
    assert proof['raw_passes']['storage']['allocated_dense_domain_cells']==0

def test_exact_global_cavities_do_not_require_component_dense_box():
    cells=connected_path();a=sparse_solid_topology(cells,max_domain_cells=2000000)
    b=sparse_solid_topology(cells,max_domain_cells=5000)
    for key in ['foreground_components','enclosed_empty_components','enclosed_empty_sizes','tunnel_dimension','euler','background_cells_accounted']:
        assert a[key]==b[key]

def test_sparse_guard_does_not_raise_resource_budget():
    with pytest.raises(ValueError,match='SPARSE_LEE_STORAGE_GUARD'):
        sparse_guarded_lee(np.indices((3,3,100)).reshape(3,-1).T,max_dense_cells=1000)
