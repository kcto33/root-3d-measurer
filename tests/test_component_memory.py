"""Storage decomposition must match the original dense estimator exactly."""
import itertools
import numpy as np
import pytest
from root_reconstruction_measurement import ReconstructedSolid
from root_tube_measurement import measure_tube_solid
from root_sparse_solid_topology import sparse_solid_topology
from root_native_axis_order import original_frame_kernel


def block(shape,offset):
    return np.indices(shape).reshape(3,-1).T+offset


@pytest.mark.parametrize('order',list(itertools.permutations(range(3))))
def test_disjoint_storage_matches_full_dense_graph(order):
    cells=np.concatenate([block((3,4,15),(0,0,0)),block((4,3,11),(32,25,28)),
                          np.array([[55,10,15]])])
    solid=ReconstructedSolid(cells,.01,{})
    dense=original_frame_kernel(solid,order,2000000)
    bounded=original_frame_kernel(solid,order,5000)
    for key in ['source_skeleton_cells','points','adjacency','cubical_boundary_radii']:
        assert dense[key]==bounded[key]
    assert bounded['thinning']['global_dense_allocation_avoided']
    assert bounded['thinning']['peak_component_dense_cells']<=5000
    assert np.array_equal(solid.cells,cells)


def test_cavity_and_cycle_additivity_matches_global_background():
    shell=block((5,5,5),(0,0,0));shell=shell[np.any((shell==0)|(shell==4),axis=1)]
    ring=block((9,9,3),(0,0,0));ring=ring[~((ring[:,0]>=3)&(ring[:,0]<=5)&(ring[:,1]>=3)&(ring[:,1]<=5))]
    cells=np.concatenate([shell,ring+[30,40,50]])
    dense=sparse_solid_topology(cells,max_domain_cells=2000000)
    bounded=sparse_solid_topology(cells,max_domain_cells=5000)
    for key in ['foreground_components','enclosed_empty_components','tunnel_dimension','euler','enclosed_empty_sizes']:
        assert bounded[key]==dense[key]
    assert bounded['enclosed_empty_components']==1 and bounded['tunnel_dimension']==1
    assert bounded['allocated_dense_cells']==0


def test_single_component_still_obeys_actual_storage_guard():
    solid=ReconstructedSolid(block((3,3,100),(0,0,0)),1.,{})
    with pytest.raises(ValueError,match='SPARSE_LEE_STORAGE_GUARD'):
        measure_tube_solid(solid,max_dense_cells=1000)


def test_nested_isolated_material_is_subtracted_from_global_cavity():
    shell=block((5,5,5),(0,0,0));shell=shell[np.any((shell==0)|(shell==4),axis=1)]
    cells=np.concatenate([shell,[[2,2,2],[80,90,100]]])
    reference=sparse_solid_topology(cells,max_domain_cells=2000000)
    bounded=sparse_solid_topology(cells,max_domain_cells=5000)
    assert bounded['enclosed_empty_sizes']==reference['enclosed_empty_sizes']==[26]
    for key in ['foreground_components','enclosed_empty_components','tunnel_dimension','euler']:
        assert reference[key]==bounded[key]
    assert bounded['background_cells_accounted']==bounded['domain_cells']-len(cells)
