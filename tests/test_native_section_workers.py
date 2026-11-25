"""Local support and spawned query ordering preserve the serial geometry."""
from itertools import product
import numpy as np
from root_reconstruction_measurement import ReconstructedSolid,cubical_ray_exit
from root_section_centroid import _covered,covered_segment
from root_native_section_workers import sections
from root_native_entity_measurement import indexed_body_measurement
from root_branch_section_measurement import branch_section_centroid_measurement


def solid():
    return ReconstructedSolid(np.array(list(product(range(-2,3),range(-2,3),range(12)))),.1,{})


def test_local_support_matches_full_scan_on_faces_gaps_and_long_segments():
    s=solid()
    rng=np.random.default_rng(5205)
    pairs=[(rng.uniform(-.3,1.3,3),rng.uniform(-.3,1.3,3)) for _ in range(80)]
    pairs += [(np.array([-.2,0,0]),np.array([.3,0,1.2])),(np.zeros(3),np.array([20.,0,0]))]
    for a,b in pairs:assert covered_segment(s,a,b)==_covered(s.cells,s.pitch,a,b)
    occupied=s.occupied
    cubical_ray_exit(s,[0,0,.5],[0,0,1])
    assert s.occupied is occupied


def test_spawned_planes_preserve_serial_order_and_values():
    s=solid()
    queries=[([0.,0.,z],n) for z in [.25,.45,.65] for n in [[0,0,1],[.2,.1,1]]]
    assert sections(s,queries,2)==sections(s,queries,1)


def test_spawned_body_matches_serial_node_ids_and_geometry():
    s=solid();points=[[.05,.05,.15+i*.1] for i in range(9)]
    adj=[[j for j in [i-1,i+1] if 0<=j<len(points)] for i in range(len(points))]
    one=indexed_body_measurement(s,points,adj,section_workers=1)
    two=indexed_body_measurement(s,points,adj,section_workers=2)
    assert one==two


def test_body_fast_path_keeps_legacy_lattice_boundary_refusals():
    s=solid()
    for x,y,z in [(0.,0.,.15),(.05,.05,.2)]:
        points=[[x,y,z+i*.1] for i in range(9)]
        adj=[[j for j in [i-1,i+1] if 0<=j<len(points)] for i in range(len(points))]
        legacy=branch_section_centroid_measurement(s,points,adj)
        fast=indexed_body_measurement(s,points,adj,section_workers=2)
        assert legacy['points']==fast['points']
        assert legacy['section_statuses']==fast['section_statuses']
        assert legacy['area_equivalent_radii']==fast['area_equivalent_radii']


