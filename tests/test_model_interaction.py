import math
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
import pytest
from root_model_interaction import ModelInteraction
from root_reconstruction_measurement import ReconstructedSolid


def model(scale=None):
    return dict(candidate_centerline=dict(points=[[0,0,0],[0,0,1],[0,0,2],[1,0,1],[-1,0,1]],
                adjacency=[[1],[0,2,3,4],[1],[1],[1]]),scale_mm_per_model_unit=scale)


def test_continuous_positions_and_user_order_no_route_search():
    engine=ModelInteraction(model(10))
    # Explicitly chooses two connected intervals through a junction.
    value=engine.length([dict(path_id='path00000',start=.123,end=1),
                         dict(path_id='path00002',start=0,end=.456)])
    assert value['length']==pytest.approx((1-.123+.456)*10)
    assert not value['graph_route_searched']
    assert value['polyline'][0]==pytest.approx([0,0,.123])


def test_disconnected_or_overlapping_user_strokes_are_rejected():
    engine=ModelInteraction(model())
    with pytest.raises(ValueError,match='DISCONNECTED'):
        engine.length([dict(path_id='path00000',start=0,end=.3),dict(path_id='path00001',start=.1,end=.5)])
    with pytest.raises(ValueError,match='REPEATED'):
        engine.length([dict(path_id='path00000',start=0,end=.8),dict(path_id='path00000',start=.8,end=.2)])
    with pytest.raises(ValueError,match='ARC_OUTSIDE'):
        engine.length([dict(path_id='path00000',start=-.01,end=1)])


def test_reverse_curved_polyline_length_is_arc_not_chord():
    report=dict(candidate_centerline=dict(points=[[0,0,0],[1,0,0],[1,1,0]],adjacency=[[1],[0,2],[1]]))
    value=ModelInteraction(report).length([dict(path_id='path00000',start=1.8,end=.2)])
    assert value['length']==pytest.approx(1.6)
    assert len(value['polyline'])==3


def test_fresh_arbitrary_section_and_actual_outline(monkeypatch):
    from root_cubical_disk_section import interior_boundary_section_class
    from root_native_entity_measurement import NativeSectionIndex
    from root_branch_section_measurement import BranchOwnership
    engine=ModelInteraction(dict(candidate_centerline=dict(points=[[.5,.5,0],[.5,.5,4]],adjacency=[[1],[0]]),
                                 entity_binding=dict(cells_sha256='test-original')))
    engine.solid=ReconstructedSolid(np.array([[0,0,i] for i in range(4)]),1.,{})
    engine.engine=NativeSectionIndex(engine.solid,interior_boundary_section_class())
    engine.owner=BranchOwnership(engine.points,engine.adjacency)
    value=engine.radius(dict(path_id='path00000',arc=1.23456789))
    assert value['radius']==pytest.approx(math.sqrt(1/math.pi))
    assert value['position'][2]==pytest.approx(1.23456789)
    assert value['outline'][0]==value['outline'][-1]
    assert len(value['outline'])>=5
    assert not value['stored_samples_interpolated']
    value2=engine.radius(dict(path_id='path00000',position=[.9,.5,2.34567]))
    assert value2['arc']==pytest.approx(2.34567)
    assert value2['requested_to_axis_distance']==pytest.approx(.4)
    with pytest.raises(ValueError,match='NONZERO'):
        engine.radius(dict(path_id='path00000',arc=1.,normal=[0,0,0]))


def test_connected_section_radius_keeps_hole_area_and_labels():
    from root_cubical_disk_section import interior_boundary_section_class
    from root_native_entity_measurement import NativeSectionIndex
    from root_branch_section_measurement import BranchOwnership
    engine=ModelInteraction(dict(candidate_centerline=dict(points=[[.5,.5,0],[.5,.5,4]],adjacency=[[1],[0]]),
                                 entity_binding=dict(cells_sha256='ring-section')))
    cells=np.array([[x,y,z] for x in range(3) for y in range(3) for z in range(4) if (x,y)!=(1,1)])
    engine.solid=ReconstructedSolid(cells,1.,{})
    engine.engine=NativeSectionIndex(engine.solid,interior_boundary_section_class())
    engine.owner=BranchOwnership(engine.points,engine.adjacency)
    result=engine.radius(dict(path_id='path00000',arc=1.234))
    assert result['radius']==pytest.approx(math.sqrt(8/math.pi))
    assert result['area']==8
    assert result['single_branch_radius'] is None
    assert result['section_classification']=='section_with_holes'
    assert len(result['outlines'])==2
    assert result['holes_filled'] is False


def test_mixed_section_keeps_complete_original_area_without_branch_radius():
    from root_cubical_disk_section import interior_boundary_section_class
    from root_native_entity_measurement import NativeSectionIndex
    from root_branch_section_measurement import BranchOwnership
    engine=ModelInteraction(dict(candidate_centerline=dict(
        points=[[.5,.5,0],[.5,.5,2],[.5,.5,4],[1.5,.5,2]],adjacency=[[1],[0,2,3],[1],[1]]),
        entity_binding=dict(cells_sha256='mixed-section')))
    engine.solid=ReconstructedSolid(np.array([[x,0,z] for x in range(2) for z in range(4)]),1.,{})
    engine.engine=NativeSectionIndex(engine.solid,interior_boundary_section_class())
    engine.owner=BranchOwnership(engine.points,engine.adjacency)
    result=engine.radius(dict(path_id='path00000',arc=1.234))
    assert result['area']==2
    assert result['radius']==pytest.approx(math.sqrt(2/math.pi))
    assert result['single_branch_radius'] is None
    assert result['section_classification']=='mixed_branch'


def test_point_contact_does_not_become_a_certified_arbitrary_radius():
    from root_cubical_disk_section import interior_boundary_section_class
    from root_native_entity_measurement import NativeSectionIndex
    from root_branch_section_measurement import BranchOwnership
    engine=ModelInteraction(dict(candidate_centerline=dict(points=[[.5,.5,0],[.5,.5,4]],adjacency=[[1],[0]]),
                                 entity_binding=dict(cells_sha256='point-contact-cut')))
    engine.solid=ReconstructedSolid(np.array([[x,x,z] for x in range(2) for z in range(4)]),1.,{})
    engine.engine=NativeSectionIndex(engine.solid,interior_boundary_section_class())
    engine.owner=BranchOwnership(engine.points,engine.adjacency)
    result=engine.radius(dict(path_id='path00000',arc=1.234))
    assert result['radius'] is None
    assert result['section_classification']=='unqualified_section'
    assert 'CONTACT' in result['status']
