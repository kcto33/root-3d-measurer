"""Point-contact qualification is consistent in oblique and axial sections."""
import itertools
import numpy as np
from root_measurer.runtime import activate
activate()
from root_reconstruction_measurement import ReconstructedSolid
from root_native_entity_measurement import NativeSectionIndex,axis_aligned_section
from root_cubical_disk_section import interior_boundary_section_class


def section(cells,point,normal):
    solid=ReconstructedSolid(np.asarray(cells,np.int64),1.,{'units':'model'})
    result=axis_aligned_section(solid,point,normal)
    return result if result is not None else NativeSectionIndex(solid,interior_boundary_section_class()).section(point,normal)


def test_oblique_two_cube_point_contact_is_not_qualified_disk():
    value=section([[0,0,0],[1,1,0]],[.5,.5,.5],[.2,.1,1.])
    assert value['status']=='NONMANIFOLD_OR_CONTACT_SECTION'
    assert value['radius'] is None
    assert value['point_contact_with_unselected_polygons']


def test_contact_rejection_is_stable_under_coordinate_permutations():
    cells=np.array([[0,0,0],[1,1,0]])
    point=np.array([.5,.5,.5]);normal=np.array([.2,.1,1.])
    for order in itertools.permutations(range(3)):
        value=section(cells[:,order],point[list(order)],normal[list(order)])
        assert value['status']=='NONMANIFOLD_OR_CONTACT_SECTION'


def test_unrelated_separate_component_does_not_invalidate_query():
    value=section([[0,0,0],[4,4,0]],[.5,.5,.5],[.2,.1,1.])
    assert value['radius'] is not None
    assert value['status']=='CLOSED_DISK_CUBICAL_SECTION'


def test_face_connected_solid_disk_remains_measurable():
    value=section([[0,0,0],[1,0,0]],[.5,.5,.5],[.2,.1,1.])
    assert value['radius'] is not None
    assert not value['point_contact_with_unselected_polygons']
