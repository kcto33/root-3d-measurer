"""Reusable occupied-solid measurement; source graphs and model assumptions explicit."""
import hashlib
import itertools
import math
import numpy as np
from root_reconstruction_measurement import ReconstructedSolid
from root_tube_measurement import measure_tube_solid
from root_graph_measurement_ledger import graph_measurement_ledger
from root_observed_volume_measurement import dyadic_cell_volume
from root_section_centroid import _covered
from root_axis_radius_profile import axis_radius_profile
from root_branch_section_measurement import branch_section_centroid_measurement
from root_flat_cap_measurement import flat_cap_graph_measurement
from root_cap_area_profile import cap_area_profile_measurement
from root_measurement_api import _units,_metric


def measure_occupied_solid(solid,*,coordinate_units,kernel_graph=None,cap_model=None,
                           include_area_radius_profile=False,max_dense_cells=350000000):
    """No GT inputs. Optional precomputed graphs avoid repeated skeletonization.

    A supplied kernel graph must use its original voxel-center prefix and
    source-supported edges. Its representation receipt is preserved, not used
    to manufacture a one-dimensional interpretation of a higher complex.
    """
    if not isinstance(solid,ReconstructedSolid):raise ValueError('QUALIFIED_OCCUPIED_SOLID_REQUIRED')
    units=_units(coordinate_units)
    if cap_model not in (None,'flat','round'):raise ValueError('EXPLICIT_SUPPORTED_CAP_MODEL_REQUIRED')
    original=hashlib.sha256(solid.cells.tobytes()).hexdigest()
    if kernel_graph is None:kernel_graph=measure_tube_solid(solid,max_dense_cells=max_dense_cells)
    p=np.asarray(kernel_graph['points'],float).reshape(-1,3);adj=kernel_graph['adjacency']
    ledger=graph_measurement_ledger(p,adj,coordinate_units=coordinate_units)
    cells=np.asarray(kernel_graph['source_skeleton_cells'],np.int64).reshape(-1,3)
    if len(cells)>len(p) or not np.array_equal(p[:len(cells)],(cells+.5)*solid.pitch):
        raise ValueError('SOURCE_KERNEL_VOXEL_CENTER_PREFIX_REQUIRED')
    occupied=set(map(tuple,solid.cells))
    if any(tuple(row) not in occupied for row in cells):raise ValueError('SOURCE_KERNEL_NODE_OUTSIDE_SOLID')
    for a,b in ledger['source_edge_ids']:
        low=np.floor(np.minimum(p[a],p[b])/solid.pitch).astype(np.int64)-1
        high=np.floor(np.maximum(p[a],p[b])/solid.pitch).astype(np.int64)
        # Every cube touching this segment must intersect its bounding box.
        if math.prod(map(int,high-low+1))<=4096:
            boxes=np.asarray([idx for idx in itertools.product(*(range(int(x),int(y)+1) for x,y in zip(low,high)))
                              if idx in occupied],np.int64).reshape(-1,3)
        else:boxes=solid.cells
        if not len(boxes) or not _covered(boxes,solid.pitch,p[a],p[b]):
            raise ValueError('SOURCE_KERNEL_EDGE_OUTSIDE_SOLID')
    representation=kernel_graph['representation']
    admitted=representation['status']=='PASS_SUPPORTED_GEOMETRIC_GRAPH_TOPOLOGY'
    if admitted:
        if representation['graph_component_count']!=ledger['source_components'] or representation['graph_cycle_rank']!=ledger['source_cycle_rank']:
            raise ValueError('SOURCE_REPRESENTATION_LEDGER_MISMATCH')
        if representation['skeleton_voxel_H2']!=0 or representation['skeleton_voxel_H1']!=ledger['source_cycle_rank']:
            raise ValueError('HIGHER_DIMENSIONAL_COMPLEX_NOT_A_GRAPH')
    candidate=None;body=None;profile=None
    if admitted and cap_model is not None:
        body=branch_section_centroid_measurement(solid,p,adj)
        if cap_model=='flat':candidate=flat_cap_graph_measurement(solid,p,adj,flat_caps_assumed=True,body_points=body['points'])
        else:candidate=cap_area_profile_measurement(solid,p,adj,round_caps_assumed=True,body_points=body['points'])
    profile_points=p if candidate is None else np.asarray(candidate['points'],float)
    if admitted and include_area_radius_profile:profile=axis_radius_profile(solid,profile_points,adj)
    complete_cap=candidate is not None and candidate['available_cap_centers']==len(ledger['endpoint_node_ids'])
    fitted_length=candidate['conditional_graph_ledger']['measurement_graph_length'] if complete_cap else None
    volume=dyadic_cell_volume(len(solid.cells),solid.pitch)
    assert hashlib.sha256(solid.cells.tobytes()).hexdigest()==original
    return dict(schema_version='root-occupied-solid-measurement-v001',representation='immutable_occupied_cubical_solid',units=units,
        source_occupied_cells=len(solid.cells),source_cells_sha256=original,pitch=solid.pitch,
        metrics=dict(
            occupied_geometric_volume=dict(volume,definition='union of all supplied canonical occupied cubes',coordinate_unit_power=3),
            source_graph_edge_sum=_metric(ledger['measurement_graph_length'],'diagnostic sum of source graph edges; may represent a higher voxel complex',1),
            supported_geometric_graph_length=_metric(ledger['measurement_graph_length'] if admitted else None,
                'edge sum of a supported 1D representation; centerline accuracy still requires GT validation',1,
                'SOURCE_SKELETON_COMPLEX_NOT_ADMITTED_AS_1D_GRAPH'),
            conditional_cap_graph_length=_metric(fitted_length,'graph length under caller-declared cap/taper geometry and complete endpoint fits',1,
                'MODEL_NOT_SUPPLIED_OR_GRAPH_UNSUPPORTED_OR_ENDPOINT_FITS_MISSING'),
            physical_root_length=_metric(None,'anatomical root centerline length',1,'PHYSICAL_CENTERLINE_AND_ENDPOINTS_NOT_QUALIFIED'),
            physical_cross_section_radius=_metric(None,'anatomical root radius',1,'GEOMETRIC_AREA_RADIUS_NOT_PHYSICAL_QUALIFICATION'),
            physical_root_tips=_metric(None,'biological root tip identities',0,'SOURCE_DEGREE1_IS_NOT_BIOLOGICAL_IDENTITY')),
        source_graph_ledger=ledger,source_representation_receipt=representation,original_kernel_graph=kernel_graph,
        declared_experimental_cap_model=cap_model,body_candidate=body,conditional_cap_candidate=candidate,
        area_equivalent_radius_profile=profile,radius_profile_requested=include_area_radius_profile,
        missing_radius_reason=None if profile is not None else 'NOT_REQUESTED_OR_SOURCE_GRAPH_UNSUPPORTED',
        original_source_geometry_and_IDs_unchanged=True,source_unknowns_filled=0,
        production_cap_metric_adopted=False,GT_reference_accepted_by_inference=False,
        physical_measurements_admitted=False)
