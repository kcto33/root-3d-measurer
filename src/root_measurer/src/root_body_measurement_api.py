"""End-fixed centerline-body measurements, distinct from full root length."""
import math
import numpy as np
from root_solid_measurement_api import measure_occupied_solid
from root_branch_section_measurement import branch_section_centroid_measurement
from root_axis_radius_profile import axis_radius_profile
from root_graph_measurement_ledger import graph_measurement_ledger
from root_measurement_api import _metric
from root_section_centroid import _covered


def measure_centerline_body(solid,*,coordinate_units,kernel_graph,include_radius_profile=False,
                            frozen_body=None,frozen_radius_profile=None):
    """Measure covered body geometry; never infer GT span or extend endpoints.

    Optional frozen results avoid repeating qualified measurements. Their
    source ID geometry, pinned points and supported shifted edges are checked.
    """
    source=measure_occupied_solid(solid,coordinate_units=coordinate_units,kernel_graph=kernel_graph)
    supported=source['metrics']['supported_geometric_graph_length']['available']
    if frozen_radius_profile is not None and not include_radius_profile:raise ValueError('RADIUS_PROFILE_CACHE_WITHOUT_REQUEST')
    body=None;profile=None;ledger=None
    if supported:
        p=np.asarray(kernel_graph['points'],float);adj=kernel_graph['adjacency']
        body=frozen_body if frozen_body is not None else branch_section_centroid_measurement(solid,p,adj)
        if body['status']!='COMPLETE_EXPERIMENTAL_BRANCH_ASSOCIATED_SECTION_GRAPH':raise ValueError('QUALIFIED_BODY_METHOD_RESULT_REQUIRED')
        if not np.array_equal(np.asarray(body['source_points'],float),p) or body['adjacency']!=[sorted(map(int,a)) for a in adj]:
            raise ValueError('BODY_CACHE_SOURCE_GRAPH_MISMATCH')
        q=np.asarray(body['points'],float)
        if q.shape!=p.shape or not np.isfinite(q).all():raise ValueError('BODY_ID_COORDINATE_SHAPE_MISMATCH')
        pinned=[i for i,a in enumerate(adj) if len(a)!=2]
        if not np.array_equal(q[pinned],p[pinned]):raise ValueError('BODY_PINNED_ENDPOINT_OR_JUNCTION_CHANGED')
        ledger=graph_measurement_ledger(q,adj,coordinate_units=coordinate_units)
        if ledger['source_edge_ids']!=source['source_graph_ledger']['source_edge_ids']:raise ValueError('BODY_SOURCE_EDGE_IDS_CHANGED')
        if ledger['measurement_graph_length']!=body['centroid_graph_ledger']['measurement_graph_length']:raise ValueError('BODY_CACHE_LEDGER_MISMATCH')
        for a,b in ledger['source_edge_ids']:
            if np.array_equal(q[a],p[a]) and np.array_equal(q[b],p[b]):continue
            if not _covered(solid.cells,solid.pitch,q[a],q[b]):raise ValueError('BODY_CHANGED_EDGE_OUTSIDE_SOURCE')
        if include_radius_profile:
            profile=frozen_radius_profile if frozen_radius_profile is not None else axis_radius_profile(solid,q,adj)
            if not np.array_equal(np.asarray(profile['points'],float),q) or profile['source_graph_ledger']['source_edge_ids']!=ledger['source_edge_ids']:
                raise ValueError('RADIUS_CACHE_BODY_GRAPH_MISMATCH')
    return dict(schema_version='root-centerline-body-measurement-v001',units=source['units'],
        source_solid_measurement=source,body_geometry=body,area_equivalent_radius_profile=profile,
        metrics=dict(
            original_source_graph_length=source['metrics']['source_graph_edge_sum'],
            end_fixed_body_graph_length=_metric(None if ledger is None else ledger['measurement_graph_length'],
                'one-pass branch-associated section-centroid graph between unchanged original endpoints/junctions; not full root length',1,
                'SOURCE_COMPLEX_NOT_ADMITTED_AS_SUPPORTED_1D_GRAPH'),
            full_root_centerline_length=_metric(None,'complete anatomical root centerline including true termini',1,'TRUE_ENDPOINT_COVERAGE_AND_ANATOMICAL_CENTERLINE_UNQUALIFIED'),
            physical_cross_section_radius=_metric(None,'biological root cross-section radius',1,'GEOMETRIC_AREA_RADIUS_NOT_BIOLOGICAL_QUALIFICATION'),
            physical_tip_count=_metric(None,'biological root tips',0,'NO_BIOLOGICAL_ENDPOINT_LABELS')),
        original_source_IDs_edges_and_endpoint_coordinates_preserved=True,
        radius_missing_reason=None if profile is not None else ('NOT_REQUESTED' if not include_radius_profile else 'SOURCE_GRAPH_UNSUPPORTED'),
        frozen_body_reused=frozen_body is not None and supported,
        source_geometry_modified=False,GT_reference_or_endpoint_span_accepted=False,
        cap_model_assumed=False,physical_tip_extension=0.,physical_measurements_admitted=False)
