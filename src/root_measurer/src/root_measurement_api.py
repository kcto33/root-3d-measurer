"""Source-preserving reconstruction measurements, with per-metric availability."""
import math
from root_triangle_measurement_input import TriangleInput
from root_triangle_area_binary_exact import binary_exact_triangle_area_sum
from root_triangle_topology import indexed_surface_topology
from root_graph_measurement_ledger import graph_measurement_ledger
from root_endpoint_measurement import measure_graph_endpoints
from root_cubical_section_normal import RangeSafeCubicalSectionMeasurer


def _units(label):
    if not isinstance(label,str) or not label.strip():raise ValueError('EXPLICIT_COORDINATE_UNIT_LABEL_REQUIRED')
    return dict(coordinate_unit_label=label,unit_label_is_caller_declaration=True,
                physical_scale_independently_verified=False)


def _metric(value,definition,power,reason=None):
    return dict(value=value,definition=definition,coordinate_unit_power=power,
                available=value is not None,missing_reason=reason if value is None else None,
                physical_or_biological_qualified=False)


def measure_triangle_input(mesh,*,coordinate_units):
    """Measure existing triangles; never invent an occupied solid or centerline."""
    units=_units(coordinate_units)
    if not isinstance(mesh,TriangleInput):raise ValueError('QUALIFIED_TRIANGLE_INPUT_REQUIRED')
    area=binary_exact_triangle_area_sum(mesh.vertices,mesh.faces)
    topology=indexed_surface_topology(mesh)
    reasons=['MATERIAL_UNION_AND_EMBEDDING_NOT_VERIFIED']
    if topology['boundary_edges']:reasons.append('OPEN_INDEXED_SURFACE')
    if topology['invalid_link_vertex_indices']:reasons.append('INVALID_VERTEX_LINKS')
    if area['zero_area_faces']:reasons.append('DEGENERATE_GEOMETRIC_FACES')
    return dict(schema_version='root-reconstruction-measurement-v001',
        representation='source_triangle_surface',source_sha256=mesh.source_sha256,units=units,
        metrics=dict(
            source_triangle_area=_metric(area['triangle_area_sum'],'sum of all source triangle areas, including overlap/internal faces',2),
            material_union_outer_area=_metric(None,'area of the exterior of the occupied material union',2,'MATERIAL_UNION_OUTER_SURFACE_NOT_SUPPLIED'),
            material_volume=_metric(None,'volume of occupied material union',3,';'.join(reasons)),
            physical_root_length=_metric(None,'anatomical root centerline length',1,'ADMITTED_ANATOMICAL_CENTERLINE_NOT_SUPPLIED'),
            physical_cross_section_radius=_metric(None,'anatomical cross-section radius',1,'ADMITTED_RADIUS_EVIDENCE_NOT_SUPPLIED'),
            physical_tip_count=_metric(None,'biological root tip count',0,'BIOLOGICAL_TIP_IDENTITIES_NOT_SUPPLIED')),
        triangle_area_diagnostics=area,indexed_surface_topology=topology,
        material_topology=None,geometry_or_input_modified=False,physical_measurements_admitted=False)


def measure_source_graph(points,adjacency,pitch,*,coordinate_units,node_radii=None,
                         radius_definition=None,display_points=None):
    """Measure a supplied ID graph; no skeleton inference or biological labeling."""
    units=_units(coordinate_units)
    if node_radii is not None and (not isinstance(radius_definition,str) or not radius_definition.strip()):
        raise ValueError('SUPPLIED_RADIUS_DEFINITION_REQUIRED')
    # Validate strict IDs first: this rejects duplicate edges before endpoint traversal.
    ledger=graph_measurement_ledger(points,adjacency,coordinate_units=coordinate_units,
        radii=node_radii,display_points=display_points)
    endpoints=measure_graph_endpoints(points,adjacency,pitch,coordinate_units=coordinate_units)
    if not math.isclose(ledger['measurement_graph_length'],endpoints['source_graph_length'],rel_tol=1e-12,abs_tol=0.):
        raise AssertionError('SOURCE_GRAPH_LENGTH_RECEIPT_DISAGREEMENT')
    assert ledger['source_nodes']==endpoints['source_nodes'] and ledger['source_edges']==endpoints['source_edges']
    return dict(schema_version='root-reconstruction-measurement-v001',representation='source_ID_graph',units=units,
        metrics=dict(
            source_graph_length=_metric(ledger['measurement_graph_length'],'sum of original unique undirected graph edges',1),
            display_curve_length=_metric(ledger['display_curve_length'],'optional display curve length, separate from source metric',1,'DISPLAY_POINTS_NOT_SUPPLIED'),
            geometric_degree1_count=_metric(endpoints['graph_endpoint_count'],'number of degree1 source graph vertices, including untyped collar',0),
            physical_root_length=_metric(None,'anatomical root centerline length',1,'CENTERLINE_GEOMETRY_AND_PHYSICAL_ENDPOINTS_NOT_ADMITTED'),
            physical_tip_count=_metric(None,'biological root tip count',0,'DEGREE1_IDS_ARE_NOT_BIOLOGICAL_TIP_LABELS'),
            physical_cross_section_radius=_metric(None,'anatomical cross-section radius',1,'SUPPLIED_GRAPH_RADII_ARE_NOT_BIOLOGICAL_RADIUS_QUALIFICATION')),
        supplied_radius_definition=radius_definition,supplied_radii_present=node_radii is not None,
        source_ID_ledger=ledger,geometric_endpoints=endpoints,
        graph_or_coordinates_or_endpoint_positions_modified=False,physical_measurements_admitted=False)


class ReconstructionSectionMeasurement:
    """Reuse one immutable cubical-solid section engine over many queries."""
    def __init__(self,solid,*,coordinate_units):
        self.units=_units(coordinate_units)
        self.measurer=RangeSafeCubicalSectionMeasurer(solid)

    def section(self,point,normal):
        result=self.measurer.section(point,normal)
        area=result['area'];radius=result['radius'];reason=result['status']
        if area is not None and (not math.isfinite(area) or area<=0):
            area=None;radius=None;reason='SECTION_AREA_FLOAT_RANGE'
        if radius is not None and (not math.isfinite(radius) or radius<=0):
            radius=None;reason='SECTION_RADIUS_FLOAT_RANGE'
        return dict(schema_version='root-reconstruction-measurement-v001',representation='cubical_solid_section',
            units=self.units,metrics=dict(
                geometric_section_area=_metric(area,'selected component of source cubical solid plane section',2,reason),
                geometric_area_equivalent_radius=_metric(radius,'sqrt(area/pi) for one admitted geometric disk section',1,reason),
                physical_cross_section_radius=_metric(None,'anatomical cross-section radius',1,'GEOMETRIC_SECTION_NOT_BIOLOGICAL_QUALIFICATION')),
            original_section_diagnostic=result,unknown_cells_filled=0,source_geometry_modified=False,
            observation_support_available=False,physical_measurements_admitted=False)
