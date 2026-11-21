"""Independent arc-conservation, coordinate and missing-support controls."""
import copy
import json
import math

import pytest

from root_trait_distributions import (build_trait_distributions, distribution_settings,
                                      flatten_distribution_rows, is_qualified_radius_sample)


def report(points, edges, *, scale=None, bounds=None):
    adjacency = [[] for _ in points]
    for a, b in edges:
        adjacency[a].append(b)
        adjacency[b].append(a)
    result = dict(status='NATIVE_OUTPUT_WITH_EXPLICIT_GAPS', coordinate_units='model_units',
                  scale_mm_per_model_unit=scale,
                  candidate_centerline=dict(points=points, adjacency=adjacency),
                  geometric_endpoints=[dict(node_id=i, status='geometric_surface_estimate')
                                       for i, row in enumerate(adjacency) if len(row) == 1])
    if bounds is not None:
        result['model_bounds'] = bounds
    return result


def sample(edge, radius, weight, position, **extra):
    return dict(edge_nodes=list(edge), radius=radius, arc_weight=weight, position=position,
                available=True, status='CLOSED_DISK_CUBICAL_SECTION', **extra)


def lengths(value, family):
    return [r['axis_length'] for r in value['distributions'][family]['bins']]


def test_diameter_uses_unique_arc_weights_including_missing_not_sample_count():
    r = report([[0, 0, 0], [4, 0, 0], [10, 0, 0]], [[0, 1], [1, 2]])
    samples = [sample((0, 1), .2, 4, [2, 0, 0]), sample((1, 2), .75, 3, [4, 0, 0]),
               sample((1, 2), 8, 3, [10, 0, 0], connected_section_classification='mixed_branch')]
    result = build_trait_distributions(r, samples=samples, options={'diameter_edges': [.5, 1.]})
    assert lengths(result, 'diameter') == pytest.approx([4, 0, 3])
    hist = result['distributions']['diameter']
    assert hist['covered_arc_fraction'] == pytest.approx(.7)
    assert hist['unassigned_axis_length'] == pytest.approx(3)
    assert sum(x['fraction_of_total_axis'] for x in hist['bins']) == pytest.approx(.7)
    assert sum(x['fraction_of_assigned_axis'] for x in hist['bins']) == pytest.approx(1)
    assert hist['status'] == 'PARTIAL_SAMPLE_SUPPORTED_ESTIMATE'


def test_actual_measurement_tables_midpoint_and_node_halfedge_contracts():
    r = report([[0, 0, 0], [0, 0, 4]], [[0, 1]])
    r['radius_curve'] = dict(sampling='edge_midpoint', rows=[
        dict(node_ids=[0, 1], radius=.5, available=True, status='CLOSED_DISK_CUBICAL_SECTION')])
    midpoint = build_trait_distributions(r)
    assert lengths(midpoint, 'diameter') == [0, 0, 4, 0]
    r['radius_curve'] = dict(sampling='node_halfedge', rows=[
        dict(node_id=0, radius=.25, available=True, status='CLOSED_DISK_CUBICAL_SECTION'),
        dict(node_id=1, radius=None, available=False, status='MISSING_MIXED_BRANCH_SECTION')])
    node = build_trait_distributions(r)
    assert lengths(node, 'diameter') == [0, 2, 0, 0]
    assert node['sample_coverage']['qualified_radius_arc_fraction'] == .5


def test_dimensional_options_use_output_units_once_and_unscaled_units_are_explicit():
    r = report([[0, 0, 0], [3, 0, 4]], [[0, 1]], scale=10)
    samples = [sample((0, 1), .1, 5, [1.5, 0, 2])]
    a = build_trait_distributions(r, samples=samples,
                                  options={'diameter_edges': [1, 3], 'origin': [0, 0, 0], 'radial_edges': [10, 20]})
    assert a['unit'] == 'mm' and a['total_axis_length'] == 50
    assert lengths(a, 'diameter') == [0, 50, 0]
    assert lengths(a, 'radial') == pytest.approx([50 / 3] * 3)
    r['scale_mm_per_model_unit'] = None
    b = build_trait_distributions(r, samples=samples, options={'diameter_edges': [1, 3]})
    assert b['unit'] == 'model_units' and lengths(b, 'diameter') == [5, 0, 0]


def test_inclination_and_axial_azimuth_invariant_to_edge_reversal():
    p = [[0, 0, 0], [2, 0, 0], [0, 2, 2], [0, 0, 3]]
    r = report(p, [[0, 1], [0, 2], [0, 3]])
    reverse = report(list(reversed(p)), [[3, 2], [3, 1], [3, 0]])
    options = dict(inclination_edges_deg=[0, 30, 60, 90], azimuth_edges_deg=[0, 60, 120, 180])
    a = build_trait_distributions(r, samples=[], options=options)
    b = build_trait_distributions(reverse, samples=[], options=options)
    assert lengths(a, 'inclination') == pytest.approx([3, math.sqrt(8), 2])
    assert lengths(a, 'azimuth') == pytest.approx([2, math.sqrt(8), 0])
    assert lengths(a, 'azimuth') == pytest.approx(lengths(b, 'azimuth'))
    assert lengths(a, 'inclination') == pytest.approx(lengths(b, 'inclination'))
    assert a['distributions']['azimuth']['unassigned_axis_length'] == pytest.approx(3)


def test_vertical_azimuth_is_unavailable_not_a_zero_degree_direction():
    r = report([[0, 0, 0], [0, 0, 8]], [[0, 1]])
    result = build_trait_distributions(r, samples=[])
    assert all(x is None for x in lengths(result, 'azimuth'))
    assert result['distributions']['azimuth']['unassigned_axis_length'] == 8
    assert result['distributions']['azimuth']['covered_arc_fraction'] == 0
    assert sum(lengths(result, 'inclination')) == 8


def test_depth_direction_splits_edges_exactly_but_diameter_depth_is_sample_supported():
    r = report([[0, 0, 0], [6, 0, 8]], [[0, 1]],
               bounds={'minimum': [0, 0, 0], 'maximum': [6, 1, 8]})
    result = build_trait_distributions(r, samples=[sample((0, 1), .25, 10, [3, 0, 4])],
        options={'depth_edges': [0, 2, 4, 6, 8], 'azimuth_edges_deg': [0, 90, 180], 'diameter_edges': [1]})
    assert lengths(result, 'depth_azimuth') == pytest.approx([2.5, 0] * 4)
    # At the internal boundary 4, recorded midpoint support belongs to [4,6).
    assert lengths(result, 'diameter_depth') == [0, 0, 10, 0, 0, 0, 0, 0]
    assert 'sample' in result['distributions']['diameter_depth']['definition']
    r['model_bounds'] = {'minimum': [0, 0, 2], 'maximum': [6, 1, 6]}
    limited = build_trait_distributions(r, samples=[], options={'depth_edges': [0, 2, 4]})
    assert limited['distributions']['depth_azimuth']['assigned_axis_length'] == pytest.approx(5)
    assert limited['distributions']['depth_azimuth']['unassigned_axis_length'] == pytest.approx(5)


def test_depth_frame_axis_sign_and_fallback_are_explicit():
    r = report([[0, 0, 0], [8, 6, 0]], [[0, 1]])
    a = build_trait_distributions(r, samples=[], options={'vertical_axis': 'x', 'down_sign': 1, 'depth_bins': 2})
    assert a['depth_reference']['source'] == 'EXTRACTED_AXIS_BOUNDS'
    assert a['depth_reference']['top_coordinate'] == 0
    assert a['horizontal_axes'] == ['y', 'z']
    assert a['distributions']['depth_azimuth']['assigned_axis_length'] == pytest.approx(10)
    b = build_trait_distributions(r, samples=[], bounds={'minimum': [-2, 0, 0], 'maximum': [10, 6, 0]},
                                 options={'vertical_axis': 'x', 'down_sign': -1})
    assert b['depth_reference']['source'] == 'SOURCE_MATERIAL_BOUNDS'
    assert b['depth_reference']['top_coordinate'] == 10


def test_radial_annuli_crossed_twice_and_constant_boundary_are_exact():
    r = report([[-3, 0, 0], [3, 0, 0], [2, 0, 0], [2, 0, 4]], [[0, 1], [2, 3]])
    result = build_trait_distributions(r, samples=[], options={'origin': [0, 0, 0], 'radial_edges': [1, 2]})
    # Horizontal diameter contributes 2 in each annulus, vertical at r=2 adds 4 to the upper class.
    assert lengths(result, 'radial') == pytest.approx([2, 2, 6])
    assert result['distributions']['radial']['covered_arc_fraction'] == 1
    subdivided = report([[-3, 0, 0], [0, 0, 0], [3, 0, 0], [2, 0, 0], [2, 0, 4]],
                        [[0, 1], [1, 2], [3, 4]])
    other = build_trait_distributions(subdivided, samples=[], options={'origin': [0, 0, 0], 'radial_edges': [1, 2]})
    assert lengths(other, 'radial') == pytest.approx(lengths(result, 'radial'))


def test_missing_radius_origin_and_axis_stay_unavailable_not_measured_zero():
    result = build_trait_distributions(report([[0, 0, 0], [2, 0, 0]], [[0, 1]]), samples=[])
    assert all(x is None for x in lengths(result, 'diameter'))
    assert result['distributions']['diameter']['covered_arc_fraction'] == 0
    assert all(x is None for x in lengths(result, 'radial'))
    assert result['distributions']['radial']['status'] == 'UNAVAILABLE_EXPLICIT_RADIAL_ORIGIN_REQUIRED'
    empty = build_trait_distributions(report([], []), samples=[])
    assert empty['total_axis_length'] is None
    assert all(h['status'] == 'UNAVAILABLE_NO_MEASURABLE_AXIS' for h in empty['distributions'].values())
    assert any(r['family'] == 'diameter_depth' and r['length'] is None for r in flatten_distribution_rows(empty))


def test_partial_status_json_safety_and_no_mutation():
    r = report([[0, 0, 0], [3, 4, 0]], [[0, 1]])
    samples = [sample((0, 1), .3, 5, [1.5, 2, 0])]
    options = {'origin': [1, 2, 3]}
    before = copy.deepcopy((r, samples, options))
    result = build_trait_distributions(r, samples=samples, options=options, partial_axis=True)
    assert result['distributions']['diameter']['status'] == 'PARTIAL_SAMPLE_SUPPORTED_ESTIMATE'
    assert result['distributions']['inclination']['status'] == 'PARTIAL_OBSERVED_AXIS_DESCRIPTOR'
    json.dumps(result, allow_nan=False)
    assert (r, samples, options) == before
    rows = flatten_distribution_rows(result)
    assert len({tuple(sorted(row)) for row in rows}) == 1
    assert next(row for row in rows if row['family'] == 'diameter')['unit'] == 'model_units'


@pytest.mark.parametrize('change', [dict(available=False), dict(radius=None), dict(radius=-1),
    dict(radius=float('nan')), dict(radius=True), dict(status='MISSING_MIXED_BRANCH_SECTION'),
    dict(status='NONDISK_SECTION_HOLES_NOT_FILLED'), dict(connected_section_classification='mixed_branch'),
    dict(connected_section_classification='section_with_holes')])
def test_radius_qualification_never_falls_back_to_whole_connected_section(change):
    value = dict(radius=.2, available=True, status='CLOSED_DISK_CUBICAL_SECTION', connected_section_radius=5)
    value.update(change)
    assert not is_qualified_radius_sample(value)


@pytest.mark.parametrize('option', [dict(diameter_edges=[1, 1]), dict(diameter_edges=[0]),
    dict(radial_edges=[-1]), dict(origin=[0, 1]), dict(origin=[0, False, 1]),
    dict(depth_edges=[1, 0]), dict(depth_bins=True), dict(vertical_axis='a'),
    dict(down_sign=0), dict(azimuth_edges_deg=[0, 360]), dict(inclination_edges_deg=[0, 180]),
    dict(unknown=1)])
def test_invalid_configuration_rejected(option):
    with pytest.raises(ValueError):
        distribution_settings(option)


def test_duplicate_sample_support_cannot_exceed_a_unique_source_edge():
    r = report([[0, 0, 0], [2, 0, 0]], [[0, 1]])
    item = sample((0, 1), .2, 2, [1, 0, 0])
    with pytest.raises(ValueError, match='SAMPLE_ARC_EXCEEDS'):
        build_trait_distributions(r, samples=[item, item])


def test_position_bearing_requires_origin_and_distinguishes_off_axis_vertical_roots():
    r = report([[0, 0, 0], [0, 0, 4], [1, 1, 0], [1, 1, 3]], [[0, 1], [2, 3]])
    missing = build_trait_distributions(r, samples=[])
    for family in ('bearing', 'depth_bearing'):
        assert missing['distributions'][family]['status'] == 'UNAVAILABLE_EXPLICIT_BEARING_ORIGIN_REQUIRED'
        assert all(value is None for value in lengths(missing, family))
    value = build_trait_distributions(r, samples=[], options={'origin': [0, 0, 0]})
    assert lengths(value, 'bearing') == pytest.approx([0, 3, 0, 0, 0, 0, 0, 0])
    histogram = value['distributions']['bearing']
    assert histogram['covered_arc_fraction'] == pytest.approx(3 / 7)
    assert histogram['unassigned_reasons']['reference_axis_bearing_undefined'] == 4
    assert histogram['unassigned_axis_length'] == 4
    # Position bearing is defined for an off-axis vertical root even though its
    # projected tangent has no azimuth direction at all.
    assert value['distributions']['azimuth']['covered_arc_fraction'] == 0


def test_bearing_crosses_zero_and_multiple_rays_with_exact_lengths():
    value = build_trait_distributions(report([[1, -3, 0], [1, 3, 0]], [[0, 1]]), samples=[],
                                     options={'origin': [0, 0, 0]})
    assert lengths(value, 'bearing') == pytest.approx([1, 2, 0, 0, 0, 0, 2, 1])
    assert value['distributions']['bearing']['covered_arc_fraction'] == pytest.approx(1)
    # Non-cardinal boundaries have independently known y=tan(theta) crossings.
    nonuniform = build_trait_distributions(report([[1, 0, 0], [1, math.sqrt(3), 0]], [[0, 1]]),
        samples=[], options={'origin': [0, 0, 0], 'bearing_edges_deg': [0, 30, 60, 360]})
    assert lengths(nonuniform, 'bearing') == pytest.approx([1 / math.sqrt(3), 2 / math.sqrt(3), 0])


@pytest.mark.parametrize('points,expected', [
    ([[-3, 0, 0], [3, 0, 0]], [3, 0, 0, 0, 3, 0, 0, 0]),
    ([[-2, -2, 0], [2, 2, 0]], [0, math.sqrt(8), 0, 0, 0, math.sqrt(8), 0, 0]),
])
def test_bearing_axis_crossing_has_zero_undefined_length(points, expected):
    value = build_trait_distributions(report(points, [[0, 1]]), samples=[], options={'origin': [0, 0, 0]})
    assert lengths(value, 'bearing') == pytest.approx(expected)
    assert value['distributions']['bearing']['unassigned_axis_length'] == pytest.approx(0)
    assert value['distributions']['bearing']['unassigned_reasons']['reference_axis_bearing_undefined'] == 0


def test_depth_bearing_splits_both_rays_and_depth_planes_and_reports_outside_arc():
    r = report([[1, -1, 0], [1, 3, 4]], [[0, 1]])
    length = math.sqrt(32)
    value = build_trait_distributions(r, samples=[], options={'origin': [0, 0, 0], 'depth_edges': [0, 1, 2, 3, 4]})
    expected = [0.] * 32
    for index in (1, 8 + 1, 16, 24 + 7):
        expected[index] = length / 4
    assert lengths(value, 'depth_bearing') == pytest.approx(expected)
    assert sum(lengths(value, 'depth_bearing')) == pytest.approx(length)
    limited = build_trait_distributions(r, samples=[], options={'origin': [0, 0, 0], 'depth_edges': [1, 3]})
    assert lengths(limited, 'depth_bearing') == pytest.approx([length / 4, length / 4, 0, 0, 0, 0, 0, 0])
    assert limited['distributions']['depth_bearing']['unassigned_reasons']['outside_depth_bins'] == pytest.approx(length / 2)


@pytest.mark.parametrize('vertical_axis', ['x', 'y', 'z'])
@pytest.mark.parametrize('down_sign', [-1, 1])
def test_bearing_respects_selected_vertical_frame_and_depth_sign(vertical_axis, down_sign):
    axis = 'xyz'.index(vertical_axis)
    horizontal = [i for i in range(3) if i != axis]
    points = [[0., 0., 0.], [0., 0., 0.]]
    for point, depth, second in zip(points, [0, 4], [-1, 3]):
        point[axis] = depth
        point[horizontal[0]] = 1
        point[horizontal[1]] = second
    value = build_trait_distributions(report(points, [[0, 1]]), samples=[], options=dict(
        origin=[0, 0, 0], vertical_axis=vertical_axis, down_sign=down_sign, depth_edges=[0, 1, 2, 3, 4]))
    assert value['horizontal_axes'] == ['xyz'[i] for i in horizontal]
    assert lengths(value, 'bearing') == pytest.approx([math.sqrt(2), math.sqrt(8), 0, 0, 0, 0, 0, math.sqrt(2)])
    layers = [lengths(value, 'depth_bearing')[i:i + 8] for i in range(0, 32, 8)]
    selected_sectors = [7, 0, 1, 1] if down_sign == 1 else [1, 1, 0, 7]
    for layer, sector in zip(layers, selected_sectors):
        expected = [0.] * 8; expected[sector] = math.sqrt(2)
        assert layer == pytest.approx(expected)


def test_position_bearing_subdivision_reversal_scale_and_translation_invariance():
    p = [[1, -1, 0], [1, 3, 4]]
    options = dict(origin=[0, 0, 0], depth_edges=[0, 2, 4, 6, 8])
    reference = build_trait_distributions(report(p, [[0, 1]], scale=2), samples=[], options=options)
    times = [0, .13, .4, .7, 1]
    offset = [32, -16, 8]
    # Changed model units and origin yield the same physical geometry, translated.
    q = [[(a + t * (b - a)) * 4 + shift for a, b, shift in zip(p[0], p[1], offset)] for t in times]
    transformed = build_trait_distributions(report(q, [[i, i + 1] for i in range(4)], scale=.5), samples=[],
        options=dict(origin=[c * .5 for c in offset], depth_edges=options['depth_edges']))
    reverse = build_trait_distributions(report(list(reversed(p)), [[0, 1]], scale=2), samples=[], options=options)
    for family in ('bearing', 'depth_bearing'):
        assert lengths(transformed, family) == pytest.approx(lengths(reference, family))
        assert lengths(reverse, family) == pytest.approx(lengths(reference, family))
    json.dumps(transformed, allow_nan=False)
    row = next(r for r in flatten_distribution_rows(transformed) if r['family'] == 'bearing')
    assert row['bearing_deg_start'] == 0 and row['bearing_deg_end'] == 45


@pytest.mark.parametrize('edges', [[0, 180], [1, 360], [0, 45, 45, 360], [0, 361]])
def test_position_bearing_requires_a_complete_ordered_circle(edges):
    with pytest.raises(ValueError):
        distribution_settings({'bearing_edges_deg': edges})


def test_position_bearing_joint_bin_budget_is_enforced():
    with pytest.raises(ValueError, match='BIN_BUDGET'):
        distribution_settings({'bearing_edges_deg': list(range(361)), 'depth_bins': 200})
