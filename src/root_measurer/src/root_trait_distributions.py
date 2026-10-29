"""Length distributions on unchanged measured edges and qualified radius samples.

Diameter classes describe sample-supported arc weights, not a continuous radius
field. Direction, position-bearing and radial distributions integrate each unique
graph edge once. Position bearing requires an explicit horizontal reference point.
No parent direction, biological root type, soil volume or missing radius is inferred.
"""
import bisect
import math

import numpy as np


METHOD = 'root-trait-distributions-v001'
_DEFAULTS = dict(vertical_axis='z', down_sign=-1, diameter_edges=[.5, 1., 2.],
                 depth_bins=10, depth_edges=None,
                 inclination_edges_deg=[0., 15., 30., 45., 60., 75., 90.],
                 azimuth_edges_deg=list(range(0, 181, 15)),
                 bearing_edges_deg=list(range(0, 361, 45)),
                 radial_edges=[1., 2., 5., 10.], origin=None)
_MAX_BINS = 65536


def is_qualified_radius_sample(sample):
    """Keep the existing single-branch radius contract; never substitute a union radius.

``ordinary_section`` is the current native connected-section label. It does not
replace the main ``radius``/``available`` qualification. Explicit mixed, holed,
missing or otherwise unqualified sections remain excluded even if they happen to
carry a number in a radius field.
    """
    available = sample.get('available')
    radius = sample.get('radius')
    if not isinstance(available, (bool, np.bool_)) or not available:
        return False
    if isinstance(radius, (bool, np.bool_)) or radius is None:
        return False
    try:
        if not math.isfinite(radius) or radius < 0:
            return False
    except (TypeError, ValueError):
        return False
    classification = str(sample.get('connected_section_classification', '')).lower()
    if classification in ('mixed_branch', 'section_with_holes', 'unqualified_section'):
        return False
    status = str(sample.get('status', '')).upper()
    if status.startswith(('MISSING', 'UNAVAILABLE', 'UNQUALIFIED', 'NOT_MEASURED',
                          'BLOCKED', 'INVALID', 'NONDISK')):
        return False
    return not any(word in status for word in ('MIXED_BRANCH', 'WITH_HOLES', 'SECTION_HOLES'))


def _numbers(value, name, *, positive=False, at_least=0):
    if not isinstance(value, (list, tuple)) or not at_least <= len(value) <= 512:
        raise ValueError(name + '_FINITE_ORDERED_LIST_REQUIRED')
    if any(isinstance(x, (bool, np.bool_)) or not isinstance(x, (int, float, np.number))
           or not math.isfinite(x) for x in value):
        raise ValueError(name + '_FINITE_ORDERED_LIST_REQUIRED')
    result = [float(x) for x in value]
    if positive and any(x <= 0 for x in result):
        raise ValueError(name + '_POSITIVE_THRESHOLDS_REQUIRED')
    if any(b <= a for a, b in zip(result, result[1:])):
        raise ValueError(name + '_STRICTLY_INCREASING_REQUIRED')
    return result


def distribution_settings(options=None):
    """All dimensional options use the output unit: mm if calibrated, else model units.

``diameter_edges`` and ``radial_edges`` are positive internal thresholds; zero
and an unbounded upper class are added automatically. Other edges are complete
boundaries. ``origin`` is an explicit point in output coordinates, never an
automatically inferred root base.
    """
    options = {} if options is None else options
    if not isinstance(options, dict) or set(options) - set(_DEFAULTS):
        raise ValueError('UNKNOWN_TRAIT_DISTRIBUTION_OPTION')
    result = {k: list(v) if isinstance(v, list) else v for k, v in _DEFAULTS.items()}
    result.update(options)
    if result['vertical_axis'] not in ('x', 'y', 'z'):
        raise ValueError('VERTICAL_AXIS_X_Y_Z_REQUIRED')
    if isinstance(result['down_sign'], (bool, np.bool_)) or result['down_sign'] not in (-1, 1):
        raise ValueError('DOWN_SIGN_MINUS_OR_PLUS_ONE_REQUIRED')
    count = result['depth_bins']
    if isinstance(count, (bool, np.bool_)) or not isinstance(count, (int, np.integer)) or not 1 <= count <= 512:
        raise ValueError('DEPTH_BINS_INTEGER_1_TO_512_REQUIRED')
    result['depth_bins'] = int(count)
    for name in ('diameter_edges', 'radial_edges'):
        result[name] = _numbers(result[name], name.upper(), positive=True)
    for name, maximum in (('inclination_edges_deg', 90.), ('azimuth_edges_deg', 180.),
                          ('bearing_edges_deg', 360.)):
        result[name] = _numbers(result[name], name.upper(), at_least=2)
        if result[name][0] != 0 or result[name][-1] != maximum:
            raise ValueError(name.upper() + '_FULL_ANGULAR_RANGE_REQUIRED')
    if result['depth_edges'] is not None:
        result['depth_edges'] = _numbers(result['depth_edges'], 'DEPTH_EDGES', at_least=2)
    if result['origin'] is not None:
        origin = result['origin']
        if not isinstance(origin, (list, tuple)) or len(origin) != 3:
            raise ValueError('FINITE_XYZ_ORIGIN_IN_OUTPUT_UNITS_REQUIRED')
        if any(isinstance(x, (bool, np.bool_)) or not isinstance(x, (int, float, np.number))
               or not math.isfinite(x) for x in origin):
            raise ValueError('FINITE_XYZ_ORIGIN_IN_OUTPUT_UNITS_REQUIRED')
        result['origin'] = [float(x) for x in origin]
    depth_count = len(result['depth_edges']) - 1 if result['depth_edges'] else count
    if max((len(result['diameter_edges']) + 1) * depth_count,
           (len(result['azimuth_edges_deg']) - 1) * depth_count,
           (len(result['bearing_edges_deg']) - 1) * depth_count) > _MAX_BINS:
        raise ValueError('TRAIT_DISTRIBUTION_BIN_BUDGET_EXCEEDED')
    return result


def _unit(report):
    scale = report.get('scale_mm_per_model_unit')
    if scale is not None and (isinstance(scale, (bool, np.bool_)) or
                              not isinstance(scale, (int, float, np.number)) or
                              not math.isfinite(scale) or scale <= 0):
        raise ValueError('POSITIVE_FINITE_SCALE_REQUIRED')
    unit = 'mm' if scale is not None else report.get('coordinate_units', 'model_units')
    if not isinstance(unit, str) or not unit.strip():
        raise ValueError('EXPLICIT_UNITS_REQUIRED')
    return unit, 1. if scale is None else float(scale)


def _bin_rows(edges, name):
    return [dict(bin=i, **{name + '_start': a, name + '_end': b})
            for i, (a, b) in enumerate(zip(edges, edges[1:]))]


def _index(value, edges):
    """Half-open intervals; only a finite final endpoint belongs to the last bin."""
    if not math.isfinite(value) or value < edges[0] or value > edges[-1]:
        return None
    return min(bisect.bisect_right(edges, value) - 1, len(edges) - 2)


def _histogram(family, rows, weights, total, *, unit, partial, method,
               reasons=None, unavailable=None):
    assigned = float(math.fsum(map(float, weights)))
    if total is not None and assigned > total + 1e-10 * max(total, np.finfo(float).tiny):
        raise AssertionError('DISTRIBUTION_LENGTH_EXCEEDS_UNIQUE_AXIS')
    assigned = min(assigned, total) if total is not None else 0.
    if total is None:
        status = 'UNAVAILABLE_NO_MEASURABLE_AXIS'
    elif unavailable:
        status = 'UNAVAILABLE_' + unavailable
    elif not assigned:
        status = 'UNAVAILABLE_NO_ASSIGNED_SUPPORT'
    else:
        limited = partial or assigned < total - 1e-12 * max(total, np.finfo(float).tiny)
        status = ('PARTIAL_' if limited else '') + method
    qualified = assigned > 0 and total is not None and unavailable is None
    for row, weight in zip(rows, weights):
        row.update(axis_length=float(weight) if qualified else None,
                   fraction_of_total_axis=float(weight) / total if qualified else None,
                   fraction_of_assigned_axis=float(weight) / assigned if qualified else None)
    return dict(family=family, status=status, unit=unit, method=method,
                total_axis_length=total,
                assigned_axis_length=assigned if total is not None else None,
                unassigned_axis_length=max(0., total - assigned) if total is not None else None,
                covered_arc_fraction=assigned / total if total is not None else None,
                unassigned_reasons=reasons or {}, bins=rows,
                bin_convention='[start,end); finite final endpoint included; null upper bound means +infinity',
                scope='observed unique measured axis only; extraction completeness and accuracy are separate')


def _depth_reference(report, bounds, points, factor, axis, sign, config):
    source = 'SOURCE_MATERIAL_BOUNDS'
    value = bounds if bounds is not None else report.get('model_bounds')
    if value is None and len(points):
        lo, hi = points.min(axis=0), points.max(axis=0)
        source = 'EXTRACTED_AXIS_BOUNDS'
    elif value is not None:
        if not isinstance(value, dict):
            raise ValueError('MODEL_UNIT_MINIMUM_MAXIMUM_BOUNDS_REQUIRED')
        lo = np.asarray(value.get('minimum'), dtype=float)
        hi = np.asarray(value.get('maximum'), dtype=float)
        if lo.shape != (3,) or hi.shape != (3,) or not np.isfinite([lo, hi]).all() or np.any(hi < lo):
            raise ValueError('FINITE_ORDERED_XYZ_BOUNDS_REQUIRED')
        lo, hi = lo * factor, hi * factor
    else:
        return dict(source='UNAVAILABLE', status='DEPTH_REFERENCE_REQUIRED',
                    top_coordinate=None, depth_span=None, edges=[])
    if not np.isfinite([lo, hi]).all():
        raise ValueError('SCALED_BOUNDS_FLOAT_RANGE')
    top = float(lo[axis] if sign == 1 else hi[axis])
    span = float(hi[axis] - lo[axis])
    explicit = config['depth_edges']
    edges = explicit if explicit is not None else (
        np.linspace(0., span, config['depth_bins'] + 1).tolist() if span > 0 else [])
    return dict(source=source, status='AVAILABLE' if edges else 'ZERO_DEPTH_SPAN',
                top_coordinate=top, depth_span=span, edges=edges,
                definition='down_sign * (coordinate - upper bound); not depth below a measured soil surface')


def _depth_weights(pairs, lengths, top, axis, sign, edges):
    """Yield exact in-layer arc fractions for straight edges, without dense resampling."""
    depth = sign * (pairs[:, :, axis] - top)
    low, high = depth.min(axis=1), depth.max(axis=1)
    span = high - low
    denominator = np.where(span > 0, span, 1.)
    for i, (a, b) in enumerate(zip(edges, edges[1:])):
        weight = np.clip((b - low) / denominator, 0., 1.) - np.clip((a - low) / denominator, 0., 1.)
        constant = (low >= a) & ((low <= b) if i == len(edges) - 2 else (low < b))
        yield np.where(span > 0, weight, constant) * lengths


def _radial_weights(pairs, lengths, origin, horizontal, thresholds):
    """Integrate line segments through horizontal annuli, including inward/outward crossings."""
    q = pairs[:, 0, horizontal] - origin[horizontal]
    v = pairs[:, 1, horizontal] - pairs[:, 0, horizontal]
    a = np.sum(v * v, axis=1)
    b = 2. * np.sum(q * v, axis=1)
    q2 = np.sum(q * q, axis=1)
    denominator = np.where(a > 0, 2. * a, 1.)
    previous = np.zeros(len(lengths))
    for radius in thresholds:
        discriminant = b * b - 4. * a * (q2 - radius * radius)
        root = np.sqrt(np.maximum(0., discriminant))
        entry, leave = (-b - root) / denominator, (-b + root) / denominator
        fraction = np.maximum(0., np.minimum(leave, 1.) - np.maximum(entry, 0.))
        fraction = np.where(discriminant > 0, fraction, 0.)
        # Constant-radius edges exactly on a boundary belong to the upper bin.
        fraction = np.where(a > 0, fraction, q2 < radius * radius)
        fraction = np.maximum(previous, np.minimum(1., fraction))
        yield float(np.dot(fraction - previous, lengths))
        previous = fraction
    yield float(np.dot(1. - previous, lengths))


def _bearing_weights(pairs, lengths, origin, horizontal, bearing_edges,
                     *, axis, sign, top, depth_edges):
    """Integrate position bearing after cuts at rays and depth planes.

    A midpoint is used only to label each already cut open interval. Its angular
    sector and depth layer are then constant on the interval; no whole-edge
    midpoint approximation is made. An extra closest-point cut resolves passages
    through the reference axis, including when a segment lies on a boundary ray.
    """
    count = len(bearing_edges) - 1
    depth_count = max(0, len(depth_edges) - 1)
    bearing = np.zeros(count)
    joint = np.zeros((depth_count, count))
    rays = []
    cardinal = {0.: (1., 0.), 90.: (0., 1.), 180.: (-1., 0.), 270.: (0., -1.)}
    for angle in bearing_edges[:-1]:
        radians = math.radians(angle)
        rays.append(cardinal.get(angle, (math.cos(radians), math.sin(radians))))
    rays = np.asarray(rays, dtype=float)
    on_axis = unresolved = outside_depth = 0.
    horizontal_origin = origin[horizontal]
    for pair, length in zip(pairs, lengths):
        if length <= 0:
            continue
        q = pair[0, horizontal] - horizontal_origin
        vector = pair[1, horizontal] - pair[0, horizontal]
        if not np.isfinite(q).all() or not np.isfinite(vector).all():
            raise ValueError('BEARING_COORDINATE_FLOAT_RANGE')
        magnitude = float(max(np.max(np.abs(q)), np.max(np.abs(vector))))
        if magnitude == 0:
            on_axis += float(length)
            continue
        # Scaling the intersection equations avoids coordinate-unit overflow and
        # makes tiny but nonzero offsets from the reference axis remain defined.
        q, vector = q / magnitude, vector / magnitude
        cuts = [0., 1.]
        denominators = vector[0] * rays[:, 1] - vector[1] * rays[:, 0]
        numerators = -(q[0] * rays[:, 1] - q[1] * rays[:, 0])
        usable = denominators != 0
        with np.errstate(over='ignore', invalid='ignore', divide='ignore'):
            times = numerators[usable] / denominators[usable]
        usable_rays = rays[usable]
        interior = np.isfinite(times) & (times > 0) & (times < 1)
        times, usable_rays = times[interior], usable_rays[interior]
        intersections = q + times[:, None] * vector
        forward = np.sum(intersections * usable_rays, axis=1) >= 0
        cuts.extend(float(t) for t in times[forward])
        squared = float(np.dot(vector, vector))
        if squared:
            nearest = -float(np.dot(q, vector)) / squared
            if 0 < nearest < 1:
                cuts.append(nearest)
        if depth_edges:
            start_depth = sign * (float(pair[0, axis]) - top)
            depth_change = sign * float(pair[1, axis] - pair[0, axis])
            if depth_change:
                cuts.extend(t for boundary in depth_edges
                            if 0 < (t := (boundary - start_depth) / depth_change) < 1)
        else:
            start_depth = depth_change = None
        cuts = sorted(set(cuts))
        for start, end in zip(cuts, cuts[1:]):
            weight = float(length) * (end - start)
            if weight <= 0:
                continue
            middle = start + (end - start) / 2
            position = q + middle * vector
            if not np.any(position):
                # A single axis crossing has zero measure. Only floating-point
                # collapse can leave an open subinterval centred exactly there.
                middle = start + (end - start) / 4
                position = q + middle * vector
            if not np.any(position):
                unresolved += weight
                continue
            angle = math.degrees(math.atan2(float(position[1]), float(position[0]))) % 360.
            sector = _index(angle, bearing_edges)
            if sector is None:
                raise AssertionError('BEARING_FULL_CIRCLE_CLASSIFICATION')
            bearing[sector] += weight
            if depth_edges:
                layer = _index(start_depth + middle * depth_change, depth_edges)
                if layer is None:
                    outside_depth += weight
                else:
                    joint[layer, sector] += weight
    return bearing, joint.reshape(-1), on_axis, unresolved, outside_depth


def build_trait_distributions(report, *, options=None, samples=None, bounds=None, partial_axis=None):
    """Build JSON-safe, conservative distributions without changing the source graph.

``bounds`` uses model coordinates, like ``report['model_bounds']``. Dimensional
options, including an explicit radial ``origin``, use the reported output unit.
``samples`` may reuse ``measurement_tables(report)``; partial sample sets retain
their missing arc in the denominator. Overlapping/overweighted sets are rejected.
    """
    from root_graph_measurement_ledger import graph_measurement_ledger
    from root_model_delivery_validation import complete_length_state

    config = distribution_settings(options)
    unit, factor = _unit(report)
    axis = 'xyz'.index(config['vertical_axis'])
    horizontal = [i for i in range(3) if i != axis]
    graph = report.get('candidate_centerline') or dict(points=[], adjacency=[])
    original_points = np.asarray(graph['points'], dtype=float).reshape(-1, 3)
    ledger = graph_measurement_ledger(original_points, graph['adjacency'],
                                     coordinate_units=report.get('coordinate_units', 'model_units')) if len(original_points) else None
    points = original_points * factor
    if not np.isfinite(points).all():
        raise ValueError('SCALED_COORDINATES_FLOAT_RANGE')
    edges = np.asarray(ledger['source_edge_ids'] if ledger else [], dtype=int).reshape(-1, 2)
    pairs = points[edges]
    delta = pairs[:, 1] - pairs[:, 0]
    lengths = np.asarray([math.hypot(*d) for d in delta])
    if not np.isfinite(lengths).all():
        raise ValueError('SCALED_LENGTH_FLOAT_RANGE')
    measured_total = math.fsum(map(float, lengths))
    total = measured_total if measured_total > 0 else None
    if partial_axis is None:
        state = complete_length_state(report)
        partial_axis = (state['complete_length_status'].startswith('BLOCKED') or
                        bool(state['isolated_component_node_ids']) or 'BLOCKED' in report.get('status', ''))
    elif not isinstance(partial_axis, (bool, np.bool_)):
        raise ValueError('BOOLEAN_PARTIAL_AXIS_REQUIRED')
    partial_axis = bool(partial_axis)
    if samples is None:
        from root_model_delivery_exports import measurement_tables
        samples = measurement_tables(report)[1] if len(edges) else []
    samples = list(samples)
    edge_lengths = {tuple(e): float(length) for e, length in zip(edges, lengths)}
    edge_weights = {key: 0. for key in edge_lengths}
    valid = []
    sample_total = 0.
    for sample in samples:
        nodes = sample.get('edge_nodes')
        if not isinstance(nodes, (list, tuple)) or len(nodes) != 2:
            raise ValueError('SAMPLE_SOURCE_EDGE_REQUIRED')
        if any(isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer)) for n in nodes):
            raise ValueError('SAMPLE_INTEGER_EDGE_IDS_REQUIRED')
        key = tuple(sorted(map(int, nodes)))
        if key not in edge_lengths:
            raise ValueError('SAMPLE_EDGE_NOT_IN_MEASUREMENT_GRAPH')
        weight = sample.get('arc_weight')
        if isinstance(weight, (bool, np.bool_)) or not isinstance(weight, (int, float, np.number)) or not math.isfinite(weight) or weight < 0:
            raise ValueError('FINITE_NONNEGATIVE_SAMPLE_ARC_REQUIRED')
        weight = float(weight) * factor
        if not math.isfinite(weight):
            raise ValueError('SCALED_SAMPLE_ARC_FLOAT_RANGE')
        edge_weights[key] += weight
        sample_total += weight
        if edge_weights[key] > edge_lengths[key] + 1e-10 * max(edge_lengths[key], np.finfo(float).tiny):
            raise ValueError('SAMPLE_ARC_EXCEEDS_UNIQUE_SOURCE_EDGE')
        if is_qualified_radius_sample(sample) and weight > 0:
            diameter = 2. * float(sample['radius']) * factor
            if not math.isfinite(diameter):
                raise ValueError('SCALED_DIAMETER_FLOAT_RANGE')
            valid.append((sample, diameter, weight))
    qualified_total = math.fsum(weight for _, _, weight in valid)
    diameter_bounds = [0., *config['diameter_edges'], None]
    diameter_rows = _bin_rows(diameter_bounds, 'diameter')
    diameter_weights = [0.] * len(diameter_rows)
    for _, diameter, weight in valid:
        diameter_weights[bisect.bisect_right(config['diameter_edges'], diameter)] += weight
    missing_radius = max(0., measured_total - qualified_total)
    histograms = {}
    histograms['diameter'] = _histogram('diameter', diameter_rows, diameter_weights, total,
        unit=unit, partial=partial_axis, method='SAMPLE_SUPPORTED_ESTIMATE',
        reasons={'radius_unqualified_or_not_sampled': missing_radius},
        unavailable='NO_QUALIFIED_RADIUS' if total is not None and not qualified_total else None)
    histograms['diameter']['definition'] = 'length classes from 2 * qualified original single-branch area-equivalent radius and sample arc_weight; no biological diameter threshold'

    depth = _depth_reference(report, bounds, points, factor, axis, config['down_sign'], config)
    depth_edges = depth['edges']
    depth_count = max(0, len(depth_edges) - 1)
    joint_rows = [dict(bin=f'd{d}:z{z}', diameter_start=da, diameter_end=db,
                       depth_start=za, depth_end=zb)
                  for d, (da, db) in enumerate(zip(diameter_bounds, diameter_bounds[1:]))
                  for z, (za, zb) in enumerate(zip(depth_edges, depth_edges[1:]))]
    joint_weights = [0.] * len(joint_rows)
    invalid_position = outside_depth = 0.
    if depth_edges:
        for sample, diameter, weight in valid:
            position = sample.get('position')
            try:
                position = np.asarray(position, dtype=float)
                good = position.shape == (3,) and np.isfinite(position).all()
            except (TypeError, ValueError):
                good = False
            if not good:
                invalid_position += weight
                continue
            location = config['down_sign'] * (float(position[axis]) * factor - depth['top_coordinate'])
            z = _index(location, depth_edges)
            if z is None:
                outside_depth += weight
            else:
                d = bisect.bisect_right(config['diameter_edges'], diameter)
                joint_weights[d * depth_count + z] += weight
    histograms['diameter_depth'] = _histogram('diameter_depth', joint_rows, joint_weights, total,
        unit=unit, partial=partial_axis, method='SAMPLE_SUPPORTED_ESTIMATE',
        reasons=dict(radius_unqualified_or_not_sampled=missing_radius,
                     invalid_sample_position=invalid_position, outside_depth_bins=outside_depth,
                     depth_reference_unavailable=qualified_total if not depth_edges else 0.),
        unavailable=depth['status'] if not depth_edges else 'NO_QUALIFIED_RADIUS' if not qualified_total else None)
    histograms['diameter_depth']['definition'] = 'sample arc_weight assigned at its recorded sample position; not exact continuous diameter-by-depth integration; radius gaps are not interpolated'

    positive = lengths > 0
    angles = np.zeros(len(lengths))
    angles[positive] = np.degrees(np.arccos(np.clip(np.abs(delta[positive, axis]) / lengths[positive], 0., 1.)))
    inclination_edges = config['inclination_edges_deg']
    inclination_weights = np.histogram(angles[positive], bins=inclination_edges, weights=lengths[positive])[0]
    histograms['inclination'] = _histogram('inclination', _bin_rows(inclination_edges, 'inclination_deg'),
        inclination_weights, total, unit=unit, partial=partial_axis, method='OBSERVED_AXIS_DESCRIPTOR')
    histograms['inclination']['definition'] = 'undirected segment inclination: 0 degrees parallel to the vertical axis, 90 horizontal; not a parent-child emergence angle'

    horizontal_norm = np.hypot(delta[:, horizontal[0]], delta[:, horizontal[1]])
    direction_defined = positive & (horizontal_norm > 64 * np.finfo(float).eps * lengths)
    azimuth = np.zeros(len(lengths))
    azimuth[direction_defined] = np.degrees(np.arctan2(delta[direction_defined, horizontal[1]],
                                                      delta[direction_defined, horizontal[0]])) % 180.
    azimuth_edges = config['azimuth_edges_deg']
    azimuth_weights = np.histogram(azimuth[direction_defined], bins=azimuth_edges, weights=lengths[direction_defined])[0]
    vertical_length = float(lengths[~direction_defined].sum())
    histograms['azimuth'] = _histogram('azimuth', _bin_rows(azimuth_edges, 'azimuth_deg'),
        azimuth_weights, total, unit=unit, partial=partial_axis, method='OBSERVED_AXIS_DESCRIPTOR',
        reasons=dict(vertical_azimuth_undefined=vertical_length))
    histograms['azimuth']['definition'] = 'undirected axial orientation modulo 180 degrees in the horizontal coordinate plane; vertical edges have undefined azimuth; not root growth direction or spatial bearing from an origin'

    azimuth_count = len(azimuth_edges) - 1
    direction_rows = [dict(bin=f'z{z}:a{a}', depth_start=za, depth_end=zb,
                           azimuth_deg_start=aa, azimuth_deg_end=ab)
                      for z, (za, zb) in enumerate(zip(depth_edges, depth_edges[1:]))
                      for a, (aa, ab) in enumerate(zip(azimuth_edges, azimuth_edges[1:]))]
    direction_weights = np.zeros(depth_count * azimuth_count)
    inside_defined = 0.
    if depth_edges and len(pairs):
        for z, weights in enumerate(_depth_weights(pairs, lengths, depth['top_coordinate'], axis,
                                                   config['down_sign'], depth_edges)):
            classified = np.histogram(azimuth[direction_defined], bins=azimuth_edges,
                                      weights=weights[direction_defined])[0]
            direction_weights[z * azimuth_count:(z + 1) * azimuth_count] = classified
            inside_defined += float(classified.sum())
    defined_length = float(lengths[direction_defined].sum())
    histograms['depth_azimuth'] = _histogram('depth_azimuth', direction_rows, direction_weights, total,
        unit=unit, partial=partial_axis, method='OBSERVED_AXIS_DESCRIPTOR',
        reasons=dict(vertical_azimuth_undefined=vertical_length,
                     outside_depth_bins=max(0., defined_length - inside_defined) if depth_edges else 0.,
                     depth_reference_unavailable=defined_length if not depth_edges else 0.),
        unavailable=depth['status'] if not depth_edges else None)
    histograms['depth_azimuth']['definition'] = 'exact straight-segment arc split at depth boundaries, classified by undirected horizontal axial orientation'

    bearing_edges = config['bearing_edges_deg']
    bearing_count = len(bearing_edges) - 1
    bearing_rows = _bin_rows(bearing_edges, 'bearing_deg')
    bearing_depth_rows = [dict(bin=f'z{z}:b{b}', depth_start=za, depth_end=zb,
                              bearing_deg_start=ba, bearing_deg_end=bb)
                         for z, (za, zb) in enumerate(zip(depth_edges, depth_edges[1:]))
                         for b, (ba, bb) in enumerate(zip(bearing_edges, bearing_edges[1:]))]
    bearing_weights = np.zeros(bearing_count)
    bearing_depth_weights = np.zeros(depth_count * bearing_count)
    on_reference_axis = bearing_unresolved = bearing_outside_depth = 0.
    if config['origin'] is not None:
        bearing_weights, bearing_depth_weights, on_reference_axis, bearing_unresolved, bearing_outside_depth = _bearing_weights(
            pairs, lengths, np.asarray(config['origin']), horizontal, bearing_edges,
            axis=axis, sign=config['down_sign'], top=depth['top_coordinate'], depth_edges=depth_edges)
    bearing_reasons = dict(origin_unavailable=measured_total if config['origin'] is None else 0.,
        reference_axis_bearing_undefined=on_reference_axis,
        numerical_bearing_unresolved=bearing_unresolved)
    histograms['bearing'] = _histogram('bearing', bearing_rows, bearing_weights, total,
        unit=unit, partial=partial_axis, method='OBSERVED_AXIS_DESCRIPTOR', reasons=bearing_reasons,
        unavailable='EXPLICIT_BEARING_ORIGIN_REQUIRED' if config['origin'] is None else None)
    histograms['bearing']['definition'] = 'position bearing from explicit origin in the horizontal plane, 0–360 degrees; exact straight-segment cuts at angular boundary rays; not tangent direction'
    histograms['bearing']['origin'] = config['origin']
    histograms['depth_bearing'] = _histogram('depth_bearing', bearing_depth_rows, bearing_depth_weights, total,
        unit=unit, partial=partial_axis, method='OBSERVED_AXIS_DESCRIPTOR',
        reasons=dict(bearing_reasons, outside_depth_bins=bearing_outside_depth,
            depth_reference_unavailable=float(np.sum(bearing_weights)) if not depth_edges else 0.),
        unavailable='EXPLICIT_BEARING_ORIGIN_REQUIRED' if config['origin'] is None else
                    depth['status'] if not depth_edges else None)
    histograms['depth_bearing']['definition'] = 'position bearing by depth: exact straight-segment arc split at angular boundary rays and depth planes; reference-axis segments have undefined bearing'
    histograms['depth_bearing']['origin'] = config['origin']

    radial_rows = _bin_rows([0., *config['radial_edges'], None], 'radial_distance')
    radial_weights = [0.] * len(radial_rows)
    if config['origin'] is not None:
        radial_weights = list(_radial_weights(pairs, lengths, np.asarray(config['origin']),
                                              horizontal, config['radial_edges']))
    histograms['radial'] = _histogram('radial', radial_rows, radial_weights, total,
        unit=unit, partial=partial_axis, method='OBSERVED_AXIS_DESCRIPTOR',
        reasons=dict(origin_unavailable=measured_total if config['origin'] is None else 0.),
        unavailable='EXPLICIT_RADIAL_ORIGIN_REQUIRED' if config['origin'] is None else None)
    histograms['radial']['definition'] = 'unique arc in horizontal annuli around the line through explicit origin parallel to the vertical axis; full straight-segment intersections, not vertex counts'
    histograms['radial']['origin'] = config['origin']

    return dict(schema_version=METHOD, unit=unit, settings=config, partial_axis=partial_axis,
                total_axis_length=total, source_geometry_unchanged=True, reference_precision_separate=True,
                depth_reference=depth, horizontal_axes=['xyz'[i] for i in horizontal],
                azimuth_reference=f"+{'xyz'[horizontal[0]]} toward +{'xyz'[horizontal[1]]}; modulo 180 degrees",
                bearing_reference=f"position from explicit origin: +{'xyz'[horizontal[0]]} toward +{'xyz'[horizontal[1]]}; 0 to 360 degrees",
                zero_length_edge_count=int(np.count_nonzero(~positive)),
                sample_coverage=dict(sample_count=len(samples), qualified_radius_sample_count=len(valid),
                    sampled_axis_length=sample_total if total is not None else None,
                    unsampled_axis_length=max(0., measured_total - sample_total) if total is not None else None,
                    qualified_radius_axis_length=qualified_total if total is not None else None,
                    qualified_radius_arc_fraction=qualified_total / total if total is not None else None,
                    unqualified_or_unsampled_axis_length=missing_radius if total is not None else None,
                    definition='sample-supported observed arc; not continuous radius coverage or reference completeness'),
                distributions=histograms,
                qualifications=['Diameter thresholds are numeric classes, not absorptive/coarse root identities.',
                                'Length normalization uses the observed unique axis, including unassigned arc.',
                                'A complete histogram does not certify a complete or biologically correct root system.'])


def flatten_distribution_rows(result):
    """Return plain records for a unified CSV; retain unavailable families and metadata."""
    dimensions = ('diameter_start', 'diameter_end', 'depth_start', 'depth_end',
                  'inclination_deg_start', 'inclination_deg_end', 'azimuth_deg_start',
                  'azimuth_deg_end', 'bearing_deg_start', 'bearing_deg_end',
                  'radial_distance_start', 'radial_distance_end')
    rows = []
    for family, histogram in result['distributions'].items():
        for value in histogram['bins'] or [dict(bin=None, axis_length=None,
                                                fraction_of_total_axis=None, fraction_of_assigned_axis=None)]:
            row = dict(family=family, bin=value['bin'], length=value['axis_length'],
                       fraction=value['fraction_of_total_axis'],
                       fraction_of_assigned_axis=value['fraction_of_assigned_axis'],
                       status=histogram['status'], unit=result['unit'],
                       covered_arc_fraction=histogram['covered_arc_fraction'],
                       total_axis_length=histogram['total_axis_length'],
                       unassigned_axis_length=histogram['unassigned_axis_length'],
                       method=result['schema_version'],
                       vertical_axis=result['settings']['vertical_axis'], down_sign=result['settings']['down_sign'],
                       depth_origin_source=result['depth_reference']['source'],
                       bin_convention=histogram['bin_convention'])
            row.update({name: value.get(name) for name in dimensions})
            rows.append(row)
    return rows
