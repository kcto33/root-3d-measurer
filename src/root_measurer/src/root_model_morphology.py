"""Finite-resolution root geometry descriptors, separate from accuracy acceptance.

All graph edges are counted once. Solid traits use the original occupied cubes;
no radius tube, biological hierarchy, soil volume or dry mass is inferred.
"""
import csv
import io
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull, QhullError

METHOD = 'root-morphology-v001'
DEFAULTS = Path(__file__).resolve().parents[1] / 'configs/morphology_defaults_v001.json'


def settings(options=None):
    value = json.loads(DEFAULTS.read_text(encoding='utf-8'))
    options = options or {}
    if set(options) - {'vertical_axis', 'down_sign'}:
        raise ValueError('UNKNOWN_MORPHOLOGY_OPTION')
    value.update(options)
    if value['vertical_axis'] not in ('x', 'y', 'z'):
        raise ValueError('VERTICAL_AXIS_X_Y_Z_REQUIRED')
    if isinstance(value['down_sign'], bool) or value['down_sign'] not in (-1, 1):
        raise ValueError('DOWN_SIGN_MINUS_OR_PLUS_ONE_REQUIRED')
    return value


def box_count(points, edges, epsilon, phase, visit_budget=2000000):
    """Count half-open boxes met by complete continuous segments (not vertices).

    Interior interval samples plus exact grid crossing points implement floor's
    boundary convention; subdivision cannot remove occupied boxes. No dense grid.
    """
    p = np.asarray(points, dtype=float)
    p = (p - p.min(axis=0)) / epsilon - np.asarray(phase)
    e = np.asarray(edges, dtype=int).reshape(-1, 2)
    a, b = p[e[:, 0]], p[e[:, 1]]
    fa, fb = np.floor(a), np.floor(b)
    visits = int(np.abs(fb-fa).sum()) + 2*len(e) + len(p)
    if visits > visit_budget:
        raise ValueError('BOX_VISIT_BUDGET_EXCEEDED')
    chunks = [np.floor(p).astype(np.int64)]
    short = np.all(np.abs(fb-fa) <= 1, axis=1)
    aa, bb = a[short], b[short]
    if len(aa):
        delta = bb-aa
        planes = np.where(delta > 0, np.floor(aa)+1, np.floor(aa))
        with np.errstate(divide='ignore', invalid='ignore'):
            t = (planes-aa)/delta
        t = np.where((t > 0) & (t < 1), t, 1.)
        times = np.sort(np.column_stack([np.zeros(len(aa)), t, np.ones(len(aa))]), axis=1)
        middle = (times[:, :-1]+times[:, 1:])/2
        chunks.append(np.floor(aa[:, None, :]+middle[:, :, None]*delta[:, None, :]).astype(np.int64).reshape(-1, 3))
        # Snap the crossing axis to its integer grid plane, avoiding roundoff.
        for axis in range(3):
            q = aa+t[:, axis, None]*delta
            mask = t[:, axis] < 1
            q[mask, axis] = planes[mask, axis]
            chunks.append(np.floor(q).astype(np.int64))
    for start, end in zip(a[~short], b[~short]):
        delta = end-start
        events = []
        for axis in range(3):
            if delta[axis] == 0:
                continue
            for plane in range(math.floor(min(start[axis], end[axis]))+1,
                               math.ceil(max(start[axis], end[axis]))):
                t = (plane-start[axis])/delta[axis]
                if 0 < t < 1:
                    events.append((t, axis, plane))
        times = np.array([0., *sorted(set(t for t, _, _ in events)), 1.])
        chunks.append(np.floor(start+(times[:-1]+times[1:])[:, None]/2*delta).astype(np.int64))
        # Merge simultaneous crossings before rounding.
        for t in sorted(set(t for t, _, _ in events)):
            q = start+t*delta
            for et, axis, plane in events:
                if et == t:
                    q[axis] = plane
            chunks.append(np.floor(q).astype(np.int64).reshape(1, 3))
    return len(np.unique(np.concatenate(chunks), axis=0))


def fractal_dimension(points, edges, resolution, config):
    result = dict(method='continuous-centerline 3D box count; eight half-box grid phases',
                  status='INSUFFICIENT_SCALE_RANGE', dimension=None, r_squared=None,
                  scales=[], definition='finite-scale descriptor, not an asymptotic fractal dimension',
                  boundary_convention='half-open floor boxes; translation anchor = graph coordinate minimum',
                  scale_selection='all dyadic scales h*2^k from h through maximum axis extent/4; no fitted window selection')
    if resolution is None or not math.isfinite(resolution) or resolution <= 0:
        result['status'] = 'INPUT_RESOLUTION_REQUIRED'
        return result
    if len(points) < 2 or not len(edges):
        result['status'] = 'AXIS_REQUIRED'
        return result
    extent = float(np.ptp(points, axis=0).max())
    epsilon = float(resolution)
    phases = list(itertools.product((0., .5), repeat=3))
    try:
        while epsilon <= extent/4*(1+1e-12):
            counts = [box_count(points, edges, epsilon, phase, config['box_visit_budget']) for phase in phases]
            result['scales'].append(dict(box_size_model_units=epsilon, phase_counts=counts, mean_count=float(np.mean(counts))))
            epsilon *= 2
    except ValueError as exc:
        result['status'] = str(exc)
        return result
    if len(result['scales']) < config['minimum_scales']:
        return result
    sizes = np.array([s['box_size_model_units'] for s in result['scales']])
    x = np.log(extent/sizes)
    def fit(y):
        slope, intercept = np.polyfit(x, y, 1)
        residual = float(np.sum((y-(slope*x+intercept))**2))
        total = float(np.sum((y-y.mean())**2))
        return float(slope), None if total == 0 else 1-residual/total
    slope, r2 = fit(np.log([s['mean_count'] for s in result['scales']]))
    phase_slopes = [fit(np.log([s['phase_counts'][i] for s in result['scales']]))[0] for i in range(8)]
    result.update(dimension=slope, r_squared=r2, phase_dimensions=phase_slopes,
                  scale_ratio=float(sizes[-1]/sizes[0]), scale_count=len(sizes),
                  status='FINITE_SCALE_ESTIMATE' if r2 is not None and r2 >= .98 else 'FINITE_SCALE_ESTIMATE_LOW_LINEARITY',
                  warning='R² describes this scale range only; grid phase, rotation, resolution and finite size affect the estimate.')
    return result


def solid_hull(cells):
    """Exact convex envelope of equal cubes, with rank-aware center hull reduction."""
    base = cells.min(axis=0)
    centers = np.asarray(cells-base, dtype=float)
    shifted = centers-centers[0]
    _, singular, basis = np.linalg.svd(shifted, full_matrices=False)
    rank = int(np.sum(singular > max(shifted.shape)*np.finfo(float).eps*(singular[0] if len(singular) else 0)))
    if rank >= 2:
        projected = shifted @ basis[:rank].T
        extremes = centers[ConvexHull(projected).vertices]
    elif rank == 1:
        projected = shifted @ basis[0]
        extremes = centers[[np.argmin(projected), np.argmax(projected)]]
    else:
        extremes = centers[:1]
    corners = (extremes[:, None, :]+np.array(list(itertools.product((0., 1.), repeat=3)))).reshape(-1, 3)
    return corners, ConvexHull(corners)


def depth_profile(pairs, lengths, top, axis, sign, depth, cells=None, pitch=None, bins=10):
    """Integrate linear segments and original cubes across equal depth intervals."""
    d = sign*(pairs[:, :, axis]-top)
    low, high = d.min(axis=1), d.max(axis=1)
    span = high-low
    def cdf(x):
        weight = np.where(span > 0, np.clip((x-low)/np.where(span > 0, span, 1), 0, 1), x >= low)
        return float(np.dot(weight, lengths))
    total = float(lengths.sum())
    def quantile(q):
        if total == 0:
            return None
        lo, hi = float(low.min()), float(high.max())
        for _ in range(70):
            middle = (lo+hi)/2
            if cdf(middle) >= q*total:
                hi = middle
            else:
                lo = middle
        return hi
    rows = []
    before_zero = cdf(float(np.nextafter(0., -np.inf)))
    previous = before_zero
    layers = counts = None
    if cells is not None:
        layers, counts = np.unique(cells[:, axis], return_counts=True)
        ends = sign*(np.column_stack([layers, layers+1])*pitch-top)
        vl, vh = ends.min(axis=1), ends.max(axis=1)
    for i in range(bins):
        a, b = depth*i/bins, depth*(i+1)/bins
        cumulative = cdf(b) if i == bins-1 else cdf(float(np.nextafter(b, -np.inf)))
        volume = None if layers is None else float(np.dot(np.maximum(0., np.minimum(vh, b)-np.maximum(vl, a))/pitch, counts)*pitch**3)
        rows.append(dict(depth_start=a, depth_end=b, axis_length=cumulative-previous,
                         axis_length_fraction=None if not total else (cumulative-previous)/total, solid_volume=volume))
        previous = cumulative
    outside = max(0., total-(cdf(depth)-before_zero))
    return rows, quantile(.5), quantile(.95), outside


def build_morphology(report, *, solid=None, mesh=None, options=None, geometry=None):
    from root_graph_measurement_ledger import graph_measurement_ledger
    from root_model_delivery_validation import complete_length_state
    config = settings(options)
    axis = 'xyz'.index(config['vertical_axis'])
    horizontal = [i for i in range(3) if i != axis]
    sign = config['down_sign']
    scale = report.get('scale_mm_per_model_unit')
    factor = 1. if scale is None else float(scale)
    if not math.isfinite(factor) or factor <= 0:
        raise ValueError('POSITIVE_FINITE_SCALE_REQUIRED')
    unit = 'mm' if scale is not None else report.get('coordinate_units', 'model_unit')
    graph = report.get('candidate_centerline') or dict(points=[], adjacency=[])
    p = np.asarray(graph['points'], dtype=float).reshape(-1, 3)
    ledger = graph_measurement_ledger(p, graph['adjacency'], coordinate_units=unit) if len(p) else None
    edges = np.asarray(ledger['source_edge_ids'] if ledger else [], dtype=int).reshape(-1, 2)
    pairs = p[edges]
    lengths = np.linalg.norm(pairs[:, 1]-pairs[:, 0], axis=1)
    total = float(lengths.sum())
    state = complete_length_state(report)
    partial = state['complete_length_status'].startswith('BLOCKED') or bool(state['isolated_component_node_ids']) or 'BLOCKED' in report.get('status', '')
    graph_status = 'PARTIAL_AXIS_DESCRIPTOR' if partial else 'MEASURED_AXIS_DESCRIPTOR_ACCURACY_SEPARATE'
    rows = []
    def add(key, label, value, power, definition, status='ORIGINAL_INPUT_GEOMETRY', coverage=None, custom_unit=None):
        rows.append(dict(parameter=key, label=label, value=None if value is None else float(value)*factor**power,
                         unit=custom_unit or ('dimensionless' if power == 0 else unit if power == 1 else unit+'^'+str(power)),
                         status=status if value is not None else 'UNAVAILABLE_'+status,
                         covered_arc_fraction=coverage, definition=definition))
    bounds = report.get('model_bounds')
    if solid is not None:
        lo, hi = solid.cells.min(axis=0)*solid.pitch, (solid.cells.max(axis=0)+1)*solid.pitch
    elif mesh is not None:
        lo, hi = mesh.vertices.min(axis=0), mesh.vertices.max(axis=0)
    elif bounds:
        lo, hi = np.array(bounds['minimum']), np.array(bounds['maximum'])
    else:
        lo = hi = None
    depth = None if lo is None else float(hi[axis]-lo[axis])
    top = None if lo is None else float(lo[axis] if sign == 1 else hi[axis])
    add('root_depth', '根系深度（实体垂直跨度）', depth, 1, 'original material vertical max-min; depth zero at upper material boundary, not soil surface')
    hull_volume = hull_area = spread = footprint = None
    hull_status = 'ORIGINAL_INPUT_GEOMETRY'
    volume = len(solid.cells)*solid.pitch**3 if solid is not None else None
    surface = (geometry or {}).get('surface_area')
    if solid is not None and surface is None:
        from root_reconstruction_measurement import boundary_faces
        surface = len(boundary_faces(solid)[0])*solid.pitch**2
    try:
        if solid is not None:
            corners, hull = solid_hull(solid.cells)
            corners *= solid.pitch
            hull_volume, hull_area = float(hull.volume*solid.pitch**3), float(hull.area*solid.pitch**2)
            footprint = len(np.unique(solid.cells[:, horizontal], axis=0))*solid.pitch**2
        elif mesh is not None:
            corners = np.asarray(mesh.vertices)-mesh.vertices[0]
            hull = ConvexHull(corners)
            hull_volume, hull_area = float(hull.volume), float(hull.area)
        else:
            corners = None
        if corners is not None:
            hp = np.unique(corners[:, horizontal], axis=0)
            projected_hull = ConvexHull(hp)
            hp = hp[projected_hull.vertices]
            projected_area = float(projected_hull.volume)
            # Bounded-memory exact diameter over projected hull vertices.
            spread = math.sqrt(max(float(np.sum((hp[i]-hp)**2, axis=1).max()) for i in range(len(hp))))
        else:
            projected_area = None
    except (QhullError, ValueError) as exc:
        projected_area = None
        hull_status = 'HULL_UNQUALIFIED_'+type(exc).__name__
    add('horizontal_spread', '最大水平展幅', spread, 1, 'maximum pairwise distance of original material horizontal projection hull', hull_status)
    add('width_depth_ratio', '展幅／深度', None if not depth or spread is None else spread/depth, 0, 'maximum horizontal spread divided by vertical material span', hull_status)
    add('convex_hull_volume', '根系凸包体积', hull_volume, 3, 'convex envelope of original material, including spaces between roots', hull_status)
    add('volume_solidity', '体积充实度', None if not hull_volume or volume is None else volume/hull_volume, 0, 'original material volume / convex hull volume; no soil core volume inferred', hull_status)
    add('projected_area', '水平材料投影面积', footprint, 2, 'union of original occupied cube columns; mesh silhouette not inferred')
    add('projected_hull_area', '水平凸包面积', projected_area, 2, 'convex hull area of horizontal original material projection', hull_status)
    add('projected_solidity', '水平投影充实度', None if not projected_area or footprint is None else footprint/projected_area, 0, 'projected material union area / projected convex hull area', hull_status)
    add('hull_length_density', '凸包内根长密度', None if partial or not hull_volume else total/hull_volume, -2, 'unique complete measured axis length / material convex hull volume; not soil-core root length density', graph_status)
    add('branching_density', '几何分叉密度', None if not total else len(ledger['junction_node_ids'])/total, -1, 'graph nodes of degree >=3 / unique measured axis length; adjacent junction nodes are not merged', graph_status)
    tortuosities = []
    if ledger:
        for branch in ledger['branches']:
            ids = branch['source_node_ids']
            chord = float(np.linalg.norm(p[ids[-1]]-p[ids[0]]))
            valid = not branch['closed'] and chord > 1e-12*max(branch['source_length'], 1e-300)
            tortuosities.append(dict(branch_id=branch['branch_id'], closed=branch['closed'], length=branch['source_length']*factor,
                                     chord=chord*factor, tortuosity=branch['source_length']/chord if valid else None))
    supported = sum(t['length'] for t in tortuosities if t['tortuosity'] is not None)
    add('tortuosity_mean', '路径曲折度（弧长加权）', None if not supported else sum(t['length']*t['tortuosity'] for t in tortuosities if t['tortuosity'] is not None)/supported,
        0, 'maximal degree-2 path arc / endpoint chord, arc weighted; closed/zero-chord paths excluded', graph_status,
        supported/(total*factor) if total else None)
    angles = np.degrees(np.arccos(np.clip(np.abs(pairs[:, 1, axis]-pairs[:, 0, axis])/np.where(lengths > 0, lengths, 1), 0, 1)))
    add('inclination_mean', '根段与竖直夹角（弧长加权）', None if not total else float(np.dot(angles, lengths)/total),
        0, 'undirected edge inclination: 0 degrees vertical, 90 horizontal; not biological lateral-root emergence angle', graph_status, custom_unit='degree')
    add('surface_volume_ratio', '比表面积（面积／体积）', None if not volume or surface is None else surface/volume, -1, 'original material boundary area / volume, including cavity boundary')
    add('sphericity', '球形度', None if not volume or not surface else math.pi**(1/3)*(6*volume)**(2/3)/surface, 0, 'equal-volume sphere area / original boundary area; cube staircasing affects this descriptor')
    profile = []; d50 = d95 = None; outside = 0.
    if top is not None and depth:
        profile, d50, d95, outside = depth_profile(pairs, lengths, top, axis, sign, depth,
            solid.cells if solid is not None else None, solid.pitch if solid is not None else None, config['depth_bins'])
    depth_status = graph_status if not outside else 'AXIS_OUTSIDE_MATERIAL_DEPTH_BOUNDS'
    for key, label, value in [('axis_depth_d50', '根长深度 D50', d50), ('axis_depth_d95', '根长深度 D95', d95)]:
        add(key, label, value, 1, 'depth below upper original-material boundary containing 50%/95% of unique measured axis arc; exact linear segment CDF', depth_status)
    for row in profile:
        for key in ('depth_start', 'depth_end', 'axis_length'):
            row[key] *= factor
        if row['solid_volume'] is not None:
            row['solid_volume'] *= factor**3
    resolution = solid.pitch if solid is not None else report.get('stages', {}).get('declared_target_pitch')
    fd = fractal_dimension(p, edges, resolution, config)
    add('fractal_dimension_3d', '三维中心线盒计数维数（有限尺度）', fd['dimension'], 0,
        'slope log(mean occupied boxes over 8 phases) vs log(max axis extent / box size); complete continuous segments; all fixed eligible dyadic scales',
        fd['status']+('_PARTIAL_AXIS' if partial else ''))
    for row in fd['scales']:
        row['box_size'] = row['box_size_model_units']*factor
    return dict(schema_version=METHOD, settings=config, unit=unit, metrics=rows, fractal=fd,
                depth_profile=profile, path_tortuosity=tortuosities, partial_axis=partial,
                outside_depth_axis_arc=outside*factor, geometry_available=solid is not None or mesh is not None,
                source_geometry_unchanged=True, reference_precision_separate=True,
                depth_origin='upper original-material boundary; soil surface unknown',
                qualifications=['No parent/root-order/dry-mass inference.', 'Fractal dimension is a finite-scale centerline descriptor.',
                                'Axis traits inherit extraction limitations; original-solid traits do not certify axis completeness.'])


def morphology_csv(value, profile=False):
    records = value['depth_profile'] if profile else value['metrics']
    fields = ['depth_start', 'depth_end', 'axis_length', 'axis_length_fraction', 'solid_volume'] if profile else ['parameter', 'label', 'value', 'unit', 'status', 'covered_arc_fraction', 'definition']
    if profile:
        records = [dict(r, unit=value['unit'], solid_volume_unit=value['unit']+'^3') for r in records]
        fields += ['unit', 'solid_volume_unit']
    metadata = dict(method=value['schema_version'], measurement_sha256=value.get('measurement',{}).get('sha256'),
                    vertical_axis=value['settings']['vertical_axis'], down_sign=value['settings']['down_sign'])
    records = [dict(r, **metadata) for r in records]
    fields += list(metadata)
    out = io.StringIO(newline='')
    writer = csv.DictWriter(out, fieldnames=fields, lineterminator='\r\n')
    writer.writeheader(); writer.writerows(records)
    return out.getvalue().encode('utf-8-sig')
