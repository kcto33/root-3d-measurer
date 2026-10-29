"""Root-axis traits with explicit direction, annotation and coverage semantics.

Axes are supplied by ``root_axis_architecture``; this module never infers a
biological identity from a graph edge or supplies missing radius measurements.
All lengths, windows and offsets use the declared output coordinate unit.
"""
from collections import defaultdict
import math

import numpy as np


METHOD = 'root-axis-traits-v001'
DEFAULTS = dict(tangent_window=1.0, diameter_offset=1.0, diameter_window=1.0,
                minimum_window_fraction=0.5, minimum_radius_coverage=0.5,
                vertical_axis='z', down_sign=-1)


def axis_trait_settings(options=None):
    """Lengths are absolute output-unit distances, never a number of vertices."""
    options = options or {}
    if not isinstance(options,dict) or set(options) - set(DEFAULTS):
        raise ValueError('UNKNOWN_AXIS_TRAIT_OPTION')
    result = dict(DEFAULTS, **options)
    for key in ('tangent_window', 'diameter_window', 'diameter_offset'):
        value = result[key]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValueError('FINITE_AXIS_TRAIT_WINDOW_REQUIRED')
        if value < 0 or (value == 0 and key != 'diameter_offset'):
            raise ValueError('POSITIVE_AXIS_TRAIT_WINDOW_REQUIRED')
    for key in ('minimum_window_fraction', 'minimum_radius_coverage'):
        value = result[key]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 < value <= 1:
            raise ValueError('AXIS_TRAIT_FRACTION_IN_ZERO_ONE_REQUIRED')
    if result['vertical_axis'] not in ('x', 'y', 'z'):
        raise ValueError('VERTICAL_AXIS_X_Y_Z_REQUIRED')
    if isinstance(result['down_sign'], bool) or result['down_sign'] not in (-1, 1):
        raise ValueError('DOWN_SIGN_MINUS_OR_PLUS_ONE_REQUIRED')
    return result


def _result(value=None, status='AVAILABLE', **details):
    return dict(value=value, status=status, **details)


def _distribution(values, weights=None):
    """Weighted empirical inverse-CDF quantiles; population SD, not sample SD."""
    values = np.asarray(values, dtype=float)
    weights = np.ones(len(values)) if weights is None else np.asarray(weights, dtype=float)
    use = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    values, weights = values[use], weights[use]
    empty = dict(count=0, weight=0.0, mean=None, minimum=None, maximum=None,
                 median=None, p05=None, p25=None, p75=None, p95=None,
                 standard_deviation=None, cv=None)
    if not len(values):
        return empty
    order = np.argsort(values, kind='stable')
    values, weights = values[order], weights[order]
    total = float(weights.sum())
    cumulative = np.cumsum(weights)
    mean = float(np.dot(values, weights) / total)
    deviation = math.sqrt(max(0.0, float(np.dot((values - mean) ** 2, weights) / total)))
    def quantile(q):
        return float(values[min(int(np.searchsorted(cumulative, q * total, side='left')), len(values) - 1)])
    return dict(count=len(values), weight=total, mean=mean,
                minimum=float(values[0]), maximum=float(values[-1]),
                median=quantile(.5), p05=quantile(.05), p25=quantile(.25),
                p75=quantile(.75), p95=quantile(.95),
                standard_deviation=deviation, cv=deviation / mean if mean > 0 else None)


def _subpolyline(points, arc, start, end):
    # Retain zero-length ID edges in the source graph, but interpolation must
    # operate on strictly increasing physical arc stations.
    keep=np.r_[True,np.diff(arc)>0]
    points,arc=points[keep],arc[keep]
    start = min(float(arc[-1]), max(0.0, float(start)))
    end = min(float(arc[-1]), max(0.0, float(end)))
    if end <= start:
        return np.empty((0, 3)), [start, start]
    def at(s):
        i = min(int(np.searchsorted(arc, s, side='right')) - 1, len(arc) - 2)
        i = max(0, i)
        return points[i] + (points[i + 1] - points[i]) * ((s - arc[i]) / (arc[i + 1] - arc[i]))
    first = int(np.searchsorted(arc,start,side='right'))
    last = int(np.searchsorted(arc,end,side='left'))
    middle = points[first:last]
    return np.vstack((at(start), middle, at(end))), [start, end]


def _continuous_pca(points):
    """Exact continuous line moments prevent dependence on vertex subdivision."""
    if len(points) < 2:
        return _result(status='INSUFFICIENT_AXIS_LENGTH', vector=None)
    p = np.asarray(points, dtype=float) - points[0]
    a, b = p[:-1], p[1:]
    lengths = np.linalg.norm(b - a, axis=1)
    total = float(lengths.sum())
    if total <= 0:
        return _result(status='INSUFFICIENT_AXIS_LENGTH', vector=None)
    mean = np.sum(lengths[:, None] * (a + b) / 2, axis=0) / total
    a, b = a - mean, b - mean
    moments = (np.einsum('ni,nj->nij', a, a) + np.einsum('ni,nj->nij', b, b)) / 3
    moments += (np.einsum('ni,nj->nij', a, b) + np.einsum('ni,nj->nij', b, a)) / 6
    covariance = np.sum(lengths[:, None, None] * moments, axis=0) / total
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    eigenvalues = np.maximum(eigenvalues, 0)
    if eigenvalues[-1] <= 0 or eigenvalues[-1] - eigenvalues[-2] <= 1e-8 * eigenvalues[-1]:
        return _result(status='AMBIGUOUS_PCA_DIRECTION', vector=None,
                       eigenvalues=eigenvalues.tolist())
    vector = eigenvectors[:, -1]
    chord = p[-1] - p[0]
    projection = float(np.dot(vector, chord))
    if np.linalg.norm(chord) <= 1e-12 * total or abs(projection) <= 1e-10 * np.linalg.norm(chord):
        return _result(status='PCA_SIGN_UNRESOLVED', vector=None,
                       eigenvalues=eigenvalues.tolist())
    if projection < 0:
        vector = -vector
    return _result(vector.tolist(), vector=vector.tolist(), eigenvalues=eigenvalues.tolist(),
                   explained_variance_fraction=float(eigenvalues[-1] / eigenvalues.sum()),
                   method='continuous arc-weighted PCA; sign aligned with base-to-apex chord')


def _direction_window(axis, start, end, requested, config):
    points, interval = _subpolyline(axis['_points'], axis['_arc'], start, end)
    effective = max(0.0, interval[1] - interval[0])
    metadata = dict(requested_window=requested, effective_window=effective,
                    arc_interval=interval, window_shortened=effective < requested * (1 - 1e-12))
    if effective + 1e-12*requested < requested * config['minimum_window_fraction'] or effective <= 0:
        return _result(status='INSUFFICIENT_TANGENT_WINDOW', vector=None, **metadata)
    result = _continuous_pca(points)
    result.update(metadata)
    return result


def _angle(a, b):
    return math.degrees(math.acos(float(np.clip(np.dot(a, b), -1.0, 1.0))))


def _profile(axis, edge_samples, factor):
    """Map original edge-supported samples to oriented axis arc without filling gaps."""
    from root_trait_distributions import is_qualified_radius_sample
    records = []
    ids, arc = axis['node_ids'], axis['_arc']
    for i, (a, b) in enumerate(zip(ids, ids[1:])):
        lo, hi = float(arc[i]), float(arc[i + 1])
        middle = (lo + hi) / 2
        intervals = []
        for sample in edge_samples.get(tuple(sorted((a, b))), []):
            node = sample.get('node_id')
            sampling = sample.get('sampling', 'edge_midpoint')
            if sampling == 'node_halfedge':
                if node not in (a, b):
                    raise ValueError('HALFEDGE_SAMPLE_NODE_NOT_ON_EDGE')
                left, right, location = (lo, middle, lo) if node == a else (middle, hi, hi)
            elif sampling == 'edge_midpoint' and node is None:
                left, right, location = lo, hi, middle
            else:
                raise ValueError('UNSUPPORTED_AXIS_RADIUS_SAMPLING')
            supplied_weight = sample.get('arc_weight')
            if supplied_weight is not None and (not math.isfinite(float(supplied_weight)) or
                    not math.isclose(float(supplied_weight) * factor, right - left, rel_tol=1e-8, abs_tol=1e-12)):
                raise ValueError('AXIS_RADIUS_SAMPLE_ARC_WEIGHT_MISMATCH')
            qualified = is_qualified_radius_sample(sample)
            radius = float(sample['radius']) * factor if qualified else None
            intervals.append(dict(s=location, support_start=left, support_end=right,
                arc_weight=right-left, radius=radius, diameter=None if radius is None else 2*radius,
                available=qualified, status=sample.get('status', 'NOT_MEASURED'),
                edge_nodes=[a, b], node_id=node, source_path_id=sample.get('path_id'),
                sampling=sampling))
        intervals.sort(key=lambda row: (row['support_start'], row['support_end']))
        cursor = lo
        for row in intervals:
            if row['support_start'] < cursor - max(1e-12, (hi-lo)*1e-10):
                raise ValueError('OVERLAPPING_AXIS_RADIUS_SAMPLES')
            if row['support_start'] > cursor:
                records.append(_missing_profile_interval(cursor, row['support_start'], a, b))
            records.append(row)
            cursor = row['support_end']
        if cursor < hi:
            records.append(_missing_profile_interval(cursor, hi, a, b))
    return records


def _missing_profile_interval(start, end, a, b):
    return dict(s=(start+end)/2, support_start=start, support_end=end,
                arc_weight=end-start, radius=None, diameter=None, available=False,
                status='NOT_MEASURED', edge_nodes=[a,b], node_id=None,
                source_path_id=None, sampling='missing_interval')


def _diameter_distribution(profile, total):
    valid = [r for r in profile if r['available'] and r['arc_weight'] > 0]
    result = _distribution([r['diameter'] for r in valid], [r['arc_weight'] for r in valid])
    covered = math.fsum(r['arc_weight'] for r in valid)
    result.update(valid_arc_length=covered, total_arc_length=total,
                  valid_arc_fraction=covered/total if total > 0 else None,
                  missing_arc_length=max(0.0, total-covered),
                  status='NO_QUALIFIED_SINGLE_ROOT_RADII' if not valid else
                         'AVAILABLE' if math.isclose(covered,total,rel_tol=1e-10,abs_tol=1e-12) else 'PARTIAL_VALID_ARC',
                  weighting='qualified single-root sample support arc', missing_interpolated=False)
    return result


def _taper(axis):
    metadata = dict(definition='weighted regression d(s)=intercept+slope*s; positive taper_rate means thinning',
                    diameter_definition='twice qualified single-root area-equivalent radius',
                    scope='valid measured sample-supported arc only', missing_interpolated=False,
                    slope=None, taper_rate=None, intercept=None, r_squared=None,
                    normalized_arc_slope=None)
    if axis.get('confirmation_status')=='UNRESOLVED':
        return dict(metadata,status='UNRESOLVED_AXIS_CONTINUATION')
    if not axis['_directed'] or axis.get('closed'):
        return dict(metadata, status='ORIENTED_OPEN_AXIS_REQUIRED')
    valid = [r for r in axis['diameter_profile'] if r['available'] and r['arc_weight'] > 0]
    if len({r['s'] for r in valid}) < 3:
        return dict(metadata, status='AT_LEAST_THREE_DISTINCT_RADIUS_STATIONS_REQUIRED',
                    distinct_station_count=len({r['s'] for r in valid}))
    x = np.array([r['s'] for r in valid]); y = np.array([r['diameter'] for r in valid])
    w = np.array([r['arc_weight'] for r in valid]); weight = float(w.sum())
    xm, ym = float(np.dot(w,x)/weight), float(np.dot(w,y)/weight)
    xx = float(np.dot(w,(x-xm)**2))
    if xx <= 0:
        return dict(metadata, status='INSUFFICIENT_RADIUS_ARC_RANGE')
    slope = float(np.dot(w,(x-xm)*(y-ym))/xx)
    intercept = ym-slope*xm
    residual = float(np.dot(w,(y-(intercept+slope*x))**2))
    variance = float(np.dot(w,(y-ym)**2))
    metadata.update(slope=slope,taper_rate=-slope,intercept=intercept,
                    r_squared=None if variance <= 0 else max(0.0,1-residual/variance),
                    normalized_arc_slope=slope*axis['length'],
                    valid_arc_length=weight,valid_arc_fraction=weight/axis['length'],
                    measured_station_range=[float(x.min()),float(x.max())],
                    distinct_station_count=len(set(x)),
                    status='AVAILABLE' if math.isclose(weight,axis['length'],rel_tol=1e-10,abs_tol=1e-12) else 'PARTIAL_VALID_ARC')
    return metadata


def _caliber_window(axis, start, end, requested, config):
    lo, hi = min(axis['length'],max(0.0,start)), max(0.0,min(axis['length'],end))
    effective = max(0.0,hi-lo)
    base = dict(requested_window=requested,effective_window=effective,
                requested_arc_interval=[start,end],arc_interval=[lo,max(lo,hi)],
                window_shortened=effective < requested*(1-1e-12),
                valid_arc_length=0.0,valid_arc_fraction=None,
                definition='arc-weighted single-root diameter; mixed/hole sections excluded')
    if effective <= 0 or effective + 1e-12*requested < requested*config['minimum_window_fraction']:
        return _result(status='INSUFFICIENT_DIAMETER_WINDOW',**base)
    first=int(np.searchsorted(axis['_profile_ends'],lo,side='right'))
    last=int(np.searchsorted(axis['_profile_starts'],hi,side='left'))
    selected = [(r['diameter'],max(0.0,min(hi,r['support_end'])-max(lo,r['support_start'])))
                for r in axis['diameter_profile'][first:last] if r['available']]
    weight = math.fsum(w for _,w in selected)
    base.update(valid_arc_length=weight,valid_arc_fraction=weight/effective)
    if weight <= 0:
        return _result(status='NO_QUALIFIED_SINGLE_ROOT_RADII',**base)
    if weight/effective < config['minimum_radius_coverage']-1e-12:
        return _result(status='INSUFFICIENT_RADIUS_COVERAGE',**base)
    value = math.fsum(d*w for d,w in selected)/weight
    return _result(value,'AVAILABLE' if math.isclose(weight,effective,rel_tol=1e-10,abs_tol=1e-12) else 'PARTIAL_VALID_ARC',**base)


def _branching(axis, children):
    base = dict(lateral_axis_ids=[a['axis_id'] for a in children],
                definition='child axes attached to this oriented mother axis; graph paths are not counted as roots',
                biological_accuracy_validated=False)
    if not axis['_directed'] or axis.get('closed') or axis.get('confirmation_status')=='UNRESOLVED':
        status='UNRESOLVED_AXIS_CONTINUATION' if axis.get('confirmation_status')=='UNRESOLVED' else 'ORIENTED_OPEN_AXIS_REQUIRED'
        missing = _result(status=status)
        base['lateral_axis_ids']=[]
        return dict(base,lateral_count=dict(missing),lrbd_whole_axis=dict(missing),lrbd_branched_zone=dict(missing),
                    interbranch_distances=[],interbranch_distance_statistics=_distribution([]),
                    attachment_arc_positions=[],basal_unbranched_length=dict(missing),
                    apical_unbranched_length=dict(missing),observed_apical_unbranched_length=dict(missing),
                    observed_basal_unbranched_length=dict(missing))
    positions = sorted(float(a['_attachment']) for a in children)
    count = len(positions)
    confirmed = axis.get('confirmation_status') == 'CONFIRMED' and all(
        a.get('confirmation_status') == 'CONFIRMED' for a in children)
    auto = axis.get('confirmation_status') == 'AUTO_INFERRED' and all(
        a.get('confirmation_status') in ('AUTO_INFERRED', 'CONFIRMED') for a in children)
    status = ('ANNOTATION_CONFIRMED' if confirmed else
              'AUTO_INFERRED_AXIS_RELATIONSHIPS' if auto else 'CANDIDATE_AXIS_RELATIONSHIPS')
    whole = _result(count/axis['length'] if axis['length'] > 0 else None,
                    status if axis['length'] > 0 else 'ZERO_MOTHER_AXIS_LENGTH',
                    numerator=count,denominator=axis['length'],denominator_name='whole_observed_mother_axis_length',
                    formula='N/L_mother')
    fragment_density = None
    if axis.get('axis_scope') == 'ROOTED_BRIDGE_FRAGMENT':
        base['scope'] = 'local bridge-connected fragment; cyclic cores and cross-core relationships excluded'
        fragment_density = dict(whole, status='LOCAL_FRAGMENT_CANDIDATE_RELATIONSHIPS',
                                denominator_name='observed_bridge_fragment_length', formula='N_local/L_fragment')
        if axis.get('cyclic_boundary_node_ids'):
            whole = _result(status='CYCLIC_CORE_INTERRUPTS_COMPLETE_MOTHER_AXIS',
                            numerator=count, denominator=None, formula='N/L_complete_mother_unavailable')
    zone = positions[-1]-positions[0] if count else None
    zone_result = _result(count/zone if zone is not None and zone > 0 else None,
                          status if zone is not None and zone > 0 else 'NO_POSITIVE_BRANCHED_ZONE',
                          numerator=count,denominator=zone,denominator_name='last_minus_first_lateral_attachment_arc',
                          formula='N/(s_last-s_first); both end attachments included')
    intervals = [b-a for a,b in zip(positions,positions[1:])]
    no_branch = _result(status='NO_LATERAL_ROOTS')
    if count:
        basal_observed = _result(positions[0],status)
        apical_observed = _result(axis['length']-positions[-1],status,
                                  endpoint_type=axis.get('apical_endpoint_type','unknown'),
                                  scope='last observed lateral to observed axis end')
        known_base = axis.get('basal_endpoint_type') == 'root_base' or axis.get('parent_axis_id') is not None
        basal = basal_observed if known_base and axis.get('basal_endpoint_type') != 'cut_end' else _result(status='NATURAL_AXIS_BASE_NOT_IDENTIFIED')
        apical = _result(apical_observed['value'],status) if axis.get('apical_endpoint_type') == 'root_tip' else _result(status='TRUE_ROOT_TIP_NOT_IDENTIFIED')
    else:
        basal = apical = basal_observed = apical_observed = no_branch
    return dict(base,lateral_count=_result(count,status),lrbd_whole_axis=whole,lrbd_branched_zone=zone_result,
                lrbd_observed_fragment=fragment_density,
                attachment_arc_positions=positions,interbranch_distances=intervals,
                interbranch_distance_statistics=dict(_distribution(intervals),
                    status='AVAILABLE' if intervals else 'AT_LEAST_TWO_LATERALS_REQUIRED',
                    cv_status='AVAILABLE' if intervals and np.mean(intervals)>0 else 'NO_POSITIVE_MEAN_INTERVAL'),
                basal_unbranched_length=basal,apical_unbranched_length=apical,
                observed_basal_unbranched_length=basal_observed,observed_apical_unbranched_length=apical_observed)


def _group_summary(axes):
    total = math.fsum(a['length'] for a in axes)
    biological_axes=all(a['_directed'] and a.get('axis_scope')!='ROOTED_BRIDGE_FRAGMENT' and a.get('confirmation_status') not in ('UNRESOLVED','GEOMETRY_ONLY') for a in axes)
    return dict(axis_count=len(axes) if biological_axes else None,
                geometric_partition_count=len(axes),total_length=total,
                count_scope='oriented root-axis interpretations' if biological_axes else 'geometry partitions only; root-axis count unresolved',
                axis_length_distribution=_distribution([a['length'] for a in axes]),
                diameter_distribution=_diameter_distribution([r for a in axes for r in a['diameter_profile']],total),
                biological_accuracy_validated=False)


def _summaries(axes):
    types,orders,confirmations = defaultdict(list),defaultdict(list),defaultdict(list)
    strahler = defaultdict(lambda:dict(edges=[],profiles=[],axis_ids=set()))
    for axis in axes:
        confirmation = axis.get('confirmation_status','GEOMETRY_ONLY')
        if axis['_directed'] and confirmation not in ('UNRESOLVED','GEOMETRY_ONLY'):
            types[(axis.get('root_type','unknown'),confirmation)].append(axis)
        confirmations[confirmation].append(axis)
        order = axis.get('developmental_order')
        order_status = ('ANNOTATION_CONFIRMED' if order is not None and confirmation == 'CONFIRMED' else
                        'AUTO_INFERRED' if order is not None and confirmation == 'AUTO_INFERRED' else 'CANDIDATE')
        if order is None:
            order = axis.get('candidate_developmental_order')
        if order is not None and confirmation!='UNRESOLVED':
            orders[(order,order_status)].append(axis)
        edge_orders = axis.get('strahler_orders')
        edge_order_status = 'AUTO_INFERRED' if confirmation == 'AUTO_INFERRED' else 'ANNOTATION_CONFIRMED'
        if not edge_orders:
            edge_orders = axis.get('candidate_strahler_orders') or []
            edge_order_status = 'CANDIDATE'
        profiles_by_edge=defaultdict(list)
        for row in axis['diameter_profile']:
            profiles_by_edge[tuple(sorted(row['edge_nodes']))].append(row)
        for i,order in enumerate(edge_orders):
            if order is None or i+1 >= len(axis['node_ids']):
                continue
            item = strahler[(int(order),edge_order_status)]
            edge = tuple(sorted(axis['node_ids'][i:i+2]))
            item['edges'].append(float(axis['_arc'][i+1]-axis['_arc'][i]))
            item['profiles'].extend(profiles_by_edge[edge])
            item['axis_ids'].add(axis['axis_id'])
    return dict(
        by_root_type=[dict(root_type=key[0],confirmation_status=key[1],**_group_summary(items))
                      for key,items in sorted(types.items())],
        by_developmental_order=[dict(order=key[0],order_status=key[1],
            convention='base-origin developmental order: primary axis=0, lateral=parent+1',**_group_summary(items))
            for key,items in sorted(orders.items())],
        by_confirmation_status=[dict(confirmation_status=key,**_group_summary(items))
                                for key,items in sorted(confirmations.items())],
        by_strahler_edge_order=[dict(order=key[0],order_status=key[1],
            convention='terminal edge=1; Strahler order may change along one biological axis',
            geometric_segment_count=len(item['edges']),contributing_axis_count=len(item['axis_ids']),
            total_length=math.fsum(item['edges']),
            diameter_distribution=_diameter_distribution(item['profiles'],math.fsum(item['edges'])),
            biological_accuracy_validated=False)
            for key,item in sorted(strahler.items())])


def build_axis_traits(graph, architecture, samples, *, options=None, scale_mm_per_model_unit=None):
    """Return JSON-safe axes, branch relations and confirmation-separated summaries.

    ``graph`` uses original node coordinates. ``samples`` is the unchanged
    ``measurement_tables(report)[1]`` table. Architecture node IDs are ordered
    base-to-apex only when orientation_status is ROOTED_BASE_TO_APEX.
    """
    config = axis_trait_settings(options)
    factor = 1.0 if scale_mm_per_model_unit is None else float(scale_mm_per_model_unit)
    if isinstance(scale_mm_per_model_unit,bool) or not math.isfinite(factor) or factor<=0:
        raise ValueError('POSITIVE_FINITE_SCALE_REQUIRED')
    unit = 'mm' if scale_mm_per_model_unit is not None else graph.get('coordinate_units',architecture.get('coordinate_units','model_units'))
    points = np.asarray(graph.get('points',[]),dtype=float).reshape(-1,3)*factor
    if not np.isfinite(points).all():
        raise ValueError('FINITE_AXIS_COORDINATES_REQUIRED')
    edge_samples = defaultdict(list)
    for sample in samples:
        edge = sample.get('edge_nodes')
        if edge is None or len(edge)!=2:
            raise ValueError('RADIUS_SAMPLE_EDGE_REQUIRED')
        edge_samples[tuple(sorted(edge))].append(sample)
    axes,by_id = [],{}
    retained = ('axis_id','component_id','node_ids','edge_ids','parent_axis_id','attachment_node_id',
                'children_axis_ids','root_type','orientation_status','confirmation_status','relationship_status',
                'developmental_order','candidate_developmental_order','strahler_orders','candidate_strahler_orders',
                'strahler_order_base','candidate_strahler_order_base','basal_endpoint_type','apical_endpoint_type','closed',
                'analysis_region_id','axis_scope','cyclic_boundary_node_ids','global_hierarchy_status',
                'local_candidate_developmental_order','local_candidate_strahler_orders')
    owned_edges = set()
    for original in architecture.get('axes',[]):
        axis = {k:original[k] for k in retained if k in original}
        axis_id = axis['axis_id']
        if axis_id in by_id:
            raise ValueError('DUPLICATE_ROOT_AXIS_ID')
        ids = axis.get('node_ids',[])
        if len(ids)<2 or any(isinstance(i,bool) or not isinstance(i,(int,np.integer)) or i<0 or i>=len(points) for i in ids):
            raise ValueError('ROOT_AXIS_REQUIRES_VALID_NODE_SEQUENCE')
        p = points[ids]
        lengths = np.linalg.norm(np.diff(p,axis=0),axis=1)
        edges = [tuple(sorted(pair)) for pair in zip(ids,ids[1:])]
        if len(set(edges))!=len(edges) or owned_edges.intersection(edges):
            raise ValueError('ROOT_AXIS_EDGE_OWNERSHIP_NOT_UNIQUE')
        owned_edges.update(edges)
        arc = np.r_[0.,np.cumsum(lengths)]
        axis.update(length=float(arc[-1]),unit=unit,_points=p,_arc=arc,
                    _node_arc={node:float(s) for node,s in zip(ids,arc)},
                    _unique_nodes=len(set(ids))==len(ids),
                    _directed=axis.get('orientation_status')=='ROOTED_BASE_TO_APEX',
                    zero_length_segment_count=int(np.count_nonzero(lengths==0)),
                    biological_accuracy_validated=False,
                    length_scope='observed extracted unique edges; anatomical completeness unverified')
        axis['diameter_profile'] = _profile(axis,edge_samples,factor)
        axis['_profile_starts']=np.asarray([r['support_start'] for r in axis['diameter_profile']])
        axis['_profile_ends']=np.asarray([r['support_end'] for r in axis['diameter_profile']])
        axis['diameter_distribution'] = _diameter_distribution(axis['diameter_profile'],axis['length'])
        axis['arc_origin'] = ('oriented_geometric_fragment_start' if axis.get('confirmation_status')=='UNRESOLVED' else
                              'axis_base' if axis['_directed'] else 'arbitrary_geometric_path_start')
        axis['taper'] = _taper(axis)
        if axis.get('axis_scope') == 'ROOTED_BRIDGE_FRAGMENT':
            axis['arc_origin'] = 'root_relative_bridge_fragment_start'
            axis['length_scope'] = 'observed bridge-connected fragment; excludes unresolved cyclic-core intervals'
            axis['taper']['scope'] = 'valid measured support on this bridge fragment only; not whole-root taper across a cyclic core'
        if axis['_directed'] and not axis.get('closed') and axis.get('confirmation_status')!='UNRESOLVED':
            direction = _continuous_pca(p)
            down = np.zeros(3);down['xyz'.index(config['vertical_axis'])]=config['down_sign']
            direction['growth_angle_from_down_deg'] = None if direction.get('vector') is None else _angle(direction['vector'],down)
            direction['angle_definition'] = 'directed whole-axis PCA against down: 0=down, 90=horizontal, 180=up'
        else:
            direction = _result(status='UNRESOLVED_AXIS_CONTINUATION' if axis.get('confirmation_status')=='UNRESOLVED' else 'ORIENTED_OPEN_AXIS_REQUIRED',
                                vector=None,growth_angle_from_down_deg=None)
        axis['whole_axis_direction'] = direction
        axes.append(axis);by_id[axis_id]=axis
    children = defaultdict(list)
    for axis in axes:
        parent_id = axis.get('parent_axis_id')
        if parent_id is None:
            continue
        if parent_id not in by_id or parent_id==axis['axis_id']:
            raise ValueError('ROOT_AXIS_PARENT_NOT_FOUND')
        parent = by_id[parent_id]
        node = axis.get('attachment_node_id')
        if node not in parent['_node_arc'] or node!=axis['node_ids'][0]:
            raise ValueError('ROOT_AXIS_ATTACHMENT_NOT_ON_PARENT')
        if not parent['_unique_nodes']:
            raise ValueError('ROOT_AXIS_PARENT_ATTACHMENT_AMBIGUOUS')
        axis['_attachment'] = parent['_node_arc'][node]
        axis['attachment_arc'] = axis['_attachment']
        children[parent_id].append(axis)
    relations=[]
    for axis in axes:
        axis['branching'] = _branching(axis,children[axis['axis_id']])
        parent_id = axis.get('parent_axis_id')
        if parent_id is None:
            continue
        parent = by_id[parent_id]
        relation = dict(parent_axis_id=parent_id,child_axis_id=axis['axis_id'],
                        attachment_node_id=axis['attachment_node_id'],attachment_arc=axis['_attachment'],
                        confirmation_status=('ANNOTATION_CONFIRMED' if axis.get('confirmation_status')=='CONFIRMED' and parent.get('confirmation_status')=='CONFIRMED' else
                            'AUTO_INFERRED' if axis.get('confirmation_status') in ('AUTO_INFERRED','CONFIRMED') and parent.get('confirmation_status') in ('AUTO_INFERRED','CONFIRMED') else 'CANDIDATE'),
                        biological_accuracy_validated=False,unit=unit)
        if axis.get('axis_scope') == 'ROOTED_BRIDGE_FRAGMENT':
            relation.update(scope='within one uniquely root-relative bridge fragment; no cyclic-core crossing',
                            analysis_region_id=axis.get('analysis_region_id'),global_parent_identity_validated=False)
        unresolved=axis.get('confirmation_status')=='UNRESOLVED' or parent.get('confirmation_status')=='UNRESOLVED'
        if not parent['_directed'] or not axis['_directed'] or parent.get('closed') or axis.get('closed') or unresolved:
            if unresolved:
                relation['confirmation_status']='UNRESOLVED'
            missing=_result(status='UNRESOLVED_AXIS_CONTINUATION' if unresolved else 'ORIENTED_OPEN_PARENT_AND_CHILD_REQUIRED')
            relation.update(local_insertion_angle_deg=dict(missing),whole_axis_pca_angle_deg=dict(missing),
                            parent_diameter=dict(missing),child_diameter=dict(missing),child_parent_diameter_ratio=dict(missing))
        else:
            s=axis['_attachment'];window=config['tangent_window']
            incoming=_direction_window(parent,s-window,s,window,config)
            outgoing=_direction_window(axis,0,window,window,config)
            angle=None if incoming.get('vector') is None or outgoing.get('vector') is None else _angle(incoming['vector'],outgoing['vector'])
            relation.update(parent_tangent=incoming,child_tangent=outgoing,
                local_insertion_angle_deg=_result(angle,'AVAILABLE' if angle is not None else 'INSUFFICIENT_OR_AMBIGUOUS_TANGENT',
                    definition='incoming mother-axis tangent versus outgoing child-axis tangent; continuous arc-weighted PCA'))
            pv,cv=parent['whole_axis_direction'].get('vector'),axis['whole_axis_direction'].get('vector')
            global_angle=None if pv is None or cv is None else _angle(pv,cv)
            relation['whole_axis_pca_angle_deg']=_result(global_angle,'AVAILABLE' if global_angle is not None else 'WHOLE_AXIS_DIRECTION_UNAVAILABLE',
                definition='angle between directed whole-axis PCA vectors; distinct from local insertion angle')
            offset,width=config['diameter_offset'],config['diameter_window']
            mother=_caliber_window(parent,s-offset-width,s-offset,width,config)
            daughter=_caliber_window(axis,offset,offset+width,width,config)
            relation.update(parent_diameter=mother,child_diameter=daughter,diameter_offset=offset)
            ratio=None if mother['value'] is None or mother['value']<=0 or daughter['value'] is None else daughter['value']/mother['value']
            relation['child_parent_diameter_ratio']=_result(ratio,'AVAILABLE' if ratio is not None else 'QUALIFIED_NONZERO_PARENT_AND_CHILD_DIAMETER_REQUIRED',
                definition='child window mean diameter / proximal mother window mean diameter')
        relations.append(relation)
    summaries=_summaries(axes)
    for group in summaries.values():
        for row in group:
            row['unit']=unit
            row['diameter_distribution']['unit']=unit
            if 'axis_length_distribution' in row:
                row['axis_length_distribution']['unit']=unit
    for relation in relations:
        for key in ('local_insertion_angle_deg','whole_axis_pca_angle_deg'):
            relation[key]['unit']='degree'
        for key in ('parent_diameter','child_diameter'):
            relation[key]['unit']=unit
        relation['child_parent_diameter_ratio']['unit']='dimensionless'
    for axis in axes:
        axis['diameter_distribution']['unit']=unit
        axis['taper']['units']=dict(slope='dimensionless',taper_rate='dimensionless',
                                    intercept=unit,normalized_arc_slope=unit)
        axis['whole_axis_direction']['unit']='dimensionless_vector; growth angle in degree'
        branching=axis['branching']
        branching['lateral_count']['unit']='count'
        for key in ('lrbd_whole_axis','lrbd_branched_zone'):
            branching[key]['unit']=unit+'^-1'
            branching[key]['denominator_unit']=unit
        for key in ('basal_unbranched_length','apical_unbranched_length',
                    'observed_basal_unbranched_length','observed_apical_unbranched_length'):
            branching[key]['unit']=unit
        branching['interbranch_distance_statistics']['unit']=unit
        for key in [k for k in axis if k.startswith('_')]:
            del axis[key]
    return dict(schema_version=METHOD,method=METHOD,unit=unit,settings=config,
                physical_scale_supplied=scale_mm_per_model_unit is not None,
                scale_mm_per_model_unit=scale_mm_per_model_unit,
                window_units=unit,axes=axes,branch_relations=relations,summaries=summaries,
                radius_scope='qualified single-root radius only; no connected-section substitution',
                statistics_convention='inverse-CDF empirical quantiles; arc-weighted population SD/CV for diameters, equal-axis weighting for axis lengths',
                annotation_confirmation_is_not_biological_accuracy_validation=True,
                source_geometry_unchanged=True,missing_interpolated=False)
