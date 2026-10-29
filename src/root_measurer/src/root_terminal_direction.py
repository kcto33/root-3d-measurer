"""Scale-invariant geometric terminal direction; no physical-tip inference."""
import math
import numpy as np


def terminal_direction(points, pitch):
    """Same 21-sample/10-pitch linear OLS definition, dimensionless arithmetic.

    Returned numerical bounds concern weighted floating arithmetic only. They
    are not statistical uncertainty or a physical tangent-error certificate.
    """
    pts=np.asarray(points,dtype=float)
    if pts.ndim!=2 or pts.shape[1]!=3 or len(pts)<2 or not np.isfinite(pts).all():
        raise ValueError('FINITE_OPEN_3D_BRANCH_AT_LEAST_TWO_POINTS_REQUIRED')
    if isinstance(pitch,bool) or not math.isfinite(pitch) or pitch<=0:
        raise ValueError('FINITE_POSITIVE_PITCH_REQUIRED')
    if np.array_equal(pts[0],pts[-1]):raise ValueError('CLOSED_BRANCH_HAS_NO_TERMINAL_DIRECTION')
    with np.errstate(over='ignore',invalid='ignore'):
        delta=np.diff(pts,axis=0)
        relative=pts-pts[0]
    if not np.isfinite(delta).all() or not np.isfinite(relative).all():
        raise ValueError('COORDINATE_DIFFERENCE_EXCEEDS_FLOAT_RANGE')
    steps=np.array([math.hypot(*row) for row in delta])
    if not np.isfinite(steps).all() or np.any(steps<=0):
        raise ValueError('POSITIVE_FINITE_SOURCE_EDGE_LENGTHS_REQUIRED_NO_SILENT_DELETION')
    scale=float(steps.max())
    arc=np.r_[0.,np.cumsum(steps/scale)]
    if np.any(np.diff(arc)<=0):raise ValueError('ARC_INCREMENT_BELOW_FLOAT_RESOLUTION')
    scaled_pitch=float(pitch)/scale
    support=min(10*scaled_pitch,float(arc[-1]))
    if not support>0 or not math.isfinite(support):
        raise ValueError('SUPPORT_EXCEEDS_FLOAT_RESOLUTION')
    locations=np.linspace(0,support,21)
    local=relative/scale
    q=np.column_stack([np.interp(locations,arc,local[:,axis]) for axis in range(3)])/support
    if not np.isfinite(q).all():raise ValueError('LOCAL_COORDINATES_EXCEED_FLOAT_RANGE')
    # For t=i/20: sum(t-.5)^2=77/40; derivative weights=2(i-10)/77.
    weights=[2*(i-10)/77 for i in range(21)]
    terms=[[weights[i]*float(q[i,axis]) for i in range(21)] for axis in range(3)]
    derivative=np.array([math.fsum(row) for row in terms])
    epsilon=np.finfo(float).eps;gamma=(5*21*epsilon)/(1-5*21*epsilon)
    bounds=[gamma*math.fsum(abs(v) for v in row) for row in terms]
    norm=math.hypot(*derivative);resolution=math.hypot(*bounds)
    base=dict(definition='linear OLS on 21 uniform arc samples over min(10pitch,branchlength)',
        samples=21,source_point_count=len(pts),window_pitch_multiple=10,
        endpoint_moved=False,source_geometry_changed=False,adopted_extension=0.,
        physical_tip=None,physical_direction_error_bound=None,
        statistical_confidence_interval=None,
        weighted_arithmetic_resolution_only=True,
        normalized_derivative=derivative.tolist(),weighted_arithmetic_resolution=bounds,
        interpolation_and_input_quantization_error_certified=False)
    if norm<=resolution:
        return dict(base,status='MISSING_NUMERICALLY_UNRESOLVED_WEIGHTED_DERIVATIVE',direction=None)
    direction=(derivative/norm).tolist()
    return dict(base,status='COMPLETE_SCALE_INVARIANT_GEOMETRIC_TERMINAL_DIRECTION',direction=direction)
