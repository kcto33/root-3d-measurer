"""Fixed linear/quadratic terminal direction sensitivity, never model selection."""
from fractions import Fraction
import math
import numpy as np
from root_terminal_direction import terminal_direction


def _weights(degree,coefficient):
    t=[Fraction(i,20) for i in range(21)]
    matrix=[[sum(x**(i+j) for x in t) for j in range(degree+1)]+
            [x**i for x in t] for i in range(degree+1)]
    for i in range(degree+1):
        pivot=matrix[i][i];matrix[i]=[v/pivot for v in matrix[i]]
        for j in range(degree+1):
            if i==j:continue
            multiple=matrix[j][i]
            matrix[j]=[a-multiple*b for a,b in zip(matrix[j],matrix[i])]
    return tuple(matrix[coefficient][degree+1:])


WEIGHTS=[_weights(2,k) for k in range(3)]
LINEAR=[Fraction(2*(i-10),77) for i in range(21)]


def terminal_model_diagnostic(points,pitch):
    linear=terminal_direction(points,pitch)
    base=dict(linear=linear,source_endpoint_changed=False,source_geometry_changed=False,
        physical_tip=None,physical_tangent=None,adopted_extension=0.,
        model_adopted=False,model_selected_by_fit_residual=False,
        sample_count=21,window_pitch_multiple=10,quadratic_degree=2,
        noise_gain_is_IID_sample_assumption_not_actual_interpolated_noise=True)
    if linear['direction'] is None:
        return dict(base,status='MISSING_NUMERICALLY_UNRESOLVED_ORIGINAL_GEOMETRIC_DIRECTION',quadratic_direction=None)
    p=np.asarray(points,float);delta=np.diff(p,axis=0)
    steps=np.array([math.hypot(*row) for row in delta]);scale=float(steps.max())
    arc=np.r_[0.,np.cumsum(steps/scale)]
    support=min(10*(pitch/scale),float(arc[-1]));local=(p-p[0])/scale
    locations=np.linspace(0,support,21)
    q=np.column_stack([np.interp(locations,arc,local[:,i]) for i in range(3)])/support
    coefficients=np.array([[math.fsum(float(w)*float(x) for w,x in zip(weights,q[:,axis]))
        for axis in range(3)] for weights in WEIGHTS])
    derivative=coefficients[1];norm=math.hypot(*derivative)
    eps=np.finfo(float).eps;gamma=105*eps/(1-105*eps)
    bounds=np.array([gamma*math.fsum(abs(float(w)*float(x)) for w,x in zip(WEIGHTS[1],q[:,axis])) for axis in range(3)])
    if norm<=math.hypot(*bounds):
        return dict(base,status='MISSING_NUMERICALLY_UNRESOLVED_QUADRATIC_DERIVATIVE',quadratic_direction=None)
    direction=derivative/norm
    t=np.linspace(0,1,21)
    quadratic=np.c_[np.ones(21),t,t*t]@coefficients
    linear_slope=np.array(linear['normalized_derivative'])
    intercept=q.mean(axis=0)-.5*linear_slope
    linear_fit=intercept+t[:,None]*linear_slope
    old=np.asarray(linear['direction'])
    separation=math.degrees(math.atan2(math.hypot(*np.cross(old,direction)),float(old@direction)))
    gain=math.sqrt(float(sum(w*w for w in WEIGHTS[1])/sum(w*w for w in LINEAR)))
    return dict(base,status='COMPLETE_TERMINAL_CURVATURE_MODEL_SENSITIVITY_NOT_PHYSICAL_QUALIFICATION',
        quadratic_direction=direction.tolist(),model_separation_degrees=separation,
        linear_residual_RMS=float(np.sqrt(np.mean((q-linear_fit)**2))),
        quadratic_residual_RMS=float(np.sqrt(np.mean((q-quadratic)**2))),
        quadratic_derivative_IID_noise_gain_vs_linear=gain,
        normalized_quadratic_coefficients=coefficients.tolist(),
        weighted_arithmetic_resolution=bounds.tolist(),interpolation_uncertainty_certified=False)
