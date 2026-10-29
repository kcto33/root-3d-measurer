"""Numerical range binding for plane normals; original section kernel intact."""
import math
import numpy as np
from root_cubical_sections import CubicalSectionMeasurer


def bind_section_normal(normal):
    array=np.asarray(normal,dtype=float)
    if array.shape!=(3,) or not np.isfinite(array).all():
        raise ValueError('FINITE_SECTION_PLANE_NORMAL_REQUIRED')
    largest=float(np.max(np.abs(array)))
    if largest==0:raise ValueError('NONZERO_SECTION_NORMAL_REQUIRED')
    exponent=math.frexp(largest)[1]
    # max component >= 2**(e-1); sum of 3 squares < 3*2**(2e).
    # These bounds keep that squared norm normal and below float64 overflow.
    minimum=(np.finfo(float).minexp+2)//2
    maximum=(np.finfo(float).maxexp-2)//2
    safe=minimum<=exponent<=maximum
    if safe:
        effective=normal
    else:
        scaled=array/largest
        effective=scaled/math.hypot(*scaled)
    return effective,dict(source_normal_passed_through=safe,normal_rescaled_for_numeric_range=not safe,
        source_max_component_frexp_exponent=exponent,square_sum_safe_exponents=[int(minimum),int(maximum)],
        physical_or_geometry_gate_changed=False)


class RangeSafeCubicalSectionMeasurer:
    def __init__(self,solid):
        self.original=CubicalSectionMeasurer(solid)

    def section(self,point,normal):
        effective,receipt=bind_section_normal(normal)
        result=self.original.section(point,effective)
        return dict(result,normal_numeric_range_receipt=receipt)
