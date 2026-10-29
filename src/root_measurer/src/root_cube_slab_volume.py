"""Exact volume of immutable canonical cubes between parallel binary planes."""
from fractions import Fraction
import itertools
import math
import numpy as np


def _integers(values):
    pairs=[float(v).as_integer_ratio() for v in values]
    power=max(d.bit_length()-1 for n,d in pairs)
    return [n << (power-(d.bit_length()-1)) for n,d in pairs],power


def cube_union_slab_volume(solid,normal,first_plane_point,second_plane_point):
    n=np.asarray(normal,float);a=np.asarray(first_plane_point,float);b=np.asarray(second_plane_point,float)
    if any(v.shape!=(3,) or not np.isfinite(v).all() for v in (n,a,b)) or not np.any(n):
        raise ValueError('FINITE_PARALLEL_SLAB_PLANES_REQUIRED')
    ns,npower=_integers(n);cs,cpower=_integers([solid.pitch,*a,*b]);h=cs[0];pa=cs[1:4];pb=cs[4:7]
    levels=sorted((sum(x*y for x,y in zip(ns,pa)),sum(x*y for x,y in zip(ns,pb))))
    weights=[abs(x)*h for x in ns if x];dimensions=len(weights);span=sum(weights)
    denominator=math.factorial(dimensions)*math.prod(weights)
    subsets=[(sum(s),(-1)**len(s)) for size in range(dimensions+1) for s in itertools.combinations(weights,size)]
    def below(t):
        if t<=0:return Fraction(0)
        if t>=span:return Fraction(1)
        return Fraction(sum(sign*max(t-shift,0)**dimensions for shift,sign in subsets),denominator)
    fraction=Fraction(0);cut=[]
    for i,row in enumerate(solid.cells):
        base=sum(ni*(int(index)*h if ni>0 else (int(index)+1)*h) for ni,index in zip(ns,row))
        ratio=below(levels[1]-base)-below(levels[0]-base)
        assert 0<=ratio<=1
        if ratio:cut.append(i);fraction+=ratio
    exact=fraction*Fraction(float(solid.pitch))**3
    return dict(status='COMPLETE_EXACT_CUBE_SLAB_VOLUME',volume=float(exact),
        exact_volume_numerator=str(exact.numerator),exact_volume_denominator=str(exact.denominator),
        intersecting_source_cell_indices=cut,source_cells_removed_or_modified=False,
        formula='uniform-sum CDF inclusion-exclusion; exact dyadic binary planes and pitch',
        active_normal_dimensions=dimensions,normal_integer_power=npower,coordinate_integer_power=cpower)
