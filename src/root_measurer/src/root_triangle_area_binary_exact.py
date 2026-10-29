"""Exact input-binary cross products with rounded square root and area sum.

This is triangle-soup arithmetic only: no welding, repair, unit inference,
material-union exterior identification, or biological surface qualification.
"""
from fractions import Fraction
import math
import numpy as np


def _sum_scaled(pieces):
    if not pieces: return None
    exponent = max(e for _, e in pieces)
    mantissa = math.fsum(math.ldexp(m, e-exponent) for m, e in pieces)
    return mantissa, exponent


def binary_exact_triangle_area_sum(vertices, faces, *, chunk_size=4096):
    """Preserve float64 input geometry and compute exact dyadic cross products.

    Large/small or nearly parallel edges are never classified by a threshold.
    Zero area means exact collinearity of the supplied binary coordinates.
    Irrational roots and positive summation remain rounded float64 arithmetic.
    """
    raw = np.asarray(vertices)
    if raw.ndim != 2 or raw.shape[1:] != (3,) or raw.dtype.kind not in 'fiu':
        raise ValueError('FINITE_FLOAT64_COMPATIBLE_N_BY_3_VERTICES')
    if raw.dtype.kind == 'f' and raw.dtype.itemsize > 8:
        raise ValueError('FLOAT64_COORDINATE_EXACTNESS')
    if raw.dtype.kind in 'iu' and any(abs(int(v)) > 2**53 for v in raw.flat):
        raise ValueError('FLOAT64_COORDINATE_EXACTNESS')
    points = np.asarray(raw, dtype=np.float64)
    if not np.isfinite(points).all(): raise ValueError('FINITE_FLOAT64_COMPATIBLE_N_BY_3_VERTICES')
    indices = np.asarray(faces)
    if indices.ndim != 2 or indices.shape[1:] != (3,) or indices.dtype.kind not in 'iu':
        raise ValueError('INTEGER_M_BY_3_FACES')
    if len(indices) and (indices.min() < 0 or indices.max() >= len(points)):
        raise ValueError('FACE_INDEX_RANGE')
    if isinstance(chunk_size, bool) or not isinstance(chunk_size, int) or not 1 <= chunk_size <= 1000000:
        raise ValueError('CHUNK_SIZE_RANGE')
    encoded = {}; chunks = []; pending = []; zero = 0

    def vertex(index):
        index = int(index)
        if index not in encoded:
            row = []
            for value in points[index]:
                numerator, denominator = float(value).as_integer_ratio()
                row.append((numerator, denominator.bit_length()-1))
            encoded[index] = row
        return encoded[index]

    for face in indices:
        source = [vertex(index) for index in face]
        power = max(e for row in source for _, e in row)
        a, b, c = [[n << (power-e) for n, e in row] for row in source]
        u = [b[i]-a[i] for i in range(3)]; v = [c[i]-a[i] for i in range(3)]
        cross = [u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]]
        squared = sum(value*value for value in cross)
        if squared == 0:
            zero += 1
        else:
            root_power = (squared.bit_length()-1)//2
            normalized = float(Fraction(squared, 1 << (2*root_power)))
            pending.append((math.sqrt(normalized)/2, root_power-2*power))
        if len(pending) >= chunk_size:
            chunks.append(_sum_scaled(pending)); pending = []
    if pending: chunks.append(_sum_scaled(pending))
    result = _sum_scaled(chunks)
    if result is None:
        area = 0.
    else:
        try: area = math.ldexp(*result)
        except OverflowError as error: raise ValueError('AREA_SUM_FLOAT_RANGE') from error
        if not math.isfinite(area) or area == 0: raise ValueError('AREA_SUM_FLOAT_RANGE')
    return dict(status='COMPLETE_EXACT_BINARY_CROSS_PRODUCT_TRIANGLE_AREA_SUM',
        triangle_area_sum=area, source_faces=len(indices), zero_area_faces=zero,
        exact_binary_coordinate_differences_and_cross_products=True,
        square_root_and_sum_rounded_float64=True, coordinates_or_faces_modified=False,
        material_union_outer_area=None, physical_unit=None, biological_surface_area_qualified=False)
