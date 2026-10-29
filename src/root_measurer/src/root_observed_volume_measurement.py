"""Range-aware volume of known-negative cells, without filling unknown space."""
import math
from fractions import Fraction
import numpy as np
from root_observed_measurement_input import ObservedMeasurementField


def dyadic_cell_volume(count,pitch):
    """Round the exact product of an integer count and the supplied binary pitch."""
    if isinstance(count,(bool,np.bool_)) or not isinstance(count,(int,np.integer)) or count<0:
        raise ValueError('NONNEGATIVE_INTEGER_CELL_COUNT_REQUIRED')
    if isinstance(pitch,(bool,np.bool_)) or not isinstance(pitch,(int,float,np.integer,np.floating)):
        raise ValueError('FINITE_POSITIVE_BINARY64_PITCH_REQUIRED')
    h=float(pitch)
    if not math.isfinite(h) or h<=0:raise ValueError('FINITE_POSITIVE_BINARY64_PITCH_REQUIRED')
    product=Fraction.from_float(h)**3*int(count)
    denominator=product.denominator
    assert denominator>0 and denominator & (denominator-1)==0
    coefficient=product.numerator;exponent=-(denominator.bit_length()-1)
    if coefficient:
        shift=(coefficient & -coefficient).bit_length()-1
        coefficient >>= shift;exponent+=shift
    try:rounded=float(product)
    except OverflowError:rounded=math.inf
    reason=None
    if not math.isfinite(rounded):reason='POSITIVE_VOLUME_BINARY64_OVERFLOW'
    elif coefficient and rounded==0:reason='POSITIVE_VOLUME_BINARY64_UNDERFLOW_TO_ZERO'
    return dict(value=None if reason else rounded,available=reason is None,missing_reason=reason,
        source_cell_count=int(count),source_binary64_pitch=h,
        exact_dyadic_volume=dict(nonnegative_integer_coefficient_hex=hex(coefficient),power_of_two_exponent=exponent),
        arithmetic='exact integer_count times binary64_pitch cubed, then one binary64 rounding',
        volume_is_exact_in_real_arithmetic=coefficient==0,
        rounded_binary64_value=reason is None and coefficient!=0,
        physical_or_biological_qualified=False)


def measure_observed_negative_volume(field,*,coordinate_units):
    """Public volume-only measurement of an existing complete observed field."""
    if not isinstance(field,ObservedMeasurementField):raise ValueError('QUALIFIED_OBSERVED_FIELD_REQUIRED')
    if not isinstance(coordinate_units,str) or not coordinate_units.strip():
        raise ValueError('EXPLICIT_COORDINATE_UNIT_LABEL_REQUIRED')
    negative=int(np.count_nonzero(field.tsdf<0));nonnegative=len(field.indices)-negative
    metric=dyadic_cell_volume(negative,field.pitch)
    metric.update(definition='volume of union of known-negative native lattice cells only',coordinate_unit_power=3)
    return dict(schema_version='root-observed-volume-measurement-v001',representation='native_observed_field',
        units=dict(coordinate_unit_label=coordinate_units,unit_label_is_caller_declaration=True,
                   physical_scale_independently_verified=False),
        observed_cells=len(field.indices),known_negative_cells=negative,known_nonnegative_cells=nonnegative,
        metrics=dict(known_negative_geometric_volume=metric,
            semantic_material_volume=dict(value=None,available=False,missing_reason='NEGATIVE_NATIVE_DOMAIN_NOT_SEMANTIC_MATERIAL'),
            physical_root_volume=dict(value=None,available=False,missing_reason='PHYSICAL_MATERIAL_AND_SCALE_NOT_QUALIFIED')),
        source_provenance=field.provenance,unknown_cells_filled=0,unknown_is_not_observed_empty=True,
        original_source_geometry_TSDF_and_weights_modified=False,physical_measurements_admitted=False)
