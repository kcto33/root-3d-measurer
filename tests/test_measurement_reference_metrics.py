import math
import numpy as np
import pytest
from root_measurement_reference_metrics import (nearest_segment_attributes,
    native_segments_from_vtp,weighted_error_summary,endpoint_diagnostics)


@pytest.mark.parametrize('seed',range(64))
def test_axis_segment_projection_exact_oracle(seed):
    q=np.random.default_rng(seed).uniform([-2,-1,-1],[2,1,2],(17,3))
    a=np.array([[-1,0,0],[0,0,0],[1,0,0]],float);b=a+[0,0,1];r=np.array([.1,.2,.3])
    got=nearest_segment_attributes(q,a,b,r,batch=3)
    expected=np.sqrt((q[:,None,0]-a[None,:,0])**2+q[:,None,1]**2+
                     (q[:,None,2]-np.clip(q[:,None,2],0,1))**2)
    ix=np.argmin(expected,axis=1)
    np.testing.assert_array_equal(got['segment_indices'],ix)
    np.testing.assert_allclose(got['distances'],expected[np.arange(len(q)),ix],atol=1e-14)
    np.testing.assert_array_equal(got['native_radius'],r[ix])


def test_one_to_one_endpoints_do_not_double_count_reference():
    q=endpoint_diagnostics([[0,0,0],[.01,0,0]],[[0,0,0]],[.1])
    assert q['one_to_one_diagnostics'][0]['one_to_one_matches']==1
    assert q['one_to_one_diagnostics'][0]['precision']==.5
    q=endpoint_diagnostics([[0,0,0],[1,0,0]],[[0,0,0],[1,0,0]],[.1])
    assert q['one_to_one_diagnostics'][0]['recall']==1
    assert endpoint_diagnostics([],[[0,0,0]],[.1])['one_to_one_diagnostics'][0]['one_to_one_matches']==0


def test_native_segment_parser_scale_and_endpoint_ids():
    xml='''<VTKFile><PolyData><Piece><Points><DataArray Name="Coordinates">0 0 0  2 0 0  2 2 0</DataArray></Points><Lines><DataArray Name="connectivity">0 1 1 2</DataArray><DataArray Name="offsets">2 4</DataArray></Lines><CellData><DataArray Name="radius">.2 .4</DataArray></CellData></Piece></PolyData></VTKFile>'''
    d=native_segments_from_vtp(xml,dict(original_bounds=[0,2,0,2,0,2],normalized_scale=.5))
    np.testing.assert_allclose(d['lengths'],[1,1]);np.testing.assert_allclose(d['radii'],[.1,.2])
    np.testing.assert_array_equal(d['endpoint_ids'],[0,2]);assert d['total_length']==2
    assert d['arc_weighted_native_radius']==pytest.approx(.15)


def test_weighted_radius_error_and_empty():
    d=weighted_error_summary([2,4],[1,2],[1,3])
    assert d['signed_bias']==pytest.approx(1.75) and d['MARE']==1
    assert d['RMSE']==pytest.approx(math.sqrt(3.25))
    assert weighted_error_summary([],[],[])['status']=='NO_VALID_WEIGHT'
