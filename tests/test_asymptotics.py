import numpy as np
from cstr_convergence import (
    series_pfr, series_cstr, series_cstr_exact, series_asymptotic_error,
    series_asymptotic_sign, logarithmic_mean, power_law_asymptotic_coefficient,
    first_order_asymptotic_diagnostic, nth_order_cstr, nth_order_pfr,
    __version__,
)

def test_version():
    assert __version__ == '0.4.0'

def test_closed_form_matches_independent_recursion():
    for a,b in [(0.5,3),(2,1),(6,0.4),(800,1),(1,800),(2,2),(0,2),(2,0)]:
        for n in (1,10,100):
            np.testing.assert_allclose(series_cstr_exact(a,b,n),series_cstr(a,b,n),atol=2e-12,rtol=2e-12)

def test_logarithmic_mean():
    assert abs(logarithmic_mean(2,2)-2)<1e-14
    assert abs(logarithmic_mean(.5,3)-((3-.5)/np.log(3/.5)))<1e-13
    assert abs(logarithmic_mean(2,2+1e-12)-(2+0.5e-12))<1e-12

def test_asymptotic_series():
    for a,b in [(2,1),(6,.4),(.5,3),(2,2)]:
        n=20000
        delta = n*(series_cstr_exact(a,b,n)-series_pfr(a,b))
        np.testing.assert_allclose(delta,series_asymptotic_error(a,b),atol=0.0004,rtol=0.004)

def test_sign():
    assert series_asymptotic_sign(.5,3) == -1
    assert series_asymptotic_sign(6,.4) == 1
    assert series_asymptotic_sign(2,2) == 0

def test_power_law():
    for order, da in ((.5,.7),(1,2),(2,2),(3,1)):
        n=5000
        k=power_law_asymptotic_coefficient(da,order)
        e=n*(nth_order_pfr(da,order)-nth_order_cstr(da,order,n))
        assert abs(k-e)<0.001
    assert np.isnan(power_law_asymptotic_coefficient(3,.5))

def test_diagnostic():
    d=first_order_asymptotic_diagnostic(10,.01)
    assert d['estimated_stages'] < 1
    assert d['second_to_first']>1
    assert d['validated'] is False
