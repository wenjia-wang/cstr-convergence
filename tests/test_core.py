import numpy as np
import pytest
from cstr_convergence import (
    first_order_pfr, first_order_cstr, first_order_minimum_stages,
    nth_order_pfr, nth_order_cstr, series_pfr, series_cstr, minimum_stages
)

def test_analytic_first_order():
    assert first_order_cstr(2., 1) == pytest.approx(2/3)
    assert first_order_pfr(2.) == pytest.approx(1-np.exp(-2))

def test_nth_order_first_order_agreement():
    for n in (1, 4, 20):
        assert nth_order_cstr(2., 1., n) == pytest.approx(first_order_cstr(2., n))
    assert nth_order_pfr(2., 1.) == pytest.approx(first_order_pfr(2.))

def test_pfr_limit():
    assert first_order_cstr(2., 100000) == pytest.approx(first_order_pfr(2.), abs=1e-5)
    assert nth_order_cstr(2., 2., 100) == pytest.approx(nth_order_pfr(2., 2.), abs=0.005)

def test_minimality():
    da, eps = 2., 0.01
    n = first_order_minimum_stages(da, eps)
    assert n is not None
    assert abs(first_order_cstr(da, n)-first_order_pfr(da)) <= eps
    assert abs(first_order_cstr(da, n-1)-first_order_pfr(da)) > eps

def test_series_mass_conservation():
    for fn in (lambda: series_pfr(2, 1), lambda: series_cstr(2, 1, 8), lambda: series_pfr(2, 2)):
        c = fn()
        assert np.sum(c) == pytest.approx(1.)
        assert np.all(c >= -1e-12)

def test_series_pfr_limit():
    assert series_cstr(2, 1, 5000) == pytest.approx(series_pfr(2, 1), abs=2e-4)

def test_zero_da_and_inputs():
    assert first_order_cstr(0., 4) == 0
    assert np.array_equal(series_cstr(0, 0, 5), [1, 0, 0])
    with pytest.raises(ValueError): first_order_cstr(1, 0)
    with pytest.raises(ValueError): first_order_cstr(-1, 3)
    with pytest.raises(ValueError): first_order_minimum_stages(2, -0.1)

def test_no_solution_under_cap():
    assert minimum_stages(lambda n: 1/n, 0.001, max_stages=10) is None
