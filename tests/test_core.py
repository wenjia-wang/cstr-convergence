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


def test_zero_order_depletion():
    from cstr_convergence import nth_order_minimum_stages
    for da in (0., 0.3, 1., 2., 10.):
        for n in (1, 2, 13):
            assert nth_order_cstr(da, 0, n) == pytest.approx(min(da, 1.))
        assert nth_order_minimum_stages(da, 0, 1e-12, max_stages=10) == 1


def test_fractional_order_and_pfr_limit():
    for order in (0.5, 1.5, 2.):
        assert nth_order_cstr(2, order, 2500) == pytest.approx(
            nth_order_pfr(2, order), abs=5e-4)


def test_stage_cap_and_cache():
    assert minimum_stages(lambda n: 1/n, .6, max_stages=1) is None
    observed = []
    def error(n):
        observed.append(n)
        return 1/n
    assert minimum_stages(error, .01, max_stages=128) == 100
    assert len(observed) == len(set(observed))
    assert max(observed) <= 128


def test_nonmonotone_search():
    from cstr_convergence import minimum_stages_scan
    error = lambda n: {1: .001, 2: .1, 3: .001, 4: .001}.get(n, 0.)
    assert minimum_stages_scan(error, .01, max_stages=4, criterion='first') == 1
    assert minimum_stages_scan(error, .01, max_stages=4, criterion='persistent') == 3
    assert minimum_stages_scan(error, .01, max_stages=2, criterion='persistent') is None


def test_series_scan_against_bruteforce():
    from cstr_convergence import series_minimum_stages
    a, b, tol, cap = 6., .4, 5e-5, 200
    target = series_pfr(a, b)[1]
    passing = [n for n in range(1, cap+1)
               if abs(series_cstr(a, b, n)[1]-target) <= tol]
    assert series_minimum_stages(a, b, tol, cap, 'first') == (min(passing) if passing else None)
    failures = [n for n in range(1, cap+1) if n not in passing]
    expected = None if cap in failures else max([0]+failures)+1
    assert series_minimum_stages(a, b, tol, cap, 'persistent') == expected


def test_nearly_equal_rate_constants():
    d1 = 2.
    for delta in (0, 1e-12, -1e-12, 1e-9, -1e-9):
        p = series_pfr(d1, d1+delta)
        assert p[1] == pytest.approx(d1*np.exp(-d1), rel=2e-9)
        assert sum(p) == pytest.approx(1.)


def test_series_pfr_extreme_rate_separation():
    # Closed-form reference: with very fast A consumption, B tends to exp(-Da2).
    for d1, d2 in ((800., 1.), (1000., 0.5)):
        values = series_pfr(d1, d2)
        assert np.all(np.isfinite(values))
        assert values[1] == pytest.approx(
            d1 * np.exp(-d2) / (d1-d2), rel=1e-12)
        assert np.sum(values) == pytest.approx(1.)
    for d1, d2 in ((0.5, 1000.), (1., 800.)):
        values = series_pfr(d1, d2)
        assert np.all(np.isfinite(values))
        assert values[1] == pytest.approx(
            d1 * np.exp(-d1) / (d2-d1), rel=1e-12)
        assert np.sum(values) == pytest.approx(1.)


def test_series_pfr_against_independent_ode():
    from scipy.integrate import solve_ivp
    for d1, d2 in ((2., 0.5), (0.5, 2.), (2., 2.), (6., 0.4)):
        def rhs(_, y):
            a, b, c = y
            return [-d1*a, d1*a-d2*b, d2*b]
        ref = solve_ivp(rhs, (0, 1), [1, 0, 0], method='DOP853',
                        rtol=1e-12, atol=1e-13).y[:, -1]
        assert series_pfr(d1, d2) == pytest.approx(ref, abs=2e-12)
