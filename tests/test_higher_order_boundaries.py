import pytest
from cstr_convergence import (
    XyloseModel, optimized_staging_coefficients,
    second_order_stage_estimate, retention_margin,
    find_concentration_boundary, find_parameter_boundary,
)

@pytest.mark.parametrize('c,phi,psi,n95,n_second', [
    (0.033, 0.488891, -0.579825, 9, 9),
    (0.176, 0.596736, -1.344407, 11, 9),
    (0.340, 0.555396, -1.354615, 10, 8),
])
def test_direct_psi_and_integer_stage(c, phi, psi, n95, n_second):
    m = XyloseModel()
    coeff = optimized_staging_coefficients(m, c)
    assert coeff.phi == pytest.approx(phi, abs=3e-6)
    assert coeff.psi == pytest.approx(psi, abs=3e-6)
    assert coeff.phi == pytest.approx(m.phi(c), abs=1e-8)
    assert m.stage_requirement(c, .95, 15) == n95
    assert second_order_stage_estimate(coeff, .95) == n_second

@pytest.mark.parametrize('bracket,expected', [
    ((0.06, 0.09), 0.0726479195),
    ((0.28, 0.33), 0.3082568016),
])
def test_exact_ten_stage_concentration_boundaries(bracket, expected):
    m = XyloseModel()
    concentration = find_concentration_boundary(m, 10, bracket)
    assert concentration == pytest.approx(expected, abs=2e-7)
    assert abs(retention_margin(m, concentration, 10, .95)) < 1e-9

def test_exact_cross_degradation_boundary():
    m = XyloseModel()
    factor = find_parameter_boundary(m, .0769, 10, (.90, 1.0))
    assert factor == pytest.approx(0.944706366, abs=2e-6)

def test_optimization_shift_nonnegative():
    coeff = optimized_staging_coefficients(XyloseModel(), .033)
    assert coeff.residence_optimization_correction >= 0


def test_direct_psi_against_large_n_full_cascade():
    """Independent finite-CSTR optimization corroborates direct Psi."""
    m = XyloseModel()
    c = 0.176
    coeff = optimized_staging_coefficients(m, c)
    n = 400
    _, optimized_n_yield = m.optimum(c, n)
    exact_loss = (coeff.y_pfr_star - optimized_n_yield) / coeff.y_pfr_star
    direct_second_order_loss = coeff.phi / n + coeff.psi / n**2
    assert abs(exact_loss - direct_second_order_loss) < 5e-7
