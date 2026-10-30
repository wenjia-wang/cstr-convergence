"""Validate a decision boundary independently of the printed example."""

from cstr_convergence import first_order_cstr, first_order_minimum_stages, first_order_pfr


def test_broad_da_variation_preserves_stage_choice():
    for da in [1.8, 1.9, 2.0, 2.1, 2.2]:
        assert first_order_minimum_stages(da, 0.01) == 27


def test_small_da_perturbation_crosses_integer_boundary():
    lo, hi = 2.39 * 0.99, 2.39 * 1.01
    assert first_order_minimum_stages(lo, 0.01) == 27
    assert first_order_minimum_stages(hi, 0.01) == 26
    for da in (lo, hi):
        n = first_order_minimum_stages(da, 0.01)
        assert abs(first_order_pfr(da) - first_order_cstr(da, n)) <= 0.01
        assert abs(first_order_pfr(da) - first_order_cstr(da, n - 1)) > 0.01
