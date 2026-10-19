"""Steady-state equal-volume CSTR staging with fixed total residence time."""
from __future__ import annotations
import numpy as np
from scipy.optimize import brentq


def first_order_pfr(da: float) -> float:
    _positive(da, 'Da', allow_zero=True)
    return float(-np.expm1(-da))


def first_order_cstr(da: float, n: int) -> float:
    _positive(da, 'Da', allow_zero=True)
    _stages(n)
    return float(-np.expm1(-n * np.log1p(da / n)))


def nth_order_cstr(da: float, order: float, n: int) -> float:
    """Conversion for r=k*C^order, Da=k*C0^(order-1)*tau."""
    _positive(da, 'Da', allow_zero=True)
    _positive(order, 'order', allow_zero=True)
    _stages(n)
    c = 1.0
    for _ in range(n):
        cin = c
        if da == 0:
            continue
        c = brentq(lambda x: x + da / n * x ** order - cin, 0.0, cin, xtol=1e-14)
    return float(1.0 - c)


def nth_order_pfr(da: float, order: float) -> float:
    _positive(da, 'Da', allow_zero=True)
    _positive(order, 'order', allow_zero=True)
    if order == 1:
        return first_order_pfr(da)
    if order == 0:
        return min(da, 1.0)
    if order < 1:
        return float(1.0 - max(0.0, 1.0 + (order-1.0)*da) ** (1.0/(1.0-order)))
    return float(1.0 - (1.0 + (order-1.0)*da) ** (-1.0/(order-1.0)))


def series_pfr(da1: float, da2: float) -> np.ndarray:
    """A -> B -> C, all first order; dimensionless outlet concentrations."""
    _positive(da1, 'Da1', allow_zero=True)
    _positive(da2, 'Da2', allow_zero=True)
    a = np.exp(-da1)
    b = da1 * np.exp(-da1) if np.isclose(da1, da2, rtol=1e-10, atol=1e-12) else da1 * (np.exp(-da1)-np.exp(-da2))/(da2-da1)
    return np.array([a, b, 1.0-a-b], dtype=float)


def series_cstr(da1: float, da2: float, n: int) -> np.ndarray:
    _positive(da1, 'Da1', allow_zero=True)
    _positive(da2, 'Da2', allow_zero=True)
    _stages(n)
    a, b, c = 1., 0., 0.
    d1, d2 = da1/n, da2/n
    for _ in range(n):
        old_a, old_b = a, b
        a = old_a/(1.+d1)
        b = (old_b + d1*a)/(1.+d2)
        c = 1.-a-b
    return np.array([a, b, c], dtype=float)


def minimum_stages(error_fn, tolerance: float, max_stages: int = 100000) -> int | None:
    """Find minimal stage count for a MONOTONIC error measure; None if unresolved."""
    _positive(tolerance, 'tolerance')
    _stages(max_stages)
    if error_fn(1) <= tolerance:
        return 1
    lo, hi = 1, 2
    while hi < max_stages and error_fn(hi) > tolerance:
        lo, hi = hi, min(2*hi, max_stages)
    if error_fn(hi) > tolerance:
        return None
    while lo+1 < hi:
        mid = (lo+hi)//2
        if error_fn(mid) <= tolerance:
            hi = mid
        else:
            lo = mid
    return hi


def first_order_minimum_stages(da: float, tolerance: float, max_stages: int = 100000) -> int | None:
    target = first_order_pfr(da)
    return minimum_stages(lambda n: abs(first_order_cstr(da, n)-target), tolerance, max_stages)


def _positive(x, label, allow_zero=False):
    if not np.isfinite(x) or (x < 0 if allow_zero else x <= 0):
        raise ValueError(f'{label} must be finite and {"nonnegative" if allow_zero else "positive"}')


def _stages(n):
    if isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer)) or n < 1:
        raise ValueError('number of stages must be a positive integer')
