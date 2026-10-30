"""Isothermal equal-volume CSTR staging at fixed total residence time."""
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
    """Conversion for -r=k*C^order, Da=k*C0^(order-1)*tau.

    Zero-order kinetics stop at reactant depletion; subsequent tanks remain at zero.
    """
    _positive(da, 'Da', allow_zero=True)
    _positive(order, 'order', allow_zero=True)
    _stages(n)
    if order == 0:
        return float(min(da, 1.0))
    c = 1.0
    for _ in range(n):
        if da == 0 or c == 0:
            break
        cin = c
        c = brentq(lambda x: x + da / n * x ** order - cin, 0.0, cin, xtol=1e-14)
    return float(1.0 - c)


def nth_order_pfr(da: float, order: float) -> float:
    _positive(da, 'Da', allow_zero=True)
    _positive(order, 'order', allow_zero=True)
    if order == 1:
        return first_order_pfr(da)
    if order == 0:
        return float(min(da, 1.0))
    if order < 1:
        return float(1.0 - max(0.0, 1.0 + (order-1.0)*da) ** (1.0/(1.0-order)))
    return float(1.0 - (1.0 + (order-1.0)*da) ** (-1.0/(order-1.0)))


def series_pfr(da1: float, da2: float) -> np.ndarray:
    """A -> B -> C, all first order, independently stable B evaluation."""
    _positive(da1, 'Da1', allow_zero=True)
    _positive(da2, 'Da2', allow_zero=True)
    a = np.exp(-da1)
    gap = abs(da2 - da1)
    ratio = -np.expm1(-gap) / gap if gap != 0 else 1.0
    b = da1 * np.exp(-min(da1, da2)) * ratio
    return np.array([a, b, 1.0-a-b], dtype=float)


def series_cstr(da1: float, da2: float, n: int) -> np.ndarray:
    _positive(da1, 'Da1', allow_zero=True)
    _positive(da2, 'Da2', allow_zero=True)
    _stages(n)
    a, b, c = 1.0, 0.0, 0.0
    d1, d2 = da1/n, da2/n
    for _ in range(n):
        old_a, old_b, old_c = a, b, c
        a = old_a/(1.0+d1)
        b = (old_b+d1*a)/(1.0+d2)
        c = old_c+d2*b  # independent product balance
    return np.array([a, b, c], dtype=float)


def minimum_stages(error_fn, tolerance: float, max_stages: int = 100000) -> int | None:
    """Fast binary search. Valid ONLY for nonincreasing, nonnegative error(n)."""
    _positive(tolerance, 'tolerance', allow_zero=True)
    _stages(max_stages)
    cache = {}
    def err(n):
        if n not in cache:
            v = float(error_fn(n))
            if not np.isfinite(v) or v < 0:
                raise ValueError('error_fn must return finite nonnegative values')
            cache[n] = v
        return cache[n]
    if err(1) <= tolerance:
        return 1
    if max_stages == 1:
        return None
    lo, hi = 1, min(2, max_stages)
    while err(hi) > tolerance:
        if hi == max_stages:
            return None
        lo, hi = hi, min(2*hi, max_stages)
    while lo + 1 < hi:
        mid = (lo+hi)//2
        if err(mid) <= tolerance:
            hi = mid
        else:
            lo = mid
    return hi


def minimum_stages_scan(error_fn, tolerance: float, max_stages: int = 1000,
                        criterion: str = 'first') -> int | None:
    """Exhaustively scan 1..max_stages for possibly nonmonotonic errors.

    criterion='first': first n meeting tolerance.
    criterion='persistent': earliest n for which *all sampled integer stages*
    through max_stages meet tolerance. This is NOT a guarantee beyond the cap.
    """
    _positive(tolerance, 'tolerance', allow_zero=True)
    _stages(max_stages)
    if criterion not in ('first', 'persistent'):
        raise ValueError("criterion must be 'first' or 'persistent'")
    last_failure = 0
    for n in range(1, max_stages+1):
        v = float(error_fn(n))
        if not np.isfinite(v) or v < 0:
            raise ValueError('error_fn must return finite nonnegative values')
        if v <= tolerance:
            if criterion == 'first':
                return n
        else:
            last_failure = n
    if criterion == 'first':
        return None
    if last_failure == max_stages:
        return None
    return max(1, last_failure+1)


def nth_order_minimum_stages(da: float, order: float, tolerance: float,
                             max_stages: int = 1000) -> int | None:
    """First stage meeting absolute conversion tolerance, exhaustive within cap."""
    target = nth_order_pfr(da, order)
    return minimum_stages_scan(lambda n: abs(nth_order_cstr(da, order, n)-target),
                               tolerance, max_stages, 'first')


def series_minimum_stages(da1: float, da2: float, tolerance: float,
                          max_stages: int = 1000, criterion: str = 'persistent',
                          metric: str = 'B') -> int | None:
    """Search an explicitly bounded range for B yield or full-vector error.

    metric='B': absolute B yield error; metric='vector': max absolute species error.
    """
    target = series_pfr(da1, da2)
    if metric not in ('B', 'vector'):
        raise ValueError("metric must be 'B' or 'vector'")
    def error(n):
        difference = np.abs(series_cstr_exact(da1, da2, n)-target)
        return difference[1] if metric == 'B' else np.max(difference)
    return minimum_stages_scan(error, tolerance, max_stages, criterion)


def first_order_minimum_stages(da: float, tolerance: float,
                               max_stages: int = 100000) -> int | None:
    target = first_order_pfr(da)
    return minimum_stages(lambda n: abs(first_order_cstr(da, n)-target), tolerance, max_stages)



def series_cstr_exact(da1: float, da2: float, n: int) -> np.ndarray:
    """Closed-form N-CSTR A->B->C, stable for equal and disparate rates."""
    _positive(da1, 'Da1', allow_zero=True)
    _positive(da2, 'Da2', allow_zero=True)
    _stages(n)
    lo, gap = min(da1, da2), abs(da2-da1)
    a = float(np.exp(-n*np.log1p(da1/n)))
    if gap == 0:
        b = float(da1*np.exp(-(n+1)*np.log1p(da1/n)))
    else:
        w = n*np.log1p(gap/(n+lo))
        b = float(da1*np.exp(-n*np.log1p(lo/n))*(-np.expm1(-w))/gap)
    return np.array([a, b, 1.0-a-b], dtype=float)


def logarithmic_mean(a: float, b: float) -> float:
    """Positive logarithmic mean with stable near-equal limit."""
    _positive(a, 'a')
    _positive(b, 'b')
    lo, hi = min(a,b), max(a,b)
    q = (hi-lo)/lo
    return float(lo * q / np.log1p(q)) if q != 0 else float(lo)


def series_asymptotic_error(da1: float, da2: float) -> np.ndarray:
    """Leading coefficient K in y_N-y_PFR = K/N + O(N^-2)."""
    _positive(da1, 'Da1', allow_zero=True)
    _positive(da2, 'Da2', allow_zero=True)
    p = series_pfr(da1, da2)
    # y''(1) = M^2 y(1); the dimensionless rate matrix M is constant.
    a, b, _ = p
    k_a = 0.5 * da1**2 * a
    dA = -da1*a
    dB = da1*a-da2*b
    ddB = da1*dA-da2*dB
    k_b = 0.5*ddB
    return np.array([k_a, k_b, -k_a-k_b], dtype=float)


def power_law_asymptotic_coefficient(da: float, order: float) -> float:
    """K in X_PFR-X_N = K/N + O(N^-2), for positive PFR exit c.

    Returns nan at finite-time depletion for 0<order<1, since
    standard smooth asymptotics do not apply there.
    """
    _positive(da, 'Da', allow_zero=True)
    _positive(order, 'order', allow_zero=True)
    if order == 0 or da == 0:
        return 0.0
    c = 1.0-nth_order_pfr(da, order)
    if c <= 0:
        return float('nan')
    return float(0.5*order*da*c**order*(-np.log(c)))


def first_order_asymptotic_diagnostic(da: float, tolerance: float) -> dict:
    """Leading estimate and second/first term ratio (NOT a certified bound)."""
    _positive(da, 'Da', allow_zero=True)
    _positive(tolerance, 'tolerance')
    k = 0.5*da*da*np.exp(-da)
    estimated_n = k/tolerance
    ratio = abs((da*da/4-2*da/3)/estimated_n) if estimated_n>0 else float('nan')
    return dict(estimated_stages=float(estimated_n), second_to_first=float(ratio),
                validated=False)


def series_asymptotic_sign(da1: float, da2: float) -> int:
    """Sign of large-N B_N-B_PFR; 0 denotes zero leading coefficient."""
    _positive(da1, 'Da1')
    _positive(da2, 'Da2')
    mean = logarithmic_mean(da1, da2)
    return 1 if mean>2 else (-1 if mean<2 else 0)

def _positive(x, label, allow_zero=False):
    if not np.isfinite(x) or (x < 0 if allow_zero else x <= 0):
        raise ValueError(f'{label} must be finite and {"nonnegative" if allow_zero else "positive"}')


def _stages(n):
    if isinstance(n, (bool, np.bool_)) or not isinstance(n, (int, np.integer)) or n < 1:
        raise ValueError('number of stages must be a positive integer')
