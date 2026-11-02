"""Exact stage-count configuration boundaries with optimized reactor yields."""
from dataclasses import replace
import numpy as np
from scipy.optimize import brentq
from .xylose import XyloseModel


def retention_margin(model: XyloseModel, concentration: float,
                     stages: int, retention: float = 0.95):
    """g_N = Y_N*(C) - retention * Y_PFR*(C); nonnegative means feasible."""
    if stages < 1 or int(stages) != stages:
        raise ValueError('stages must be a positive integer')
    if not 0 < retention <= 1:
        raise ValueError('retention must be within (0,1]')
    if concentration <= 0:
        raise ValueError('concentration must be positive')
    _, y_n = model.optimum(concentration, n=stages)
    _, y_pfr = model.optimum(concentration)
    return y_n - retention * y_pfr


def find_concentration_boundary(model: XyloseModel, stages: int,
                                bracket: tuple[float,float],
                                retention: float = 0.95, xtol: float = 1e-12):
    """Find a bracketed zero of the exact optimized-yield margin."""
    lo, hi = bracket
    if not 0 < lo < hi:
        raise ValueError('bracket requires 0 < lo < hi')
    return brentq(lambda c: retention_margin(model, c, stages, retention),
                  lo, hi, xtol=xtol, rtol=1e-13)


def find_parameter_boundary(model: XyloseModel, concentration: float,
                            stages: int, bracket: tuple[float,float],
                            retention: float = 0.95,
                            parameter: str = 'gamma_per_molar'):
    """Find the multiplicative parameter factor at an exact stage boundary."""
    if parameter not in ('alpha','beta','gamma_per_molar'):
        raise ValueError('select alpha, beta, or gamma_per_molar')
    low, high = bracket
    if not 0 < low < high:
        raise ValueError('factor bracket must be positive and ordered')
    value = getattr(model, parameter)

    def margin(factor):
        modified = replace(model, **{parameter: value * factor})
        return retention_margin(modified, concentration, stages, retention)
    return brentq(margin, low, high, xtol=1e-12, rtol=1e-13)
