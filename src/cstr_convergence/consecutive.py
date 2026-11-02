"""Objective-specific, optimized finite-CSTR staging for first-order A -> B -> C."""
import numpy as np
from scipy.optimize import minimize_scalar
from .core import series_pfr, series_cstr_exact

def optimized_series(k1,k2,n=None):
    if k1<=0 or k2<=0: raise ValueError('k1,k2 must be positive')
    def yy(t):
        a,b=k1*t,k2*t
        return (series_pfr(a,b) if n is None else series_cstr_exact(a,b,n))[1]
    r=minimize_scalar(lambda t:-yy(t),bounds=(1e-8,30/min(k1,k2)),method='bounded',options={'xatol':1e-10})
    return float(r.x),float(-r.fun)

def optimized_phi_series(k1,k2):
    """Relative loss coefficient Phi=lim N*(1-Y_N*/Y_PFR*)."""
    t,y=optimized_series(k1,k2)
    a,b=k1*t,k2*t
    # y''(t)/2 for the intermediate yield; k units cancel in the normalized ratio.
    # BE e_B = (t**2/2)*d²B/dt²; since optimum y'=0, y'' = -k1*k2*y.
    return float(0.5*(k1*k2*t*t))

def optimized_minimum_stages(k1,k2,retention=0.95,max_stages=1000):
    if not 0<retention<=1: raise ValueError('retention outside (0,1]')
    _,ref=optimized_series(k1,k2)
    for n in range(1,max_stages+1):
        if optimized_series(k1,k2,n)[1]>=retention*ref: return n
    return None
