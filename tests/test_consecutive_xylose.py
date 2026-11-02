import numpy as np
from cstr_convergence import XyloseModel, optimized_phi_series, optimized_series, optimized_minimum_stages

def test_series_optimized_phi():
    assert abs(optimized_phi_series(1.,1.)-.5)<1e-6
    assert abs(optimized_phi_series(1.,.01)-.1082)<0.002
    assert optimized_minimum_stages(1.,1.,.95,20)==10

def test_xylose_reference_s10():
    m=XyloseModel()
    for c,phi in [(0.033,.4889),(.176,.5967),(.340,.5554)]:
        assert abs(m.phi(c)-phi)<0.009,(c,m.phi(c),phi)

def test_xylose_cstr_converges():
    m=XyloseModel();c=.176;da=m.optimum(c)[0]
    for n in (50,200):
        err=m.cstr(da,c,n)[2]-m.pfr(da,c)[2]
        assert err<0
        assert abs(n*err/m.pfr(da,c)[2]+m.phi(c))<0.04

def test_x_yield_and_nonnegative():
    m=XyloseModel()
    for n in (1,4,20):
        z=m.cstr(2,.176,n)
        assert np.min(z)>=0 and np.max(z)<=1


def test_stage_decision_crossing():
    from dataclasses import replace
    m=XyloseModel()
    c=.0769
    assert replace(m,gamma_per_molar=m.gamma_per_molar*.7).stage_requirement(c,.95,20)==10
    assert m.stage_requirement(c,.95,20)==11
