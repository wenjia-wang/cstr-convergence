"""A reproducible parameter uncertainty / integer-stage decision demonstration."""
from dataclasses import replace
from cstr_convergence import XyloseModel
m=XyloseModel(); c=0.0769
for factor in (0.7,1.0,1.3):
    x=replace(m,gamma_per_molar=m.gamma_per_molar*factor)
    print(f'gamma factor={factor:4.2f}, Phi={x.phi(c):.4f}, N95={x.stage_requirement(c,0.95,60)}')
