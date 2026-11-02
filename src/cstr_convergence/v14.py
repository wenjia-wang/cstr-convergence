"""Reproducible V14 benchmarks: optimized PFR vs finite equal-CSTR cascade.

Numerical targets are from the nonlinear xylose section of manuscript V14,
sharing the SI kinetic coefficients used by XyloseModel.
"""
from .xylose import XyloseModel

V14_CONCENTRATIONS = (0.033, 0.0769, 0.176, 0.340)


def v14_benchmark(model=None, concentrations=V14_CONCENTRATIONS,
                  retention=0.95, max_stages=30):
    """Compute independently optimized PFR and CSTR 95%-retention stages."""
    if model is None:
        model = XyloseModel()
    rows = []
    for c in concentrations:
        da_star, yield_star = model.optimum(c)
        rows.append({
            'C_X0_M': float(c),
            'gamma': model.gamma_per_molar * c,
            'Da_PFR_star': da_star,
            'Y_F_PFR_star': yield_star,
            'tau_PFR_star_min': da_star / model.k_x_per_second / 60.0,
            'Phi': model.phi(c),
            'N95_exact': model.stage_requirement(c, retention, max_stages),
        })
    return rows
