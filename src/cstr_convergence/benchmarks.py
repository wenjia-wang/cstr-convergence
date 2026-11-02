"""Numerical benchmarks for optimized isothermal xylose reactor staging.

Includes PFR reference performance, exact finite-stage requirements,
and the objective-specific leading finite-staging coefficient Phi.
"""

from .xylose import XyloseModel

REFERENCE_CONCENTRATIONS = (0.033, 0.0769, 0.176, 0.340)


def xylose_benchmark(model=None, concentrations=REFERENCE_CONCENTRATIONS,
                     retention=0.95, max_stages=30):
    """Return reference metrics for independently optimized PFR/CSTR models."""
    if model is None:
        model = XyloseModel()
    rows = []
    for concentration in concentrations:
        da_star, yield_star = model.optimum(concentration)
        rows.append({
            "C_X0_M": float(concentration),
            "gamma": model.gamma_per_molar * concentration,
            "Da_PFR_star": da_star,
            "Y_F_PFR_star": yield_star,
            "tau_PFR_star_min": da_star / model.k_x_per_second / 60.0,
            "Phi": model.phi(concentration),
            "N95_exact": model.stage_requirement(
                concentration, retention, max_stages
            ),
        })
    return rows
