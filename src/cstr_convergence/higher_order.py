"""Direct first- and second-order equal-CSTR staging corrections.

For a PFR d y/d u=f(y), finite stages use backward Euler with h=Da/N.
The modified equation is
  dy/du=f + h J f/2 + h² [J² f/3 + H[f,f]/12] + O(h³).
The ODEs below integrate the first and second state corrections directly;
no coefficient is obtained by fitting a finite-N CSTR sweep.

For optimized yield Y*(N), define
  (Y_PFR* - Y_N*)/Y_PFR* = Phi/N + Psi/N² + O(N^-3).
"""
from dataclasses import dataclass
from math import ceil, sqrt
import numpy as np
from scipy.integrate import solve_ivp
from .xylose import XyloseModel


@dataclass(frozen=True)
class StagingCoefficients:
    concentration_molar: float
    da_star: float
    y_pfr_star: float
    phi: float
    psi: float
    first_order_yield: float
    second_order_yield_fixed_da: float
    residence_optimization_correction: float
    first_order_yield_da_derivative: float
    pfr_yield_second_da_derivative: float


def optimized_staging_coefficients(model: XyloseModel, concentration: float):
    """Compute Phi and Psi for independently optimized PFR/CSTR yields.

    Analytical Jacobian and bilinear Hessian terms are integrated in the
    augmented PFR ODE. Formula includes the displacement of optimal Da
    for finite N in the second-order optimized objective.
    """
    if concentration <= 0:
        raise ValueError('Inlet concentration must be positive')
    da, yopt = model.optimum(concentration)
    gamma = model.gamma_per_molar * concentration

    def hessian_bilinear(a, b):
        # Only the I*F loss rate has nonzero Hessian entries.
        cross = -gamma * (a[1] * b[2] + a[2] * b[1])
        return np.array([0., cross, cross])

    def rhs(u, state):
        y = state[:3]
        first = state[3:6]
        second = state[6:9]
        f = model.rates(y, concentration)
        j = model.jacobian(y, concentration)
        jf = j @ f
        first_derivative = j @ first + 0.5 * jf
        second_derivative = (
            j @ second
            + 0.5 * hessian_bilinear(first, first)
            + 0.5 * (hessian_bilinear(first, f) + j @ (j @ first))
            + (j @ jf) / 3.0
            + hessian_bilinear(f, f) / 12.0
        )
        return np.concatenate((f, first_derivative, second_derivative))

    initial = np.r_[1., 0., 0., np.zeros(6)]
    sol = solve_ivp(rhs, (0., da), initial,
                    method='DOP853', rtol=1e-11, atol=1e-13)
    if not sol.success:
        raise RuntimeError(sol.message)
    y = sol.y[:3, -1]
    first = sol.y[3:6, -1]
    second = sol.y[6:9, -1]
    f = model.rates(y, concentration)
    j = model.jacobian(y, concentration)

    # h = Da/N converts the integrated state corrections to 1/N terms.
    y1 = da * first[2]
    y2_fixed_da = da**2 * second[2]
    y1_prime = first[2] + da * (j @ first + 0.5 * (j @ f))[2]
    y0_second = (j @ f)[2]
    if y0_second >= 0:
        raise RuntimeError('PFR optimum has no negative yield curvature')
    optimization_correction = -y1_prime**2 / (2 * y0_second)
    phi = -y1 / yopt
    psi = -(y2_fixed_da + optimization_correction) / yopt
    return StagingCoefficients(
        concentration_molar=float(concentration), da_star=da,
        y_pfr_star=yopt, phi=float(phi), psi=float(psi),
        first_order_yield=float(y1),
        second_order_yield_fixed_da=float(y2_fixed_da),
        residence_optimization_correction=float(optimization_correction),
        first_order_yield_da_derivative=float(y1_prime),
        pfr_yield_second_da_derivative=float(y0_second),
    )


def second_order_stage_estimate(coefficients: StagingCoefficients,
                                retention: float = 0.95):
    """Return the upper positive root of the quadratic loss estimate.

    This is a large-N estimate, not a certified finite-stage minimum.
    Call XyloseModel.stage_requirement for the exact integer design.
    """
    if not 0 < retention < 1:
        raise ValueError('retention must be between zero and one')
    epsilon = 1 - retention
    phi, psi = coefficients.phi, coefficients.psi
    discriminant = phi**2 + 4 * epsilon * psi
    if discriminant < 0:
        return None
    n = (phi + sqrt(discriminant)) / (2 * epsilon)
    return max(1, ceil(n))
