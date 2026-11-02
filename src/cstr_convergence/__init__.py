"""Finite CSTR staging and convergence calculations."""
from .core import (
    first_order_pfr, first_order_cstr, first_order_minimum_stages,
    nth_order_pfr, nth_order_cstr, nth_order_minimum_stages,
    series_pfr, series_cstr, series_cstr_exact, series_minimum_stages,
    minimum_stages, minimum_stages_scan, logarithmic_mean,
    series_asymptotic_error, series_asymptotic_sign,
    power_law_asymptotic_coefficient, first_order_asymptotic_diagnostic,
)
__version__ = "0.5.0"

from .xylose import XyloseModel
from .consecutive import optimized_series, optimized_phi_series, optimized_minimum_stages

from .benchmarks import xylose_benchmark, REFERENCE_CONCENTRATIONS

from .higher_order import (StagingCoefficients, optimized_staging_coefficients, second_order_stage_estimate)
from .boundaries import retention_margin, find_concentration_boundary, find_parameter_boundary
