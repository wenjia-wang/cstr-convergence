"""Isothermal CSTR staging and finite-range convergence calculations."""
from .core import (
    first_order_pfr, first_order_cstr, first_order_minimum_stages,
    nth_order_pfr, nth_order_cstr, nth_order_minimum_stages,
    series_pfr, series_cstr, series_minimum_stages,
    minimum_stages, minimum_stages_scan,
)
