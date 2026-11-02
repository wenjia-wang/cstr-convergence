# Numerical verification and reference cases

This document gives reproducible targets for the finite-CSTR calculations. All results are computed from the same four-reaction xylose kinetic model, without fitting to finite-stage output.

## Objective and finite-stage coefficients

For a PFR and an equal-volume CSTR cascade, independently optimized with respect to residence time, define the normalized objective loss

`epsilon_N = (Y_PFR* - Y_N*) / Y_PFR* = Phi/N + Psi/N² + O(N^-3)`.

`higher_order.py` integrates the first two modified-equation state corrections, including the residence-time optimization contribution to `Psi`. `boundaries.py` calculates the exact optimized-yield margin `g_N(C) = Y_N*(C) - p Y_PFR*(C)` and finds roots using bracketed scalar solvers.

## Reference results

| Inlet xylose (M) | Phi | Psi | Exact N95 | Second-order estimate |
| ---: | ---: | ---: | ---: | ---: |
| 0.033 | 0.488891 | -0.579825 | 9 | 9 |
| 0.176 | 0.596736 | -1.344407 | 11 | 9 |
| 0.340 | 0.555396 | -1.354615 | 10 | 8 |

At a 95% retention target, the exact 10-stage yield margin has concentration roots of **0.07264792** and **0.30825680 mol/L**. Between these roots, the 10-stage cascade falls below the target. At 0.0769 mol/L, the I+F cross-degradation multiplier that makes the margin zero is **0.94470637** relative to the reference value.

The difference between exact N95 and second-order estimates at two concentrations demonstrates why integer-stage design is evaluated using finite-CSTR optimization.

## Run the verification

```bash
python -m pip install -e ".[dev,plot]"
python -m pytest -q
python examples/xylose_benchmarks.py
python examples/second_order_boundaries.py --plot
```

`tests/test_benchmarks.py` and `tests/test_higher_order_boundaries.py` verify the coefficients, reference stage requirements, boundary roots and a large-stage convergence check.

See [KINETIC_MODEL.md](KINETIC_MODEL.md) for parameter provenance, model definitions and the underlying research citation.
