# Manuscript alignment

## v0.4.2a: V14 research baseline

This release retains the verified v0.4.2 implementation and adds an
independently executable V14 benchmark table (`examples/v14_reproduction.py`).
Its numeric targets are the optimized xylose concentration-dependent Phi
results and exact stage requirements for 95% retention:

| Inlet xylose (M) | Phi | Exact N95 |
| ---: | ---: | ---: |
| 0.033 | 0.488891 | 9 |
| 0.176 | 0.596736 | 11 |
| 0.340 | 0.555396 | 10 |

Kinetic coefficients are from the previously established SI Table S9
(160 °C, model hydronium concentration 0.10 M) and are implemented in
`src/cstr_convergence/xylose.py`. The optimum residence time is obtained
separately for PFR and each finite-CSTR cascade.

`v0.4.2a` is the manuscript-alignment label. Python distribution metadata
uses `0.4.2.post1`, a valid PEP 440 version later than `0.4.2`.

The V14 release includes first-order finite-staging sensitivity Phi,
exact finite-stage optimization, and an existing kinetic-perturbation
decision example. Direct second-order Psi and exact concentration-boundary
root finding are introduced in v0.5.0 aligned with V15.

Main associated reference: Wang, *Kinetic Criteria for Finite CSTR
Staging toward Plug-Flow Performance*, ChemRxiv (2026),
https://doi.org/10.26434/chemrxiv.15007934/v1.
