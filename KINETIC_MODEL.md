# Xylose reaction network and reproducibility

## Kinetic model

The isothermal network contains four reactions: `X -> I`, `I -> F`, `F -> degradation` and `I + F -> degradation`. Concentrations are normalized by the inlet xylose concentration `C_X0`.

Use `Da = k_X t`, `alpha = k_I/k_X`, `beta = k_F/k_X`, and `gamma = k_FI C_X0/k_X`, with furfural yield `F` as the objective. The implemented reference case uses **160 °C** and model hydronium concentration **0.10 M**.

The PFR solves `dy/dDa = f(y)`. Each ideal stage satisfies `y_j = y_(j-1) + (Da/N) f(y_j)`. Total residence time is independently optimized for each reactor architecture.

## Finite-stage coefficients

At the PFR optimum, the objective-normalized leading coefficient is `Phi = lim_(N -> infinity) N (Y_PFR* - Y_N*)/Y_PFR*`. A directly integrated second-order perturbation yields `Psi`; the optimized correction includes the change in the residence-time optimum. Both coefficients characterize large-N behavior, while exact stage requirements use direct finite-CSTR calculations.

## Reference values

| Inlet xylose (M) | Phi |
| ---: | ---: |
| 0.033 | 0.4889 |
| 0.176 | 0.5967 |
| 0.340 | 0.5554 |

The numerical coefficients are derived from the kinetic study associated with: W. Wang, [*Kinetic Criteria for Finite CSTR Staging toward Plug-Flow Performance*](https://doi.org/10.26434/chemrxiv.15007934/v1), ChemRxiv (2026), with the underlying xylose kinetics attributed to Krzelj et al., *Industrial & Engineering Chemistry Research* (2020), DOI: [10.1021/acs.iecr.0c01197](https://doi.org/10.1021/acs.iecr.0c01197). The program uses the reported effective rate coefficients; it does not refit the kinetic constants.

## Run

```bash
python -m pip install -e ".[dev,plot]"
python -m pytest -q
python examples/xylose_phi_sweep.py
python examples/decision_sensitivity.py
python examples/second_order_boundaries.py --plot
```

Additional reference outputs and numerical tolerances are listed in [NUMERICAL_VERIFICATION.md](NUMERICAL_VERIFICATION.md).
