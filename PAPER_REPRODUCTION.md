# Reproducing finite-CSTR staging kinetics

**Reference:** Wenjia Wang, [*Kinetic Criteria for Finite CSTR Staging toward Plug-Flow Performance*](https://doi.org/10.26434/chemrxiv.15007934/v1), ChemRxiv (2026).

## Kinetic model and parameter provenance

The four-rate xylose reaction network and numerical coefficients come from the manuscript's Supporting Information, Sections S7–S9 and Tables S8–S13. This implementation uses the **effective** coefficients reported for 160 °C and model hydronium concentration 0.10 M (Table S9). No refitting of Arrhenius parameters is performed.

The reaction network is `X -> I`, `I -> F`, `F -> degradation`, and `I + F -> degradation`. Concentrations `X`, `I`, and `F` are normalized by inlet xylose concentration `C_X0`. Define `Da = k_X * t`, `alpha = k_I/k_X`, `beta = k_F/k_X`, and `gamma = k_FI * C_X0/k_X`; the optimized objective is the furfural yield `F`.

The PFR satisfies `dz/du = f(z)` and each of `N` equal CSTR stages satisfies `z_j = z_(j-1) + (Da/N) f(z_j)`. For the leading finite-staging correction, the first-order state perturbation obeys `e' = J e + (Da/2) J f`, with `e(0) = 0`, evaluated over the normalized residence coordinate. The manuscript's objective-normalized coefficient is `Phi = -e_F(u*)/F_PFR(u*)`, where `u*` maximizes the PFR furfural yield. The PFR residence-time optimum is stationary, so its first-order shift does not contribute to Φ.

## Reproduce

```bash
python -m pip install -e ".[dev,plot]"
python -m pytest -q
python examples/paper_reproduction.py
python examples/decision_sensitivity.py
```

The concentration sweep writes `results/xylose_phi.csv` and `results/xylose_phi.png`. Stage-count calculations optimize the finite-CSTR cascade independently for each integer stage number.

## Verification targets

Reference Φ values from Supporting Information Tables S10–S11:

| Inlet xylose, M | Φ |
| ---: | ---: |
| 0.033 | 0.4889 |
| 0.176 | 0.5967 |
| 0.340 | 0.5554 |

The sensitivity example uses `C_X0 = 0.0769 M`, a concentration included in Table S10. The ±30% perturbations of the `I + F` cross-degradation coefficient are illustrative parameter changes. Tests also check the matched-rate consecutive-reaction benchmark (`Phi = 0.5`, `N95 = 10`) and direct finite-CSTR predictions.

The archived supporting-information tables provide the numerical regression targets used here. The linked ChemRxiv record identifies the associated research publication.
