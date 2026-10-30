# CSTR Convergence

**Decision-oriented kinetic criteria for finite reactor staging and process-design decisions.**

This Python package connects reaction-network kinetics to the number of equal-volume continuous stirred-tank reactors (CSTRs) needed to retain a specified fraction of plug-flow reactor (PFR) performance. It implements the finite-staging coefficient **Φ**, reproduces nonlinear xylose-dehydration results from a research manuscript, and illustrates when kinetic changes alter a discrete reactor-design choice.

**Research reference:** W. Wang, [*Kinetic Criteria for Finite CSTR Staging toward Plug-Flow Performance*](https://doi.org/10.26434/chemrxiv.15007934/v1), ChemRxiv (2026).

## Main result: concentration-dependent staging sensitivity

![Optimized finite-staging coefficient Phi versus inlet xylose concentration](results/xylose_phi.png)

At the published kinetic condition (160 °C, model hydronium concentration 0.10 M), the four-pathway xylose network gives an **interior maximum** in the optimized, objective-normalized coefficient Φ.

| Inlet xylose (M) | Φ |
| ---: | ---: |
| 0.033 | 0.4889 |
| 0.176 | **0.5967** |
| 0.340 | 0.5554 |

Φ connects reaction-network structure and the chosen performance objective to the leading fractional loss from finite staging:

\[
\frac{J_{\mathrm{PFR}}^*-J_N^*}{J_{\mathrm{PFR}}^*}
=\frac{\Phi}{N}+O(N^{-2}).
\]

Here each architecture independently optimizes total residence time for the target furfural yield. The exact integer stage requirement is obtained by optimizing the finite-CSTR cascade, rather than rounding the leading asymptotic estimate.

## Decision example: kinetic sensitivity versus stage-count sensitivity

At **0.0769 M inlet xylose**, vary the cross-degradation coefficient for `I + F -> degradation` by ±30%, keeping all other model coefficients fixed:

| Cross-degradation coefficient | Φ | Minimum stages for 95% optimized PFR yield |
| --- | ---: | ---: |
| 0.70 × baseline | 0.5287 | **10** |
| 1.00 × baseline | 0.5577 | **11** |
| 1.30 × baseline | 0.5764 | **11** |

The continuous sensitivity of Φ does not translate into a stage-count change in both directions. This is a **controlled kinetic-perturbation demonstration**; the ±30% range is chosen for illustration.

```bash
python examples/decision_sensitivity.py
```

This example provides a computational basis for examining which kinetic uncertainties are consequential to a reactor design decision. The 0.0769 M point appears in the manuscript Supporting Information concentration sweep (Table S10).

## Install and reproduce

```bash
python -m pip install -e ".[dev,plot]"
python -m pytest -q
python examples/paper_reproduction.py
python examples/decision_sensitivity.py
```

`paper_reproduction.py` writes `results/xylose_phi.csv` and `results/xylose_phi.png`. Detailed model equations, kinetic parameters, and verification targets are in [PAPER_REPRODUCTION.md](PAPER_REPRODUCTION.md).

## Available models and numerical methods

- **Nonlinear xylose dehydration:** `X -> I -> F`, first-order furfural degradation, and bimolecular `I + F` degradation; PFR, finite CSTR, optimized furfural yield, Φ, and local kinetic sensitivity.
- **Consecutive first-order reactions:** `A -> B -> C`, including analytical PFR and finite-CSTR benchmarks, optimized yield, asymptotic coefficients, and exact stage-count calculations.
- **Single power-law reaction:** PFR and finite-CSTR conversion for order `n >= 0`, including zero-order depletion.
- **Error and stage-count criteria:** absolute conversion error; intermediate-yield and composition errors; first-passing and bounded persistent stage scans.

The mathematical connection is the backward-Euler equivalence of an equal-CSTR cascade and the PFR reaction ODE. The reaction network and performance objective determine the finite-stage error coefficient; equal-CSTR hydrodynamics supply the `1/N` scale.

## Validation

```bash
python -m pytest -q
```

Tests cover analytical first-order solutions, zero- and fractional-order kinetics, CSTR mass balances, exact versus recursive consecutive-reaction solutions, extreme Damköhler-number numerical stability, independent PFR ODE comparisons, Φ values from manuscript Supporting Information Tables S10–S11, optimized consecutive-reaction benchmarks, and the discrete decision example. GitHub Actions runs the test suite on pushes and pull requests.

## Scope

Models assume steady isothermal operation, constant volumetric flow, and equal-volume ideal CSTR stages at fixed total residence time for each comparison. `minimum_stages_scan(..., criterion="persistent")` checks persistence through the supplied `max_stages`; its finite search interval is specified in the function documentation. Asymptotic coefficients describe the large-`N` limit and are complemented by explicit finite-stage calculations.

## Further development

- Extend the objective-specific Φ interface to additional nonlinear chemical reaction networks.
- Propagate kinetic and chemical-state uncertainty into stage-selection boundaries.
- Connect reactor-performance objectives to downstream separation recovery and process economics.
- Integrate nonisothermal reaction histories and experimentally characterized residence-time distributions.

See [CHANGELOG.md](CHANGELOG.md) for version history and [LICENSE](LICENSE) for the MIT license.
