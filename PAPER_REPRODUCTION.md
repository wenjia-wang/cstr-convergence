# Research reproduction: finite-CSTR staging kinetics

Source: Wenjia Wang, *Kinetic Criteria for Finite CSTR Staging toward Plug-Flow Performance*, ChemRxiv (2026), DOI [10.26434/chemrxiv.15007934/v1](https://doi.org/10.26434/chemrxiv.15007934/v1).

## Scope and provenance

The normalized four-rate xylose kinetic system is transcribed from the manuscript Supporting Information, Sections S7–S9, Tables S8–S13. This v0.4 implementation uses **effective** 160 °C, model hydronium 0.10 M coefficients **as reported in Table S9**, with no Arrhenius refitting. It reproduces the key interior maximum in the **optimized, objective-normalized** finite-stage coefficient Phi. This is distinct from an absolute fixed-residence-time error coefficient.

Network: X -> I (k_X X); I -> F (k_I I); F -> degradation (k_F F); I+F -> degradation (k_FI I F). X,I,F are normalized by inlet xylose C_X0, with Da=k_X*t, alpha=k_I/k_X, beta=k_F/k_X, gamma=k_FI*C_X0/k_X. F is the optimized target yield.

The governing balance is dz/du=f(z). Each CSTR satisfies z_j=z_(j-1)+(Da/N)f(z_j). The perturbation e to the PFR state solves e'=J e + (Da/2) J f with e(0)=0. The manuscript's normalized coefficient is Phi=-e_F(u*)/F_PFR(u*), where u* independently maximizes PFR F. The first-order residence-time optimum shift does not enter Phi because the PFR objective is stationary.

## Run

```bash
python -m pip install -e '.[dev,plot]'
python -m pytest -q
python examples/paper_reproduction.py
python examples/decision_sensitivity.py
```

The examples are intended as a computational reproduction and a **new exploratory decision example**. Sensitivity values from the gamma perturbation are numerical derivatives, not independently published or measured experimental estimates. Explicit 95% stage counts use independent optimization of each finite cascade. Stage-count crossings should be checked under uncertainty rather than inferred from Phi alone.

## Verification targets

Supporting Information Table S10 / S11: at C_X0=0.033, 0.176 and 0.340 M, Phi approximately 0.4889, 0.5967 and 0.5554, respectively. Tests enforce tolerance to account for rounded Table S9 kinetic coefficients. Table S1: matched-rate consecutive network Phi=0.5 and N95=10. For publication-grade figure-by-figure replication, compare against the accepted manuscript and current SI version; the archived SI used here was V7 and the main manuscript was a later V11 draft.
