# cstr-convergence

A compact, reproducible chemical reaction engineering tool for quantifying how a series of equal-volume continuous stirred-tank reactors approaches ideal plug-flow behavior at the same total residence time.

## Models

- Single irreversible reaction with power-law kinetics: `-r_A = k C_A^n` (order n >= 0).
- First-order A -> B -> C network, with individual rate constants.
- The first-order reaction has an exact CSTR formula and an exact PFR reference.
- Calculates the minimum integer stage count for a given **absolute conversion error** for first-order and power-law reactions, and finite-range stage scans for series-reaction intermediate yields.

For the first-order case, at `Da = k * tau`:

`X_PFR = 1 - exp(-Da)`

`X_N = 1 - (1 + Da/N)^(-N)`

`error = abs(X_N - X_PFR)`

`N_min = min {N >= 1 : error <= tolerance}`

## Installation

```bash
python -m pip install -e ".[dev,plot]"
```

## Example

```python
from cstr_convergence import first_order_minimum_stages, series_cstr

print(first_order_minimum_stages(da=2.0, tolerance=0.01))
print(series_cstr(da1=2.0, da2=1.0, n=10))
```

Run a graphical example:

```bash
python examples/demo.py
```

It creates `outputs/convergence.png`.

## Verification

```bash
python -m pytest -q
```

Tests check the exact first-order solution, PFR convergence, conservation for A -> B -> C, and minimum-stage search.

## Assumptions and scope

1. Constant fluid density and volumetric flow; steady state; isothermal reactors.
2. Identical perfectly mixed tanks with equal residence time `tau/N`. Total residence time is held fixed as N varies.
3. Power-law reaction with concentration normalized by its inlet value: `Da = k C_A0^(n-1) tau`. For a zero-order model, complete depletion is treated as a physical concentration floor.
4. For A -> B -> C, both steps are first order and dimensionless `Da1 = k1*tau`, `Da2 = k2*tau`.
5. `minimum_stages` uses binary search and requires a nonincreasing error. Do not use it for intermediate-species yield errors. `minimum_stages_scan` exhaustively evaluates integer stages up to `max_stages` and provides `first` or `persistent` criteria. The latter guarantees the criterion **within the inspected finite range only**, not for all larger N.
6. These functions demonstrate established reactor design relationships. No new mathematical convergence theorem is claimed.

## Additional stage-count functions

```python
from cstr_convergence import nth_order_minimum_stages, series_minimum_stages
print(nth_order_minimum_stages(da=2, order=2, tolerance=.01))
print(series_minimum_stages(da1=6, da2=.4, tolerance=5e-5,
                            max_stages=200, criterion="persistent", metric="B"))
```

For intermediate yields, `criterion="first"` returns the first passing integer stage. `criterion="persistent"` returns the earliest stage at which all remaining **checked** stages through the inclusive `max_stages` pass. No claim is made for stages beyond this cap. `metric="vector"` uses the largest absolute error among A, B and C. These finite-range searches may be computationally expensive.

The exact first-order result is also a model-verification benchmark. As N becomes large, the leading absolute conversion error is Da² exp(-Da)/(2N); this is an asymptotic approximation and should not be used as an exact stage-count criterion.

## Roadmap

- Dimensionless criterion maps for multiple kinetic orders.
- Explicit error definition for desired intermediate product yield and full composition vectors.
- Optional nonisothermal rate constants and experimental RTD coupling.


### v0.2.1 numerical correction

The PFR intermediate-product expression evaluates the exponential difference using
`expm1` and factors out `exp(-min(Da1, Da2))`, avoiding overflow/underflow
for strongly separated nonnegative Damköhler numbers. Regression tests include
Da values up to 1000 and independent ODE integration for moderate values.
The power-law stage count remains an exhaustive bounded scan; acceleration
by binary search requires a monotonicity justification.


## Version 0.3.0: asymptotic analysis and exact series evaluation

- `series_cstr_exact(da1, da2, n)` evaluates the closed-form A->B->C cascade independently of the recursive solver.
- `logarithmic_mean(da1, da2)` and `series_asymptotic_sign` classify the **large-N** B-yield error using L(da1,da2) compared with 2. This is **not** a finite-N monotonicity theorem.
- `series_asymptotic_error` gives the leading **vector** coefficient of `(y_N-y_PFR)`, evaluated at dimensionless outlet time 1.
- `power_law_asymptotic_coefficient` returns the first term of `X_PFR-X_N` where the PFR exit concentration is positive; returns NaN for finite-time depletion.
- `first_order_asymptotic_diagnostic` reports an estimate and relative second-order correction. The estimate is **not** an exact minimum stage count nor a validated error bound; use `first_order_minimum_stages` for the actual design choice.
- The `persistent` criterion in `series_minimum_stages` remains certified only over the enumerated finite range 1..max_stages.

Run `python examples/convergence_v03.py` after installing the optional plot dependencies to reproduce B-yield error crossing zero.

## Research connection

This repository provides reproducible isothermal staging examples relevant to dimensionless finite-stage performance analysis. It does not yet implement the manuscript-specific generalized coefficient Phi, composition-dependent kinetics, uncertainty propagation, or decision-regret optimization. These remain separate proposed extensions.
