# cstr-convergence

A compact, reproducible chemical reaction engineering tool for quantifying how a series of equal-volume continuous stirred-tank reactors approaches ideal plug-flow behavior at the same total residence time.

## Models

- Single irreversible reaction with power-law kinetics: `-r_A = k C_A^n` (order n >= 0).
- First-order A -> B -> C network, with individual rate constants.
- The first-order reaction has an exact CSTR formula and an exact PFR reference.
- Calculates the minimum integer stage count for a given **absolute conversion error** in the single first-order reaction.

For the first-order case, at `Da = k * tau`:

`X_PFR = 1 - exp(-Da)`

`X_N = 1 - (1 + Da/N)^(-N)`

`error = abs(X_N - X_PFR)`

`N_min = min {N >= 1 : error <= tolerance}`

## Installation

```bash
python -m pip install -e ".[dev]"
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
5. The generic `minimum_stages` binary search assumes the chosen error decreases monotonically with N. For an intermediate-species yield in a network, monotonicity must be verified before this search can be used.
6. These functions demonstrate established reactor design relationships. No new mathematical convergence theorem is claimed.

## Roadmap

- Dimensionless criterion maps for multiple kinetic orders.
- Explicit error definition for desired intermediate product yield and full composition vectors.
- Optional nonisothermal rate constants and experimental RTD coupling.
