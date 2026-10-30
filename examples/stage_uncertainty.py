"""Contrast parameter uncertainty with discrete first-order reactor stage decisions.

The two cases use the same absolute conversion tolerance and a uniform,
illustrative range of Damkohler numbers. Changing Da represents uncertainty
in k at fixed feed concentration and total residence time.
"""

from collections import Counter

import numpy as np

from cstr_convergence import first_order_minimum_stages


def evaluate_case(label: str, center_da: float, fractional_range: float) -> None:
    tolerance = 0.01
    da_samples = np.linspace(
        center_da * (1.0 - fractional_range),
        center_da * (1.0 + fractional_range),
        201,
    )
    stages = [first_order_minimum_stages(float(da), tolerance) for da in da_samples]
    counts = Counter(stages)
    print(f"{label}: Da={center_da:g} +/- {fractional_range:.0%}")
    print(f"  N_min range: {min(stages)} to {max(stages)}")
    print(f"  Sampled stage counts: {dict(sorted(counts.items()))}")
    print(f"  Central estimate: {first_order_minimum_stages(center_da, tolerance)}")
    print()


def main() -> None:
    print("First-order A -> products; absolute conversion tolerance = 0.01")
    evaluate_case("Broad kinetic variation, stable stage decision", 2.0, 0.10)
    evaluate_case("Small kinetic variation, stage boundary", 2.39, 0.01)


if __name__ == "__main__":
    main()
