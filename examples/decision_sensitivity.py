"""Demonstrate how kinetic sensitivity affects an integer reactor-stage decision.

Reproduces the baseline xylose model at C_X0 = 0.0769 M and varies
only the I + F cross-degradation coefficient by a specified factor.
"""

from dataclasses import replace

from cstr_convergence import XyloseModel


def main():
    baseline = XyloseModel()
    inlet_xylose_molar = 0.0769
    retention_target = 0.95
    max_stages = 60

    print("Kinetic sensitivity and minimum CSTR stage requirement")
    print(f"Inlet xylose: {inlet_xylose_molar:.4f} M")
    print(f"Optimized PFR yield retention: {retention_target:.0%}")
    print()

    for factor in (0.70, 1.00, 1.30):
        # Hold all other kinetic parameters fixed while perturbing I + F loss.
        model = replace(
            baseline,
            gamma_per_molar=baseline.gamma_per_molar * factor,
        )

        phi = model.phi(inlet_xylose_molar)
        minimum_n95 = model.stage_requirement(
            inlet_xylose_molar,
            retention_target,
            max_stages,
        )

        print(
            f"gamma factor={factor:.2f}, "
            f"Phi={phi:.4f}, "
            f"N95={minimum_n95}"
        )


if __name__ == "__main__":
    main()
