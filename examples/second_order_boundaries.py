"""Compute Phi, direct Psi, stage estimates, and exact design boundaries.

Use --plot to generate the two publication-scale figures in results/.
"""
import argparse
import csv
from pathlib import Path
import numpy as np
from cstr_convergence import (
    XyloseModel, optimized_staging_coefficients,
    second_order_stage_estimate, find_concentration_boundary,
    find_parameter_boundary, retention_margin,
)

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / 'results'


def save_csv(filename, rows):
    with (RESULTS / filename).open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main(make_plots=False):
    RESULTS.mkdir(exist_ok=True)
    model = XyloseModel()
    rows = []
    # Compare large-stage coefficients against exact optimized reactor performance.
    for c in (0.033, 0.176, 0.340):
        coeff = optimized_staging_coefficients(model, c)
        exact_n95 = model.stage_requirement(c, .95, 15)
        n95_second = second_order_stage_estimate(coeff, .95)
        rows.append(dict(C_X0_M=c, Phi=coeff.phi, Psi=coeff.psi,
                         Da_PFR_star=coeff.da_star,
                         Y_F_PFR_star=coeff.y_pfr_star,
                         N95_exact=exact_n95, N95_second_order=n95_second))
        print(f'C={c:.3f} M; Phi={coeff.phi:.6f}; Psi={coeff.psi:.6f}; '
              f'N95 exact={exact_n95}; second-order={n95_second}')
    save_csv('staging_coefficients.csv', rows)

    # Find both concentration roots of the exact 10-stage retention margin.
    left = find_concentration_boundary(model, 10, (.06, .09))
    right = find_concentration_boundary(model, 10, (.28, .33))
    gamma_factor = find_parameter_boundary(model, .0769, 10, (.9, 1.0))
    boundaries = [
        dict(boundary='lower concentration', value=left, units='M'),
        dict(boundary='upper concentration', value=right, units='M'),
        dict(boundary='gamma multiplier', value=gamma_factor, units='1'),
    ]
    save_csv('configuration_boundaries.csv', boundaries)
    print(f'10-stage concentration boundaries: {left:.8f}, {right:.8f} M')
    print(f'At 0.0769 M, gamma multiplier at 10-stage boundary: {gamma_factor:.8f}')

    if make_plots:
        import matplotlib.pyplot as plt
        plt.rcParams.update({
            'font.size': 15,
            'axes.labelsize': 18,
            'axes.labelweight': 'bold',
            'axes.linewidth': 2.3,
            'xtick.labelsize': 14,
            'ytick.labelsize': 14,
        })
        b = '#800020'
        blue = '#002FA7'
        green = '#00754A'
        concentrations = np.linspace(.033, .34, 17)
        margin = [retention_margin(model, float(c), 10, .95)
                  for c in concentrations]
        fig, ax = plt.subplots(figsize=(10, 8), facecolor='white')
        ax.plot(concentrations, margin, color=b, linewidth=3.1)
        ax.axhline(0, color=blue, linewidth=2.2, linestyle='--')
        for c in (left, right):
            ax.plot(c, 0, 'o', color=green, markersize=11)
        ax.set_xlabel('Inlet xylose concentration (mol/L)')
        ax.set_ylabel(r'$g_{10}=Y^*_{10}-0.95Y^*_{PFR}$')
        ax.tick_params(which='both', width=2, length=6)
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight('bold')
        ax.grid(False)
        for s in ax.spines.values(): s.set_visible(True); s.set_linewidth(2.3)
        fig.tight_layout()
        fig.savefig(RESULTS / 'configuration_boundaries.png', dpi=190, facecolor='white')
        plt.close(fig)

        cc = np.linspace(.033, .34, 15)
        coefficients = [optimized_staging_coefficients(model, float(c)) for c in cc]
        fig, ax = plt.subplots(figsize=(10, 8), facecolor='white')
        ax.plot(cc, [r.phi for r in coefficients], color=b, linewidth=3,
                label=r'$\Phi$')
        ax.plot(cc, [r.psi for r in coefficients], color=blue, linewidth=3,
                label=r'$\Psi$')
        ax.set_xlabel('Inlet xylose concentration (mol/L)')
        ax.set_ylabel('Optimized finite-stage coefficients')
        ax.tick_params(which='both', width=2, length=6)
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight('bold')
        ax.legend(frameon=False)
        ax.grid(False)
        for s in ax.spines.values(): s.set_visible(True); s.set_linewidth(2.3)
        fig.tight_layout()
        fig.savefig(RESULTS / 'staging_coefficients.png', dpi=190, facecolor='white')
        plt.close(fig)
        print('Saved coefficient and boundary figures to results/')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plot', action='store_true', help='generate coefficient and boundary figures')
    args = parser.parse_args()
    main(args.plot)
