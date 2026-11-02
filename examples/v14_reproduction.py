"""Reproduce V14 optimized xylose staging metrics and write a CSV table."""
import csv
from pathlib import Path
from cstr_convergence import v14_benchmark

def main():
    rows = v14_benchmark()
    output = Path(__file__).resolve().parents[1] / 'results' / 'v14_metrics.csv'
    output.parent.mkdir(exist_ok=True)
    with output.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(f"C={row['C_X0_M']:.4f} M; Phi={row['Phi']:.6f}; "
              f"N95={row['N95_exact']}; Da*={row['Da_PFR_star']:.6f}")
    print(f'Saved {output}')

if __name__ == '__main__':
    main()
