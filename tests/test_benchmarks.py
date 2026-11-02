import pytest
from cstr_convergence import xylose_benchmark, XyloseModel

@pytest.mark.parametrize('c,phi,n95', [
    (0.033, 0.488890744, 9),
    (0.176, 0.596736335, 11),
    (0.340, 0.555395636, 10),
])
def test_optimized_phi_and_exact_n95(c, phi, n95):
    m = XyloseModel()
    assert m.phi(c) == pytest.approx(phi, abs=2e-6)
    assert m.stage_requirement(c, 0.95, 20) == n95

def test_xylose_benchmark_physical_scaling():
    row = xylose_benchmark(concentrations=(0.033,))[0]
    assert row['Da_PFR_star'] == pytest.approx(2.1998, abs=2e-4)
    assert row['tau_PFR_star_min'] == pytest.approx(141.38, abs=0.02)
    assert row['Y_F_PFR_star'] == pytest.approx(0.6118, abs=2e-4)
