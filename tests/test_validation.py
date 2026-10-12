import numpy as np

from core.poisson import bulk_electron_density, fermi_level_from_density
import pytest
from core.validation import (current_diagnostics, central_second_difference,
                             refinement_agrees, require_plot_data)


def test_conservation_is_required_even_when_solver_converges():
    d = current_diagnostics([-1e-9, -1e-9, 0., np.nan], [1.02e-9, 1e-9, 1e-16, 1.],
                            [True, False, True, True])
    np.testing.assert_array_equal(d["current_valid"], [False, False, True, False])


def test_fixed_doping_changes_fermi_level_across_temperature():
    mass, density = .28*9.1093837015e-31, 5.03636495047155e24
    levels = []
    for temperature in (10, 77, 150, 300):
        kt = .02585*temperature/300
        ef = fermi_level_from_density(density, kt, mass)
        levels.append(ef)
        np.testing.assert_allclose(bulk_electron_density(ef, kt, mass), density, rtol=1e-9)
    np.testing.assert_allclose(levels[-1], .02, atol=1e-8)
    assert levels[0] > .038
    assert np.all(np.diff(levels) < 0)


def test_spectral_refinement_detects_error_hidden_in_current():
    fine = np.array([1., 2., 3.]) * 1e-9
    coarse = fine + [0, 1e-12, 0]  # only 0.05% current error, but spurious curvature
    assert np.all(refinement_agrees(coarse, fine, rtol=.01, atol=1e-14))
    assert not refinement_agrees(central_second_difference(coarse, .016),
                                central_second_difference(fine, .016), .05, 1e-10)


def test_plot_gate_rejects_uncomputed_device_current_and_unvalidated_spectrum():
    d = dict(V=np.array([0., .1, .2]), I_L=np.array([0., -1e-9, -2e-9]),
             I_R=np.array([0., 1e-9, 2e-9]), I_device=np.zeros(3),
             converged=np.ones(3, dtype=bool), transverse=False)
    with pytest.raises(ValueError, match="not calculated"):
        require_plot_data(d, "I_device", allow_unvalidated=True)
    with pytest.raises(ValueError, match="grid validation"):
        require_plot_data(d, "I_R")
    require_plot_data(d, "I_R", allow_unvalidated=True)


def test_sweep_saves_failure_diagnostics_and_preserves_existing_results(tmp_path):
    from run.run_rank1_sweep import run_sweep
    kw = dict(molecule="Baseline", device="ZnO_MgZnO_symmetric", V_points=3,
              V_max=.032, a_nm=.5, dE=.01, T_K=77, donor_density=5.03636495047155e24,
              scba_max_iter=1, fast_solver=True, verbose=False, out_dir=str(tmp_path))
    with pytest.raises(RuntimeError, match="failed solver"):
        run_sweep(**kw)
    path = next(tmp_path.glob("*.npz"))
    before = path.read_bytes()
    with np.load(path) as d:
        assert str(d["temperature_policy"]) == "fixed_doping"
        assert float(d["Ef"]) > .037
        assert not np.all(d["numerical_valid"])
        assert not bool(d["spectrum_grid_validated"])
        assert np.all(np.isnan(d["I_device"]))
    with pytest.raises(FileExistsError):
        run_sweep(**kw)
    assert path.read_bytes() == before
