"""The chunked low-memory SCBA must reproduce the fast solver exactly."""
import numpy as np
import pytest

from core.scba_rank1_keldysh import (run_rank1_keldysh_single_bias_fast,
                                     run_rank1_keldysh_single_bias_lowmem)
from run.run_rank1_sweep import (build_phonon_modes, build_stack, emitter_barrier_center_nm,
                                 gaussian_chi, linear_bias_profile)
from core.validation import current_diagnostics


def _molecule_inputs(a_m=.5e-9, dE=.004, V=.4):
    H, UB, t0, Np, z, b = build_stack("ZnO_MgZnO_symmetric", a_m)
    chi = gaussian_chi(z, emitter_barrier_center_nm(b), .3)
    D, hi, nb, ch, chl = build_phonon_modes("Mol_A", dE, .02585, chi, 1, 1, UB=UB)
    return dict(V=V, E_grid=np.arange(-.25, .5+dE/2, dE), H_z=H, UB=UB,
                bias_profile=linear_bias_profile(V, Np, 1, 1), t0=t0, Ef=.02, kT=.02585,
                chi_diag=ch, D0_sq_per_mode=D, hnu_idx_per_mode=hi,
                N_bose_per_mode=nb, chi_per_mode=chl, max_iter=500, tol=1e-7, mix=.3)


def test_anderson_rejects_the_observed_negative_rate_branch():
    # Previously depth=8 returned converged=True with min Sigma_out=-0.865 eV.
    kw = _molecule_inputs(.2e-9, .0005, .56)
    linear = run_rank1_keldysh_single_bias_lowmem(method="rgf", **kw)
    guarded = run_rank1_keldysh_single_bias_lowmem(method="rgf", anderson_depth=8, **kw)
    assert guarded.mixer_fallbacks > 0
    for r in (linear, guarded):
        assert r.converged and r.physical_valid and r.self_energy_residual < kw["tol"]
        assert r.sigma_in_ph.real.min() >= 0 and r.sigma_out_ph.real.min() >= 0
        assert r.n_diag.min() >= 0 and r.p_diag.min() >= 0
    np.testing.assert_allclose(guarded.I_right, linear.I_right, rtol=2e-6)


def test_negative_warm_start_cannot_enter_scba():
    kw = _molecule_inputs()
    sigma = np.zeros((kw["E_grid"].size, kw["H_z"].shape[0]), complex)
    sigma[0, 0] = -.1
    with pytest.raises(ValueError, match="nonnegative rates"):
        run_rank1_keldysh_single_bias_lowmem(method="rgf", sigma_init=(sigma, sigma), **kw)


@pytest.mark.parametrize("method", ["rgf", "dense"])
def test_current_check_continues_beyond_self_energy_tolerance(method):
    kw = _molecule_inputs()
    kw["tol"] = .01
    loose = run_rank1_keldysh_single_bias_lowmem(method=method, **kw)
    tight = run_rank1_keldysh_single_bias_lowmem(method=method, current_rtol=2e-3, **kw)
    assert loose.converged
    assert not current_diagnostics(loose.I_left, loose.I_right, True, rtol=2e-3)["current_valid"]
    assert tight.iters_used > loose.iters_used
    assert current_diagnostics(tight.I_left, tight.I_right, tight.converged, rtol=2e-3)["current_valid"]


def test_lowmem_matches_fast():
    a_m, dE, kT = 0.5e-9, 0.004, 0.02585
    H, UB, t0, Np, z, b = build_stack("ZnO_MgZnO_symmetric", a_m)
    chi = gaussian_chi(z, emitter_barrier_center_nm(b), 0.3)
    D0, hi, Nb, chd, chl = build_phonon_modes("Mol_A", dE, kT, chi, 1, 1, UB=UB)
    E = np.arange(-0.1, 0.5, dE)
    kw = dict(V=0.4, E_grid=E, H_z=H, UB=UB, bias_profile=linear_bias_profile(0.4, Np, 1, 1),
              t0=t0, Ef=0.02, kT=kT, chi_diag=chd, D0_sq_per_mode=D0, hnu_idx_per_mode=hi,
              N_bose_per_mode=Nb, chi_per_mode=chl, max_iter=300, tol=1e-6, mix=0.3)
    a = run_rank1_keldysh_single_bias_fast(**kw)
    lo = run_rank1_keldysh_single_bias_lowmem(chunk=17, **kw)   # chunk not dividing NE
    assert lo.iters_used == a.iters_used and lo.converged == a.converged
    n_fast = np.real(np.einsum("eii->ei", a.G_lesser))
    assert np.abs(lo.n_diag - n_fast).max() <= 1e-12 * np.abs(n_fast).max()
    assert abs(lo.I_right / a.I_right - 1) < 1e-12
    assert abs(lo.I_left / a.I_left - 1) < 1e-12


def test_lowmem_warm_start_reaches_same_answer():
    a_m, dE, kT = 0.5e-9, 0.004, 0.02585
    H, UB, t0, Np, z, b = build_stack("ZnO_MgZnO_symmetric", a_m)
    chi = gaussian_chi(z, emitter_barrier_center_nm(b), 0.3)
    D0, hi, Nb, chd, chl = build_phonon_modes("Mol_A", dE, kT, chi, 1, 1, UB=UB)
    kw = dict(V=0.4, E_grid=np.arange(-0.1, 0.5, dE), H_z=H, UB=UB,
              bias_profile=linear_bias_profile(0.4, Np, 1, 1), t0=t0, Ef=0.02, kT=kT,
              chi_diag=chd, D0_sq_per_mode=D0, hnu_idx_per_mode=hi, N_bose_per_mode=Nb,
              chi_per_mode=chl, max_iter=300, tol=1e-6, mix=0.3)
    cold = run_rank1_keldysh_single_bias_lowmem(**kw)
    warm = run_rank1_keldysh_single_bias_lowmem(
        sigma_init=(cold.sigma_in_ph, cold.sigma_out_ph), **kw)
    assert cold.converged and warm.converged
    assert warm.iters_used <= 2 < cold.iters_used
    assert abs(warm.I_right / cold.I_right - 1) < 1e-6


def test_rgf_matches_dense():
    """Recursive Green's functions give the dense solver's diagonals."""
    a_m, dE, kT = 0.5e-9, 0.004, 0.02585
    H, UB, t0, Np, z, b = build_stack("ZnO_MgZnO_symmetric", a_m)
    chi = gaussian_chi(z, emitter_barrier_center_nm(b), 0.3)
    D0, hi, Nb, chd, chl = build_phonon_modes("Mol_A", dE, kT, chi, 1, 1, UB=UB)
    kw = dict(V=0.4, E_grid=np.arange(-0.1, 0.5, dE), H_z=H, UB=UB,
              bias_profile=linear_bias_profile(0.4, Np, 1, 1), t0=t0, Ef=0.02, kT=kT,
              chi_diag=chd, D0_sq_per_mode=D0, hnu_idx_per_mode=hi, N_bose_per_mode=Nb,
              chi_per_mode=chl, max_iter=300, tol=1e-6, mix=0.3, anderson_depth=0)
    d = run_rank1_keldysh_single_bias_lowmem(method="dense", **kw)
    r = run_rank1_keldysh_single_bias_lowmem(method="rgf", **kw)
    assert r.iters_used == d.iters_used
    assert np.abs(r.n_diag - d.n_diag).max() <= 1e-10 * np.abs(d.n_diag).max()
    assert abs(r.I_right / d.I_right - 1) < 1e-10
    assert abs(r.I_left / d.I_left - 1) < 1e-10


def test_mixer_accepts_roundoff_and_rejects_real_negatives():
    from core.scba_rank1_keldysh import _ScbaMixer
    def proposal(target):
        # residuals +1 and -1 give DIIS weights (1/2, 1/2); with mix = 1 the
        # proposal is the mean of the two outputs, placed exactly at `target`
        m = _ScbaMixer(mix=1., depth=8, max_iter=10, tol=1e-6)
        f0 = np.ones((1, 3), complex); f1 = 2*target - f0
        m.update(f0 - 1, f0 - 1, f0, f0)
        return m, m.update(f1 + 1, f1 + 1, f1, f1)[0]
    m, si = proposal(-1e-16)                 # roundoff-level negative: accepted
    assert m.fallbacks == 0 and -1e-15 < si.real.min() < 0
    m, si = proposal(-1e-3)                  # physically negative: rejected
    assert m.fallbacks == 1 and not np.allclose(si, -1e-3)


def test_legacy_reproduction_cannot_require_current_conservation():
    kw = _molecule_inputs()
    from core.scba_rank1_keldysh import run_rank1_keldysh_single_bias
    with pytest.raises(ValueError, match="legacy_reproduction"):
        run_rank1_keldysh_single_bias(legacy_reproduction=True, current_rtol=2e-3, **kw)
