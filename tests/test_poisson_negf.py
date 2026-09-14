"""Integration smoke tests for the self-consistent Poisson–NEGF loop.

These run the *real* rank-1 solver on the paper device at deliberately coarse
resolution (large grid spacing, coarse energy grid) so they finish in a few
seconds. They check that the machinery is wired correctly and physically sane —
not that the production sweep is converged (that is an open validation item in
``docs/POISSON_INTEGRATION.md``).
"""
import numpy as np
import pytest

from core.poisson_negf import (
    build_electrostatics,
    compute_density_prefactor,
    deep_contact_sites,
    run_self_consistent_bias,
)
from run.run_rank1_sweep import (
    build_phonon_modes,
    build_stack,
    emitter_barrier_center_nm,
    gaussian_chi,
)

DEVICE = "ZnO_MgZnO_symmetric"


@pytest.fixture(scope="module")
def coarse_setup():
    """Coarse-grid setup on the paper device (fast: Np~54, NE~75)."""
    a_nm = 0.5
    a_m = a_nm * 1e-9
    dE = 0.01
    E_grid = np.arange(-0.25, 0.5 + 0.5 * dE, dE)
    kT = 0.02585
    Ef = 0.02
    H_z, UB, t0, Np, z_nm, bounds = build_stack(DEVICE, a_m)
    z0 = emitter_barrier_center_nm(bounds)
    chi_mol = gaussian_chi(z_nm, z0, 0.3)
    D0_sq, hnu_idx, N_bose, chi_def, chi_list = build_phonon_modes(
        "Baseline", dE, kT, chi_mol, 1, 1, UB=UB)
    eps_r, N_D, cmask = build_electrostatics(DEVICE, a_m)
    return dict(a_m=a_m, dE=dE, E_grid=E_grid, kT=kT, Ef=Ef, H_z=H_z, UB=UB,
                t0=t0, Np=Np, chi_def=chi_def, D0_sq=D0_sq, hnu_idx=hnu_idx,
                N_bose=N_bose, chi_list=chi_list, eps_r=eps_r, N_D=N_D,
                cmask=cmask)


def test_electrostatics_grid_aligns_with_hamiltonian(coarse_setup):
    s = coarse_setup
    assert s["eps_r"].size == s["Np"]
    assert s["N_D"].size == s["Np"]
    # Contacts doped, active region undoped.
    assert s["N_D"][s["cmask"]].min() > 0
    assert np.all(s["N_D"][~s["cmask"]] == 0)
    # Both contact blocks present (emitter + collector).
    assert s["cmask"][0] and s["cmask"][-1]


def test_density_prefactor_is_positive_and_finite(coarse_setup):
    s = coarse_setup
    C = compute_density_prefactor(
        E_grid=s["E_grid"], H_z=s["H_z"], UB=s["UB"], t0=s["t0"], Ef=s["Ef"],
        kT=s["kT"], chi_diag=s["chi_def"], D0_sq_per_mode=s["D0_sq"],
        hnu_idx_per_mode=s["hnu_idx"], N_bose_per_mode=s["N_bose"],
        chi_per_mode=s["chi_list"], N_D=s["N_D"], contact_mask=s["cmask"],
        scba_max_iter=20, scba_mix=0.3, scba_tol=1e-4)
    assert np.isfinite(C) and C > 0


def test_equilibrium_self_consistent_is_neutral_and_flat(coarse_setup):
    """At V=0 the converged potential is ~flat and contacts are charge-neutral."""
    s = coarse_setup
    C = compute_density_prefactor(
        E_grid=s["E_grid"], H_z=s["H_z"], UB=s["UB"], t0=s["t0"], Ef=s["Ef"],
        kT=s["kT"], chi_diag=s["chi_def"], D0_sq_per_mode=s["D0_sq"],
        hnu_idx_per_mode=s["hnu_idx"], N_bose_per_mode=s["N_bose"],
        chi_per_mode=s["chi_list"], N_D=s["N_D"], contact_mask=s["cmask"],
        scba_max_iter=20, scba_mix=0.3, scba_tol=1e-4)
    sc = run_self_consistent_bias(
        V=0.0, E_grid=s["E_grid"], H_z=s["H_z"], UB=s["UB"], t0=s["t0"],
        Ef=s["Ef"], kT=s["kT"], chi_diag=s["chi_def"], D0_sq_per_mode=s["D0_sq"],
        hnu_idx_per_mode=s["hnu_idx"], N_bose_per_mode=s["N_bose"],
        chi_per_mode=s["chi_list"], eps_r=s["eps_r"], N_D=s["N_D"],
        a_m=s["a_m"], density_prefactor=C, contact_mask=s["cmask"],
        scba_max_iter=20, scba_mix=0.3,
        scba_tol=1e-4, poisson_max_iter=20, poisson_tol=1e-3)

    # Boundary conditions held (V=0 → both ends grounded).
    assert sc.U[0] == pytest.approx(0.0, abs=1e-9)
    assert sc.U[-1] == pytest.approx(0.0, abs=1e-9)
    # Potential bounded and physical (no runaway).
    assert np.all(np.isfinite(sc.U))
    assert np.abs(sc.U).max() < 0.5  # eV — equilibrium band bending is small
    # Density non-negative everywhere; contacts near charge neutrality.
    assert np.all(sc.n_e >= -1e-3 * s["N_D"].max())
    deep = deep_contact_sites(s["cmask"])
    n_contact = float(np.mean(sc.n_e[deep]))
    N_D_contact = float(np.mean(s["N_D"][deep]))
    assert n_contact == pytest.approx(N_D_contact, rel=0.5)


def test_finite_bias_no_spurious_contact_bending(coarse_setup):
    """Thesis-style reservoir BC (Akkala Fig 3.1/3.3): deep flat-band terminals
    carry semiclassical TF/Boltzmann charge → Thomas-Fermi screening pins them
    at the rail (small drop). Guards against the earlier ~1 eV spurious-bending
    regression *and* against the over-aggressive whole-contact flat-clamp that
    suppressed the emitter accumulation layer."""
    from core.poisson_negf import terminal_masks
    s = coarse_setup
    C = compute_density_prefactor(
        E_grid=s["E_grid"], H_z=s["H_z"], UB=s["UB"], t0=s["t0"], Ef=s["Ef"],
        kT=s["kT"], chi_diag=s["chi_def"], D0_sq_per_mode=s["D0_sq"],
        hnu_idx_per_mode=s["hnu_idx"], N_bose_per_mode=s["N_bose"],
        chi_per_mode=s["chi_list"], N_D=s["N_D"], contact_mask=s["cmask"],
        scba_max_iter=20, scba_mix=0.3, scba_tol=1e-4)
    V = 0.4
    sc = run_self_consistent_bias(
        V=V, E_grid=s["E_grid"], H_z=s["H_z"], UB=s["UB"], t0=s["t0"],
        Ef=s["Ef"], kT=s["kT"], chi_diag=s["chi_def"], D0_sq_per_mode=s["D0_sq"],
        hnu_idx_per_mode=s["hnu_idx"], N_bose_per_mode=s["N_bose"],
        chi_per_mode=s["chi_list"], eps_r=s["eps_r"], N_D=s["N_D"],
        a_m=s["a_m"], density_prefactor=C, contact_mask=s["cmask"],
        scba_max_iter=20, scba_mix=0.3, scba_tol=1e-4,
        poisson_max_iter=10, poisson_tol=5e-3)
    # Endpoints Dirichlet-exact.
    assert sc.U[0] == pytest.approx(+V / 2.0, abs=1e-9)
    assert sc.U[-1] == pytest.approx(-V / 2.0, abs=1e-9)
    # No ~1 eV runaway anywhere.
    assert np.abs(sc.U).max() < V / 2 + 0.20
    # Deep terminals stay near the rails (TF screening pins them — far from a
    # whole-V/2 drop or any runaway).
    em_term, col_term = terminal_masks(s["cmask"])
    if em_term.sum() > 0:
        assert np.max(np.abs(sc.U[em_term] - V / 2.0)) < 0.10
    if col_term.sum() > 0:
        assert np.max(np.abs(sc.U[col_term] - (-V / 2.0))) < 0.10
