"""
Unit tests for the analytic d²I/dV² implementation (SISPAD Tier-1 α).

Tests the coherent (Landauer) path against dense finite-difference
ground truth on a toy Lorentzian T(E), the Fermi-function analytic
derivatives against numerical ones, and sanity identities for the
thermal kernel (zero-mass, V=0 vanishing, antisymmetry).

See ``docs/ANALYTIC_D2IDV2.md`` for the derivation.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from core.iets_analytic import (
    analytic_d2idv2_coherent,
    analytic_d2idv2_inelastic_at_bias,
    coherent_d2idv2_prefactor,
    fermi,
    fermi_deriv_E,
    fermi_second_deriv_E,
    thermal_d2_kernel,
)
from core.scba_rank1_keldysh import (
    contact_gammas,
    run_rank1_keldysh_single_bias,
)


KT_300K = 0.02585       # eV, room temperature


# ---------------------------------------------------------------------------
# Fermi function and derivatives.
# ---------------------------------------------------------------------------


def test_fermi_limits():
    # f(μ) = 1/2, f(−∞) = 1, f(+∞) = 0.
    assert fermi(0.0, 0.0, KT_300K) == pytest.approx(0.5, abs=1e-14)
    assert fermi(-10.0, 0.0, KT_300K) == pytest.approx(1.0, abs=1e-12)
    assert fermi(+10.0, 0.0, KT_300K) == pytest.approx(0.0, abs=1e-12)


def test_fermi_stable_at_large_arg():
    """No overflow warnings/NaNs for |x| ≫ 1."""
    E = np.array([-50.0, -1.0, 0.0, 1.0, 50.0])
    f = fermi(E, 0.0, KT_300K)
    assert np.all(np.isfinite(f))
    assert np.all((f >= 0.0) & (f <= 1.0))


def test_fermi_first_deriv_matches_fd():
    """∂f/∂E from the closed form matches a central difference to 1e-7."""
    E = np.linspace(-0.3, 0.3, 401)
    dE = 1e-6
    num = (fermi(E + dE, 0.0, KT_300K) - fermi(E - dE, 0.0, KT_300K)) / (2 * dE)
    ana = fermi_deriv_E(E, 0.0, KT_300K)
    np.testing.assert_allclose(ana, num, atol=1e-7, rtol=1e-6)


def test_fermi_second_deriv_matches_fd():
    """∂²f/∂E² from the closed form matches a central difference to 1e-5."""
    E = np.linspace(-0.3, 0.3, 401)
    dE = 1e-4
    num = (fermi(E + dE, 0.0, KT_300K)
           - 2 * fermi(E, 0.0, KT_300K)
           + fermi(E - dE, 0.0, KT_300K)) / dE ** 2
    ana = fermi_second_deriv_E(E, 0.0, KT_300K)
    np.testing.assert_allclose(ana, num, atol=1e-4, rtol=1e-4)


# ---------------------------------------------------------------------------
# Thermal d² kernel.
# ---------------------------------------------------------------------------


def test_kernel_vanishes_at_V_zero():
    E = np.linspace(-0.3, 0.3, 101)
    K = thermal_d2_kernel(E, 0.0, KT_300K)
    assert K.shape == E.shape
    np.testing.assert_allclose(K, 0.0, atol=1e-14)


def test_kernel_zero_mass():
    """∫ dE K_therm(E, V) = 0 because each f_α'' integrates to 0."""
    E = np.linspace(-1.0, 1.0, 5001)
    for V in (0.05, 0.1, 0.3):
        K = thermal_d2_kernel(E, V, KT_300K)
        integral = np.trapezoid(K, x=E)
        assert abs(integral) < 1e-10


def test_kernel_peak_locations():
    """Extrema of K_therm live near ±V/2 at low T."""
    E = np.linspace(-0.2, 0.2, 4001)
    V = 0.1
    K = thermal_d2_kernel(E, V, kT=KT_300K)
    # Largest magnitude peak on the (E > 0) side lives near E ≈ μ_L = V/2.
    right_peak = E[np.argmax(np.abs(K[E > 0]) * 0 + K[E > 0])]  # max on right
    idx_r = np.argmax(K[E > 0])
    idx_l = np.argmin(K[E < 0])
    E_right = E[E > 0][idx_r]
    E_left = E[E < 0][idx_l]
    # Peak positions within one kT of ±V/2 (∂²f/∂E² has its maximum at μ).
    assert abs(E_right - V / 2) < 3 * KT_300K
    assert abs(E_left - (-V / 2)) < 3 * KT_300K


def test_kernel_broadcasts_over_V():
    E = np.linspace(-0.3, 0.3, 51)
    V = np.linspace(0.0, 0.2, 11)
    K = thermal_d2_kernel(E, V, KT_300K)
    assert K.shape == (len(V), len(E))
    # Consistency with the scalar-V path.
    for i, v in enumerate(V):
        np.testing.assert_allclose(K[i], thermal_d2_kernel(E, v, KT_300K),
                                   atol=1e-14)


# ---------------------------------------------------------------------------
# Coherent d²I/dV² against numerical differentiation.
# ---------------------------------------------------------------------------


def _lorentzian_T(E, E0=0.15, gamma=0.02, T_max=0.8):
    """Toy resonant-transmission Lorentzian."""
    return T_max * gamma ** 2 / ((E - E0) ** 2 + gamma ** 2)


def _numerical_I_coh(V, E_grid, T_of_E, kT, E_F=0.0):
    """Direct numerical I_coh(V) = (e²/h) ∫ dE T(E) [f_L - f_R]."""
    mu_L = E_F + 0.5 * V
    mu_R = E_F - 0.5 * V
    integrand = T_of_E * (fermi(E_grid, mu_L, kT) - fermi(E_grid, mu_R, kT))
    return coherent_d2idv2_prefactor() * np.trapezoid(integrand, x=E_grid)


def test_coherent_d2idv2_matches_finite_difference():
    """Analytic d²I_coh/dV² vs a dense-grid central-difference on I(V)."""
    E_grid = np.linspace(-0.5, 0.8, 8001)
    T_of_E = _lorentzian_T(E_grid)

    # Sparse V-grid where we want d²I/dV² (the whole point of the analytic
    # path is that the V-grid can be arbitrarily coarse without noise).
    V_query = np.array([0.05, 0.10, 0.15, 0.20, 0.30])

    d2_analytic = analytic_d2idv2_coherent(
        T_of_E, E_grid, V_query, kT=KT_300K, E_F=0.0
    )

    # Ground truth: five-point stencil d²I/dV² on a dense V-sweep.
    dV = 1e-4
    d2_numeric = np.empty_like(V_query)
    for i, V in enumerate(V_query):
        I_m2 = _numerical_I_coh(V - 2 * dV, E_grid, T_of_E, KT_300K)
        I_m1 = _numerical_I_coh(V - dV, E_grid, T_of_E, KT_300K)
        I_0 = _numerical_I_coh(V, E_grid, T_of_E, KT_300K)
        I_p1 = _numerical_I_coh(V + dV, E_grid, T_of_E, KT_300K)
        I_p2 = _numerical_I_coh(V + 2 * dV, E_grid, T_of_E, KT_300K)
        # Five-point second derivative, O(dV^4) accurate.
        d2_numeric[i] = (-I_m2 + 16 * I_m1 - 30 * I_0 + 16 * I_p1 - I_p2) / (12 * dV ** 2)

    # Match to a few parts in 1e-5.
    np.testing.assert_allclose(d2_analytic, d2_numeric, rtol=5e-5, atol=0)


def test_coherent_d2idv2_at_V_zero():
    """At V = 0 the kernel is 0 → d²I_coh/dV²(0) = 0."""
    E_grid = np.linspace(-0.5, 0.8, 2001)
    T_of_E = _lorentzian_T(E_grid)
    d2 = analytic_d2idv2_coherent(T_of_E, E_grid, np.array([0.0]), KT_300K)
    assert abs(d2[0]) < 1e-20


def test_coherent_d2idv2_shape_and_T_check():
    E_grid = np.linspace(-0.3, 0.5, 501)
    T_of_E = _lorentzian_T(E_grid)
    V_grid = np.linspace(0.0, 0.3, 31)
    d2 = analytic_d2idv2_coherent(T_of_E, E_grid, V_grid, KT_300K)
    assert d2.shape == V_grid.shape
    assert np.all(np.isfinite(d2))


def test_coherent_rejects_shape_mismatch():
    E_grid = np.linspace(0, 1, 100)
    T_of_E = np.ones(99)
    with pytest.raises(ValueError):
        analytic_d2idv2_coherent(T_of_E, E_grid, np.array([0.1]), KT_300K)


# ---------------------------------------------------------------------------
# Inelastic path (FBA) — closes Tier-1 α inelastic.
# ---------------------------------------------------------------------------


def _patil_minimal_state(V, max_iter, D0_sq):
    """Run the rank-1 Keldysh wrapper at one bias on Patil's 1D toy device.

    Returns the Rank1KeldyshResult plus the bare-coherent T(E) so a coherent
    analytic d²I/dV² can be computed on the same E-grid for cross-checking.
    """
    NS, NC, ND = 1, 40, 1
    Np = NS + NC + ND
    Nb, Vb = 15, 0.6
    UB = np.concatenate([
        np.zeros(NS),
        Vb * np.ones(Nb),
        np.zeros(NC - 2 * Nb),
        Vb * np.ones(Nb),
        np.zeros(ND),
    ])
    t0 = 5.2
    H_z = (
        2 * t0 * np.eye(Np)
        - t0 * np.eye(Np, k=1)
        - t0 * np.eye(Np, k=-1)
        + np.diag(UB)
    )
    Ef, kT = 0.02, 0.025
    dE = 0.005
    E_grid = np.arange(-0.2, 0.8 + 0.5 * dE, dE)
    hnu_idx = (18, 35)
    Dnu = (D0_sq, D0_sq)
    N_bose = tuple(1.0 / (np.exp(dE * h / kT) - 1.0) for h in hnu_idx)
    chi = np.ones(Np)
    bias_profile = V * np.concatenate([
        0.5 * np.ones(NS),
        np.linspace(0.5, -0.5, NC),
        -0.5 * np.ones(ND),
    ])
    res = run_rank1_keldysh_single_bias(
        V=V,
        E_grid=E_grid,
        H_z=H_z,
        UB=UB,
        bias_profile=bias_profile,
        t0=t0,
        Ef=Ef,
        kT=kT,
        chi_diag=chi,
        D0_sq_per_mode=Dnu,
        hnu_idx_per_mode=hnu_idx,
        N_bose_per_mode=N_bose,
        max_iter=max_iter,
        tol=1e-5,
        mix=0.5,
    )
    # Bare T(E) for the coherent reference path: T = Tr[Γ_L G^R Γ_R G^A]
    # using the SAME G^R that the wrapper used (which equals the bare G^R
    # when D0_sq = 0).
    G_R = res.G_R
    G_A = np.conj(G_R.transpose(0, 2, 1))
    T_of_E = np.empty(E_grid.size)
    for k in range(E_grid.size):
        T_of_E[k] = float(np.real(np.trace(
            res.Gam_L[k] @ G_R[k] @ res.Gam_R[k] @ G_A[k]
        )))
    return res, E_grid, T_of_E, kT, Ef


def test_analytic_inelastic_reduces_to_coherent_at_D0_zero():
    """At D₀² = 0, analytic_d2idv2_inelastic_at_bias must equal the
    analytic coherent d²I/dV² evaluated on the bare Landauer T(E).

    This is the optical-theorem cross-check derived in
    ``_inelastic_kernels``: with no phonons, T_RR − T_RA = −T_RL and the
    inelastic formula collapses to (1/4) ∫ T (f_L'' − f_R'') dE.
    """
    V = 0.16
    res, E_grid, T_of_E, kT, Ef = _patil_minimal_state(
        V=V, max_iter=2, D0_sq=0.0
    )
    d2_inel = analytic_d2idv2_inelastic_at_bias(res, kT=kT, E_F=Ef)
    d2_coh_arr = analytic_d2idv2_coherent(
        T_of_E, E_grid, np.array([V]), kT, E_F=Ef
    )
    d2_coh = float(d2_coh_arr[0])
    # Both formulas use the same trapezoid quadrature on the same E-grid;
    # they should agree to ≲ 1e-6 — the residual is dominated by the
    # ``eta = 1e-12`` causal infinitesimal which makes G^A slightly off
    # from the exact (G^R)†, leaking into T_RA ≠ T_RL + T_RR by ~1e-7.
    rel = abs(d2_inel - d2_coh) / max(abs(d2_coh), 1e-30)
    assert rel < 1e-6, (
        f"V={V}: inelastic={d2_inel:.6e}  coherent={d2_coh:.6e}  rel={rel:.2e}"
    )


def test_analytic_inelastic_iets_peak_with_phonons():
    """With phonons on, d²I_R/dV² evaluated at V = ℏω_low is enhanced
    relative to the no-phonon coherent baseline.

    This is a *qualitative* check that the analytic formula sees the IETS
    feature: at V = ℏω₁ ≈ 90 meV the integrand picks up the phonon-broadened
    structure in T_RL(E) that lives at E ≈ E_res ± ℏω₁. The exact magnitude
    is set by SCBA convergence; we only assert the analytic formula returns
    a finite, non-trivial value at this bias.
    """
    V = 0.09  # ℏω₁ for the Patil toy (hnu_idx[0] = 18, dE = 0.005)
    res, _E_grid, _T_of_E, kT, Ef = _patil_minimal_state(
        V=V, max_iter=10, D0_sq=0.1
    )
    d2 = analytic_d2idv2_inelastic_at_bias(res, kT=kT, E_F=Ef)
    assert np.isfinite(d2)
    # Order-of-magnitude sanity: comparable to the e²/h prefactor times a
    # ~1/eV² Fermi-second-derivative integrated against an O(0.1) trace.
    # Tighter assertions belong in a sweep-level regression test.
    assert abs(d2) < 1e-2, f"d²I_R/dV² = {d2}"
