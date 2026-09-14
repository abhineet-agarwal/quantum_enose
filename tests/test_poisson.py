"""Unit tests for the 1D Poisson solver and NEGF density extraction.

These exercise :mod:`core.poisson` in isolation (no NEGF solve): analytic
parabolic potential for uniform charge, exact linear ramp for the
applied-bias-only case, the charge-neutral null case, and the density
extraction + contact-neutrality calibration.
"""
import numpy as np
import pytest

from core.poisson import (
    _EPS0,
    _Q,
    calibrate_density_prefactor,
    extract_density_1d,
    fermi_dirac_half,
    fermi_dirac_minus_half,
    poisson_newton_update,
    solve_poisson_1d,
)


def test_fermi_dirac_half_known_values():
    # Blakemore convention: F_{1/2}(0) = (2/√π)·(1−1/√2)·Γ(3/2)·ζ(3/2) ≈ 0.7654.
    assert fermi_dirac_half(0.0) == pytest.approx(0.7654, rel=1e-2)
    # Deep Boltzmann limit: F_{1/2}(η) → e^η.
    assert fermi_dirac_half(-8.0) == pytest.approx(np.exp(-8.0), rel=1e-3)
    # Deep degenerate limit: F_{1/2}(η) → (4/(3√π)) η^(3/2).
    asy = (4.0 / (3.0 * np.sqrt(np.pi))) * 20.0**1.5
    assert fermi_dirac_half(20.0) == pytest.approx(asy, rel=1e-2)


def test_fermi_dirac_minus_half_is_derivative():
    # F_{-1/2}(η) = dF_{1/2}/dη — check at η=0.5 against a tighter FD.
    eta0 = 0.5
    h = 1e-3
    expected = (fermi_dirac_half(eta0 + h) - fermi_dirac_half(eta0 - h)) / (2 * h)
    got = fermi_dirac_minus_half(eta0)
    assert got == pytest.approx(expected, rel=1e-3)


def test_uniform_charge_gives_analytic_parabola():
    # Uniform donor charge, fully depleted (n=0), grounded contacts.
    Np = 51
    a = 0.3e-9
    eps_r = np.full(Np, 10.0)
    D0 = 1.0e24
    N_D = np.full(Np, D0)
    n_e = np.zeros(Np)

    U = solve_poisson_1d(eps_r, N_D, n_e, a, U_left=0.0, U_right=0.0)

    # Continuous solution of eps U'' = q N_D with U(0)=U(L)=0:
    #   U(z) = (q N_D / 2 eps) z (z - L)
    z = np.arange(Np) * a
    L = (Np - 1) * a
    eps = 10.0 * _EPS0
    U_analytic = (_Q * D0 / (2.0 * eps)) * z * (z - L)

    # The 3-point Laplacian is exact for quadratics → machine precision.
    assert np.allclose(U, U_analytic, atol=1e-9 * np.abs(U_analytic).max())
    assert U.min() < 0.0  # positive donor charge dips the electron energy


def test_charge_neutral_is_zero_potential():
    # n_e == N_D everywhere, grounded contacts → U ≡ 0.
    Np = 40
    a = 0.2e-9
    eps_r = np.full(Np, 8.5)
    N_D = np.full(Np, 1e25)
    n_e = N_D.copy()
    U = solve_poisson_1d(eps_r, N_D, n_e, a, U_left=0.0, U_right=0.0)
    assert np.allclose(U, 0.0, atol=1e-12)


def test_applied_bias_only_is_linear_ramp():
    # No net charge, Dirichlet ±V/2 → exact linear drop across the device.
    Np = 31
    a = 0.2e-9
    eps_r = np.full(Np, 8.5)
    N_D = np.zeros(Np)
    n_e = np.zeros(Np)
    V = 0.4
    U = solve_poisson_1d(eps_r, N_D, n_e, a, U_left=+V / 2, U_right=-V / 2)
    U_linear = np.linspace(+V / 2, -V / 2, Np)
    assert np.allclose(U, U_linear, atol=1e-12)


def test_neumann_right_zero_charge_is_flat():
    # Dirichlet left = +V/2, zero-field Neumann right, no charge → Laplace with
    # one value BC and zero slope at the other end forces a *constant* profile
    # (flat band pinned at the Dirichlet value, zero field everywhere).
    Np = 31
    a = 0.2e-9
    eps_r = np.full(Np, 8.5)
    N_D = np.zeros(Np)
    n_e = np.zeros(Np)
    V = 0.4
    U = solve_poisson_1d(eps_r, N_D, n_e, a, U_left=+V / 2, U_right=0.0,
                         right_bc="neumann")
    assert np.allclose(U, +V / 2, atol=1e-12)


def test_neumann_right_uniform_charge_is_parabola_zero_slope():
    # Uniform donor charge, Dirichlet left = 0, zero-field Neumann right.
    # Continuous solution of eps U'' = q N_D with U(0)=0, U'(L)=0:
    #   U(z) = (q N_D / 2 eps)(z² − 2 L z),   U'(L) = 0 ✓
    Np = 41
    a = 0.3e-9
    eps_r = np.full(Np, 10.0)
    D0 = 1.0e24
    N_D = np.full(Np, D0)
    n_e = np.zeros(Np)

    U = solve_poisson_1d(eps_r, N_D, n_e, a, U_left=0.0, U_right=0.0,
                         right_bc="neumann")

    z = np.arange(Np) * a
    L = (Np - 1) * a
    eps = 10.0 * _EPS0
    U_analytic = (_Q * D0 / (2.0 * eps)) * (z * z - 2.0 * L * z)

    # FV half-cell Neumann is 2nd-order; quadratic solution recovered closely.
    assert np.allclose(U, U_analytic, rtol=2e-3, atol=1e-9 * np.abs(U_analytic).max())
    # Discrete zero-field at the right end: last two nodes nearly equal.
    assert abs(U[-1] - U[-2]) < 1e-3 * abs(U_analytic.min())


def test_newton_step_neumann_shares_direct_fixed_point():
    # A Newton step (Dirichlet-left + Neumann-right) started from the exact
    # direct Neumann solve on the same charge must not move it.
    Np = 41
    a = 0.2e-9
    eps_r = np.full(Np, 8.5)
    N_D = np.zeros(Np)
    N_D[:8] = 1e25
    N_D[-8:] = 1e25
    rng = np.random.default_rng(1)
    n_e = N_D * (0.9 + 0.2 * rng.random(Np))
    U_L = 0.1
    U_exact = solve_poisson_1d(eps_r, N_D, n_e, a, U_left=U_L, U_right=0.0,
                               right_bc="neumann")
    U_step = poisson_newton_update(U_exact, eps_r, N_D, n_e, a, U_L, 0.0,
                                   kT_screen=0.002, right_bc="neumann")
    assert np.allclose(U_step, U_exact, atol=1e-9)
    # Left end pinned (gauge); right end is free (zero-field face), NOT pinned
    # to the passed U_right=0 — that is what distinguishes it from Dirichlet.
    assert U_step[0] == pytest.approx(U_L)
    assert abs(U_step[-1]) > 1e-3  # floated well away from the nominal 0


def test_newton_both_ends_neumann_rejected():
    Np = 20
    a = 0.2e-9
    eps_r = np.full(Np, 8.5)
    N_D = np.zeros(Np)
    n_e = np.zeros(Np)
    U0 = np.zeros(Np)
    with pytest.raises(ValueError, match="singular"):
        poisson_newton_update(U0, eps_r, N_D, n_e, a, 0.1, -0.1,
                              kT_screen=0.002, left_bc="neumann", right_bc="neumann")


def test_both_ends_neumann_rejected():
    Np = 20
    a = 0.2e-9
    eps_r = np.full(Np, 8.5)
    N_D = np.zeros(Np)
    n_e = np.zeros(Np)
    with pytest.raises(ValueError, match="singular"):
        solve_poisson_1d(eps_r, N_D, n_e, a, 0.2, -0.2,
                         left_bc="neumann", right_bc="neumann")


def test_variable_epsilon_continuity():
    # With a permittivity step but zero charge, flux ε dU/dz is continuous,
    # so U is piecewise-linear and monotone between the Dirichlet ends.
    Np = 21
    a = 0.2e-9
    eps_r = np.full(Np, 8.5)
    eps_r[Np // 2:] = 25.0  # high-k half
    N_D = np.zeros(Np)
    n_e = np.zeros(Np)
    U = solve_poisson_1d(eps_r, N_D, n_e, a, U_left=+0.2, U_right=-0.2)
    assert U[0] == pytest.approx(0.2)
    assert U[-1] == pytest.approx(-0.2)
    # Monotone decreasing (no charge → no local extrema).
    assert np.all(np.diff(U) <= 1e-12)


def test_newton_step_fixes_the_exact_solution():
    # poisson_newton_update must share the fixed point of solve_poisson_1d:
    # starting from the exact direct solution, one Newton step changes nothing.
    Np = 41
    a = 0.2e-9
    eps_r = np.full(Np, 8.5)
    N_D = np.zeros(Np)
    N_D[:8] = 1e25
    N_D[-8:] = 1e25
    rng = np.random.default_rng(0)
    n_e = N_D * (0.9 + 0.2 * rng.random(Np))  # arbitrary charge
    U_L, U_R = 0.1, -0.1
    U_exact = solve_poisson_1d(eps_r, N_D, n_e, a, U_L, U_R)
    U_step = poisson_newton_update(U_exact, eps_r, N_D, n_e, a, U_L, U_R,
                                   kT_screen=0.002)
    assert np.allclose(U_step, U_exact, atol=1e-9)


def test_newton_step_preserves_boundary_values():
    Np = 30
    a = 0.2e-9
    eps_r = np.full(Np, 8.5)
    N_D = np.zeros(Np)
    n_e = np.full(Np, 1e24)
    U0 = np.linspace(0.2, -0.2, Np)
    U1 = poisson_newton_update(U0, eps_r, N_D, n_e, a, 0.2, -0.2,
                               kT_screen=0.002)
    assert U1[0] == pytest.approx(0.2)
    assert U1[-1] == pytest.approx(-0.2)


def test_extract_density_integral():
    # Synthetic G^< with a known real diagonal → N_raw = (1/2π)∫dE diag.
    NE, Np = 100, 5
    dE = 0.002
    G = np.zeros((NE, Np, Np), dtype=complex)
    # constant diagonal value d on every site/energy
    d = 0.3
    for z in range(Np):
        G[:, z, z] = d + 0.0j
    N_raw = extract_density_1d(G, dE)
    expected = NE * d * dE / (2.0 * np.pi)
    assert np.allclose(N_raw, expected)


def test_calibrate_prefactor_neutrality():
    Np = 20
    N_raw = np.full(Np, 0.5)
    N_D = np.zeros(Np)
    N_D[:5] = 1e25
    C = calibrate_density_prefactor(N_raw, N_D, slice(0, 5))
    assert C == pytest.approx(1e25 / 0.5)
    # Applying C to the contact occupation recovers N_D there.
    assert (C * N_raw[:5]) == pytest.approx(np.full(5, 1e25))


def test_calibrate_prefactor_rejects_zero_occupation():
    with pytest.raises(ValueError):
        calibrate_density_prefactor(np.zeros(10), np.ones(10) * 1e25, slice(0, 5))
