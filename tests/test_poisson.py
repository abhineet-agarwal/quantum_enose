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
    poisson_newton_update,
    solve_poisson_1d,
)


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
