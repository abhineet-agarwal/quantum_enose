"""
Ground-truth test for the rank-1 projected SCBA.

This test is deliberately written **before** ``core/scba_solver_rank1.py``
exists. It constructs a tiny (Nm=2, Np=6) coupled-mode problem whose phonon
self-energy is rank-1 in mode space *by construction*, solves the full
coupled Dyson equation in the combined (mode ⊗ site) basis as ground truth,
then checks that the rank-1 projected solver reproduces the same projected
Green's function ``⟨u|G^R|u⟩`` to machine precision.

Three levels are tested, each progressively exercising more of the solver:

    1. ``test_projected_dyson_single_step`` — a single Dyson solve at fixed
       σ̃^R, no SCBA. Validates the core algebra
       G̃^R = [I − G̃_0^R σ̃^R]^(-1) G̃_0^R
       with G̃_0^R = Σ_nm |u_nm|² g_nm^R.

    2. ``test_projected_scba_selfconsistent`` — the full projected SCBA
       loop (several iterations), closing σ̃^R self-consistently from
       G̃^<.

    3. ``test_parity_decoupling_symmetric_box`` — a 3×3 transverse mode
       basis with an analytic ``u_nm = sin(nπ/2) sin(mπ/2)`` centred
       molecule. The (even, *) and (*, even) modes should decouple
       *exactly* from the phonon sector, and the projected solver should
       reduce to solving only the (odd, odd) block.

Ground truth: the brute-force solve in the full (mode ⊗ site) basis, i.e.
a (Nm·Np)×(Nm·Np) linear system. That's cheap at Nm=2, Np=6.

This file lives next to the existing solver smoke tests in ``tests/``. Run
with::

    source venv/bin/activate && python -m pytest tests/test_rank1_vs_full.py

and expect machine-precision (≤ 1e-12) agreement on all three levels.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pytest

# Make the repo importable when run via `pytest` from the project root.
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

try:
    from core.scba_solver_rank1 import (
        build_projected_bare_g,
        projected_dyson_step,
        rank1_scba_loop,
    )
    RANK1_SOLVER_AVAILABLE = True
except ImportError:  # module does not exist yet
    RANK1_SOLVER_AVAILABLE = False


# ---------------------------------------------------------------------------
# Helpers: ground-truth brute-force solver in the (mode ⊗ site) basis.
# ---------------------------------------------------------------------------


def build_full_sigma_r(u, sigma_tilde):
    """Full Σ^R in the (mode ⊗ site) basis as a rank-1 tensor.

    Parameters
    ----------
    u : ndarray, shape (Nm,)
        Transverse-mode projection vector u_nm = ψ_nm(r_mol).
    sigma_tilde : ndarray, shape (Np, Np)
        Scalar-in-mode-space phonon self-energy kernel σ̃^R(z, z').

    Returns
    -------
    Sigma_full : ndarray, shape (Nm*Np, Nm*Np)
        Σ^R[(nm, z), (n'm', z')] = u_nm * u_{n'm'} * σ̃^R(z, z').
    """
    Nm = u.shape[0]
    Np = sigma_tilde.shape[0]
    # outer product in mode space, tensor with the site kernel
    uu = np.outer(u, u.conj())                  # (Nm, Nm)
    Sigma_full = np.kron(uu, sigma_tilde)       # (Nm*Np, Nm*Np)
    assert Sigma_full.shape == (Nm * Np, Nm * Np)
    return Sigma_full


def build_full_bare_g(g_bare_per_mode):
    """Block-diagonal bare G_0^R in the (mode ⊗ site) basis.

    Parameters
    ----------
    g_bare_per_mode : ndarray, shape (Nm, Np, Np)

    Returns
    -------
    G0_full : ndarray, shape (Nm*Np, Nm*Np)
    """
    Nm, Np, _ = g_bare_per_mode.shape
    G0_full = np.zeros((Nm * Np, Nm * Np), dtype=complex)
    for nm in range(Nm):
        G0_full[nm * Np:(nm + 1) * Np, nm * Np:(nm + 1) * Np] = g_bare_per_mode[nm]
    return G0_full


def full_dyson_solve(g_bare_per_mode, u, sigma_tilde):
    """Solve G^R = (I − G_0^R Σ^R)^(-1) G_0^R in the full basis.

    Returns both the full G^R and its projection ⟨u|G^R|u⟩.
    """
    Nm, Np, _ = g_bare_per_mode.shape
    G0_full = build_full_bare_g(g_bare_per_mode)
    Sigma_full = build_full_sigma_r(u, sigma_tilde)

    I = np.eye(Nm * Np, dtype=complex)
    G_full = np.linalg.solve(I - G0_full @ Sigma_full, G0_full)

    # Projection: G̃^R(z, z') = Σ_{nm, n'm'} u*_nm G^R_[(nm, z), (n'm', z')] u_{n'm'}
    G_tilde = np.zeros((Np, Np), dtype=complex)
    for nm in range(Nm):
        for nmp in range(Nm):
            block = G_full[nm * Np:(nm + 1) * Np, nmp * Np:(nmp + 1) * Np]
            G_tilde += np.conj(u[nm]) * u[nmp] * block
    return G_full, G_tilde


def projected_bare_g_reference(g_bare_per_mode, u):
    """Reference G̃_0^R = Σ_nm |u_nm|² g_nm^R."""
    Nm = g_bare_per_mode.shape[0]
    weights = np.abs(u) ** 2
    return np.einsum("m,mij->ij", weights, g_bare_per_mode)


# ---------------------------------------------------------------------------
# Fixtures: a reproducible tiny toy.
# ---------------------------------------------------------------------------


@pytest.fixture
def tiny_toy():
    """2-mode, 6-site toy problem with random complex g_nm^R and σ̃^R."""
    rng = np.random.default_rng(20260407)
    Nm, Np = 2, 6

    # Random ballistic g_nm^R per mode. Retarded GFs have negative imaginary
    # part on the diagonal; bake that in so (I − G_0^R Σ^R) is well-conditioned.
    g_bare = np.zeros((Nm, Np, Np), dtype=complex)
    for nm in range(Nm):
        A = rng.standard_normal((Np, Np)) + 1j * rng.standard_normal((Np, Np))
        A = 0.1 * (A + A.conj().T)
        A -= 1j * np.eye(Np)       # retarded pole structure
        g_bare[nm] = A

    # Random real u (matches the physical u_nm = ψ_nm(r_mol) case)
    u = rng.standard_normal(Nm)
    u /= np.linalg.norm(u)

    # Random small σ̃^R so the Born series converges
    sigma_tilde = 0.05 * (
        rng.standard_normal((Np, Np)) + 1j * rng.standard_normal((Np, Np))
    )

    return dict(g_bare=g_bare, u=u, sigma_tilde=sigma_tilde, Nm=Nm, Np=Np)


# ---------------------------------------------------------------------------
# Test 1: single-step projected Dyson solve matches the full basis.
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not RANK1_SOLVER_AVAILABLE,
    reason="core/scba_solver_rank1.py not implemented yet",
)
def test_projected_dyson_single_step(tiny_toy):
    """G̃^R from the rank-1 solver == ⟨u|G^R|u⟩ from the brute-force solve."""
    g_bare = tiny_toy["g_bare"]
    u = tiny_toy["u"]
    sigma_tilde = tiny_toy["sigma_tilde"]

    # --- Ground truth: full (mode ⊗ site) Dyson solve ---
    _, G_tilde_truth = full_dyson_solve(g_bare, u, sigma_tilde)

    # --- Rank-1 solver under test ---
    G0_tilde = build_projected_bare_g(g_bare, u)
    G_tilde_rank1 = projected_dyson_step(G0_tilde, sigma_tilde)

    # (a) G̃_0^R matches the analytic |u|²-weighted sum.
    G0_tilde_ref = projected_bare_g_reference(g_bare, u)
    np.testing.assert_allclose(
        G0_tilde, G0_tilde_ref, atol=1e-14, rtol=0,
        err_msg="build_projected_bare_g does not equal Σ_nm |u_nm|² g_nm^R",
    )

    # (b) G̃^R from the rank-1 solver matches ⟨u|G^R|u⟩ from the full solve.
    np.testing.assert_allclose(
        G_tilde_rank1, G_tilde_truth, atol=1e-12, rtol=0,
        err_msg="Rank-1 projected Dyson step disagrees with the full Dyson solve",
    )


# ---------------------------------------------------------------------------
# Test 2: full self-consistent SCBA loop agrees with the full-basis loop.
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not RANK1_SOLVER_AVAILABLE,
    reason="core/scba_solver_rank1.py not implemented yet",
)
def test_projected_scba_selfconsistent(tiny_toy):
    """After N SCBA iterations, both solvers must converge to the same G̃^R.

    We define a trivial update rule for σ̃^R that closes self-consistency:

        σ̃^R ← α · (G̃^R)  + β · σ̃_0

    with α, β small. This is not physical FBA — the point is only to exercise
    the iterative loop with a nontrivial functional dependence so that any
    sign or projection error would blow up after a few iterations.
    """
    g_bare = tiny_toy["g_bare"]
    u = tiny_toy["u"]
    sigma_seed = tiny_toy["sigma_tilde"]
    Nm, Np = tiny_toy["Nm"], tiny_toy["Np"]

    alpha, beta = 0.15, 0.60
    n_iter = 6

    def update(sigma_tilde, G_tilde):
        return alpha * G_tilde + beta * sigma_seed

    # --- Ground-truth loop in the full basis ---
    sigma_full = sigma_seed.copy()
    for _ in range(n_iter):
        _, G_tilde_full = full_dyson_solve(g_bare, u, sigma_full)
        sigma_full = update(sigma_full, G_tilde_full)

    # --- Rank-1 loop ---
    sigma_r1 = sigma_seed.copy()
    G0_tilde = build_projected_bare_g(g_bare, u)
    for _ in range(n_iter):
        G_tilde_r1 = projected_dyson_step(G0_tilde, sigma_r1)
        sigma_r1 = update(sigma_r1, G_tilde_r1)

    np.testing.assert_allclose(
        sigma_r1, sigma_full, atol=1e-11, rtol=0,
        err_msg="σ̃^R diverges between rank-1 and full solvers after SCBA iterations",
    )
    np.testing.assert_allclose(
        G_tilde_r1, G_tilde_full, atol=1e-11, rtol=0,
        err_msg="G̃^R diverges between rank-1 and full solvers after SCBA iterations",
    )


# ---------------------------------------------------------------------------
# Test 3: parity decoupling in a 3×3 symmetric box.
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not RANK1_SOLVER_AVAILABLE,
    reason="core/scba_solver_rank1.py not implemented yet",
)
def test_parity_decoupling_symmetric_box():
    """For u_nm ∝ sin(nπ/2)·sin(mπ/2), only the 4 odd-odd modes couple.

    This is the METHOD_DERIVATION.md §7 theorem: with the molecule at the
    symmetric centre of the box, u_nm vanishes whenever n or m is even, so
    5 of the 9 transverse modes decouple from the phonon sector identically.
    The projected solver must therefore give the same answer as solving only
    the 4-mode (odd-odd) subproblem and leaving the other 5 ballistic.
    """
    rng = np.random.default_rng(20260408)
    n_max = m_max = 3
    Np = 6

    # Build u_nm = sin(nπ/2) sin(mπ/2) for (n,m) ∈ {1,2,3}². Computed via
    # np.sin, even indices give ~1e-16 rather than exact zero, so we snap.
    n_idx = np.arange(1, n_max + 1)
    m_idx = np.arange(1, m_max + 1)
    U = np.sin(np.pi * n_idx / 2)[:, None] * np.sin(np.pi * m_idx / 2)[None, :]
    U[np.abs(U) < 1e-12] = 0.0
    u = U.reshape(-1)                             # (9,)
    assert np.count_nonzero(u) == 4, "centred (odd,odd) parity should give exactly 4 nonzero modes"

    odd_mask = np.abs(u) > 1e-14                  # (9,)
    u_odd = u[odd_mask]                           # (4,)

    # Random bare g_nm^R for all 9 modes, random σ̃^R
    g_bare = np.zeros((9, Np, Np), dtype=complex)
    for nm in range(9):
        A = rng.standard_normal((Np, Np)) + 1j * rng.standard_normal((Np, Np))
        A = 0.1 * (A + A.conj().T) - 1j * np.eye(Np)
        g_bare[nm] = A

    sigma_tilde = 0.05 * (
        rng.standard_normal((Np, Np)) + 1j * rng.standard_normal((Np, Np))
    )

    # (a) Full 9-mode projected solve
    G0_tilde_9 = build_projected_bare_g(g_bare, u)
    G_tilde_9 = projected_dyson_step(G0_tilde_9, sigma_tilde)

    # (b) Reduced 4-mode (odd-odd) projected solve — same molecule, same
    #     phonon, different bare-G set
    g_bare_odd = g_bare[odd_mask]                 # (4, Np, Np)
    G0_tilde_4 = build_projected_bare_g(g_bare_odd, u_odd)
    G_tilde_4 = projected_dyson_step(G0_tilde_4, sigma_tilde)

    # The projected Green's function only sees Σ_nm |u_nm|² g_nm^R, so the
    # even-mode g_nm^R contributions are all killed by their zero u weight.
    # Both builds must therefore give identical G̃_0^R and identical G̃^R.
    np.testing.assert_allclose(
        G0_tilde_9, G0_tilde_4, atol=1e-14, rtol=0,
        err_msg="G̃_0^R differs between full 9-mode and odd-only 4-mode builds — parity theorem broken",
    )
    np.testing.assert_allclose(
        G_tilde_9, G_tilde_4, atol=1e-14, rtol=0,
        err_msg="G̃^R differs between full 9-mode and odd-only 4-mode builds — parity theorem broken",
    )


# ---------------------------------------------------------------------------
# Stand-alone runner (for quick iteration without pytest)
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    if not RANK1_SOLVER_AVAILABLE:
        print("core/scba_solver_rank1.py not importable — run the test via pytest "
              "once the solver exists.")
        sys.exit(1)

    # minimal manual run — mirrors the fixture
    rng = np.random.default_rng(20260407)
    Nm, Np = 2, 6
    g_bare = np.zeros((Nm, Np, Np), dtype=complex)
    for nm in range(Nm):
        A = rng.standard_normal((Np, Np)) + 1j * rng.standard_normal((Np, Np))
        A = 0.1 * (A + A.conj().T) - 1j * np.eye(Np)
        g_bare[nm] = A
    u = rng.standard_normal(Nm); u /= np.linalg.norm(u)
    sigma_tilde = 0.05 * (rng.standard_normal((Np, Np)) + 1j * rng.standard_normal((Np, Np)))

    _, G_tilde_truth = full_dyson_solve(g_bare, u, sigma_tilde)
    G0_tilde = build_projected_bare_g(g_bare, u)
    G_tilde_r1 = projected_dyson_step(G0_tilde, sigma_tilde)
    err = np.max(np.abs(G_tilde_r1 - G_tilde_truth))
    print(f"max |G̃^R rank-1 − G̃^R full| = {err:.3e}")
    assert err < 1e-12, "Rank-1 solver diverges from full basis"
    print("OK")
