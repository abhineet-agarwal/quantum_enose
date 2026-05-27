"""
Rank-1 projected SCBA solver.

Implements the rank-1 projected Dyson/SCBA equations derived in
``docs/METHOD_DERIVATION.md`` §3-§5. The key structural result is that for a
molecular form factor g(r) with σ_mol ≪ L_⊥/n_max, the electron-vibron
coupling matrix in the transverse-mode basis factorises as

    M_{nm, n'm'}(z) = D₀ · χ(z − z_mol) · u_nm · u*_{n'm'}

with u_nm ≡ ψ_nm(r_⊥,mol). That rank-1 structure propagates to the phonon
self-energy (Σ^R = |u⟩⟨u| · σ̃^R) and collapses the full (mode ⊗ site) Dyson
equation to a single *scalar-in-mode-space*, 1D-in-longitudinal-space problem
for

    G̃^R(z, z'; E) ≡ ⟨u| G^R(z, z'; E) |u⟩.

The projected Dyson equation is

    G̃^R = [I − G̃_0^R σ̃^R]⁻¹ G̃_0^R,    G̃_0^R = Σ_nm |u_nm|² g_nm^R.

This module provides the minimum surface the test suite
(``tests/test_rank1_vs_full.py``) exercises: the projected bare G builder, a
single projected Dyson step, and a generic SCBA iteration helper. Keldysh
bookkeeping and the Meir-Wingreen current are left to the caller / follow-up
module — this file is the *algebraic core*.
"""

from __future__ import annotations

from typing import Callable

import numpy as np


def build_projected_bare_g(g_bare_per_mode: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Return G̃_0^R = Σ_nm |u_nm|² g_nm^R.

    Parameters
    ----------
    g_bare_per_mode : ndarray, shape (Nm, Np, Np)
        Bare (ballistic) retarded Green's function for each transverse mode nm,
        evaluated on the Np-site longitudinal grid at a single energy.
    u : ndarray, shape (Nm,)
        Transverse-mode projection vector u_nm = ψ_nm(r_⊥,mol).

    Returns
    -------
    G0_tilde : ndarray, shape (Np, Np)
        The |u|²-weighted mode sum of the bare per-mode Green's functions.
    """
    if g_bare_per_mode.ndim != 3:
        raise ValueError(
            f"g_bare_per_mode must have shape (Nm, Np, Np), got {g_bare_per_mode.shape}"
        )
    if u.shape[0] != g_bare_per_mode.shape[0]:
        raise ValueError(
            f"u has {u.shape[0]} modes but g_bare_per_mode has {g_bare_per_mode.shape[0]}"
        )
    weights = np.abs(u) ** 2
    return np.einsum("m,mij->ij", weights, g_bare_per_mode)


def projected_dyson_step(
    G0_tilde: np.ndarray, sigma_tilde: np.ndarray
) -> np.ndarray:
    """One projected Dyson solve: G̃^R = [I − G̃_0^R σ̃^R]⁻¹ G̃_0^R.

    Parameters
    ----------
    G0_tilde : ndarray, shape (Np, Np)
        Projected bare retarded Green's function (from
        :func:`build_projected_bare_g`).
    sigma_tilde : ndarray, shape (Np, Np)
        Scalar-in-mode-space phonon self-energy kernel σ̃^R(z, z'; E).

    Returns
    -------
    G_tilde : ndarray, shape (Np, Np)
        Projected dressed retarded Green's function.
    """
    if G0_tilde.shape != sigma_tilde.shape or G0_tilde.ndim != 2:
        raise ValueError(
            f"shape mismatch: G0_tilde {G0_tilde.shape}, sigma_tilde {sigma_tilde.shape}"
        )
    Np = G0_tilde.shape[0]
    I = np.eye(Np, dtype=complex)
    return np.linalg.solve(I - G0_tilde @ sigma_tilde, G0_tilde)


def rank1_scba_loop(
    g_bare_per_mode: np.ndarray,
    u: np.ndarray,
    sigma_initial: np.ndarray,
    update: Callable[[np.ndarray, np.ndarray], np.ndarray],
    n_iter: int,
    mix: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Iterate σ̃^R ← update(σ̃^R, G̃^R) for n_iter steps.

    A minimal SCBA driver that closes the loop between the projected Dyson
    equation and a caller-supplied ``update`` rule for σ̃^R. The physical
    Keldysh-consistent update (Σ^R from G^<, G^< from Σ^<) lives one layer up;
    this helper exists so the iteration machinery can be unit-tested against
    the brute-force full-basis solve in :mod:`tests.test_rank1_vs_full`.

    Parameters
    ----------
    g_bare_per_mode : ndarray, shape (Nm, Np, Np)
    u : ndarray, shape (Nm,)
    sigma_initial : ndarray, shape (Np, Np)
        Seed σ̃^R.
    update : callable (sigma_tilde, G_tilde) -> new_sigma_tilde
        Functional dependence of σ̃^R on the dressed G̃^R.
    n_iter : int
        Number of SCBA iterations.
    mix : float, optional
        Linear mixing factor between old and new σ̃^R each iteration
        (``σ̃ ← (1-mix)·σ̃ + mix·update(σ̃, G̃)``). Default 1.0 (no damping).

    Returns
    -------
    sigma_tilde : ndarray, shape (Np, Np)
    G_tilde : ndarray, shape (Np, Np)
    """
    G0_tilde = build_projected_bare_g(g_bare_per_mode, u)
    sigma_tilde = sigma_initial.copy()
    G_tilde = projected_dyson_step(G0_tilde, sigma_tilde)
    for _ in range(n_iter):
        G_tilde = projected_dyson_step(G0_tilde, sigma_tilde)
        sigma_new = update(sigma_tilde, G_tilde)
        sigma_tilde = (1.0 - mix) * sigma_tilde + mix * sigma_new
    return sigma_tilde, G_tilde
