"""
Keldysh-consistent SCBA driver for the rank-1 projected solver.

Wraps the algebraic core in :mod:`core.scba_solver_rank1` with the bookkeeping
needed to actually compute observables: bare contact self-energies, the
Keldysh equation closing G̃^< from σ̃^< + contact in-scattering, the FBA/SCBA
phonon update on the energy axis, and the Meir-Wingreen current at both
contacts (for δ-current-conservation checks).

Scope (2026-04-09):

* **Single-mode (Nm = 1) is fully implemented and tested** against
  ``tests/patil_reference_1d.npz`` via ``tests/test_rank1_keldysh.py``.
  This is the configuration that closes SISPAD Tier-1 γ-comparison —
  the same physics Patil's `rtd2modes_1d.m` solves, plus the I_L/I_R
  decomposition needed for δ.

* **Multimode (Nm > 1)** is left as a documented stub. The rank-1
  algebra carries through verbatim for the *phonon* piece, but the
  *contact* in-scattering is mode-diagonal-but-not-rank-1, so the
  reconstruction of G^< from G̃^R requires explicitly summing the
  rank-1 correction over modes (Eq. 11 of ``docs/METHOD_DERIVATION.md``).
  That sum is the next-turn task — see ``_multimode_not_implemented``.

The energy axis follows Patil's convention: ℏω is an integer multiple of
dE so that ``G^<(E ± ℏω)`` is an array shift, not an interpolation.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from core.scba_solver_rank1 import (
    build_projected_bare_g,
    projected_dyson_step,
)


def _robust_solve(M, I_broadcast):
    """np.linalg.solve with regularization fallback for near-singular M.

    When the phonon self-energy pushes a pole of G^R onto the real axis,
    M becomes singular at specific (E, V) points. Adding a small imaginary
    part (1e-8) regularizes the solve without perceptibly changing the
    current at other energies.
    """
    try:
        return np.linalg.solve(M, I_broadcast)
    except np.linalg.LinAlgError:
        Np = M.shape[1]
        reg = 1e-8j * np.eye(Np, dtype=complex)[None, :, :]
        return np.linalg.solve(M + reg, I_broadcast)


# ---------------------------------------------------------------------------
# Bare per-mode retarded Green's function with tight-binding contacts
# ---------------------------------------------------------------------------
def build_bare_gR_single_mode(
    E_grid: np.ndarray,
    H_z: np.ndarray,
    U_bias: np.ndarray,
    UB: np.ndarray,
    t0: float,
    eta: float = 1e-12,
) -> np.ndarray:
    """Bare retarded Green's function for a single (1D) mode at every E.

    Implements Patil's tight-binding contact self-energy
    ``Σ_{1,1} = -t₀ exp(i ka)`` with ``ka = acos(1 − (E+iη − U − UB)/(2t₀))``,
    matching ``papers/rtd2modes_1d.m`` lines 85–88. The result is the same
    object the Python and Octave references in
    ``tests/patil_reference_1d.py`` and ``papers/rtd2modes_1d_octave.m``
    compute, lifted out so the SCBA loop can be reused for non-Patil
    geometries.

    Parameters
    ----------
    E_grid : ndarray, shape (NE,)
        Energy grid in eV.
    H_z : ndarray, shape (Np, Np)
        Tight-binding longitudinal Hamiltonian (no contacts, no bias drop).
    U_bias : ndarray, shape (Np,)
        Diagonal bias potential profile (linear drop across the device).
    UB : ndarray, shape (Np,)
        Static barrier profile (the on-site potential of the bare device,
        before bias).
    t0 : float
        Hopping amplitude in eV (Patil: 5.2 for GaAs).
    eta : float, optional
        Causal infinitesimal (Patil uses ``i·1e-12``).

    Returns
    -------
    G_R : ndarray, shape (NE, Np, Np), complex
        Bare retarded Green's function at each energy.
    """
    NE = E_grid.size
    Np = H_z.shape[0]
    G_R = np.zeros((NE, Np, Np), dtype=complex)
    H_eff = H_z + np.diag(U_bias)
    z_eta = 1j * eta
    eye = np.eye(Np, dtype=complex)
    for k, E in enumerate(E_grid):
        sig1 = np.zeros((Np, Np), dtype=complex)
        sig2 = np.zeros((Np, Np), dtype=complex)
        ck = 1.0 - ((E + z_eta - U_bias[0] - UB[0]) / (2.0 * t0))
        sig1[0, 0] = -t0 * np.exp(1j * np.arccos(ck))
        ck = 1.0 - ((E + z_eta - U_bias[-1] - UB[-1]) / (2.0 * t0))
        sig2[-1, -1] = -t0 * np.exp(1j * np.arccos(ck))
        M = (E + z_eta) * eye - H_eff - sig1 - sig2
        G_R[k] = np.linalg.solve(M, eye)
    return G_R, H_eff


def contact_gammas(
    E_grid: np.ndarray,
    U_bias: np.ndarray,
    UB: np.ndarray,
    t0: float,
    Np: int,
    eta: float = 1e-12,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (Γ_L, Γ_R), each shape (NE, Np, Np)."""
    NE = E_grid.size
    Gam_L = np.zeros((NE, Np, Np), dtype=complex)
    Gam_R = np.zeros((NE, Np, Np), dtype=complex)
    z_eta = 1j * eta
    for k, E in enumerate(E_grid):
        ck = 1.0 - ((E + z_eta - U_bias[0] - UB[0]) / (2.0 * t0))
        sig1 = -t0 * np.exp(1j * np.arccos(ck))
        Gam_L[k, 0, 0] = 1j * (sig1 - np.conj(sig1))
        ck = 1.0 - ((E + z_eta - U_bias[-1] - UB[-1]) / (2.0 * t0))
        sig2 = -t0 * np.exp(1j * np.arccos(ck))
        Gam_R[k, -1, -1] = 1j * (sig2 - np.conj(sig2))
    return Gam_L, Gam_R


# ---------------------------------------------------------------------------
# FBA/SCBA phonon self-energy update on the energy axis
# ---------------------------------------------------------------------------
def fba_phonon_sigma(
    G_lesser: np.ndarray,
    G_greater: np.ndarray,
    chi_diag: np.ndarray,
    D0_sq_per_mode: Sequence[float],
    hnu_idx_per_mode: Sequence[int],
    N_bose_per_mode: Sequence[float],
    chi_per_mode: Sequence[np.ndarray] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Einstein-vibron self-energies on the energy grid (Patil convention).

    Implements Eq. (6) of ``docs/METHOD_DERIVATION.md`` in the energy-shift
    form Patil uses (``ne(nn,nn,1:NE-inu) = n(nn,nn,inu+1:NE)``):

        σ̃^<(E) = Σ_modes D₀² · χ_ν² · [(N+1) G̃^<(E−ℏω) + N G̃^<(E+ℏω)]
        σ̃^>(E) = Σ_modes D₀² · χ_ν² · [N G̃^>(E−ℏω) + (N+1) G̃^>(E+ℏω)]

    Parameters
    ----------
    G_lesser, G_greater : ndarray, shape (NE, Np, Np)
        Projected G̃^{<,>} on the energy grid.
    chi_diag : ndarray, shape (Np,)
        Default longitudinal form factor χ(z), used for modes that don't
        have a per-mode override.
    D0_sq_per_mode : sequence of float
        D₀² for each phonon mode in eV² (Patil: ``[0.1, 0.1]``).
    hnu_idx_per_mode : sequence of int
        ℏω/dE for each mode (Patil: ``[18, 35]`` at dE=0.005 eV → 0.09, 0.175 eV).
    N_bose_per_mode : sequence of float
        Bose-Einstein occupation N(ℏω) for each mode at the lattice T.
    chi_per_mode : sequence of ndarray or None
        Per-mode χ(z) vectors, each shape (Np,). If None, all modes use
        chi_diag. This allows bulk phonons (χ=1 in active region) and
        molecular vibrons (Gaussian χ) to coexist in the same SCBA.

    Returns
    -------
    sigma_in, sigma_out : ndarray, shape (NE, Np, Np)
    """
    NE, Np, _ = G_lesser.shape
    sig_in = np.zeros_like(G_lesser)
    sig_out = np.zeros_like(G_lesser)
    diag_idx = np.arange(Np)
    Gl_diag = G_lesser[:, diag_idx, diag_idx]   # (NE, Np)
    Gg_diag = G_greater[:, diag_idx, diag_idx]
    for m, (D0_sq, inu, N) in enumerate(zip(
            D0_sq_per_mode, hnu_idx_per_mode, N_bose_per_mode)):
        if inu >= NE:
            continue
        # Per-mode χ² — allows different spatial profiles per phonon mode
        if chi_per_mode is not None:
            chi_m = chi_per_mode[m]
        else:
            chi_m = chi_diag
        chi_sq = (chi_m * chi_m).astype(float)
        ne_l = np.zeros_like(Gl_diag)
        na_l = np.zeros_like(Gl_diag)
        ne_g = np.zeros_like(Gg_diag)
        na_g = np.zeros_like(Gg_diag)
        ne_l[: NE - inu] = Gl_diag[inu:]
        na_l[inu:] = Gl_diag[: NE - inu]
        ne_g[: NE - inu] = Gg_diag[inu:]
        na_g[inu:] = Gg_diag[: NE - inu]
        coef_emit = (N + 1.0) * D0_sq
        coef_abs = N * D0_sq
        sig_in_diag = chi_sq[None, :] * (coef_emit * ne_l + coef_abs * na_l)
        sig_out_diag = chi_sq[None, :] * (coef_abs * ne_g + coef_emit * na_g)
        sig_in[:, diag_idx, diag_idx] += sig_in_diag
        sig_out[:, diag_idx, diag_idx] += sig_out_diag
    return sig_in, sig_out


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
@dataclass
class Rank1KeldyshResult:
    """Container for one bias point's worth of SCBA output."""

    V: float
    E_grid: np.ndarray
    G_R: np.ndarray  # (NE, Np, Np)
    G_lesser: np.ndarray  # (NE, Np, Np)
    G_greater: np.ndarray  # (NE, Np, Np)
    sigma_in_ph: np.ndarray  # (NE, Np, Np)
    sigma_out_ph: np.ndarray
    Gam_L: np.ndarray  # (NE, Np, Np), bare contact broadening — V-independent kernels need it
    Gam_R: np.ndarray  # (NE, Np, Np)
    I_left: float  # A, computed at the left contact
    I_right: float  # A, computed at the right contact (sign-flipped to be positive for L→R)
    iters_used: int
    converged: bool


_HBAR = 1.054571817e-34
_QE = 1.602176634e-19
_IE_PREFACTOR = (_QE * _QE) / (2.0 * np.pi * _HBAR)  # Patil "IE" constant


def run_rank1_keldysh_single_bias(
    *,
    V: float,
    E_grid: np.ndarray,
    H_z: np.ndarray,
    UB: np.ndarray,
    bias_profile: np.ndarray,
    t0: float,
    Ef: float,
    kT: float,
    chi_diag: np.ndarray,
    D0_sq_per_mode: Sequence[float],
    hnu_idx_per_mode: Sequence[int],
    N_bose_per_mode: Sequence[float],
    chi_per_mode: Sequence[np.ndarray] | None = None,
    max_iter: int = 200,
    tol: float = 1e-5,
    mix: float = 0.5,
    eta: float = 1e-12,
) -> Rank1KeldyshResult:
    """Run Patil's bulk-phonon SCBA at a single bias and return all observables.

    This is the *single-mode* (Nm = 1) Keldysh-consistent SCBA. The rank-1
    projected algebra in :mod:`core.scba_solver_rank1` collapses to a no-op
    here (|u|² = 1), and the entire loop runs on the standard 1D Green's
    function. The function exists in the rank-1 module rather than alongside
    the legacy Patil port because:

    1. Its API is the one the multimode wrapper will inherit.
    2. It exposes I_L *and* I_R separately so SISPAD δ
       (current-conservation) can be computed without re-running anything.
    3. It returns σ̃^{<,>} so SISPAD α (analytic d²I/dV² inelastic part)
       can be computed downstream as a one-liner integral.

    Parameters mirror Patil's MATLAB and the Python port; see
    ``tests/patil_reference_1d.py`` for a working example of input shapes.
    """
    Np = H_z.shape[0]
    NE = E_grid.size
    U_bias = bias_profile  # alias for clarity

    mu_L = Ef + V / 2.0
    mu_R = Ef - V / 2.0
    f_L = 1.0 / (1.0 + np.exp((E_grid - mu_L) / kT))
    f_R = 1.0 / (1.0 + np.exp((E_grid - mu_R) / kT))

    # Bare retarded G^R with tight-binding contacts (per energy).
    G_R_bare, _H_eff = build_bare_gR_single_mode(
        E_grid, H_z, U_bias, UB, t0, eta=eta
    )
    Gam_L, Gam_R = contact_gammas(E_grid, U_bias, UB, t0, Np, eta=eta)

    # Contact in/out scattering (per energy).
    sig_in_L = f_L[:, None, None] * Gam_L
    sig_in_R = f_R[:, None, None] * Gam_R
    sig_out_L = (1.0 - f_L)[:, None, None] * Gam_L
    sig_out_R = (1.0 - f_R)[:, None, None] * Gam_R

    # SCBA loop. Initialize σ̃_ph = 0 (FBA seed).
    sig_in_ph = np.zeros_like(G_R_bare)
    sig_out_ph = np.zeros_like(G_R_bare)
    G_lesser = np.zeros_like(G_R_bare)
    G_greater = np.zeros_like(G_R_bare)

    # Anderson mixing history (DIIS). We work with the diagonal elements
    # only: sig_{in,out}_ph are diagonal in Patil's formulation, so we
    # flatten the diagonals into a single real vector for the mixing.
    _diag = np.arange(Np)
    _anderson_depth = 8  # number of history vectors to keep
    _x_hist: list[np.ndarray] = []  # input vectors
    _r_hist: list[np.ndarray] = []  # residual vectors

    def _pack(si, so):
        """Pack diagonal elements of sig_in/sig_out into a 1D real vector."""
        di = si[:, _diag, _diag].ravel()
        do = so[:, _diag, _diag].ravel()
        return np.concatenate([di.real, di.imag, do.real, do.imag])

    def _unpack(vec, template):
        """Unpack 1D real vector back into (NE, Np, Np) diagonal arrays."""
        n = NE * Np
        di = vec[:n] + 1j * vec[n:2*n]
        do = vec[2*n:3*n] + 1j * vec[3*n:]
        si = np.zeros_like(template)
        so = np.zeros_like(template)
        si[:, _diag, _diag] = di.reshape(NE, Np)
        so[:, _diag, _diag] = do.reshape(NE, Np)
        return si, so

    # Precompute quantities that don't change per iteration.
    eye = np.eye(Np, dtype=complex)
    z_eta = 1j * eta
    E_arr = np.asarray(E_grid, dtype=complex) + z_eta
    H_full = H_z + np.diag(U_bias)
    sig_L = np.zeros((NE, Np, Np), dtype=complex)
    sig_R = np.zeros((NE, Np, Np), dtype=complex)
    ck_L = 1.0 - ((E_arr - U_bias[0] - UB[0]) / (2.0 * t0))
    ck_R = 1.0 - ((E_arr - U_bias[-1] - UB[-1]) / (2.0 * t0))
    sig_L[:, 0, 0] = -t0 * np.exp(1j * np.arccos(ck_L))
    sig_R[:, -1, -1] = -t0 * np.exp(1j * np.arccos(ck_R))
    _eye_broadcast = np.broadcast_to(eye, (NE, Np, Np)).copy()

    iters_used = 0
    converged = False
    for it in range(max_iter):
        iters_used = it + 1

        # Build G^R with current phonon self-energy.
        gamp_ph = sig_in_ph + sig_out_ph
        M = (
            E_arr[:, None, None] * eye[None, :, :]
            - H_full[None, :, :]
            - sig_L
            - sig_R
            + 0.5j * gamp_ph
        )
        G_R = _robust_solve(M, _eye_broadcast)
        G_A = np.conj(G_R.transpose(0, 2, 1))

        # Keldysh equation for G^< and G^>.
        sig_in_total = sig_in_L + sig_in_R + sig_in_ph
        sig_out_total = sig_out_L + sig_out_R + sig_out_ph
        G_lesser = np.matmul(np.matmul(G_R, sig_in_total), G_A)
        G_greater = np.matmul(np.matmul(G_R, sig_out_total), G_A)

        sig_in_new, sig_out_new = fba_phonon_sigma(
            G_lesser,
            G_greater,
            chi_diag=chi_diag,
            D0_sq_per_mode=D0_sq_per_mode,
            hnu_idx_per_mode=hnu_idx_per_mode,
            N_bose_per_mode=N_bose_per_mode,
            chi_per_mode=chi_per_mode,
        )

        # Relative convergence criterion.
        abs_change = float(np.sum(np.abs(sig_in_new - sig_in_ph))) + float(
            np.sum(np.abs(sig_out_new - sig_out_ph))
        )
        norm_old = float(np.sum(np.abs(sig_in_ph))) + float(
            np.sum(np.abs(sig_out_ph))
        )
        change = abs_change / max(norm_old, 1e-30)

        # --- Anderson mixing (DIIS) ---
        # x = current input, f(x) = sig_new (output of SCBA map), r = f(x) - x
        x_vec = _pack(sig_in_ph, sig_out_ph)
        f_vec = _pack(sig_in_new, sig_out_new)
        r_vec = f_vec - x_vec

        _x_hist.append(x_vec)
        _r_hist.append(r_vec)
        if len(_x_hist) > _anderson_depth:
            _x_hist.pop(0)
            _r_hist.pop(0)

        m = len(_r_hist)
        if m >= 2:
            # Build the residual overlap matrix and solve for optimal coeffs.
            # min ||Σ α_i r_i||² s.t. Σ α_i = 1
            # Equivalent to solving [R^T R, 1; 1^T, 0] [α; λ] = [0; 1]
            R = np.column_stack(_r_hist)  # (n, m)
            RtR = R.T @ R
            # Tikhonov regularization for numerical stability
            RtR += 1e-12 * np.eye(m)
            # Constrained least squares via Lagrange multiplier
            A = np.zeros((m + 1, m + 1))
            A[:m, :m] = RtR
            A[:m, m] = 1.0
            A[m, :m] = 1.0
            b = np.zeros(m + 1)
            b[m] = 1.0
            try:
                coeffs = np.linalg.solve(A, b)[:m]
            except np.linalg.LinAlgError:
                coeffs = np.ones(m) / m  # fallback to equal weights

            # Anderson update: x_next = Σ α_i (x_i + β * r_i)
            x_next = sum(c * (x + mix * r) for c, x, r
                         in zip(coeffs, _x_hist, _r_hist))
            sig_in_ph, sig_out_ph = _unpack(x_next, G_R_bare)
        else:
            # Not enough history yet — use simple linear mixing.
            sig_in_ph = (1.0 - mix) * sig_in_ph + mix * sig_in_new
            sig_out_ph = (1.0 - mix) * sig_out_ph + mix * sig_out_new

        if change < tol:
            converged = True
            break

    # ── FINAL PASS: recompute G with converged/updated phonon self-energy ──
    # Without this, FBA (max_iter=1) computes phonon Σ from ballistic G but
    # never feeds it back — making all molecules produce identical currents.
    # This is the rank-1 analogue of the "FINAL PASS" fix in scba_solver_hybrid.
    gamp_ph = sig_in_ph + sig_out_ph
    M = (
        E_arr[:, None, None] * eye[None, :, :]
        - H_full[None, :, :]
        - sig_L
        - sig_R
        + 0.5j * gamp_ph
    )
    G_R = _robust_solve(M, np.broadcast_to(eye, (NE, Np, Np)).copy())
    G_A = np.conj(G_R.transpose(0, 2, 1))
    sig_in_total = sig_in_L + sig_in_R + sig_in_ph
    sig_out_total = sig_out_L + sig_out_R + sig_out_ph
    G_lesser = np.matmul(np.matmul(G_R, sig_in_total), G_A)
    G_greater = np.matmul(np.matmul(G_R, sig_out_total), G_A)

    # Meir-Wingreen current at L and R contacts.
    # I_L = (e/h) ∫ dE Tr[Σ_out_L · n − Σ_in_L · p]
    # with n = G^< / i, p = -G^>/i (Patil convention has imag absorbed; we
    # follow his exact trace formula).
    # Patil's expression (lines 126-128 of rtd2modes_1d.m) operates on n,p
    # which are *real* arrays defined as n = real(G·Σ_in·G†), p = A−n. We
    # replicate that here for byte-equivalence with the Octave reference.
    n_arr = np.real(G_lesser)  # G^< has the same trace structure as Patil's "n"
    A_spec = 1j * (G_R - G_A)
    p_arr = np.real(A_spec) - n_arr

    dE = E_grid[1] - E_grid[0]
    I1 = 0.0  # I_R contact (Patil's "I1")
    I2 = 0.0  # I_L contact (Patil's "I2")
    for k in range(NE):
        I1 += float(
            np.real(
                np.trace(sig_out_R[k] @ n_arr[k]) - np.trace(sig_in_R[k] @ p_arr[k])
            )
        )
        I2 += float(
            np.real(
                np.trace(sig_out_L[k] @ n_arr[k]) - np.trace(sig_in_L[k] @ p_arr[k])
            )
        )
    I_right = I1 * dE * _IE_PREFACTOR
    I_left = I2 * dE * _IE_PREFACTOR

    return Rank1KeldyshResult(
        V=V,
        E_grid=E_grid,
        G_R=G_R,
        G_lesser=G_lesser,
        G_greater=G_greater,
        sigma_in_ph=sig_in_ph,
        sigma_out_ph=sig_out_ph,
        Gam_L=Gam_L,
        Gam_R=Gam_R,
        I_left=I_left,
        I_right=I_right,
        iters_used=iters_used,
        converged=converged,
    )


def _multimode_not_implemented(*_args, **_kwargs):
    """Stub for the Nm > 1 driver.

    The rank-1 *phonon* algebra carries through verbatim — what is missing
    is the reconstruction of G^< from the *contact* in-scattering, which is
    mode-diagonal-but-not-rank-1. The path forward (next turn) is:

    1. Build per-mode bare g_nm^R(E) for each transverse mode (n, m), with
       the energy argument shifted by ε_nm = (ℏ²/2m*)((nπ/L_y)² + (mπ/L_z)²).
    2. Build G̃_0^R = Σ_nm |u_nm|² g_nm^R via build_projected_bare_g.
    3. Solve the projected Dyson equation for G̃^R.
    4. Reconstruct each block of the *full* G^R via Eq. (11) of
       docs/METHOD_DERIVATION.md (rank-1 correction added to the diagonal).
    5. Compute G^< for each mode via the Keldysh equation with mode-diagonal
       contact in-scattering, then re-project to G̃^<.
    6. Run the SCBA loop on G̃^<, σ̃^<.

    The cost is dominated by step 4 which is one Np×Np matmul per (mode, E),
    not per (mode, SCBA iter), so the total work stays O(Np³ · NE) per
    iteration regardless of Nm.
    """
    raise NotImplementedError(
        "Multimode rank-1 Keldysh wrapper not yet implemented — see "
        "core.scba_rank1_keldysh._multimode_not_implemented for the path forward."
    )
