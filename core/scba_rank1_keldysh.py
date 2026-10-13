"""Single-transverse-mode local NEGF-SCBA and current integration helpers.

The production solver uses standard 1D Green functions and longitudinally
diagonal phonon self-energies. Its historic rank1 name does not imply an
implemented multimode projection: the corresponding Keldysh reconstruction is
still a stub. Separate tests validate the algebra in core.scba_solver_rank1.

Contact currents are per transverse mode per spin. Transverse quadrature sums
independent channels and omits localized intermode coupling. Energy shifts for
phonons use integer multiples of the uniform energy spacing. Solver convergence
and contact conservation do not establish grid convergence of molecular spectra.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from core.validation import current_diagnostics

from core.scba_solver_rank1 import (
    build_projected_bare_g,
    projected_dyson_step,
)


def _nonnegative_finite(*arrays):
    """Diagonal rates/populations must be real and nonnegative, up to roundoff."""
    for value in arrays:
        value = np.asarray(value)
        if not np.all(np.isfinite(value)):
            return False
        floor = 1e-14 + 1e-12 * float(np.max(np.abs(value)))
        if np.min(value.real) < -floor or np.max(np.abs(value.imag)) > floor:
            return False
    return True


def _relative_sigma_residual(si, so, si_new, so_new):
    return float((np.abs(si_new-si).sum() + np.abs(so_new-so).sum())
                 / max(np.abs(si).sum() + np.abs(so).sum(), 1e-30))


class _ScbaMixer:
    """Optional DIIS, with rejection of nonphysical proposals and history reset.

    Never clip a negative extrapolated rate into an apparently converged state.
    Fall back to the convex linear update from the current physical iterate.
    """

    def __init__(self, mix, depth, max_iter, tol, guard=True):
        if not np.isfinite(mix) or not 0 < mix <= 1:
            raise ValueError("SCBA mix must be in (0, 1]")
        if depth < 0 or int(depth) != depth or max_iter < 1 or not np.isfinite(tol) or tol <= 0:
            raise ValueError("Require nonnegative integer Anderson depth, max_iter >= 1 and tol > 0")
        self.mix, self.depth, self.guard = mix, depth, guard
        self.x, self.r = [], []
        self.fallbacks = 0

    def update(self, si, so, si_new, so_new):
        linear = ((1-self.mix)*si + self.mix*si_new,
                  (1-self.mix)*so + self.mix*so_new)
        if self.depth < 2:
            return linear
        pack = lambda a, b: np.concatenate([a.real.ravel(), a.imag.ravel(),
                                            b.real.ravel(), b.imag.ravel()])
        x, f = pack(si, so), pack(si_new, so_new)
        self.x.append(x); self.r.append(f-x)
        if len(self.x) > self.depth:
            self.x.pop(0); self.r.pop(0)
        m = len(self.r)
        if m < 2:
            return linear
        R = np.column_stack(self.r)
        A = np.zeros((m+1, m+1)); A[:m, :m] = R.T @ R + 1e-12*np.eye(m)
        A[:m, m] = 1; A[m, :m] = 1
        b = np.zeros(m+1); b[m] = 1
        try:
            coeffs = np.linalg.solve(A, b)[:m]
        except np.linalg.LinAlgError:
            if self.guard:
                coeffs = None
            else:                       # historical fallback: equal weights
                coeffs = np.ones(m) / m
        if coeffs is not None:
            proposal = sum(c*(x+self.mix*r) for c, x, r in zip(coeffs, self.x, self.r))
            n = si.size
            candidate = ((proposal[:n]+1j*proposal[n:2*n]).reshape(si.shape),
                         (proposal[2*n:3*n]+1j*proposal[3*n:]).reshape(so.shape))
            # _nonnegative_finite already allows a roundoff floor; a strict
            # ">= 0" on top of it rejected roundoff-level negatives too
            if not self.guard or _nonnegative_finite(*candidate):
                return candidate
        self.x.clear(); self.r.clear()
        self.fallbacks += 1
        return linear


def _contact_current_pair(n, p, si_c, so_c, dE):
    """Outward contact currents from diagonal populations and contact rates."""
    return tuple(float(np.sum(np.real(so_c[:, j]*n[:, j] - si_c[:, j]*p[:, j])))
                 * dE * _IE_PREFACTOR for j in (0, -1))


def _current_ok(left, right, rtol, atol):
    return rtol is None or bool(current_diagnostics(left, right, True, rtol, atol)["current_valid"])


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
    physical_valid: bool = True
    self_energy_residual: float = float("inf")
    mixer_fallbacks: int = 0


_HBAR = 1.054571817e-34
_QE = 1.602176634e-19
_IE_PREFACTOR = (_QE * _QE) / (2.0 * np.pi * _HBAR)  # Patil "IE" constant


_HBAR_JS = 1.054571817e-34  # J*s


def coherent_transmission(G_R: np.ndarray, Gam_L: np.ndarray,
                          Gam_R: np.ndarray) -> np.ndarray:
    """Elastic transmission T(E) = Tr[Gam_L G^R Gam_R G^A] (Fisher-Lee/Caroli).

    With a converged SCBA G^R this is the *coherent* part of transport: the
    phonon self-energy enters through the dressed G^R, but the inelastic
    (incoherent) contribution to the current is not included. Compare against
    ``Rank1KeldyshResult.I_right`` (Meir-Wingreen, which includes both) to see
    how much of the current is inelastic.
    """
    # Gam_L is nonzero only at (0,0) and Gam_R only at (Np-1,Np-1) (see
    # contact_gammas), so the trace collapses to a single element:
    #     Tr[Gam_L G Gam_R G^dag] = Gam_L[0,0] * Gam_R[-1,-1] * |G[0,-1]|^2
    # Doing this as a 4-operand einsum instead costs O(Np^4) per energy under
    # numpy's default (un-optimised) contraction order -- minutes per call at
    # NE=376, Np=135 -- so the sparse form is used whenever it applies.
    Np = G_R.shape[1]
    off = np.arange(Np) != 0
    corner_only = (not np.any(Gam_L[:, off, :]) and not np.any(Gam_L[:, :, off])
                   and not np.any(Gam_R[:, :-1, :]) and not np.any(Gam_R[:, :, :-1]))
    if corner_only:
        return np.real(Gam_L[:, 0, 0] * Gam_R[:, -1, -1] * np.abs(G_R[:, 0, -1]) ** 2)
    G_A = np.conj(np.transpose(G_R, (0, 2, 1)))
    return np.real(np.einsum('kij,kjl,klm,kmi->k', Gam_L, G_R, Gam_R, G_A,
                             optimize=True))


def tsu_esaki_current(E_grid: np.ndarray, T_E: np.ndarray, mu_L: float,
                      mu_R: float, kT: float, m_eff_kg: float,
                      area_m2: float) -> float:
    """Transverse-integrated (device) current from an elastic transmission.

    ``run_rank1_keldysh_single_bias`` returns the current of a **single spin in
    a single transverse mode**: its prefactor is q/h, the Landauer quantum per
    mode per spin. A real pixel of area A_perp >> lambda_dB carries a continuum
    of transverse modes. For a separable parabolic transverse dispersion the
    mode sum becomes an integral that can be done analytically,

        int_0^inf dE_t [f(E+E_t-mu_L) - f(E+E_t-mu_R)]
            = kT [ ln(1+e^{(mu_L-E)/kT}) - ln(1+e^{(mu_R-E)/kT}) ],

    giving the Tsu-Esaki form (docs/STACK_DECISION.md Sec. 3,
    docs/METHOD_DERIVATION.md Sec. 9.4):

        I = A_perp * (q m* kT) / (2 pi^2 hbar^3)
              * int dE T(E) [ ln(1+e^{(mu_L-E)/kT}) - ln(1+e^{(mu_R-E)/kT}) ]

    Spin degeneracy 2 is included in the prefactor. Energies (``E_grid``,
    ``mu_*``, ``kT``) are in eV; the grid may be non-uniform.

    NOTE this is exact for elastic transport only. With inelastic scattering
    the transverse integral does not factorise, because the phonon
    self-energies depend on the occupations and therefore on E_t; an exact
    treatment needs the SCBA solved per transverse energy. Check the inelastic
    fraction (see :func:`coherent_transmission`) before relying on it.
    """
    E = np.asarray(E_grid, dtype=float)

    def _supply(mu):
        x = (mu - E) / kT
        return np.where(x > 30.0, x, np.log1p(np.exp(np.clip(x, -600.0, 30.0))))

    integrand = np.asarray(T_E, dtype=float) * (_supply(mu_L) - _supply(mu_R))
    trapz = getattr(np, "trapezoid", None) or np.trapz
    integral = float(trapz(integrand, E))          # eV
    # q^3: one q for the eV->J of kT, one for dE, one for the charge itself.
    prefac = area_m2 * m_eff_kg * kT * _QE ** 3 / (2.0 * np.pi ** 2 * _HBAR_JS ** 3)
    return prefac * integral


def landauer_current_1mode(E_grid: np.ndarray, T_E: np.ndarray, mu_L: float,
                           mu_R: float, kT: float) -> float:
    """Per-mode, per-spin elastic current in A: (q^2/h) int T(E) [f_L - f_R] dE.

    Same q/h-per-spin convention as the solver's ``I_right`` and the same
    trapezoid rule in E as :func:`tsu_esaki_current`, so its E_t integral is
    exactly ``tsu_esaki_current`` -- which is what makes it usable as the
    control in :func:`transverse_integrated_current`.
    """
    E = np.asarray(E_grid, dtype=float)
    def _f(mu):  # tanh form cannot overflow at low kT
        return 0.5 * (1.0 - np.tanh(0.5 * (E - mu) / kT))
    trapz = getattr(np, "trapezoid", None) or np.trapz
    return _IE_PREFACTOR * float(trapz(np.asarray(T_E, dtype=float) * (_f(mu_L) - _f(mu_R)), E))


def transverse_integrated_current(current_at_Ef, Ef: float, kT: float,
                                  m_eff_kg: float, area_m2: float,
                                  e_max: float, n_nodes: int = 12,
                                  control=None, control_exact=None):
    """Exact transverse-mode integration of a per-mode current.

    The transverse dispersion is separable and parabolic, and the longitudinal
    Hamiltonian does not depend on the transverse energy E_t. Running the
    longitudinal problem at E_t is therefore identical to running it with BOTH
    contact chemical potentials lowered by E_t, i.e. with ``Ef -> Ef - E_t`` at
    fixed V. Hence

        I_device = A_perp * 2 * m*/(2 pi hbar^2) * int_0^inf dE_t I_1mode(Ef - E_t)

    This integrates inelastic transport **under the independent-channel
    assumption**. It does not describe transverse-mode mixing by a localized
    molecule or general polar-optical scattering. Unlike
    :func:`tsu_esaki_current`, the phonon
    self-energies depend on the contact occupations and hence on E_t, so the
    transverse integral does not otherwise factorise. Substituting the Landauer
    form for ``I_1mode`` reduces this expression analytically to Tsu-Esaki, so
    the two agree in the elastic limit.

    Parameters
    ----------
    current_at_Ef : callable
        ``current_at_Ef(Ef_shifted) -> float`` returning the per-mode, per-spin
        current in A at the shifted Fermi level (V held fixed).
    Ef, kT, e_max : float
        Fermi level, thermal energy and upper limit of the E_t integral, in eV.
        The integrand falls off once ``Ef - E_t`` drops below the transmitting
        window; check convergence in ``e_max`` as well as in ``n_nodes``.
    n_nodes : int
        Gauss-Legendre nodes on ``[0, e_max]``. Each one costs a full solve.
    control, control_exact : callable, float, optional
        Control variate. ``control(Ef_shifted)`` is a cheap per-mode current
        with the same E_t structure as the integrand, and ``control_exact`` is
        its E_t integral done analytically, already in device units (A). The
        quadrature then only sees ``current_at_Ef - control``. The result is
        exact for ANY control; a closer one just needs fewer nodes. Use
        :func:`landauer_current_1mode` on the central solve's coherent
        transmission, whose exact integral is :func:`tsu_esaki_current`. At
        low T the per-mode current is confined to a narrow E_t range with
        kT-wide edges, which a fixed Gauss-Legendre rule misses (8 nodes read
        +29 % at 10 K); the elastic control carries that structure exactly.

    Returns
    -------
    (I_device, nodes, values) : the current in A, the E_t nodes, and the
    per-mode currents at those nodes (returned so convergence can be inspected
    without re-solving).
    """
    x, w = np.polynomial.legendre.leggauss(int(n_nodes))
    nodes = 0.5 * e_max * (x + 1.0)
    weights = 0.5 * e_max * w
    values = np.array([float(current_at_Ef(Ef - float(Et))) for Et in nodes])
    prefac = area_m2 * 2.0 * m_eff_kg / (2.0 * np.pi * _HBAR_JS ** 2) * _QE
    if control is None:
        return prefac * float(np.dot(weights, values)), nodes, values
    if control_exact is None:
        raise ValueError("control requires control_exact")
    ctrl = np.array([float(control(Ef - float(Et))) for Et in nodes])
    return (float(control_exact) + prefac * float(np.dot(weights, values - ctrl)),
            nodes, values)


def fba_phonon_sigma_diag(Gl_diag, Gg_diag, chi_diag, D0_sq_per_mode,
                          hnu_idx_per_mode, N_bose_per_mode, chi_per_mode=None):
    """Diagonal-only form of :func:`fba_phonon_sigma`.

    Identical algebra, but takes and returns the (NE, Np) diagonals rather than
    full (NE, Np, Np) arrays. The phonon self-energy is diagonal in Patil's
    formulation, so nothing is lost -- and at NE=376, Np=135 a full array is
    110 MB, which makes the SCBA loop memory-bandwidth bound.
    """
    NE, Np = Gl_diag.shape
    sig_in = np.zeros((NE, Np), dtype=complex)
    sig_out = np.zeros((NE, Np), dtype=complex)
    for m, (D0_sq, inu, N) in enumerate(zip(
            D0_sq_per_mode, hnu_idx_per_mode, N_bose_per_mode)):
        if inu >= NE:
            continue
        chi_m = chi_diag if chi_per_mode is None else chi_per_mode[m]
        chi_sq = (np.asarray(chi_m) * np.asarray(chi_m)).astype(float)
        ne_l = np.zeros_like(Gl_diag); na_l = np.zeros_like(Gl_diag)
        ne_g = np.zeros_like(Gg_diag); na_g = np.zeros_like(Gg_diag)
        ne_l[: NE - inu] = Gl_diag[inu:]; na_l[inu:] = Gl_diag[: NE - inu]
        ne_g[: NE - inu] = Gg_diag[inu:]; na_g[inu:] = Gg_diag[: NE - inu]
        ce, ca = (N + 1.0) * D0_sq, N * D0_sq
        sig_in += chi_sq[None, :] * (ce * ne_l + ca * na_l)
        sig_out += chi_sq[None, :] * (ca * ne_g + ce * na_g)
    return sig_in, sig_out


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
    anderson_depth: int = 0,
    current_rtol: float | None = None,
    legacy_reproduction: bool = False,
    current_atol: float = 1e-14,
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

    Linear mixing is the default. ``anderson_depth >= 2`` enables guarded
    Anderson/DIIS: proposals with negative/nonfinite diagonal rates are rejected
    and history is reset. The final rates, populations and self-energy residual
    are checked again before reporting convergence. ``current_rtol`` optionally
    requires contact conservation as well, continuing SCBA within ``max_iter``
    even after the self-energy residual first meets ``tol``. This is useful for
    small currents, for which a global self-energy norm is an insufficient test.
    """
    if legacy_reproduction and current_rtol is not None:
        raise ValueError("legacy_reproduction reproduces the archived stopping rule; "
                         "it cannot also require current conservation")
    mixer = _ScbaMixer(mix, anderson_depth, max_iter, tol, guard=not legacy_reproduction)
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

    _diag = np.arange(Np)
    si_c = sig_in_L[:, _diag, _diag] + sig_in_R[:, _diag, _diag]
    so_c = sig_out_L[:, _diag, _diag] + sig_out_R[:, _diag, _diag]

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

        can_stop = change < tol
        if legacy_reproduction:
            # Archived order: mix first, then stop on the pre-mix change. The
            # returned Sigma is therefore one mixing step past the tested one.
            si, so = mixer.update(sig_in_ph[:, _diag, _diag], sig_out_ph[:, _diag, _diag],
                                  sig_in_new[:, _diag, _diag], sig_out_new[:, _diag, _diag])
            sig_in_ph[:, _diag, _diag] = si
            sig_out_ph[:, _diag, _diag] = so
            if can_stop:
                converged = True
                break
            continue
        if can_stop and current_rtol is not None:
            nd = G_lesser[:, _diag, _diag].real
            gd = G_R[:, _diag, _diag]
            pd = (1j*(gd-gd.conj())).real - nd
            can_stop = _current_ok(*_contact_current_pair(nd, pd, si_c, so_c, E_grid[1]-E_grid[0]),
                                   current_rtol, current_atol)
        if can_stop:
            converged = True
            break
        si, so = mixer.update(sig_in_ph[:, _diag, _diag], sig_out_ph[:, _diag, _diag],
                              sig_in_new[:, _diag, _diag], sig_out_new[:, _diag, _diag])
        sig_in_ph[:, _diag, _diag] = si
        sig_out_ph[:, _diag, _diag] = so

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

    si_check, so_check = fba_phonon_sigma_diag(G_lesser[:, _diag, _diag], G_greater[:, _diag, _diag], chi_diag,
        D0_sq_per_mode, hnu_idx_per_mode, N_bose_per_mode, chi_per_mode)
    residual = _relative_sigma_residual(sig_in_ph[:, _diag, _diag], sig_out_ph[:, _diag, _diag], si_check, so_check)
    physical = _nonnegative_finite(sig_in_ph[:, _diag, _diag], sig_out_ph[:, _diag, _diag], n_arr[:, _diag, _diag], p_arr[:, _diag, _diag])
    if not legacy_reproduction:
        converged = bool(converged and physical and residual < tol
                         and _current_ok(I_left, I_right, current_rtol, current_atol))
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
        converged=converged, physical_valid=physical, self_energy_residual=residual,
        mixer_fallbacks=mixer.fallbacks,
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


def run_rank1_keldysh_single_bias_fast(
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
    anderson_depth: int = 0,
    current_rtol: float | None = None,
    current_atol: float = 1e-14,
) -> Rank1KeldyshResult:
    """Numerically identical to :func:`run_rank1_keldysh_single_bias`, ~5x faster.

    The original is memory-bandwidth bound, not compute bound. At NE=376,
    Np=135 one (NE, Np, Np) complex array is 110 MB, and the loop builds about
    ten of them per iteration; a line-level profile put 23% of the time in
    *adding* three such arrays together (``sig_in_L + sig_in_R + sig_in_ph``)
    and only 14% in the linear solve.

    The waste is structural: ``Gam_L`` is nonzero only at (0,0), ``Gam_R`` only
    at (Np-1, Np-1), and the phonon self-energy is diagonal, so every Sigma in
    the loop is a *diagonal* matrix stored as a full one. And G^< / G^> are used
    only through their diagonals (for Sigma_ph) and two corner elements (for the
    current). So the loop here carries (NE, Np) diagonals, using

        diag(G Sigma G^dagger)_z = sum_w |G_zw|^2 Sigma_ww    (Sigma diagonal)

    and materialises the full G^< / G^> exactly once, after convergence, so the
    returned object is unchanged for downstream consumers (extract_density_1d,
    analytic d2I/dV2, coherent_transmission).

    Two further savings: the retarded denominator ``E - H - Sigma_L - Sigma_R``
    is bias-independent, so it is built once and only its diagonal updated per
    iteration; and ``|G|^2`` is formed once and reused for both G^< and G^>.
    """
    mixer = _ScbaMixer(mix, anderson_depth, max_iter, tol)
    Np = H_z.shape[0]
    NE = E_grid.size
    U_bias = bias_profile
    idx = np.arange(Np)

    mu_L, mu_R = Ef + V / 2.0, Ef - V / 2.0
    f_L = 1.0 / (1.0 + np.exp((E_grid - mu_L) / kT))
    f_R = 1.0 / (1.0 + np.exp((E_grid - mu_R) / kT))

    Gam_L, Gam_R = contact_gammas(E_grid, U_bias, UB, t0, Np, eta=eta)
    sig_in_L = f_L[:, None, None] * Gam_L
    sig_in_R = f_R[:, None, None] * Gam_R
    sig_out_L = (1.0 - f_L)[:, None, None] * Gam_L
    sig_out_R = (1.0 - f_R)[:, None, None] * Gam_R

    # Diagonal (NE, Np) views of the contact in/out scattering.
    gamL_d = np.zeros((NE, Np), dtype=complex); gamL_d[:, 0] = Gam_L[:, 0, 0]
    gamR_d = np.zeros((NE, Np), dtype=complex); gamR_d[:, -1] = Gam_R[:, -1, -1]
    sig_in_L_d = f_L[:, None] * gamL_d
    sig_in_R_d = f_R[:, None] * gamR_d
    sig_out_L_d = (1.0 - f_L)[:, None] * gamL_d
    sig_out_R_d = (1.0 - f_R)[:, None] * gamR_d

    # Retarded denominator without the phonon term: bias-independent, built once.
    z_eta = 1j * eta
    E_arr = np.asarray(E_grid, dtype=complex) + z_eta
    H_full = H_z + np.diag(U_bias)
    ck_L = 1.0 - ((E_arr - U_bias[0] - UB[0]) / (2.0 * t0))
    ck_R = 1.0 - ((E_arr - U_bias[-1] - UB[-1]) / (2.0 * t0))
    sigL_c = -t0 * np.exp(1j * np.arccos(ck_L))
    sigR_c = -t0 * np.exp(1j * np.arccos(ck_R))
    M_base = np.broadcast_to(-H_full.astype(complex), (NE, Np, Np)).copy()
    M_base[:, idx, idx] += E_arr[:, None]
    M_base[:, 0, 0] -= sigL_c
    M_base[:, -1, -1] -= sigR_c
    eye_b = np.broadcast_to(np.eye(Np, dtype=complex), (NE, Np, Np)).copy()

    si_ph = np.zeros((NE, Np), dtype=complex)
    so_ph = np.zeros((NE, Np), dtype=complex)
    iters_used = 0
    converged = False
    G_R = None
    for it in range(max_iter):
        iters_used = it + 1
        M = M_base.copy()
        M[:, idx, idx] += 0.5j * (si_ph + so_ph)
        G_R = _robust_solve(M, eye_b)

        absG2 = np.abs(G_R)
        np.square(absG2, out=absG2)              # |G_zw|^2, reused for G< and G>
        si_tot = sig_in_L_d + sig_in_R_d + si_ph
        so_tot = sig_out_L_d + sig_out_R_d + so_ph
        Gl_d = np.einsum('kzw,kw->kz', absG2, si_tot, optimize=True)
        Gg_d = np.einsum('kzw,kw->kz', absG2, so_tot, optimize=True)

        si_new, so_new = fba_phonon_sigma_diag(
            Gl_d, Gg_d, chi_diag, D0_sq_per_mode, hnu_idx_per_mode,
            N_bose_per_mode, chi_per_mode)

        abs_change = float(np.sum(np.abs(si_new - si_ph))) + float(np.sum(np.abs(so_new - so_ph)))
        norm_old = float(np.sum(np.abs(si_ph))) + float(np.sum(np.abs(so_ph)))
        change = abs_change / max(norm_old, 1e-30)

        can_stop = change < tol
        if can_stop and current_rtol is not None:
            gd = G_R[:, idx, idx]
            pd = (1j*(gd-gd.conj())).real - Gl_d.real
            can_stop = _current_ok(*_contact_current_pair(Gl_d.real, pd,
                sig_in_L_d+sig_in_R_d, sig_out_L_d+sig_out_R_d, E_grid[1]-E_grid[0]),
                current_rtol, current_atol)
        if can_stop:
            converged = True
            break
        si_ph, so_ph = mixer.update(si_ph, so_ph, si_new, so_new)

    # FINAL PASS: rebuild G with the current Sigma and materialise the full
    # G^< / G^> once, so the returned object matches the reference exactly.
    M = M_base.copy()
    M[:, idx, idx] += 0.5j * (si_ph + so_ph)
    G_R = _robust_solve(M, eye_b)
    G_A = np.conj(G_R.transpose(0, 2, 1))
    sig_in_ph = np.zeros((NE, Np, Np), dtype=complex)
    sig_out_ph = np.zeros((NE, Np, Np), dtype=complex)
    sig_in_ph[:, idx, idx] = si_ph
    sig_out_ph[:, idx, idx] = so_ph
    sig_in_total = sig_in_L + sig_in_R + sig_in_ph
    sig_out_total = sig_out_L + sig_out_R + sig_out_ph
    G_lesser = np.matmul(np.matmul(G_R, sig_in_total), G_A)
    G_greater = np.matmul(np.matmul(G_R, sig_out_total), G_A)

    n_arr = np.real(G_lesser)
    A_spec = 1j * (G_R - G_A)
    p_arr = np.real(A_spec) - n_arr
    dE = E_grid[1] - E_grid[0]
    # Sigma_R is nonzero only at (Np-1, Np-1) and Sigma_L only at (0,0), so
    # Tr[Sigma X] = Sigma_cc X_cc -- identical to the reference's dense trace.
    I1 = float(np.sum(np.real(sig_out_R[:, -1, -1] * n_arr[:, -1, -1]
                              - sig_in_R[:, -1, -1] * p_arr[:, -1, -1])))
    I2 = float(np.sum(np.real(sig_out_L[:, 0, 0] * n_arr[:, 0, 0]
                              - sig_in_L[:, 0, 0] * p_arr[:, 0, 0])))
    si_check, so_check = fba_phonon_sigma_diag(G_lesser[:, idx, idx], G_greater[:, idx, idx], chi_diag,
        D0_sq_per_mode, hnu_idx_per_mode, N_bose_per_mode, chi_per_mode)
    residual = _relative_sigma_residual(si_ph, so_ph, si_check, so_check)
    physical = _nonnegative_finite(si_ph, so_ph, n_arr[:, idx, idx], p_arr[:, idx, idx])
    converged = bool(converged and physical and residual < tol
                     and _current_ok(I2*dE*_IE_PREFACTOR, I1*dE*_IE_PREFACTOR, current_rtol, current_atol))
    return Rank1KeldyshResult(
        V=V, E_grid=E_grid, G_R=G_R, G_lesser=G_lesser, G_greater=G_greater,
        sigma_in_ph=sig_in_ph, sigma_out_ph=sig_out_ph,
        Gam_L=Gam_L, Gam_R=Gam_R,
        I_left=I2 * dE * _IE_PREFACTOR, I_right=I1 * dE * _IE_PREFACTOR,
        iters_used=iters_used, converged=converged, physical_valid=physical, self_energy_residual=residual,
        mixer_fallbacks=mixer.fallbacks,
    )


@dataclass
class Rank1DiagResult:
    """What :func:`run_rank1_keldysh_single_bias_lowmem` keeps: diagonals only."""
    V: float
    E_grid: np.ndarray
    n_diag: np.ndarray      # (NE, Np)  Re diag G^<
    p_diag: np.ndarray      # (NE, Np)  Re diag A - n (holes), as the fast solver uses
    sigma_in_ph: np.ndarray   # (NE, Np)  converged phonon in/out-scattering, to
    sigma_out_ph: np.ndarray  #           warm-start a nearby solve (sigma_init)
    I_left: float
    I_right: float
    iters_used: int
    converged: bool
    physical_valid: bool = True
    self_energy_residual: float = float("inf")
    mixer_fallbacks: int = 0


def run_rank1_keldysh_single_bias_lowmem(
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
    anderson_depth: int = 0,
    current_rtol: float | None = None,
    current_atol: float = 1e-14,
    chunk: int = 64,
    sigma_init: tuple[np.ndarray, np.ndarray] | None = None,
    method: str = "dense",
) -> Rank1DiagResult:
    """The SCBA of :func:`run_rank1_keldysh_single_bias_fast` in O(chunk * Np^2) memory.

    The fast solver's loop already works on (NE, Np) diagonals, but it holds
    every energy's (Np, Np) matrices at once (denominator, identity, G^R, |G|^2)
    and materialises full G^< / G^> at the end: ~9.6 MB per energy at Np = 435,
    so a 0.5 meV grid over a 0.6 eV window needs > 11 GB. Here each iteration
    sweeps the energies in chunks and keeps only diagonals, and the result
    carries only what the density and current need. Same arithmetic, same
    Anderson mixing, same final pass.

    ``sigma_init`` = (sigma_in_ph, sigma_out_ph) from a previous result on the
    same energy grid starts the iteration there instead of at zero, e.g. inside
    an outer Poisson loop where U changes little between calls.

    ``method="rgf"`` gets the same diagonals by recursive Green's functions,
    O(Np) per energy instead of a dense O(Np^3) solve (~100x at Np = 435).
    It needs a tridiagonal Hamiltonian with one hopping value, which
    build_stack produces; every self-energy here is diagonal.
    """
    mixer = _ScbaMixer(mix, anderson_depth, max_iter, tol)
    Np = H_z.shape[0]
    NE = E_grid.size
    idx = np.arange(Np)
    mu_L, mu_R = Ef + V / 2.0, Ef - V / 2.0
    f_L = 1.0 / (1.0 + np.exp((E_grid - mu_L) / kT))
    f_R = 1.0 / (1.0 + np.exp((E_grid - mu_R) / kT))

    E_arr = np.asarray(E_grid, dtype=complex) + 1j * eta
    U_bias = bias_profile
    H_full = (H_z + np.diag(U_bias)).astype(complex)
    ck_L = 1.0 - ((E_arr - U_bias[0] - UB[0]) / (2.0 * t0))
    ck_R = 1.0 - ((E_arr - U_bias[-1] - UB[-1]) / (2.0 * t0))
    sigL_c = -t0 * np.exp(1j * np.arccos(ck_L))
    sigR_c = -t0 * np.exp(1j * np.arccos(ck_R))
    gamL = 1j * (sigL_c - np.conj(sigL_c))      # Gamma_L(0,0), as contact_gammas
    gamR = 1j * (sigR_c - np.conj(sigR_c))
    si_c = np.zeros((NE, Np), dtype=complex); so_c = np.zeros((NE, Np), dtype=complex)
    si_c[:, 0] += f_L * gamL; so_c[:, 0] += (1 - f_L) * gamL
    si_c[:, -1] += f_R * gamR; so_c[:, -1] += (1 - f_R) * gamR

    def sweep_dense(si_ph, so_ph):
        """Diagonals of G^<, G^> and G^R for the current phonon self-energy."""
        Gl = np.empty((NE, Np), dtype=complex)
        Gg = np.empty((NE, Np), dtype=complex)
        GRd = np.empty((NE, Np), dtype=complex)
        for a in range(0, NE, chunk):
            b = min(a + chunk, NE)
            M = np.broadcast_to(-H_full, (b - a, Np, Np)).copy()
            M[:, idx, idx] += E_arr[a:b, None] + 0.5j * (si_ph[a:b] + so_ph[a:b])
            M[:, 0, 0] -= sigL_c[a:b]
            M[:, -1, -1] -= sigR_c[a:b]
            G = _robust_solve(M, np.broadcast_to(np.eye(Np, dtype=complex), (b - a, Np, Np)))
            GRd[a:b] = G[:, idx, idx]
            G2 = np.abs(G); np.square(G2, out=G2)
            Gl[a:b] = np.einsum("kzw,kw->kz", G2, si_c[a:b] + si_ph[a:b], optimize=True)
            Gg[a:b] = np.einsum("kzw,kw->kz", G2, so_c[a:b] + so_ph[a:b], optimize=True)
        return Gl, Gg, GRd

    hop = -float(H_z[0, 1])                    # M = E - H - Sigma has off-diagonal +hop
    if method == "rgf" and not (np.allclose(np.diag(H_z, 1), -hop) and np.allclose(H_z, H_z.T)
                                and np.count_nonzero(np.triu(H_z, 2)) == 0):
        raise ValueError("method='rgf' needs a symmetric tridiagonal H with one hopping value")
    h_d = np.diag(H_full)

    def sweep_rgf(si_ph, so_ph):
        """Same diagonals as sweep_dense, by recursion along the chain.

        M = E - H - Sigma is complex symmetric and tridiagonal (hopping +hop),
        with right/left-connected g^R_k = 1/(a_k - hop^2 g^R_{k+1}) etc., so
            G_kk = 1 / (a_k - hop^2 g^L_{k-1} - hop^2 g^R_{k+1}),
            G_zw = G_zz prod_{k=z+1..w} (-hop g^R_k)   (z < w),  G_wz = G_zw.
        With q_k = |hop g^R_k|^2, sum_w |G_zw|^2 s_w splits into
            |G_zz|^2 (s_z + S_z) + T_z,
            S_z = q_{z+1} (s_{z+1} + S_{z+1}),   T_z = q_z (T_{z-1} + |G_{z-1}|^2 s_{z-1}),
        two O(Np) recursions instead of forming |G|^2.
        """
        # Recursion advances along z and reads an entire energy column each
        # step. Keep those columns contiguous on fine energy grids.
        a = np.asfortranarray(E_arr[:, None] - h_d[None, :] + 0.5j * (si_ph + so_ph))
        a[:, 0] -= sigL_c
        a[:, -1] -= sigR_c
        h2 = hop * hop
        gR = np.empty_like(a); gL = np.empty_like(a)
        gR[:, -1] = 1.0 / a[:, -1]
        for k in range(Np - 2, -1, -1):
            gR[:, k] = 1.0 / (a[:, k] - h2 * gR[:, k + 1])
        gL[:, 0] = 1.0 / a[:, 0]
        for k in range(1, Np):
            gL[:, k] = 1.0 / (a[:, k] - h2 * gL[:, k - 1])
        dl = np.zeros_like(a); dr = np.zeros_like(a)
        dl[:, 1:] = h2 * gL[:, :-1]
        dr[:, :-1] = h2 * gR[:, 1:]
        GRd = 1.0 / (a - dl - dr)
        G2 = np.abs(GRd) ** 2
        q = np.abs(hop * gR) ** 2

        def contract(src):
            src = np.asfortranarray(src)
            S = np.zeros_like(a, dtype=float); T = np.zeros_like(a, dtype=float)
            for z in range(Np - 2, -1, -1):
                S[:, z] = q[:, z + 1] * (src[:, z + 1] + S[:, z + 1])
            for z in range(1, Np):
                T[:, z] = q[:, z] * (T[:, z - 1] + G2[:, z - 1] * src[:, z - 1])
            return G2 * (src + S) + T

        # in/out-scattering rates are real (Gamma * f, D^2 * occupation)
        Gl = contract(np.real(si_c + si_ph)).astype(complex)
        Gg = contract(np.real(so_c + so_ph)).astype(complex)
        return Gl, Gg, GRd

    if method not in ("dense", "rgf"):
        raise ValueError(f"method must be 'dense' or 'rgf', not {method!r}")
    sweep = sweep_rgf if method == "rgf" else sweep_dense

    if sigma_init is None:
        si_ph = np.zeros((NE, Np), dtype=complex)
        so_ph = np.zeros((NE, Np), dtype=complex)
    else:
        si_ph = np.array(sigma_init[0], dtype=complex)
        so_ph = np.array(sigma_init[1], dtype=complex)
    if (si_ph.shape != (NE, Np) or so_ph.shape != (NE, Np)
            or not _nonnegative_finite(si_ph, so_ph)
            or np.min(si_ph.real) < 0 or np.min(so_ph.real) < 0):
        raise ValueError("SCBA warm start must contain finite nonnegative rates on this grid")

    iters_used, converged = 0, False
    for it in range(max_iter):
        iters_used = it + 1
        Gl, Gg, GRd = sweep(si_ph, so_ph)
        si_new, so_new = fba_phonon_sigma_diag(Gl, Gg, chi_diag, D0_sq_per_mode,
                                               hnu_idx_per_mode, N_bose_per_mode, chi_per_mode)
        change = ((float(np.sum(np.abs(si_new - si_ph))) + float(np.sum(np.abs(so_new - so_ph))))
                  / max(float(np.sum(np.abs(si_ph))) + float(np.sum(np.abs(so_ph))), 1e-30))
        can_stop = change < tol
        if can_stop and current_rtol is not None:
            pd = (1j*(GRd-GRd.conj())).real - Gl.real
            can_stop = _current_ok(*_contact_current_pair(Gl.real, pd, si_c, so_c, E_grid[1]-E_grid[0]),
                                   current_rtol, current_atol)
        if can_stop:
            converged = True
            break
        si_ph, so_ph = mixer.update(si_ph, so_ph, si_new, so_new)

    # FINAL PASS, as in the fast solver: observables from the final Sigma.
    Gl, Gg, GRd = sweep(si_ph, so_ph)
    n = np.real(Gl)
    p = np.real(1j * (GRd - np.conj(GRd))) - n
    dE = E_grid[1] - E_grid[0]
    I1 = float(np.sum(np.real(so_c[:, -1] * n[:, -1] - si_c[:, -1] * p[:, -1])))
    I2 = float(np.sum(np.real(so_c[:, 0] * n[:, 0] - si_c[:, 0] * p[:, 0])))
    si_check, so_check = fba_phonon_sigma_diag(Gl, Gg, chi_diag,
        D0_sq_per_mode, hnu_idx_per_mode, N_bose_per_mode, chi_per_mode)
    residual = _relative_sigma_residual(si_ph, so_ph, si_check, so_check)
    physical = _nonnegative_finite(si_ph, so_ph, n, p)
    converged = bool(converged and physical and residual < tol
                     and _current_ok(I2*dE*_IE_PREFACTOR, I1*dE*_IE_PREFACTOR, current_rtol, current_atol))
    return Rank1DiagResult(V=V, E_grid=E_grid, n_diag=n, p_diag=p,
                           sigma_in_ph=si_ph, sigma_out_ph=so_ph,
                           I_left=I2 * dE * _IE_PREFACTOR, I_right=I1 * dE * _IE_PREFACTOR,
                           iters_used=iters_used, converged=converged, physical_valid=physical, self_energy_residual=residual,
        mixer_fallbacks=mixer.fallbacks)
