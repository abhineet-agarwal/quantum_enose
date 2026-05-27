"""
Analytic d²I/dV² — SISPAD Tier-1 item α.

Replaces the numerical double ``np.gradient`` over the V-grid in
``run/run_iets.py`` with closed-form expressions derived by differentiating
the Meir-Wingreen current symbolically under the Patil symmetric-bias
convention μ_{L,R}(V) = E_F ± V/2.

The full derivation (and a TODO sketch for the inelastic path) lives in
``docs/ANALYTIC_D2IDV2.md``. This module implements the **coherent
(Landauer) path** analytically:

    d²I_coh/dV² = (e/h) · ∫ dE · T(E) · K_therm(E, V)

    K_therm(E, V) = (1/(4 (kT)²)) · [f_L(1-f_L)(1-2f_L)
                                     - f_R(1-f_R)(1-2f_R)],

with f_α = 1/(1 + exp[(E - μ_α(V))/kT]). The inelastic piece requires the
rank-1 Keldysh wrapper and is out of scope here; see §5 and §7 of
``docs/ANALYTIC_D2IDV2.md``.

All functions accept NumPy arrays and broadcast across (E, V) for vectorised
evaluation over an arbitrary V-grid.
"""

from __future__ import annotations

import numpy as np

# Physical constants (SI) for the (e/h) prefactor in d²I/dV².
_E_CHARGE = 1.602176634e-19    # C
_H_PLANCK = 6.62607015e-34     # J·s
_E_OVER_H = _E_CHARGE / _H_PLANCK  # Hz/V == A/(V²·?) -- see docstrings

# We work with E and V in eV below, so the appropriate current prefactor for
# an integral ∫ dE [eV] · K_therm [1/eV²] · T [dimensionless] is
#   (e/h) · e = e²/h  [Siemens] → current per volt squared when the
# kernel is d²/dV² with V in volts. Since V is in volts and E is in eV
# (i.e. we do the integral in eV), we use the G0 = 2e²/h quantum of
# conductance /volt² ×e = e²·e/(h·V²) combination below. In practice we
# return the unitful prefactor once and let the caller take its time.
# See `coherent_d2idv2_prefactor()`.


# ---------------------------------------------------------------------------
# Fermi function and its E-derivatives (analytic).
# ---------------------------------------------------------------------------


def fermi(E, mu, kT):
    """Fermi-Dirac distribution f(E, μ) = 1/(1 + exp((E-μ)/kT)).

    Numerically stable for |E - μ| ≫ kT via the where/expit trick.
    """
    x = np.asarray((np.asarray(E) - np.asarray(mu)) / kT, dtype=float)
    # Stable: for x ≥ 0 use e^{-x}/(1+e^{-x}); for x < 0 use 1/(1+e^x).
    out = np.empty_like(x)
    pos = x >= 0
    ex_neg = np.exp(-x[pos])
    out[pos] = ex_neg / (1.0 + ex_neg)
    out[~pos] = 1.0 / (1.0 + np.exp(x[~pos]))
    return out


def fermi_deriv_E(E, mu, kT):
    """∂f/∂E = -(1/kT) · f · (1 - f)."""
    f = fermi(E, mu, kT)
    return -(1.0 / kT) * f * (1.0 - f)


def fermi_second_deriv_E(E, mu, kT):
    """∂²f/∂E² = +(1/(kT)²) · f · (1 - f) · (1 - 2f)."""
    f = fermi(E, mu, kT)
    return (1.0 / kT ** 2) * f * (1.0 - f) * (1.0 - 2.0 * f)


# ---------------------------------------------------------------------------
# Thermal d² kernel (coherent path).
# ---------------------------------------------------------------------------


def thermal_d2_kernel(E, V, kT, E_F=0.0):
    """Closed-form thermal kernel K_therm(E, V) from Eq. (11) of
    ``docs/ANALYTIC_D2IDV2.md``:

        K_therm(E, V) = (1/4) · [f_L''(E) - f_R''(E)]
                      = (1/(4 (kT)²)) · [f_L(1-f_L)(1-2f_L)
                                         - f_R(1-f_R)(1-2f_R)],

    with μ_L = E_F + V/2, μ_R = E_F - V/2 (Patil symmetric bias).

    Parameters
    ----------
    E : ndarray, shape (NE,)
    V : float or ndarray, shape (NV,)
    kT : float
        Thermal energy in the same units as E (eV).
    E_F : float, optional
        Equilibrium Fermi level (default 0).

    Returns
    -------
    K : ndarray, shape broadcast(E, V)
        The kernel evaluated at every (E, V). For E shape (NE,) and V shape
        (NV,), returns shape (NV, NE). Scalar V returns shape (NE,).
    """
    E_arr = np.asarray(E, dtype=float)
    V_arr = np.asarray(V, dtype=float)
    scalar_V = (V_arr.ndim == 0)
    V_arr = np.atleast_1d(V_arr)

    # Broadcast: (NV, 1) vs (NE,) → (NV, NE)
    mu_L = (E_F + 0.5 * V_arr)[:, None]
    mu_R = (E_F - 0.5 * V_arr)[:, None]
    Ebr = E_arr[None, :]

    d2_fL = fermi_second_deriv_E(Ebr, mu_L, kT)
    d2_fR = fermi_second_deriv_E(Ebr, mu_R, kT)
    K = 0.25 * (d2_fL - d2_fR)

    return K[0] if scalar_V else K


# ---------------------------------------------------------------------------
# Coherent d²I/dV² (Landauer path).
# ---------------------------------------------------------------------------


def coherent_d2idv2_prefactor():
    """Return the physical prefactor e²/h in SI (S = A/V).

    The convention used here: with E and V in volts (so E_F, μ, V all in V),
    the integral ∫ dE K(E,V) T(E) is dimensionless in V, and

        d²I/dV²  [A/V²]  =  (e²/h)  ·  ∫ dE [V] · T(E) · K(E, V) [1/V²]

    where the [1/V²] in K comes from the double ∂/∂V on the Fermi step.
    Units check: [A/V] = [S] = [e²/h]; dividing by V once more for the
    second derivative gives A/V². Good.

    If E and V are in *eV* (common in this codebase — kT = 0.02585 eV at
    300 K), then E and V are numerically equal to their values in volts
    because 1 eV = 1 e·V, and the prefactor is still e²/h in SI.
    """
    return _E_CHARGE ** 2 / _H_PLANCK


def analytic_d2idv2_coherent(T_of_E, E_grid, V_grid, kT, E_F=0.0):
    """Closed-form d²I/dV² for the coherent (Landauer) current.

        d²I_coh/dV²(V) = (e²/h) · ∫ dE T(E) · K_therm(E, V)

    The integral is done by trapezoidal quadrature on ``E_grid``.

    Parameters
    ----------
    T_of_E : ndarray, shape (NE,)
        Elastic transmission T(E) = Tr[Γ_L G₀^R Γ_R G₀^A]. V-independent
        (no Poisson self-consistency assumed).
    E_grid : ndarray, shape (NE,)
        Energy grid in eV (same units as kT, V_grid, E_F).
    V_grid : ndarray, shape (NV,)
        Bias grid in V (numerically equal to eV).
    kT : float
        Thermal energy in eV.
    E_F : float, optional
        Equilibrium Fermi level in eV (default 0).

    Returns
    -------
    d2I : ndarray, shape (NV,)
        d²I_coh/dV² in A/V² (SI).

    Notes
    -----
    Cost: one K_therm evaluation per (V, E) plus a single np.trapezoid per V.
    For NV = 201 and NE = 500, this is ~10⁵ ops total — subsecond.
    """
    E_arr = np.asarray(E_grid, dtype=float)
    V_arr = np.asarray(V_grid, dtype=float)
    T_arr = np.asarray(T_of_E, dtype=float)

    if T_arr.shape != E_arr.shape:
        raise ValueError(
            f"T_of_E shape {T_arr.shape} must match E_grid shape {E_arr.shape}"
        )

    # K has shape (NV, NE); integrate over the E axis.
    K = thermal_d2_kernel(E_arr, V_arr, kT, E_F=E_F)
    integrand = K * T_arr[None, :]
    integral = np.trapezoid(integrand, x=E_arr, axis=1)  # shape (NV,)

    return coherent_d2idv2_prefactor() * integral


# ---------------------------------------------------------------------------
# Inelastic d²I/dV² (FBA path) — closes SISPAD Tier-1 α inelastic.
# ---------------------------------------------------------------------------


def _inelastic_kernels(G_R, Gam_L, Gam_R):
    """Per-energy V-independent kernels for the FBA d²I_R/dV² formula.

    Returns three real (NE,) arrays:

        T_RL(E) = Tr[Γ_R · G^R Γ_L G^A]            (transmission, real & ≥ 0)
        T_RR(E) = Tr[Γ_R · G^R Γ_R G^A]            (right-self loop, real)
        T_RA(E) = Tr[Γ_R · A]   with A = i(G^R - G^A)

    These are the kernels appearing in

        I_R(V) = ∫ dE [f_L(E,V) T_RL(E)
                       + f_R(E,V) (T_RR(E) - T_RA(E))
                       + (V-independent constant from σ^<_ph)]

    derived from Patil's Meir-Wingreen trace at the right contact:
    I_R = ∫ dE Tr[(1-f_R)Γ_R · n - f_R Γ_R · p] with n = real(G^R Σ_in G^A)
    and p = A - n. The σ^<_ph piece of Σ_in is V-independent at FBA, so it
    contributes only a constant to I_R and drops out of d²I/dV².

    Optical-theorem cross-check (no phonons): with Γ_L + Γ_R = i(Σ^R - Σ^A)
    and A = G^R(Γ_L+Γ_R)G^A, T_RA = T_RL + T_RR, hence T_RR - T_RA = -T_RL
    and the formula collapses to (1/4) ∫ dE T_RL · (f_L'' - f_R''), which is
    the coherent-path expression in Eq. (10) of docs/ANALYTIC_D2IDV2.md.
    """
    G_A = np.conj(G_R.transpose(0, 2, 1))
    A_spec = 1j * (G_R - G_A)
    NE = G_R.shape[0]
    T_RL = np.empty(NE)
    T_RR = np.empty(NE)
    T_RA = np.empty(NE)
    for k in range(NE):
        GR_GamL_GA = G_R[k] @ Gam_L[k] @ G_A[k]
        GR_GamR_GA = G_R[k] @ Gam_R[k] @ G_A[k]
        T_RL[k] = float(np.real(np.trace(Gam_R[k] @ GR_GamL_GA)))
        T_RR[k] = float(np.real(np.trace(Gam_R[k] @ GR_GamR_GA)))
        T_RA[k] = float(np.real(np.trace(Gam_R[k] @ A_spec[k])))
    return T_RL, T_RR, T_RA


def analytic_d2idv2_inelastic_at_bias(result, kT, E_F=0.0):
    """Closed-form d²I_R/dV² at a single bias from a converged Keldysh state.

    Implements the FBA inelastic path of §5 / §7 of
    ``docs/ANALYTIC_D2IDV2.md``:

        d²I_R/dV²(V) = (e²/h) · (1/4) · ∫ dE
            [ (∂²f_L/∂E²)(E,V) · T_RL(E)
            + (∂²f_R/∂E²)(E,V) · (T_RR(E) - T_RA(E)) ]

    The kernels T_RL, T_RR, T_RA are built from the converged G^R, Γ_L, Γ_R
    inside ``result``; the phonon broadening enters G^R through the
    ``0.5j · (σ^<_ph + σ^>_ph)`` lump, which is what gives the IETS satellite
    structure in T_RL(E) and ultimately the d²I/dV² peaks at V ≈ ℏω/e once
    f_L'' (peaked at E = μ_L = E_F + V/2) sweeps over those satellites.

    Parameters
    ----------
    result : Rank1KeldyshResult
        Output of ``run_rank1_keldysh_single_bias`` at the bias V of interest.
        Must carry G_R, Gam_L, Gam_R (added 2026-04-09).
    kT : float
        Thermal energy in eV.
    E_F : float, optional
        Equilibrium Fermi level in eV.

    Returns
    -------
    d2I : float
        d²I_R/dV² in A/V² (SI), evaluated at ``result.V``.

    Notes
    -----
    For a sweep over V, call this once per converged ``result`` — the cost is
    one O(NE · Np³) kernel build. The σ^<_ph term in I_R is V-independent at
    FBA and so contributes nothing to d²I/dV² — see the proof in
    ``_inelastic_kernels`` and §5 of the doc.
    """
    V = float(result.V)
    E_arr = np.asarray(result.E_grid, dtype=float)
    G_R = result.G_R
    Gam_L = result.Gam_L
    Gam_R = result.Gam_R

    T_RL, T_RR, T_RA = _inelastic_kernels(G_R, Gam_L, Gam_R)

    mu_L = E_F + 0.5 * V
    mu_R = E_F - 0.5 * V
    f_L_dd = fermi_second_deriv_E(E_arr, mu_L, kT)
    f_R_dd = fermi_second_deriv_E(E_arr, mu_R, kT)

    integrand = 0.25 * (f_L_dd * T_RL + f_R_dd * (T_RR - T_RA))
    integral = float(np.trapezoid(integrand, x=E_arr))
    return coherent_d2idv2_prefactor() * integral
