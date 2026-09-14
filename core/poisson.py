"""
1D Poisson solver and NEGF density extraction for self-consistent transport.

This module is deliberately dependency-light (numpy only) and free of any
device/solver coupling so it can be unit-tested in isolation. The
self-consistent outer loop that ties it to the rank-1 Keldysh SCBA lives in
:mod:`core.poisson_negf`.

Conventions
-----------
We solve for ``U(z)`` = the **electron potential energy** in eV — exactly the
quantity that is added to the tight-binding Hamiltonian diagonal (the
``bias_profile`` argument of :func:`core.scba_rank1_keldysh.run_rank1_keldysh_single_bias`).

The electrostatic Poisson equation in SI is ``∇·(ε ∇φ) = −ρ`` with
``ρ = q(N_D − n)`` (ionized donors positive, electrons negative). With the
electron potential energy ``U = −qφ`` (so numerically ``U[eV] = −φ[V]``) this
becomes

    d/dz [ ε(z) dU/dz ] = q ( N_D(z) − n(z) )                      (★)

which is what :func:`solve_poisson_1d` discretizes. The source-term sign makes
the self-consistency *stabilizing*: extra electrons → RHS more negative → U
bulges up locally → electrons are pushed away (negative feedback).

(The archived ``quantum_enose_poisson/core/poisson.py`` solved the same matrix
but then divided the volt-valued solution by ``q`` — a unit bug that inflated
the potential by ~1e19. That step is intentionally absent here.)
"""
from __future__ import annotations

import numpy as np

# Physical constants (SI)
_Q = 1.602176634e-19        # C
_EPS0 = 8.854187817e-12     # F/m
_HBAR = 1.054571817e-34     # J·s


def physical_transverse_density(
    G_R: np.ndarray,
    Gam_L: np.ndarray,
    Gam_R: np.ndarray,
    E_grid: np.ndarray,
    mu_L: float,
    mu_R: float,
    kT: float,
    m_eff_kg: float,
    a_m: float,
) -> np.ndarray:
    """First-principles transverse-integrated electron density (no fit constant).

    Replaces the single calibrated prefactor ``C`` with the physical transverse
    integration (Lake 1997 / Datta). Transverse motion is a parabolic 2D
    continuum: each longitudinal energy ``E`` carries a subband, and integrating
    ``∫d²k⊥/(2π)² = (m*/2πħ²)∫dE⊥`` over the contact filling gives the analytic
    supply function ``∫dE⊥ f(E+E⊥) = kT·ln(1+e^{(μ−E)/kT})``. Hence

        n(z) = [2 m* kT / (2π ħ² a)] · (1/2π) ∫dE  Σ_{c=L,R}
                   Re[Gᴿ Γ_c Gᴬ]_zz(E) · ln(1 + e^{(μ_c − E)/kT})

    The leading factor 2 is spin. Energies (``E_grid``, ``mu_*``, ``kT``) are in
    eV; ``kT`` also appears in SI (× q) in the areal-DOS prefactor.

    Parameters
    ----------
    G_R : ndarray (NE, Np, Np) complex   — retarded Green's function
    Gam_L, Gam_R : ndarray (NE, Np, Np)  — contact broadening functions
    E_grid : ndarray (NE,)               — uniform energy grid (eV)
    mu_L, mu_R : float                   — contact electrochemical potentials (eV)
    kT : float                           — thermal energy (eV)
    m_eff_kg : float                     — transverse effective mass (kg)
    a_m : float                          — grid spacing (m)

    Returns
    -------
    n : ndarray (Np,)  — electron density (m⁻³)
    """
    dE = float(E_grid[1] - E_grid[0])
    G_A = np.conj(np.transpose(G_R, (0, 2, 1)))
    dl = np.real(np.diagonal(np.matmul(np.matmul(G_R, Gam_L), G_A), axis1=1, axis2=2))
    dr = np.real(np.diagonal(np.matmul(np.matmul(G_R, Gam_R), G_A), axis1=1, axis2=2))

    def _supply(mu):
        x = (mu - E_grid) / kT
        # stable kT·ln(1+e^x) with kT folded in the prefactor → return ln(1+e^x)
        return np.where(x > 30.0, x, np.log1p(np.exp(np.clip(x, -600.0, 30.0))))

    integ = np.sum(dl * _supply(mu_L)[:, None] + dr * _supply(mu_R)[:, None],
                   axis=0) * dE / (2.0 * np.pi)          # (Np,) dimensionless
    prefac = 2.0 * m_eff_kg * (kT * _Q) / (2.0 * np.pi * _HBAR ** 2 * a_m)  # 1/m³
    return prefac * integ


# ───────────────────────────────────────────────────────────────────────────
# Fermi-Dirac integrals F_{1/2}(η) and F_{−1/2}(η)
#
# Convention: ``n = N_c · F_{1/2}((E_F − E_c)/kT)`` with the normalization
# F_{1/2}(η) → exp(η) in the Boltzmann (non-degenerate) limit η → −∞ and
# F_{1/2}(η) → (4/(3√π))·η^(3/2) in the deep-degenerate limit η → +∞ — i.e.
# the thesis (Akkala Eq 3.1) Fermi-Dirac normalization (the factor 2/√π is
# absorbed into F).
#
# Implementation: 64-point Gauss–Laguerre quadrature of the integral
# F_{1/2}(η) = (2/√π) ∫_0^∞ √u · e^u / (1+e^(u−η)) · e^(−u) du,
# with Boltzmann series for η < −10 and the leading asymptotic expansion for
# η > 30. F_{−1/2}(η) ≡ dF_{1/2}/dη is computed by central finite difference,
# which inherits the Gauss–Laguerre accuracy (≲ 10⁻⁶ relative).
# ───────────────────────────────────────────────────────────────────────────
_GL_X, _GL_W = np.polynomial.laguerre.laggauss(64)


def _F12_scalar(eta: float) -> float:
    if eta < -10.0:
        z = float(np.exp(eta))
        return z * (1.0 - z / 2.828427124746190 + (z * z) / 5.196152422706632
                    - (z ** 3) / 8.0)
    if eta > 30.0:
        e2 = eta * eta
        return (4.0 / (3.0 * np.sqrt(np.pi))) * eta ** 1.5 * (
            1.0 + np.pi ** 2 / (8.0 * e2)
            + 7.0 * np.pi ** 4 / (640.0 * e2 * e2))
    # f(u) e^(-u) with f(u) = √u e^u / (1 + e^(u-η)). The exponentials cancel
    # so the integrand stays bounded for all u ≥ 0.
    u = _GL_X
    # Use log-sum-exp style: 1/(1+e^(u-η)) = e^(η)/(e^η + e^u) → numerically stable.
    denom = np.exp(eta) + np.exp(u)
    f = np.sqrt(u) * np.exp(u + eta) / denom
    return float((2.0 / np.sqrt(np.pi)) * np.sum(_GL_W * f))


def fermi_dirac_half(eta):
    """Fermi–Dirac integral ``F_{1/2}(η)`` (Blakemore convention: ``F_{1/2}(η) →
    exp(η)`` as η→−∞; ``F_{1/2}(0) ≈ 0.7654``).

    Scalar in → scalar out; array in → same-shape array out.
    """
    arr = np.asarray(eta, dtype=float)
    if arr.ndim == 0:
        return _F12_scalar(float(arr))
    out = np.empty_like(arr)
    flat = out.ravel()
    src = arr.ravel()
    for i in range(src.size):
        flat[i] = _F12_scalar(float(src[i]))
    return out


def fermi_dirac_minus_half(eta):
    """``F_{−1/2}(η) = dF_{1/2}/dη``, central difference of F_{1/2}."""
    h = 1.0e-4
    return (fermi_dirac_half(np.asarray(eta, dtype=float) + h)
            - fermi_dirac_half(np.asarray(eta, dtype=float) - h)) / (2.0 * h)


def solve_poisson_1d(
    eps_r: np.ndarray,
    N_D: np.ndarray,
    n_e: np.ndarray,
    a_m: float,
    U_left: float,
    U_right: float,
    *,
    left_bc: str = "dirichlet",
    right_bc: str = "dirichlet",
) -> np.ndarray:
    """Solve the 1D Poisson equation (★) with per-end Dirichlet or Neumann BCs.

    Parameters
    ----------
    eps_r : ndarray, shape (Np,)
        Relative permittivity at each grid site (dimensionless).
    N_D : ndarray, shape (Np,)
        Ionized donor density at each site (m⁻³).
    n_e : ndarray, shape (Np,)
        Electron density at each site (m⁻³).
    a_m : float
        Uniform grid spacing (m).
    U_left, U_right : float
        For a Dirichlet end, the fixed value of the electron potential energy
        at that contact (eV). For an applied bias ``V`` the natural Dirichlet
        choice is ``U_left = +V/2``, ``U_right = −V/2`` (matching the solver's
        ``μ_L = E_F + V/2`` convention). For a **Neumann** end the corresponding
        argument is ignored (the zero-field condition carries no value).
    left_bc, right_bc : {"dirichlet", "neumann"}
        Boundary condition at each end. ``"neumann"`` imposes a **zero-field**
        (``dU/dz = 0``) condition via a finite-volume half-cell balance — the
        physical condition for a charge-neutral ohmic reservoir edge, where no
        electric field penetrates the neutral bulk.

        .. note::
           A **both-ends Neumann** problem is singular: the Poisson operator's
           nullspace is the constant vector, so the potential is defined only
           up to a gauge and the total drop ``U[0]−U[-1]`` is fixed by the
           charge — it *cannot* be set to an applied bias ``V``. This function
           therefore rejects ``left_bc == right_bc == "neumann"``. A swept bias
           needs at least one Dirichlet end (or an explicit ``U[0]−U[-1]=V``
           constraint row, handled by the caller).

    Returns
    -------
    U : ndarray, shape (Np,)
        Electron potential energy (eV), ready to add to the H diagonal.
    """
    eps_r = np.asarray(eps_r, dtype=float)
    N_D = np.asarray(N_D, dtype=float)
    n_e = np.asarray(n_e, dtype=float)
    Np = N_D.size
    if not (eps_r.size == n_e.size == Np):
        raise ValueError("eps_r, N_D, n_e must have the same length")
    if left_bc not in ("dirichlet", "neumann") or right_bc not in ("dirichlet", "neumann"):
        raise ValueError("left_bc/right_bc must be 'dirichlet' or 'neumann'")
    if left_bc == "neumann" and right_bc == "neumann":
        raise ValueError(
            "both-ends Neumann is singular (nullspace = const) and cannot "
            "carry an applied bias; use at least one Dirichlet end."
        )

    eps = eps_r * _EPS0          # F/m
    rho = _Q * (N_D - n_e)       # C/m³
    inv_a2 = 1.0 / (a_m * a_m)

    A = np.zeros((Np, Np))
    b = np.zeros(Np)
    for i in range(Np):
        if i == 0 and left_bc == "dirichlet":
            A[i, i] = 1.0
            b[i] = U_left
        elif i == Np - 1 and right_bc == "dirichlet":
            A[i, i] = 1.0
            b[i] = U_right
        elif i == 0:
            # Zero-field Neumann, left face: FV half-cell balance
            #   ε_R (U[1]−U[0])/a² = ρ[0]/2   (left-face flux = 0).
            eps_R = 0.5 * (eps[i] + eps[i + 1])
            A[i, i] = -eps_R * inv_a2
            A[i, i + 1] = eps_R * inv_a2
            b[i] = 0.5 * rho[i]
        elif i == Np - 1:
            # Zero-field Neumann, right face: FV half-cell balance
            #   ε_L (U[-2]−U[-1])/a² = ρ[-1]/2  (right-face flux = 0).
            eps_L = 0.5 * (eps[i - 1] + eps[i])
            A[i, i - 1] = eps_L * inv_a2
            A[i, i] = -eps_L * inv_a2
            b[i] = 0.5 * rho[i]
        else:
            eps_L = 0.5 * (eps[i - 1] + eps[i])
            eps_R = 0.5 * (eps[i] + eps[i + 1])
            A[i, i - 1] = eps_L * inv_a2
            A[i, i] = -(eps_L + eps_R) * inv_a2
            A[i, i + 1] = eps_R * inv_a2
            b[i] = rho[i]
    # Np ~ 135: a dense solve is trivial and avoids a scipy dependency.
    return np.linalg.solve(A, b)


def _dirichlet_setup(Np, U_left, U_right, dirichlet_mask, dirichlet_vals):
    """Resolve which sites are held Dirichlet and at what value.

    Default (mask None): only the two endpoints, at ``U_left``/``U_right``.
    With a mask (e.g. the charge-neutral doped-contact region), every masked
    site is held at the corresponding ``dirichlet_vals`` entry and Poisson is
    solved only on the unmasked (active) interior.
    """
    if dirichlet_mask is None:
        fixed = np.zeros(Np, dtype=bool)
        fixed[0] = fixed[-1] = True
        fvals = np.zeros(Np)
        fvals[0] = U_left
        fvals[-1] = U_right
        return fixed, fvals
    return np.asarray(dirichlet_mask, dtype=bool), np.asarray(dirichlet_vals, dtype=float)


def poisson_newton_update(
    U_old: np.ndarray,
    eps_r: np.ndarray,
    N_D: np.ndarray,
    n_e: np.ndarray,
    a_m: float,
    U_left: float,
    U_right: float,
    kT_screen: float,
    dirichlet_mask: np.ndarray | None = None,
    dirichlet_vals: np.ndarray | None = None,
    neutral_mask: np.ndarray | None = None,
    dn_dU_override: np.ndarray | None = None,
    left_bc: str = "dirichlet",
    right_bc: str = "dirichlet",
) -> np.ndarray:
    """One predictor–corrector (quasi-Newton) Poisson step.

    Plain Gummel iteration (``solve_poisson_1d`` fed the raw NEGF charge) limit-
    cycles for resonant-tunneling devices: a few-meV change in U can swing the
    resonant well charge by orders of magnitude, so the fixed-point map
    ``U → solve_poisson(n(U))`` is wildly non-contractive. The standard cure
    (Trellakis/Datta) is to fold the *linearized* charge response into the
    Poisson operator. With the Thomas–Fermi estimate ``dn/dU ≈ −n/kT`` the
    Newton system is

        [ L + diag(q·dn/dU) ] δU = −F(U_old),
        F(U_old) = L[U_old] − q(N_D − n_e),

    where ``L`` is the discrete ``d/dz[ε dU/dz]``. Because ``dn/dU < 0`` the
    added diagonal is negative, strengthening the (already negative) Laplacian
    diagonal exactly where the charge is large — self-damping the resonant well.

    The approximate Jacobian affects only the convergence *path*: at the fixed
    point ``δU → 0`` and ``F(U*) = 0`` is the exact Poisson/NEGF balance, so a
    crude (even over-damped) ``dn/dU`` is safe. ``kT_screen`` is the screening
    temperature in eV; smaller values over-damp (slower but more robust).

    Parameters
    ----------
    U_old : ndarray, shape (Np,)
        Current potential (eV); must satisfy the Dirichlet BC at the contacts.
    eps_r, N_D, n_e : ndarray, shape (Np,)
        Permittivity (rel.), donor and electron densities (m⁻³).
    a_m : float
        Grid spacing (m).
    U_left, U_right : float
        Dirichlet boundary values (eV).
    kT_screen : float
        Screening temperature (eV) for ``dn/dU ≈ −n/kT_screen``.

    Returns
    -------
    U_new : ndarray, shape (Np,)
        Updated potential (eV), satisfying the same Dirichlet BC.
    """
    eps_r = np.asarray(eps_r, dtype=float)
    N_D = np.asarray(N_D, dtype=float)
    n_e = np.asarray(n_e, dtype=float)
    U_old = np.asarray(U_old, dtype=float)
    Np = N_D.size

    eps = eps_r * _EPS0
    inv_a2 = 1.0 / (a_m * a_m)
    rho = _Q * (N_D - n_e)           # C/m³
    if dn_dU_override is not None:
        dn_dU = np.asarray(dn_dU_override, dtype=float)
    else:
        dn_dU = -n_e / kT_screen     # Boltzmann predictor (m⁻³ / V, ≤ 0)
    if left_bc not in ("dirichlet", "neumann") or right_bc not in ("dirichlet", "neumann"):
        raise ValueError("left_bc/right_bc must be 'dirichlet' or 'neumann'")
    if left_bc == "neumann" and right_bc == "neumann":
        raise ValueError(
            "both-ends Neumann is singular (nullspace = const) and cannot carry "
            "an applied bias; keep one Dirichlet end as the gauge/bias reference."
        )
    fixed, fvals = _dirichlet_setup(Np, U_left, U_right, dirichlet_mask, dirichlet_vals)
    # A Neumann endpoint is a free unknown (zero-field), not a pinned value.
    if left_bc == "neumann":
        fixed[0] = False
    if right_bc == "neumann":
        fixed[-1] = False
    _ = neutral_mask  # legacy arg; charge model is now assembled by the caller

    J = np.zeros((Np, Np))
    F = np.zeros(Np)
    for i in range(Np):
        if fixed[i]:
            J[i, i] = 1.0
            F[i] = U_old[i] - fvals[i]
        elif i == 0:
            # Zero-field Neumann, left face (half-cell): flux in from the right
            # balances the half-cell charge. dn/dU enters at half weight.
            cR = 0.5 * (eps[i] + eps[i + 1]) * inv_a2
            F[i] = cR * (U_old[i + 1] - U_old[i]) - 0.5 * rho[i]
            J[i, i] = -cR + 0.5 * _Q * dn_dU[i]
            J[i, i + 1] = cR
        elif i == Np - 1:
            # Zero-field Neumann, right face (half-cell).
            cL = 0.5 * (eps[i - 1] + eps[i]) * inv_a2
            F[i] = cL * (U_old[i - 1] - U_old[i]) - 0.5 * rho[i]
            J[i, i - 1] = cL
            J[i, i] = -cL + 0.5 * _Q * dn_dU[i]
        else:
            eps_L = 0.5 * (eps[i - 1] + eps[i])
            eps_R = 0.5 * (eps[i] + eps[i + 1])
            cL = eps_L * inv_a2
            cR = eps_R * inv_a2
            # Residual F = L[U_old] − rho.
            L_U = cL * (U_old[i - 1] - U_old[i]) + cR * (U_old[i + 1] - U_old[i])
            F[i] = L_U - rho[i]
            # Jacobian row: L stencil + q·dn/dU on the diagonal.
            J[i, i - 1] = cL
            J[i, i] = -(cL + cR) + _Q * dn_dU[i]
            J[i, i + 1] = cR
    delta = np.linalg.solve(J, -F)
    return U_old + delta


def density_response_jacobian(
    G_R: np.ndarray,
    G_lesser: np.ndarray,
    dE: float,
    prefactor: float,
) -> np.ndarray:
    """Quantum density-response Jacobian ``∂n_i/∂U_j`` for full Newton-Raphson.

    Linearizing the Keldysh density ``n = (C/2π)∫dE Re G^<`` with
    ``G^< = G^R Σⁱⁿ G^A`` in a diagonal potential perturbation
    ``δH = diag(δU)`` (so ``δG^R = G^R δH G^R``, ``δG^A = G^A δH G^A``, and the
    contact/phonon self-energies are held fixed — exact for interior sites,
    where the Dirichlet BC pins ``δU=0`` at the contacts, and a weak-coupling
    approximation for the phonon Σ) gives the closed form

        ∂n_i/∂U_j = (C/2π) ∫dE  Re[ G^R_ij G^<_ji + G^<_ij G^A_ji ].

    This reuses the Green's functions the SCBA already computed, so the full
    Jacobian costs one O(Np²·NE) assembly — no finite-difference NEGF re-solves.
    It is the response that makes the Poisson–NEGF Newton step converge
    quadratically (vs. the diagonal Thomas–Fermi guess in
    :func:`poisson_newton_update`).

    Parameters
    ----------
    G_R, G_lesser : ndarray, shape (NE, Np, Np), complex
        Retarded and lesser Green's functions on the (uniform) energy grid.
    dE : float
        Energy grid spacing (eV).
    prefactor : float
        Density prefactor ``C`` (so the result is in m⁻³ per V).

    Returns
    -------
    dn_dU : ndarray, shape (Np, Np)
        ``∂n_i/∂U_j`` in m⁻³ / V (real).
    """
    G_A = np.conj(np.transpose(G_R, (0, 2, 1)))
    Gl_T = np.transpose(G_lesser, (0, 2, 1))     # G^<_ji at [E, i, j]
    GA_T = np.transpose(G_A, (0, 2, 1))          # G^A_ji at [E, i, j]
    integrand = G_R * Gl_T + G_lesser * GA_T     # elementwise (NE, Np, Np)
    summed = np.real(integrand.sum(axis=0))      # (Np, Np)
    return prefactor * (dE / (2.0 * np.pi)) * summed


def poisson_newton_full_step(
    U_old: np.ndarray,
    eps_r: np.ndarray,
    N_D: np.ndarray,
    n_e: np.ndarray,
    dn_dU: np.ndarray,
    a_m: float,
    U_left: float,
    U_right: float,
    max_step: float | None = None,
    dirichlet_mask: np.ndarray | None = None,
    dirichlet_vals: np.ndarray | None = None,
) -> np.ndarray:
    """One full Newton-Raphson step for coupled Poisson–NEGF.

    Solves ``J δU = −F`` with the *exact* (dense) Jacobian

        F_i = L[U]_i − q(N_D_i − n_i),
        J   = L + q · ∂n/∂U                                   (interior rows)

    where ``L = d/dz[ε d/dz]`` is the discrete Poisson operator and ``∂n/∂U`` is
    :func:`density_response_jacobian`. Dirichlet rows are the identity with
    ``δU = 0``. Near a fixed point this converges quadratically.

    With ``dirichlet_mask`` the whole doped-contact region is held Dirichlet at
    ``dirichlet_vals`` (charge-neutral reservoir BC) and Poisson is solved only
    on the active barriers+well region — otherwise the contacts develop spurious
    (~1 eV) band bending at finite bias. Default (None) fixes only the endpoints.

    ``max_step`` optionally damps the Newton step (trust-region style): if
    ``max|δU| > max_step`` the step is scaled down. This guards against
    overshoot where the Jacobian is poorly conditioned (e.g. close to RTD
    bistability); the converged fixed point is unaffected.
    """
    eps_r = np.asarray(eps_r, dtype=float)
    N_D = np.asarray(N_D, dtype=float)
    n_e = np.asarray(n_e, dtype=float)
    U_old = np.asarray(U_old, dtype=float)
    Np = N_D.size

    eps = eps_r * _EPS0
    inv_a2 = 1.0 / (a_m * a_m)
    rho = _Q * (N_D - n_e)
    fixed, fvals = _dirichlet_setup(Np, U_left, U_right, dirichlet_mask, dirichlet_vals)

    J = np.zeros((Np, Np))
    F = np.zeros(Np)
    for i in range(Np):
        if fixed[i]:
            J[i, i] = 1.0
            F[i] = U_old[i] - fvals[i]
        else:
            eps_L = 0.5 * (eps[i - 1] + eps[i])
            eps_R = 0.5 * (eps[i] + eps[i + 1])
            cL = eps_L * inv_a2
            cR = eps_R * inv_a2
            F[i] = (cL * (U_old[i - 1] - U_old[i])
                    + cR * (U_old[i + 1] - U_old[i]) - rho[i])
            J[i, i - 1] += cL
            J[i, i] += -(cL + cR)
            J[i, i + 1] += cR
            J[i, :] += _Q * dn_dU[i, :]   # full quantum charge-response row
    delta = np.linalg.solve(J, -F)
    if max_step is not None:
        m = float(np.max(np.abs(delta)))
        if m > max_step:
            delta *= max_step / m
    return U_old + delta


def extract_density_1d(G_lesser: np.ndarray, dE: float) -> np.ndarray:
    """Per-site longitudinal occupation from the projected ``G^<``.

    In the rank-1 Keldysh solver, ``G_lesser = G^R Σ^in G^A`` has a real,
    non-negative diagonal (it is a sum of ``Σ_in[j] |G^R[z,j]|²``), which plays
    the role of Patil's electron correlation array ``n(z, E)``. The per-site
    occupation (per transverse mode) is the energy integral

        N_raw(z) = (1 / 2π) ∫ dE  n(z, E).

    The absolute m⁻³ density is obtained by multiplying this by a single
    geometric prefactor (transverse modes per area / site volume), which
    :mod:`core.poisson_negf` calibrates from contact charge neutrality.

    Parameters
    ----------
    G_lesser : ndarray, shape (NE, Np, Np), complex
        Lesser Green's function on the energy grid.
    dE : float
        Energy grid spacing (eV); the grid must be uniform.

    Returns
    -------
    N_raw : ndarray, shape (Np,)
        Per-site occupation (dimensionless), before the geometric prefactor.
    """
    Np = G_lesser.shape[1]
    diag = np.real(G_lesser[:, np.arange(Np), np.arange(Np)])  # (NE, Np)
    return np.sum(diag, axis=0) * dE / (2.0 * np.pi)


def calibrate_density_prefactor(
    N_raw: np.ndarray,
    N_D: np.ndarray,
    contact_sites: np.ndarray | slice,
) -> float:
    """Geometric prefactor ``C`` such that ``n_e[m⁻³] = C · N_raw``.

    Fixed by demanding charge neutrality in the doped contact at equilibrium:
    deep inside the reservoir the electron density must equal the donor
    density, so ``C = mean(N_D[contact]) / mean(N_raw[contact])``. This folds
    the transverse-mode count and site volume into one constant, consistent
    with the Datta analytic transverse-area ("Fix A") treatment used for the
    current elsewhere in the pipeline.

    Parameters
    ----------
    N_raw : ndarray, shape (Np,)
        Equilibrium (V≈0) per-site occupation from :func:`extract_density_1d`.
    N_D : ndarray, shape (Np,)
        Donor density (m⁻³).
    contact_sites : ndarray of int or slice
        Indices of deep-contact sites used for the neutrality anchor.
    """
    Nc = float(np.mean(N_raw[contact_sites]))
    Dc = float(np.mean(N_D[contact_sites]))
    if Nc <= 0.0:
        raise ValueError(
            "Non-positive contact occupation; cannot calibrate density "
            "prefactor (check energy grid covers the occupied band)."
        )
    return Dc / Nc
