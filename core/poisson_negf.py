"""
Self-consistent Poisson–NEGF outer loop for the rank-1 Keldysh SCBA.

This wraps the *unchanged* :func:`core.scba_rank1_keldysh.run_rank1_keldysh_single_bias`
in a Gummel-style outer iteration:

    1. Solve NEGF-SCBA with the current potential U(z)  →  G^<  →  n(z).
    2. Solve 1D Poisson for U_new(z) given n(z) (Dirichlet ±V/2).
    3. Mix U ← (1−β)U + βU_new, repeat until max|ΔU| < tol.

The transport is 1D along the growth direction z; the transverse plane is the
analytic Datta-area prefactor ("Fix A"), so the electron density is normalized
to m⁻³ by a single geometric prefactor calibrated from contact charge
neutrality at equilibrium (see :func:`compute_density_prefactor`).

The default (Poisson OFF) reproduction path through ``run/run_rank1_sweep.py``
is untouched: nothing here is called unless ``--poisson`` is passed.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from config.device_library import DEVICES, MATERIALS
from core.poisson import (
    calibrate_density_prefactor,
    density_response_jacobian,
    extract_density_1d,
    poisson_newton_full_step,
    poisson_newton_update,
)
from core.scba_rank1_keldysh import (
    Rank1KeldyshResult,
    run_rank1_keldysh_single_bias,
)


# ──────────────────────────────────────────────────────────────────────────────
# Device electrostatics on the transport grid
# ──────────────────────────────────────────────────────────────────────────────
def build_electrostatics(device_name: str, a_m: float):
    """Per-site permittivity, donor density, and contact mask.

    Replicates the per-layer site assignment of
    ``run.run_rank1_sweep.build_stack`` (``n_sites = max(1, round(t/a))`` per
    layer, layers in order) so the returned arrays align site-for-site with the
    Hamiltonian grid ``z_nm``.

    Returns
    -------
    eps_r : ndarray, shape (Np,)
        Relative permittivity at each site.
    N_D : ndarray, shape (Np,)
        Ionized donor density (m⁻³) at each site.
    contact_mask : ndarray of bool, shape (Np,)
        True on sites belonging to a doped (contact) layer.
    """
    dev = DEVICES[device_name]
    eps_list: list[float] = []
    nd_list: list[float] = []
    contact_list: list[bool] = []
    for layer in dev["layers"]:
        mat = layer["material"]
        thick_m = layer["thickness"]
        doping = float(layer.get("doping", 0.0))
        n_sites = max(1, int(round(thick_m / a_m)))
        eps_r = MATERIALS[mat]["epsilon_r"]
        for _ in range(n_sites):
            eps_list.append(eps_r)
            nd_list.append(doping)
            contact_list.append(doping > 0.0)
    return (
        np.asarray(eps_list, dtype=float),
        np.asarray(nd_list, dtype=float),
        np.asarray(contact_list, dtype=bool),
    )


def deep_contact_sites(contact_mask: np.ndarray, drop_edge: int = 2,
                       n_take: int = 5) -> np.ndarray:
    """Indices of stable bulk sites inside the left (emitter) doped contact.

    Skips the outermost ``drop_edge`` sites (which carry the tight-binding lead
    self-energy and Friedel oscillations) and returns the next ``n_take``,
    used to anchor the density prefactor to contact charge neutrality.
    """
    left = np.where(contact_mask)[0]
    if left.size == 0:
        raise ValueError("No doped contact sites found; cannot anchor density.")
    # First contiguous contact block (the emitter side).
    start = left[0]
    block = [start]
    for idx in left[1:]:
        if idx == block[-1] + 1:
            block.append(idx)
        else:
            break
    block = np.asarray(block)
    if block.size <= drop_edge:
        return block
    chosen = block[drop_edge: drop_edge + n_take]
    return chosen if chosen.size else block[drop_edge:]


def linear_bias_profile(V: float, Np: int, NS: int = 1, ND: int = 1) -> np.ndarray:
    """Linear potential drop: flat +V/2 in emitter, ramp, flat −V/2 in collector.

    Identical to ``run.run_rank1_sweep.linear_bias_profile`` — used here only as
    the warm-start for the Poisson outer loop.
    """
    NC = Np - NS - ND
    return V * np.concatenate([
        0.5 * np.ones(NS),
        np.linspace(0.5, -0.5, NC),
        -0.5 * np.ones(ND),
    ])


# ──────────────────────────────────────────────────────────────────────────────
# Density prefactor (equilibrium calibration)
# ──────────────────────────────────────────────────────────────────────────────
def compute_density_prefactor(
    *,
    E_grid: np.ndarray,
    H_z: np.ndarray,
    UB: np.ndarray,
    t0: float,
    Ef: float,
    kT: float,
    chi_diag: np.ndarray,
    D0_sq_per_mode,
    hnu_idx_per_mode,
    N_bose_per_mode,
    chi_per_mode,
    N_D: np.ndarray,
    contact_mask: np.ndarray,
    NS: int = 1,
    ND: int = 1,
    scba_max_iter: int = 10,
    scba_mix: float = 0.4,
    scba_tol: float = 1e-5,
    eta: float = 1e-12,
) -> float:
    """Run one equilibrium (V=0) NEGF solve and calibrate the density prefactor.

    Returns ``C`` such that ``n_e[m⁻³] = C · N_raw`` (charge-neutral contacts).
    """
    Np = H_z.shape[0]
    bp = linear_bias_profile(0.0, Np, NS, ND)  # all zeros at V=0
    res = run_rank1_keldysh_single_bias(
        V=0.0, E_grid=E_grid, H_z=H_z, UB=UB, bias_profile=bp, t0=t0,
        Ef=Ef, kT=kT, chi_diag=chi_diag, D0_sq_per_mode=D0_sq_per_mode,
        hnu_idx_per_mode=hnu_idx_per_mode, N_bose_per_mode=N_bose_per_mode,
        chi_per_mode=chi_per_mode, max_iter=scba_max_iter, tol=scba_tol,
        mix=scba_mix, eta=eta,
    )
    dE = float(E_grid[1] - E_grid[0])
    N_raw = extract_density_1d(res.G_lesser, dE)
    sites = deep_contact_sites(contact_mask)
    return calibrate_density_prefactor(N_raw, N_D, sites)


# ──────────────────────────────────────────────────────────────────────────────
# Self-consistent single-bias solve
# ──────────────────────────────────────────────────────────────────────────────
@dataclass
class SelfConsistentResult:
    """One bias point's self-consistent Poisson–NEGF output."""

    result: Rank1KeldyshResult   # final NEGF solve at the converged U
    U: np.ndarray                # (Np,) converged electron potential energy (eV)
    n_e: np.ndarray              # (Np,) converged electron density (m⁻³)
    poisson_iters: int
    poisson_converged: bool
    dU_final: float              # max|ΔU| at the last iteration (eV)


def run_self_consistent_bias(
    *,
    V: float,
    E_grid: np.ndarray,
    H_z: np.ndarray,
    UB: np.ndarray,
    t0: float,
    Ef: float,
    kT: float,
    chi_diag: np.ndarray,
    D0_sq_per_mode,
    hnu_idx_per_mode,
    N_bose_per_mode,
    chi_per_mode,
    eps_r: np.ndarray,
    N_D: np.ndarray,
    a_m: float,
    density_prefactor: float,
    NS: int = 1,
    ND: int = 1,
    scba_max_iter: int = 10,
    scba_mix: float = 0.4,
    scba_tol: float = 1e-5,
    eta: float = 1e-12,
    poisson_max_iter: int = 40,
    poisson_tol: float = 5e-3,
    kT_screen: float = 0.002,
    newton_max_step: float = 0.05,
    newton_stall_patience: int = 3,
    U_init: np.ndarray | None = None,
    verbose: bool = False,
) -> SelfConsistentResult:
    """Self-consistent Poisson–NEGF loop at a single applied bias ``V``.

    The Dirichlet boundary values are ``U[0] = +V/2``, ``U[-1] = −V/2`` so the
    applied bias is carried by Poisson; the interior potential redistributes
    self-consistently with the NEGF charge. Convergence is measured by the
    magnitude of the **actual potential update** ``max|U_new − U|`` per
    iteration: in Newton mode this is the true Newton step (well-scaled by the
    exact Jacobian; → 0 only at self-consistency, so it does not short-circuit
    before the well potential forms), in predictor mode the over-damped step.

    **Hybrid Newton/predictor stepping.** Each iteration takes a **full
    Newton-Raphson** step using the exact quantum density-response Jacobian
    ``∂n/∂U`` (:func:`core.poisson.density_response_jacobian`) — quadratically
    convergent and the standard method (off resonance, V=0 converges in ~2
    steps). But near an RTD resonance ``∂n/∂U`` is sign-indefinite (raising the
    potential can shift a resonance *into* the window, *increasing* density), so
    the Newton Jacobian is indefinite and the damped step can stall. When Newton
    fails to reduce the residual for ``newton_stall_patience`` iterations, the
    loop falls back to the **over-damped Thomas–Fermi predictor**
    (:func:`core.poisson.poisson_newton_update` with ``dn/dU ≈ −n/kT_screen``),
    whose Jacobian is always negative-definite — a guaranteed (if slow) descent
    that grinds the hard resonance points down to the inner-SCBA noise floor.
    This is Levenberg–Marquardt's logic (Newton ↔ gradient descent) as a clean
    method switch.

    ``newton_max_step`` trust-region-caps the Newton step (eV); ``kT_screen``
    (default 2 meV) sets the predictor damping. ``U_init`` (e.g. the converged U
    from the previous bias) gives a warm start; if None the linear drop is used.

    Returns
    -------
    SelfConsistentResult
    """
    Np = H_z.shape[0]
    dE = float(E_grid[1] - E_grid[0])
    U_L, U_R = +V / 2.0, -V / 2.0

    U = linear_bias_profile(V, Np, NS, ND) if U_init is None else U_init.copy()
    U[0], U[-1] = U_L, U_R  # enforce BC on the warm start

    def negf_and_density(U_profile):
        r = run_rank1_keldysh_single_bias(
            V=V, E_grid=E_grid, H_z=H_z, UB=UB, bias_profile=U_profile, t0=t0,
            Ef=Ef, kT=kT, chi_diag=chi_diag, D0_sq_per_mode=D0_sq_per_mode,
            hnu_idx_per_mode=hnu_idx_per_mode, N_bose_per_mode=N_bose_per_mode,
            chi_per_mode=chi_per_mode, max_iter=scba_max_iter, tol=scba_tol,
            mix=scba_mix, eta=eta,
        )
        return r, density_prefactor * extract_density_1d(r.G_lesser, dE)

    res: Rank1KeldyshResult | None = None
    n_e = np.zeros(Np)
    converged = False
    dU = np.inf
    best_dU = np.inf
    stall = 0
    use_predictor = False  # latches once full Newton proves it cannot progress
    it = 0
    for it in range(1, poisson_max_iter + 1):
        res, n_e = negf_and_density(U)
        # Take the step for this iteration, then measure it. The convergence
        # metric is the magnitude of the **actual update** ``max|U_new − U|``.
        # In Newton mode this is the true Newton step (well-scaled by the exact
        # Jacobian, → 0 only at self-consistency) — so it does NOT short-circuit
        # before the well potential forms. In predictor-fallback mode it is the
        # over-damped step, which floors at the hard resonance points (accepted).
        if use_predictor:
            # Over-damped Thomas–Fermi predictor (dn/dU ≈ −n/kT_screen makes the
            # Jacobian negative-definite → a guaranteed, if slow, descent).
            U_new = poisson_newton_update(U, eps_r, N_D, n_e, a_m,
                                          U_L, U_R, kT_screen)
        else:
            # Full Newton-Raphson with the exact quantum response Jacobian.
            J = density_response_jacobian(res.G_R, res.G_lesser, dE,
                                          density_prefactor)
            U_new = poisson_newton_full_step(U, eps_r, N_D, n_e, J, a_m,
                                             U_L, U_R, max_step=newton_max_step)
        dU = float(np.max(np.abs(U_new - U)))
        U = U_new
        if verbose:
            mode = "pred " if use_predictor else "newton"
            print(f"    [poisson {it:2d} {mode}] step={dU*1e3:8.3f} meV  "
                  f"I_R={res.I_right:+.3e} A  Umin={U.min():+.3f}  "
                  f"n_well_max={n_e.max():.2e} m⁻³", flush=True)
        if dU < poisson_tol:
            converged = True
            break

        # Track progress; if full Newton stalls (step stuck/capped, residual not
        # shrinking) latch to the always-stable predictor for the rest.
        if dU < best_dU - 1e-12:
            best_dU = dU
            stall = 0
        else:
            stall += 1
        if not use_predictor and stall >= newton_stall_patience:
            use_predictor = True

    # Final consistency solve so the returned observables match the returned U.
    res, n_e = negf_and_density(U)

    return SelfConsistentResult(
        result=res, U=U, n_e=n_e, poisson_iters=it,
        poisson_converged=converged, dU_final=dU,
    )
