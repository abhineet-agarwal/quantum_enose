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
    ballistic_transverse_density,
    calibrate_density_prefactor,
    extract_density_1d,
    fermi_dirac_half,
    fermi_dirac_minus_half,
    physical_transverse_density,
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


def terminal_masks(contact_mask: np.ndarray,
                   quantum_buffer_sites: int = 15):
    """Split the doped contacts into (deep flat-band emitter | collector)
    terminals and the inner quantum-buffer region (Akkala thesis Fig 3.1/3.3).

    Returns ``(em_term, col_term)``: boolean masks marking the *deep* part of
    each contact (outer ``block_size − quantum_buffer_sites`` sites) where
    semiclassical **Thomas-Fermi / Boltzmann charge** is used — providing the
    Thomas-Fermi screening that pins the contact at flat band while the inner
    ``quantum_buffer_sites`` (≈3 nm) hosts the quantum NEGF charge so the
    **emitter accumulation / triangular-well quasi-bound state** can form.
    """
    cm = np.asarray(contact_mask, dtype=bool)
    Np = cm.size
    em_term = np.zeros(Np, dtype=bool)
    col_term = np.zeros(Np, dtype=bool)
    idx = np.where(cm)[0]
    if idx.size == 0:
        return em_term, col_term
    # Emitter block: contiguous run from idx[0]
    p = idx[0]; em = []
    for k in idx:
        if k == p: em.append(k); p += 1
        else: break
    # Collector block: contiguous run back from idx[-1]
    p = idx[-1]; col = []
    for k in idx[::-1]:
        if k == p: col.append(k); p -= 1
        else: break
    em = np.asarray(em); col = np.asarray(sorted(col))
    if em.size > quantum_buffer_sites:
        em_term[em[:-quantum_buffer_sites]] = True
    if col.size > quantum_buffer_sites:
        col_term[col[quantum_buffer_sites:]] = True
    return em_term, col_term


def flat_band_terminal_mask(contact_mask: np.ndarray,
                            quantum_buffer_sites: int = 15) -> np.ndarray:
    """Combined deep flat-band terminal mask (emitter ∪ collector).

    See :func:`terminal_masks` for the partition into individual sides.
    """
    em, col = terminal_masks(contact_mask, quantum_buffer_sites)
    return em | col


def outer_clamp_masks(contact_mask: np.ndarray, n_clamp_sites: int):
    """Outermost ``n_clamp_sites`` of each doped contact, for Dirichlet clamping.

    A pragmatic complement to the FD/F_{1/2} terminal screening for *short*
    contacts (like the 10 nm ZnO in the SISPAD stack), where natural screening
    alone can't fully pin the deep contact at the rail at high bias. Clamping
    the outermost N sites of each contact at ±V/2 (Dirichlet) extends the
    "true reservoir" depth, leaving the middle of the contact for F_{1/2}
    screening and the inner quantum-buffer for the NEGF accumulation layer.
    """
    cm = np.asarray(contact_mask, dtype=bool)
    Np = cm.size
    em_clamp = np.zeros(Np, dtype=bool)
    col_clamp = np.zeros(Np, dtype=bool)
    if n_clamp_sites <= 0:
        return em_clamp, col_clamp
    idx = np.where(cm)[0]
    if idx.size == 0:
        return em_clamp, col_clamp
    p = idx[0]; em = []
    for k in idx:
        if k == p: em.append(k); p += 1
        else: break
    p = idx[-1]; col = []
    for k in idx[::-1]:
        if k == p: col.append(k); p -= 1
        else: break
    em = np.asarray(em); col = np.asarray(sorted(col))
    n = min(n_clamp_sites, em.size)
    em_clamp[em[:n]] = True
    n = min(n_clamp_sites, col.size)
    col_clamp[col[-n:]] = True
    return em_clamp, col_clamp


def contact_dirichlet(contact_mask: np.ndarray, V: float):
    """[Legacy] Whole-contact flat clamp at ±V/2.

    .. deprecated:: superseded by :func:`flat_band_terminal_mask` + ρ=0
        neutralization in the deep terminals. The whole-contact flat clamp
        suppressed the emitter accumulation layer and gave the wrong sign for
        the peak-current change vs the fixed-bias baseline (Akkala thesis
        Fig 3.8: Hartree peak should be *higher* than the no-charge baseline).
        Kept for the tests that exercised this earlier path.
    """
    Np = contact_mask.size
    vals = np.zeros(Np)
    idx = np.where(contact_mask)[0]
    if idx.size == 0:
        return contact_mask.copy(), vals
    # Emitter block: contiguous run from the first contact site.
    p = idx[0]
    emitter = []
    for k in idx:
        if k == p:
            emitter.append(k)
            p += 1
        else:
            break
    # Collector block: contiguous run back from the last contact site.
    p = idx[-1]
    collector = []
    for k in idx[::-1]:
        if k == p:
            collector.append(k)
            p -= 1
        else:
            break
    vals[np.asarray(emitter)] = +V / 2.0
    vals[np.asarray(collector)] = -V / 2.0
    return contact_mask.copy(), vals


def flat_contact_profile(V: float, contact_mask: np.ndarray,
                         dirichlet_vals: np.ndarray) -> np.ndarray:
    """Initial U: contacts flat at ±V/2, linear ramp across the active region."""
    U = dirichlet_vals.copy()
    active = ~contact_mask
    if active.any():
        ai = np.where(active)[0]
        U[active] = np.linspace(+V / 2.0, -V / 2.0, ai.size + 2)[1:-1]
    return U


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
    contact_mask: np.ndarray | None = None,
    NS: int = 1,
    ND: int = 1,
    scba_max_iter: int = 10,
    scba_mix: float = 0.4,
    scba_tol: float = 1e-5,
    eta: float = 1e-12,
    poisson_max_iter: int = 40,
    poisson_tol: float = 5e-3,
    kT_screen: float = 0.002,
    outer_clamp_sites: int = 15,
    bc_scheme: str = "clamp",
    density_mode: str = "calibrated",
    density_E_grid: np.ndarray | None = None,
    m_eff_kg: float | None = None,
    U_init: np.ndarray | None = None,
    final_scba: bool = True,
    density_offset: np.ndarray | None = None,
    density_correction=None,
    density_method: str = "banded",
    verbose: bool = False,
) -> SelfConsistentResult:
    """Self-consistent Poisson–NEGF loop at a single applied bias ``V``
    (Gummel–Hartree, following Akkala 2011 thesis Eq 3.50–3.53).

    Dirichlet BC at the two outer ends only: ``U[0] = +V/2``, ``U[Np-1] = −V/2``.
    Each iteration takes a **semiclassical Thomas–Fermi Newton-Raphson Poisson
    step** with the **NEGF quantum charge** as the source. The Jacobian uses
    ``dn/dU ≈ −n/kT_screen`` (Eq 3.53 in the Boltzmann limit) — well-conditioned
    (negative-definite → strictly stable descent) and the textbook standard for
    NEGF–Poisson RTD self-consistency.

    **Charge model (thesis Fig 3.1/3.3, set by** ``contact_mask`` **).** Each
    doped contact splits into three zones along its 10 nm depth:

    1. **Outer Dirichlet clamp** (outermost ``outer_clamp_sites``) — held at
       ±V/2 to extend the "true reservoir" depth, which is necessary for the
       short SISPAD ZnO contact where natural F_{1/2} screening over 10 nm
       can't fully pin the deep contact at high bias.
    2. **F_{1/2} screening zone** (middle of the contact) — semiclassical
       Fermi-Dirac charge ``n = Nc_FD · F_{1/2}((E_F±V/2 − U)/kT)`` (thesis
       Eq 3.1), calibrated so ``n=N_D`` at the rail. Provides Thomas-Fermi
       screening that smooths the transition into the active region.
    3. **Inner quantum buffer** (innermost ``quantum_buffer_sites``, ≈3 nm
       next to each barrier) — NEGF quantum charge, where the **emitter
       accumulation / triangular-quasi-bound state** forms. That accumulation
       is what raises the resonant peak (Akkala Fig 3.8).

    Convergence: magnitude of the actual update ``max|U_new − U|``;
    ``kT_screen`` (default 2 meV) damps the predictor; ``U_init`` warm-starts
    from the previous bias (else a linear drop).

    Returns
    -------
    SelfConsistentResult
    """
    Np = H_z.shape[0]
    dE = float(E_grid[1] - E_grid[0])
    U_L, U_R = +V / 2.0, -V / 2.0

    # Thesis (Akkala 2011) Gummel–Hartree charge model: Poisson over the WHOLE
    # domain with Dirichlet only at the two outer ends. The doped contacts split
    # into a deep **flat-band terminal** part (semiclassical Boltzmann/TF charge
    # n_TF=Nc·exp((E_F±V/2 − U)/kT), calibrated so n_TF=N_D at the rail) and an
    # inner **quantum-buffer** part nearest the barrier (NEGF charge, where the
    # emitter accumulation / triangular-well quasi-bound state forms). The TF
    # charge responds to U via Boltzmann → Thomas-Fermi screening pins the deep
    # contact at flat-band; this is the *correct* terminal treatment (the earlier
    # ρ=0 simplification removed the screening response and let bias drop freely
    # in the contacts; the flat-clamp before that killed accumulation entirely).
    if bc_scheme not in ("clamp", "neumann"):
        raise ValueError("bc_scheme must be 'clamp' or 'neumann'")
    # PI scheme (bc_scheme='neumann'): replace the pragmatic outer Dirichlet clamp
    # with a physical zero-field Neumann condition at the collector edge. Justified
    # because the field screens to zero over the sub-nm Thomas–Fermi length of the
    # degenerate n+ contact, so it has already vanished at the outer boundary. The
    # emitter endpoint stays Dirichlet at +V/2 as the gauge/bias reference; the
    # collector level is then free and *emerges* at ≈−V/2 from charge neutrality —
    # an internal validation that bias-carrying and neutrality are consistent.
    right_bc = "neumann" if bc_scheme == "neumann" else "dirichlet"

    em_term = col_term = None
    em_fd = col_fd = None
    d_mask = d_vals = None
    Nc_FD = 0.0
    if contact_mask is not None:
        em_term, col_term = terminal_masks(contact_mask)
        # Outer Dirichlet clamp on each contact (clamp scheme only): extends the
        # "true reservoir" depth; the middle of the contact keeps F_{1/2}
        # screening; the inner quantum-buffer keeps the NEGF charge for the
        # accumulation layer. The neumann scheme drops the clamp entirely.
        if bc_scheme == "clamp":
            em_clamp, col_clamp = outer_clamp_masks(contact_mask, outer_clamp_sites)
        else:
            em_clamp = col_clamp = np.zeros(Np, dtype=bool)
        em_fd = em_term & ~em_clamp
        col_fd = col_term & ~col_clamp
        if em_clamp.any() or col_clamp.any():
            d_mask = em_clamp | col_clamp
            d_vals = np.zeros(Np)
            d_vals[em_clamp] = +V / 2.0
            d_vals[col_clamp] = -V / 2.0
        # Nc_FD calibrated so n_TF = Nc_FD · F_{1/2}(Ef/kT) = N_D at the rail.
        N_D_contact = float(N_D[contact_mask].max())
        Nc_FD = N_D_contact / float(fermi_dirac_half(Ef / kT))

    if U_init is not None:
        U = U_init.copy()
    else:
        U = linear_bias_profile(V, Np, NS, ND)
    U[0] = U_L
    # Emitter endpoint is the gauge/bias reference in both schemes. The collector
    # endpoint is pinned only in the clamp scheme; under Neumann it is a free
    # (zero-field) node whose value emerges, so seed it but don't hold it.
    if right_bc == "dirichlet":
        U[-1] = U_R
    # Apply outer-clamp Dirichlet values to the warm start as well.
    if d_mask is not None:
        U[d_mask] = d_vals[d_mask]

    if density_mode not in ("calibrated", "physical"):
        raise ValueError("density_mode must be 'calibrated' or 'physical'")
    if density_mode == "physical" and m_eff_kg is None:
        raise ValueError("density_mode='physical' requires m_eff_kg")
    mu_L, mu_R = Ef + V / 2.0, Ef - V / 2.0

    def negf_and_density(U_profile, need_scba=True):
        # With a separate density grid the loop's density is ballistic and the
        # potential update never reads the SCBA result, so the SCBA is only
        # needed once, at the converged U (identical U trajectory and result).
        if density_mode == "physical" and density_E_grid is not None and not need_scba:
            return None, ballistic_transverse_density(
                density_E_grid, H_z, UB, U_profile, t0,
                mu_L, mu_R, kT, m_eff_kg, a_m, eta, method=density_method)
        r = run_rank1_keldysh_single_bias(
            V=V, E_grid=E_grid, H_z=H_z, UB=UB, bias_profile=U_profile, t0=t0,
            Ef=Ef, kT=kT, chi_diag=chi_diag, D0_sq_per_mode=D0_sq_per_mode,
            hnu_idx_per_mode=hnu_idx_per_mode, N_bose_per_mode=N_bose_per_mode,
            chi_per_mode=chi_per_mode, max_iter=scba_max_iter, tol=scba_tol,
            mix=scba_mix, eta=eta,
        )
        if density_mode == "physical":
            if density_E_grid is not None:
                # Evaluate the density on its own (finer) grid via the banded
                # ballistic route. The density integrand peaks at the band edge,
                # on the 1-D 1/sqrt(E) van Hove singularity, so it needs far
                # finer dE than the current does (which is converged to 0.3% on
                # the transport grid because f_L - f_R vanishes there). Phonons
                # shift the density by 0.003%, so dropping them here is safe.
                n = ballistic_transverse_density(
                    density_E_grid, H_z, UB, U_profile, t0,
                    mu_L, mu_R, kT, m_eff_kg, a_m, eta, method=density_method)
            else:
                n = physical_transverse_density(
                    r.G_R, r.Gam_L, r.Gam_R, E_grid, mu_L, mu_R, kT, m_eff_kg, a_m)
        else:
            n = density_prefactor * extract_density_1d(r.G_lesser, dE)
        return r, n

    res: Rank1KeldyshResult | None = None
    n_e = np.zeros(Np)
    converged = False
    dU = np.inf
    it = 0
    for it in range(1, poisson_max_iter + 1):
        res, n_q = negf_and_density(U, need_scba=False)
        # Thesis Eq 3.1 charge model assembled per-site:
        #   • deep emitter terminal: n_TF = Nc_FD · F_{1/2}((E_F+V/2 − U)/kT)
        #   • deep collector terminal: n_TF = Nc_FD · F_{1/2}((E_F−V/2 − U)/kT)
        #   • quantum region (inner contact buffer + barriers + well): NEGF n_q.
        # Jacobian dn/dU per region:
        #   • terminals: −Nc_FD · F_{−1/2}(η)/kT (exact Fermi-Dirac derivative)
        #   • quantum region: −n_q/kT_screen (Boltzmann predictor — Newton's
        #     approximation; the fixed point is unaffected).
        if density_offset is not None:
            # held fixed through this loop; a caller-side outer iteration
            # refreshes it (e.g. a scattering correction to the ballistic n)
            n_q = np.maximum(n_q + density_offset, 0.0)
        if density_correction is not None:
            # callable(U) -> additive density, re-evaluated every step, so a
            # correction defined relative to the ballistic density follows
            # the potential instead of lagging it
            n_q = np.maximum(n_q + density_correction(U), 0.0)
        n_e = n_q.copy()
        dn_dU = -n_q / kT_screen
        if em_fd is not None and em_fd.any():
            eta_em = (Ef + V / 2.0 - U[em_fd]) / kT
            f12 = fermi_dirac_half(eta_em)
            fmh = fermi_dirac_minus_half(eta_em)
            n_e[em_fd] = Nc_FD * f12
            dn_dU[em_fd] = -Nc_FD * fmh / kT
        if col_fd is not None and col_fd.any():
            eta_col = (Ef - V / 2.0 - U[col_fd]) / kT
            f12 = fermi_dirac_half(eta_col)
            fmh = fermi_dirac_minus_half(eta_col)
            n_e[col_fd] = Nc_FD * f12
            dn_dU[col_fd] = -Nc_FD * fmh / kT
        # Thesis Eq 3.50–3.53: Newton-Raphson Poisson step with the assembled
        # per-site charge and Jacobian; outer Dirichlet clamp on each contact.
        U_new = poisson_newton_update(
            U, eps_r, N_D, n_e, a_m, U_L, U_R, kT_screen,
            dn_dU_override=dn_dU,
            dirichlet_mask=d_mask, dirichlet_vals=d_vals,
            right_bc=right_bc,
        )
        dU = float(np.max(np.abs(U_new - U)))
        U = U_new
        if verbose:
            i_r = "   (no SCBA)  " if res is None else f"{res.I_right:+.3e}"
            print(f"    [poisson {it:2d}] step={dU*1e3:8.3f} meV  "
                  f"I_R={i_r} A  Umin={U.min():+.3f}  "
                  f"n_well_max={n_e.max():.2e} m⁻³", flush=True)
        if dU < poisson_tol:
            converged = True
            break

    # Final consistency solve so the returned observables match the returned U.
    # final_scba=False (density on its own grid only) skips it: result is None
    # and the caller works from U, e.g. a banded ballistic current.
    res, n_e = negf_and_density(U, need_scba=final_scba)

    return SelfConsistentResult(
        result=res, U=U, n_e=n_e, poisson_iters=it,
        poisson_converged=converged, dU_final=dU,
    )
