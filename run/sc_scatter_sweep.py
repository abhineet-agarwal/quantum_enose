"""Self-consistent Poisson-NEGF I-V of the spacer device WITH phonon scattering.

Bulk LO scattering on every undoped site (spacers, barriers, well), at a
coupling calibrated to ZnO's Froehlich rate: matching the model broadening
D^2 (N+1) A_1D(E - hw) to the polar-optical emission rate gives
D^2 ~ 0.011-0.018 eV^2/site for electrons 25-50 meV above hw (alpha ~ 1.07),
10-20x the 0.001 used so far. Default 0.015.

Per bias, the outer loop of run/scattering_sc_check.py: the scattering change
to the density (low-memory SCBA on an injection window, 0.5 meV) is added to
the fine-grid ballistic Poisson density, damped, until the correction is
self-consistent. Warm starts: U and the correction from the previous bias; the
phonon self-energy from the previous outer iteration (and, aligned to the
emitter band edge, from the previous bias).

Usage: python -u run/sc_scatter_sweep.py [D2=0.015] [V0 V1 dV = 0 1.2 0.032] [lo_margin=0.25]
"""
from __future__ import annotations

import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np

import run.bistability_check as bc
from core.poisson import ballistic_transmission, ballistic_transverse_density
from core.poisson_negf import run_self_consistent_bias
from core.scba_rank1_keldysh import (landauer_current_1mode,
                                     run_rank1_keldysh_single_bias_lowmem,
                                     tsu_esaki_current)
from run.run_rank1_sweep import build_phonon_modes, emitter_barrier_center_nm, gaussian_chi

DE_SCBA = 0.0005
OUTER_TOL_U, OUTER_TOL_RES, OUTER_MAX, OUTER_MIX = 5e-4, 5e-3, 15, 0.3

z = np.arange(bc.Np) * 0.2
REGIONS = {"emitter spacer": (z >= 30) & (z < 40),
           "well": np.isin(np.arange(bc.Np), bc.well_sites),
           "collector spacer": (z >= 47) & (z < 57)}


def phonons(D2):
    chi_mol = gaussian_chi(z, emitter_barrier_center_nm(bc.build_stack(bc.DEVICE, bc.A_M)[5]), 0.3)
    D0, hi, Nb, chd, chl = build_phonon_modes("Baseline", DE_SCBA, bc.kT, chi_mol, 1, 1,
                                              UB=bc.UB, bulk_mask=~bc.cmask)
    return (D2,) + tuple(D0[1:]), hi, Nb, chd, chl


def window(V, lo_margin):
    return np.arange(V / 2 - lo_margin, bc.EF + V / 2 + 12 * bc.kT, DE_SCBA)


def scba(V, U, E, ph, sigma_init=None, tol=1e-3, max_iter=60):
    D0, hi, Nb, chd, chl = ph
    on = D0 is not None
    return run_rank1_keldysh_single_bias_lowmem(
        V=V, E_grid=E, H_z=bc.H_z, UB=bc.UB, bias_profile=U, t0=bc.t0, Ef=bc.EF, kT=bc.kT,
        chi_diag=chd if on else np.zeros(bc.Np), D0_sq_per_mode=D0 if on else [],
        hnu_idx_per_mode=hi if on else [], N_bose_per_mode=Nb if on else [],
        # Inexact inner solves: the outer loop is the fixed point that matters,
        # and the warm start carries progress between outer iterations. At
        # D2 = 0.015 a 1e-4 tolerance ran 150 iterations without settling;
        # a 1e-3 change in Sigma moves the density correction by ~0.1 %.
        # Linear mixing: Anderson/DIIS returns wrong-sign, unconverged currents
        # at D2 >= 0.01 (finding-patil-gamma-anderson); at 0.015 it gave
        # -0.08 nA at 320 mV where linear mixing converges to +0.22 nA.
        chi_per_mode=chl if on else None, max_iter=max_iter, tol=tol, mix=0.3, eta=1e-12,
        sigma_init=sigma_init, anderson_depth=0, method="rgf")


def correction(V, U, ph, lo_margin, sigma_init):
    """Scattering change to the transverse-integrated density at fixed U."""
    E = window(V, lo_margin)
    rb = scba(V, U, E, (None,) * 5)
    rs = scba(V, U, E, ph, sigma_init)
    n_b = rb.n_diag.sum(0) * DE_SCBA / (2 * np.pi)
    n_s = rs.n_diag.sum(0) * DE_SCBA / (2 * np.pi)
    n_win = ballistic_transverse_density(np.arange(E[0], E[-1], bc.DE_DENSITY), bc.H_z, bc.UB, U,
                                         bc.t0, bc.EF + V / 2, bc.EF - V / 2, bc.kT, bc.m_eff,
                                         bc.A_M, 1e-12, method="rgf")
    ok = n_b > 1e-8 * n_b.max()
    dn = n_win * (np.where(ok, n_s / np.where(ok, n_b, 1.0), 1.0) - 1.0)
    cons = abs(rs.I_left + rs.I_right) / max(abs(rs.I_right), 1e-30)
    return dn, n_win, rs, rb.I_right, cons


def poisson(V, U0, dn):
    E_density, E_current, E_transport = bc.energy_grids(V)
    sc = run_self_consistent_bias(
        V=V, E_grid=E_transport, H_z=bc.H_z, UB=bc.UB, t0=bc.t0, Ef=bc.EF, kT=bc.kT,
        chi_diag=np.zeros(bc.Np), D0_sq_per_mode=[], hnu_idx_per_mode=[],
        N_bose_per_mode=[], chi_per_mode=None, eps_r=bc.eps_r, N_D=bc.N_D, a_m=bc.A_M,
        density_prefactor=None, contact_mask=bc.cmask, NS=1, ND=1,
        scba_max_iter=1, scba_mix=0.3, scba_tol=1e-4, eta=1e-12,
        poisson_max_iter=bc.POISSON_MAX_ITER, poisson_tol=bc.POISSON_TOL,
        kT_screen=0.002, bc_scheme="neumann", density_mode="physical",
        m_eff_kg=bc.m_eff, density_E_grid=E_density, U_init=U0, final_scba=False,
        density_offset=dn, density_method="rgf")
    T = ballistic_transmission(E_current, bc.H_z, bc.UB, sc.U, bc.t0, method="rgf")
    mu_L, mu_R = bc.EF + V / 2, bc.EF - V / 2
    return (sc, landauer_current_1mode(E_current, T, mu_L, mu_R, bc.kT),
            tsu_esaki_current(E_current, T, mu_L, mu_R, bc.kT, bc.m_eff, bc.area))


def align(sig, NE):
    """Previous bias's self-energy on this bias's grid. Both grids start at
    V/2 - lo_margin, so index k is the same energy above the emitter band edge."""
    if sig is None:
        return None
    out = []
    for a in sig:
        b = np.zeros((NE, a.shape[1]), dtype=complex)
        n = min(NE, a.shape[0]); b[:n] = a[:n]
        out.append(b)
    return tuple(out)


def solve_bias(V, U, dn, ph, lo_margin, sigma_prev):
    sigma = align(sigma_prev, window(V, lo_margin).size)
    dU, hist = np.inf, []
    for k in range(1, OUTER_MAX + 1):
        t = time.time()
        dn_new, n_win, rs, I_b_win, cons = correction(V, U, ph, lo_margin, sigma)
        sigma = (rs.sigma_in_ph, rs.sigma_out_ph)
        resid = float(np.abs(dn_new - dn).sum() / max(n_win.sum(), 1e-300))
        frac = {name: float(dn_new[m].sum() / max(n_win[m].sum(), 1e-300)) for name, m in REGIONS.items()}
        row = dict(k=k, resid=resid, I_scba_mode=rs.I_right, I_bal_win_mode=I_b_win,
                   scba_iters=rs.iters_used, scba_conv=rs.converged, scba_cons=cons, **frac)
        if resid < OUTER_TOL_RES and dU < OUTER_TOL_U:
            hist.append(row)
            break
        dn = OUTER_MIX * dn_new + (1.0 - OUTER_MIX) * dn
        sc, I_bal, I_dev_bal = poisson(V, U, dn)
        dU = float(np.max(np.abs(sc.U - U)))
        U = sc.U
        row.update(dU=dU, poisson_iters=sc.poisson_iters, poisson_conv=sc.poisson_converged,
                   U_well=float(U[bc.well_sites].mean()), I_bal_mode=I_bal, I_dev_bal=I_dev_bal,
                   secs=time.time() - t)
        hist.append(row)
        print(f"    outer {k:2d}: resid {resid:6.2%}  dU {dU*1e3:7.3f} meV  U_well {row['U_well']*1e3:+8.2f}  "
              f"I_SCBA {rs.I_right*1e9:8.3f} nA  SCBA {rs.iters_used} it cons {cons:.0e}  "
              f"Poisson {sc.poisson_iters} it{'' if sc.poisson_converged else ' NC'}  [{row['secs']:.0f}s]", flush=True)
    converged = hist[-1]["resid"] < OUTER_TOL_RES and dU < OUTER_TOL_U
    # The density correction tolerates a loose inner solve; the current does
    # not: off resonance it is a small difference of large contact flows.
    rf = scba(V, U, window(V, lo_margin), ph, sigma, tol=1e-5, max_iter=300)
    final = dict(I_R=rf.I_right, I_L=rf.I_left, I_mode=0.5 * (rf.I_right - rf.I_left),
                 I_err=0.5 * abs(rf.I_right + rf.I_left), iters=rf.iters_used, conv=rf.converged)
    return U, dn, sigma, hist, converged, I_bal, I_dev_bal, final


def main():
    D2 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.015
    V0, V1, dV = (float(x) for x in sys.argv[2:5]) if len(sys.argv) > 4 else (0.0, 1.2, 0.032)
    lo_margin = float(sys.argv[5]) if len(sys.argv) > 5 else 0.25
    out = os.path.join(os.path.dirname(bc.OUT), "..", "2026-10-05", f"sc_scatter_sweep_D2{D2:g}.npz")
    out = os.path.abspath(out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    db = {k: np.load(out, allow_pickle=True)[k].item() for k in np.load(out).files} if os.path.exists(out) else {}
    ph = phonons(D2)
    Vs = np.round(np.arange(V0, V1 + dV / 2, dV), 6)
    print(f"[setup] {bc.DEVICE} D2_bulk={D2:g} eV^2 on {int((~bc.cmask).sum())} undoped sites, "
          f"SCBA window [V/2-{lo_margin}, mu_L+12kT] @ {DE_SCBA*1e3:g} meV, V {V0}-{V1} step {dV}, "
          f"{len(db)} biases already done", flush=True)
    U = dn = sigma = None
    for V in Vs:
        key = f"{V*1e3:.0f}"
        if key in db:
            U, dn, sigma = db[key]["U"], db[key]["dn"], None
            continue
        if U is None:   # cold start from the converged ballistic up-sweep at this bias
            U = np.load(bc.OUT, allow_pickle=True)[f"up_{key}"].item()["U"]
            dn = np.zeros(bc.Np)
        t = time.time()
        print(f"  V = {V*1e3:.0f} mV", flush=True)
        U, dn, sigma, hist, conv, I_bal, I_dev_bal, final = solve_bias(float(V), U, dn, ph, lo_margin, sigma)
        last = hist[-1]
        ratio = final["I_mode"] / max(last["I_bal_win_mode"], 1e-30)
        db[key] = dict(V=float(V), U=U, dn=dn, hist=hist, converged=conv,
                       U_well=float(U[bc.well_sites].mean()), I_scba_mode=final["I_mode"],
                       I_scba_err=final["I_err"], final_scba=final,
                       I_bal_mode=I_bal, I_dev_bal=I_dev_bal, I_dev_scba_est=I_dev_bal * ratio,
                       D2=D2, lo_margin=lo_margin, secs=time.time() - t)
        np.savez(out, **{k: np.array(v, dtype=object) for k, v in db.items()})
        print(f"  V = {V*1e3:.0f} mV {'converged' if conv else 'NOT converged'} in {len(hist)} outer: "
              f"U_well {db[key]['U_well']*1e3:+.2f} meV  I_SCBA {final['I_mode']*1e9:.3f} "
              f"+- {final['I_err']*1e9:.3f} nA/mode (final solve {final['iters']} it{'' if final['conv'] else ' NC'})  "
              f"I_dev~{db[key]['I_dev_scba_est']*1e3:.2f} mA  [{db[key]['secs']/60:.1f} min]", flush=True)
    print("[done]", flush=True)


if __name__ == "__main__":
    main()
