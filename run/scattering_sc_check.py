"""Does the self-consistent bistability survive phonon scattering in the spacers?

Starts from a branch's converged ballistic potential (run/bistability_check.py)
and adds scattering to the Poisson density by an outer iteration:

  1. at the current U, solve the injection window [V/2 - 0.25 eV, mu_L + 12 kT]
     with phonons off and on (fast solver; the full grid does not fit in 16 GB);
  2. turn the per-mode change into a transverse-integrated one site by site,
     dn(z) = n_window(z) * (n_SCBA(z) / n_ballistic(z) - 1), assuming each
     transverse slice changes by the same fraction;
  3. re-converge Poisson with the fine-grid ballistic density + dn;
  4. repeat until the correction is self-consistent (5 %) and U has settled.

Bulk LO scattering acts on every undoped site (spacers, barriers, well); the
doped leads stay out (Fix 7). If the two branches converge to the same U, the
bistability was an artifact of scattering-free spacers.

Usage: python -u run/scattering_sc_check.py <V_volts> <up|down> [D2_bulk]
"""
from __future__ import annotations

import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
_argv, sys.argv = sys.argv, sys.argv[:1]          # bistability_check reads argv at import
import numpy as np

import run.bistability_check as bc
from core.poisson import ballistic_transmission, ballistic_transverse_density
from core.poisson_negf import run_self_consistent_bias
from core.scba_rank1_keldysh import landauer_current_1mode, run_rank1_keldysh_single_bias_lowmem
from run.run_rank1_sweep import build_phonon_modes, emitter_barrier_center_nm, gaussian_chi

sys.argv = _argv
V = float(sys.argv[1])
BRANCH = sys.argv[2]
D2_BULK = float(sys.argv[3]) if len(sys.argv) > 3 else 0.001
# 2 meV left a ~1 meV-wide emitter-notch state at mu_L under-resolved: one grid
# energy carried the whole correction, which swung -13.6 % -> -2.9 % for a
# 1.4 meV shift in U. The chunked solver makes a finer grid affordable.
DE_SCBA = 0.0005
# Near a bistable point a small charge change moves U a lot; undamped, the
# correction overshot (U_well 255 -> 250 -> 282 -> 260 meV). Mix it in.
OUTER_TOL, OUTER_MAX, OUTER_MIX = 5e-4, 25, 0.3
OUT = os.path.join(os.path.dirname(bc.OUT), f"scatter_sc_{BRANCH}_{V*1e3:.0f}mV_D2{D2_BULK:g}.npz")

z = np.arange(bc.Np) * 0.2
chi_mol = gaussian_chi(z, emitter_barrier_center_nm(bc.build_stack(bc.DEVICE, bc.A_M)[5]), 0.3)
undoped = ~bc.cmask
D0, hi, Nb, chd, chl = build_phonon_modes("Baseline", DE_SCBA, bc.kT, chi_mol, 1, 1,
                                          UB=bc.UB, bulk_mask=undoped)
D0 = (D2_BULK,) + tuple(D0[1:])
REGIONS = {"emitter spacer": (z >= 30) & (z < 40),
           "well": np.isin(np.arange(bc.Np), bc.well_sites),
           "collector spacer": (z >= 47) & (z < 57)}


def window_solve(U, phonons):
    E = np.arange(V / 2 - 0.25, bc.EF + V / 2 + 12 * bc.kT, DE_SCBA)
    r = run_rank1_keldysh_single_bias_lowmem(
        V=V, E_grid=E, H_z=bc.H_z, UB=bc.UB, bias_profile=U, t0=bc.t0, Ef=bc.EF, kT=bc.kT,
        chi_diag=chd if phonons else np.zeros(bc.Np),
        D0_sq_per_mode=D0 if phonons else [], hnu_idx_per_mode=hi if phonons else [],
        N_bose_per_mode=Nb if phonons else [], chi_per_mode=chl if phonons else None,
        max_iter=100, tol=1e-4, mix=0.3, eta=1e-12)
    n = r.n_diag.sum(0) * DE_SCBA / (2 * np.pi)
    cons = abs(r.I_left + r.I_right) / max(abs(r.I_right), 1e-30)
    out = (E, n, r.I_right, r.converged, r.iters_used, cons)
    del r
    return out


def correction(U):
    E, n_b, _, _, _, _ = window_solve(U, False)
    _, n_s, I_s, conv, its, cons = window_solve(U, True)
    E_fine = np.arange(E[0], E[-1], bc.DE_DENSITY)
    n_win = ballistic_transverse_density(E_fine, bc.H_z, bc.UB, U, bc.t0, bc.EF + V / 2,
                                         bc.EF - V / 2, bc.kT, bc.m_eff, bc.A_M, 1e-12)
    ok = n_b > 1e-8 * n_b.max()
    ratio = np.where(ok, n_s / np.where(ok, n_b, 1.0), 1.0)
    return n_win * (ratio - 1.0), n_win, I_s, conv, its, cons


def poisson(U0, offset):
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
        density_offset=offset)
    T = ballistic_transmission(E_current, bc.H_z, bc.UB, sc.U, bc.t0)
    I_bal = landauer_current_1mode(E_current, T, bc.EF + V / 2, bc.EF - V / 2, bc.kT)
    return sc, I_bal


if __name__ == "__main__":
    start = np.load(bc.OUT, allow_pickle=True)[f"{BRANCH}_{V*1e3:.0f}"].item()
    U = start["U"]
    print(f"[{BRANCH} {V*1e3:.0f} mV] D2_bulk={D2_BULK:g} on {int(undoped.sum())} undoped sites | "
          f"start U_well {start['U_well']*1e3:+.2f} meV, I_ballistic {start['I_mode']*1e9:.3f} nA", flush=True)
    hist = []
    dn = np.zeros(bc.Np)
    dU = np.inf
    for k in range(1, OUTER_MAX + 1):
        t = time.time()
        dn_new, n_win, I_s, conv, its, cons = correction(U)
        # Converged only when the correction is self-consistent, not merely
        # when the damped step is small (damping makes every step small).
        resid = float(np.abs(dn_new - dn).sum() / max(np.abs(dn_new).sum(), 1e-300))
        if resid < 0.05 and dU < OUTER_TOL:
            print(f"  outer {k}: correction self-consistent (residual {resid:.1%}), U settled", flush=True)
            break
        dn = OUTER_MIX * dn_new + (1.0 - OUTER_MIX) * dn
        sc, I_bal = poisson(U, dn)
        dU = float(np.max(np.abs(sc.U - U)))
        U = sc.U
        row = dict(k=k, dU=dU, resid=resid, U_well=float(U[bc.well_sites].mean()), I_scba_mode=I_s,
                   I_bal_mode=I_bal, scba_conv=conv, scba_iters=its, scba_cons=cons,
                   poisson_conv=sc.poisson_converged, poisson_iters=sc.poisson_iters,
                   **{f"dn_frac_{name}": float(dn_new[m].sum() / max(n_win[m].sum(), 1e-300))
                      for name, m in REGIONS.items()})
        hist.append(row)
        np.savez(OUT, U=U, dn=dn, hist=np.array(hist, dtype=object), V=V, branch=BRANCH,
                 D2_bulk=D2_BULK, U_start=start["U"])
        print(f"  outer {k}: resid {resid:5.1%}  dU {dU*1e3:7.3f} meV  U_well {row['U_well']*1e3:+8.2f} meV  "
              f"I_SCBA {I_s*1e9:8.3f} nA (conv {conv}, {its} it, cons {cons:.1e})  "
              f"I_bal {I_bal*1e9:8.3f} nA  Poisson {sc.poisson_iters} it{'' if sc.poisson_converged else ' NC'}  "
              f"[{time.time()-t:.0f}s]", flush=True)
        print("           scattering adds " + ", ".join(
            f"{name} {100*row['dn_frac_'+name]:+.1f}%" for name in REGIONS), flush=True)
    print(f"[done] {BRANCH} {V*1e3:.0f} mV: U_well {start['U_well']*1e3:+.2f} -> "
          f"{U[bc.well_sites].mean()*1e3:+.2f} meV", flush=True)
