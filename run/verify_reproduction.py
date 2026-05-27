"""
Fast reproduction gate for the SISPAD -9 abstract.

Re-runs the rank-1 projected SCBA at a handful of bias points in the *current*
working tree and checks the right-contact current against the committed golden
.npz data. Confirms in <1 min that the simulation path (config + solver +
Anderson mixing) still reproduces the published numbers, without paying for a
full multi-hour sweep.

Authoritative result: SISPAD_2026_Abhineet-9.pdf
  - ZnO/Mg0.3Zn0.7O symmetric RTD, peak |I_R| ~ 69 nA at 560 mV (300 K).
  - Production params: 0-800 mV, mix=0.3, max_iter=100, tol=1e-4, Ef=20 meV,
    a=0.2 nm (Np=135), dE=2 meV, D2_bulk=0.001 eV2, molecular D0^2=0.1 eV2.

Usage:  python run/verify_reproduction.py
Exit code 0 if every checked point matches golden within 0.5 % (rel), else 1.
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import run.run_rank1_sweep as R
from core.scba_rank1_keldysh import run_rank1_keldysh_single_bias

GOLDEN_DIR = "results/sispad_scba_2026-04-14"
DEVICE = "ZnO_MgZnO_symmetric"
# (molecule, golden run-date tag, bias points to spot-check in V)
CHECKS = [
    ("Baseline", "2026-04-14", [0.32, 0.56]),  # 0.56 = NDR peak (~69 nA)
    ("Mol_A", "2026-04-14", [0.56]),
]
TOL_REL = 5e-3  # 0.5 %

# Production SCBA parameters (match the golden npz metadata).
A_NM = 0.2
T_K = 300.0
EF = 0.02
DE = 0.002
E_MIN, E_MAX = -0.25, 0.5
SCBA_MIX = 0.3
SCBA_MAX_ITER = 100
SCBA_TOL = 1e-4
SIGMA_MOL_NM = 0.3


def _setup(mol):
    a_m = A_NM * 1e-9
    kT = 0.02585 * (T_K / 300.0)
    H_z, UB, t0, Np, z_nm, bounds = R.build_stack(DEVICE, a_m)
    z0 = R.emitter_barrier_center_nm(bounds)
    chi_mol = R.gaussian_chi(z_nm, z0, SIGMA_MOL_NM)
    E_grid = np.arange(E_MIN, E_MAX + 0.5 * DE, DE)
    D0_sq, hnu_idx, N_bose, chi_def, chi_list = R.build_phonon_modes(
        mol, DE, kT, chi_mol, 1, 1, UB=UB)
    return dict(H_z=H_z, UB=UB, t0=t0, Np=Np, E_grid=E_grid, kT=kT,
                chi_def=chi_def, D0_sq=D0_sq, hnu_idx=hnu_idx,
                N_bose=N_bose, chi_list=chi_list)


def main() -> int:
    print(f"Reproduction gate: device={DEVICE}  tol={TOL_REL:.1%}\n")
    all_ok = True
    for mol, date, biases in CHECKS:
        gpath = os.path.join(
            GOLDEN_DIR,
            f"iets_{DEVICE}_{mol}_0-800mV_300K_rank1scba_{date}.npz")
        if not os.path.exists(gpath):
            print(f"[FAIL] golden missing: {gpath}")
            all_ok = False
            continue
        g = np.load(gpath, allow_pickle=True)
        Vref, IRref = g["V"], g["I_R"]
        s = _setup(mol)
        for V in biases:
            j = int(np.argmin(np.abs(Vref - V)))
            bp = R.linear_bias_profile(V, s["Np"], 1, 1)
            t = time.time()
            res = run_rank1_keldysh_single_bias(
                V=float(V), E_grid=s["E_grid"], H_z=s["H_z"], UB=s["UB"],
                bias_profile=bp, t0=s["t0"], Ef=EF, kT=s["kT"],
                chi_diag=s["chi_def"], D0_sq_per_mode=s["D0_sq"],
                hnu_idx_per_mode=s["hnu_idx"], N_bose_per_mode=s["N_bose"],
                chi_per_mode=s["chi_list"], max_iter=SCBA_MAX_ITER,
                tol=SCBA_TOL, mix=SCBA_MIX, eta=1e-12)
            rel = abs(res.I_right - IRref[j]) / max(abs(IRref[j]), 1e-30)
            ok = rel <= TOL_REL
            all_ok &= ok
            print(f"  [{'OK ' if ok else 'FAIL'}] {mol:9s} V={V:.3f}  "
                  f"I_R={res.I_right*1e9:+8.3f} nA  golden={IRref[j]*1e9:+8.3f} nA  "
                  f"rel={rel:.2e}  ({time.time()-t:.0f}s)")
    print("\n" + ("PASS — current tree reproduces the SISPAD -9 numbers."
                  if all_ok else "FAIL — see mismatches above."))
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
