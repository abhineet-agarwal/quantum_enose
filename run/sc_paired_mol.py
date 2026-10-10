"""Molecule vs Baseline in the coupled scattering-Poisson model, as paired solves.

At 960 mV the molecular change in current is ~1 % (Mol_A: +0.134 nA on
14.35 nA), while two separately converged sweeps at the default tolerances
differ by ~2 % (Poisson step 0.1 meV on a steep I-V) and the loose final
SCBA left +-0.08 nA of left/right mismatch. So at each bias both systems are
solved from the SAME starting potential (the Baseline's, from the previous
bias) with tight tolerances (Poisson 0.01 meV, inner SCBA 1e-5, final SCBA
1e-7): start-point error cancels in the difference and the current is
conserved to ~1e-4 nA.

Up-sweep only, matching the Baseline up-branch; inside the bistable window
both start from the Baseline's charged state. Grid: 32 mV, with 16 mV over
FINE_LO-FINE_HI around the main resonance. Resumable.

Usage: python -u run/sc_paired_mol.py [Mol_A] [D2=0.015]
"""
from __future__ import annotations

import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np

import run.sc_scatter_sweep as sw

MOL = sys.argv[1] if len(sys.argv) > 1 else "Mol_A"
D2 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.015
TIGHT = dict(scba_tol=1e-5, poisson_tol=1e-5, final_tol=1e-7, final_max_iter=3000)
FINE_LO, FINE_HI = 1.040, 1.376
LO_MARGIN = 0.25
OUTDIR = os.path.abspath(os.path.join(os.path.dirname(sw.bc.OUT), "..", "2026-10-07"))


def grid():
    a = np.arange(0.0, FINE_LO - 0.008, 0.032)
    b = np.arange(FINE_LO, FINE_HI + 0.008, 0.016)
    c = np.arange(b[-1] + 0.032, 1.6 + 0.016, 0.032)
    return np.round(np.concatenate([a, b, c]), 6)


def main():
    out = os.path.join(OUTDIR, f"sc_paired_D2{D2:g}_{MOL}.npz")
    db = ({k: np.load(out, allow_pickle=True)[k].item() for k in np.load(out).files}
          if os.path.exists(out) else {})
    ph = {"Baseline": sw.phonons(D2, "Baseline"), MOL: sw.phonons(D2, MOL)}
    Vs = grid()
    print(f"[setup] paired {MOL} vs Baseline, D2={D2:g}, {len(Vs)} biases "
          f"(16 mV over {FINE_LO}-{FINE_HI} V), tolerances {TIGHT}, {len(db)} done", flush=True)
    U = np.load(os.path.join(OUTDIR, f"sc_coupled_D2{D2:g}_up.npz"), allow_pickle=True)["0"].item()["U"]
    sig = {"Baseline": None, MOL: None}
    for V in Vs:
        key = f"{V*1e3:.0f}"
        if key in db:
            U = db[key]["Baseline"]["U"]
            sig = {"Baseline": None, MOL: None}
            continue
        t = time.time()
        row = {}
        for m in ("Baseline", MOL):
            Um, _, s, info, fin, Ib, Idb = sw.solve_bias_coupled(
                float(V), U, ph[m], LO_MARGIN, sig[m], **TIGHT)
            sig[m] = s
            row[m] = dict(U=Um, U_well=float(Um[sw.bc.well_sites].mean()), I=fin["I_mode"],
                          I_err=fin["I_err"], final=fin, info=info, I_bal=Ib, I_dev_bal=Idb,
                          converged=bool(info["poisson_conv"]) and bool(fin["conv"]))
        dI = row[MOL]["I"] - row["Baseline"]["I"]
        dIe = float(np.hypot(row[MOL]["I_err"], row["Baseline"]["I_err"]))
        db[key] = dict(V=float(V), Baseline=row["Baseline"], mol=row[MOL], molecule=MOL,
                       dI=dI, dI_err=dIe, D2=D2, tolerances=TIGHT, secs=time.time() - t)
        U = row["Baseline"]["U"]
        np.savez(out, **{k: np.array(v, dtype=object) for k, v in db.items()})
        ok = "" if (row["Baseline"]["converged"] and row[MOL]["converged"]) else "  NOT CONVERGED"
        print(f"  V = {V*1e3:5.0f} mV  I_BL {row['Baseline']['I']*1e9:9.4f}  I_{MOL} {row[MOL]['I']*1e9:9.4f} nA  "
              f"dI {dI*1e9:+8.4f} +- {dIe*1e9:.4f} nA ({100*dI/max(abs(row['Baseline']['I']),1e-30):+6.2f} %)  "
              f"dU_well {(row[MOL]['U_well']-row['Baseline']['U_well'])*1e3:+6.2f} meV{ok}  "
              f"[{db[key]['secs']/60:.1f} min]", flush=True)
    print("[done]", flush=True)


if __name__ == "__main__":
    main()
