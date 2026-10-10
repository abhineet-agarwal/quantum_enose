"""Molecule's direct transport fingerprint at a fixed self-consistent potential.

At each bias, Baseline and the molecule are solved by SCBA (RGF, tight final
tolerance) on the SAME potential: the converged coupled Baseline up-sweep's U.
No Poisson, so this omits the molecule's electrostatic feedback (at 960 mV a
tight paired solve gave dU_well = +0.19 meV) but it isolates where in bias the
IETS signal lives, at a precision the Poisson loop cannot yet reach (~0.3 meV
in U, ~1 % in I, comparable to the ~1 % molecular signal).

Usage: python -u run/fixedU_fingerprint.py [Mol_A] [D2=0.015]
"""
from __future__ import annotations
import os, sys, time
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np
import run.sc_scatter_sweep as sw

MOL = sys.argv[1] if len(sys.argv) > 1 else "Mol_A"
D2 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.015
D = os.path.abspath(os.path.join(os.path.dirname(sw.bc.OUT), "..", "2026-10-07"))
SRC = os.path.join(D, f"sc_coupled_D2{D2:g}_up.npz")
OUT = os.path.join(D, f"fixedU_fingerprint_D2{D2:g}_{MOL}.npz")

if __name__ == "__main__":
    src = np.load(SRC, allow_pickle=True)
    db = {k: np.load(OUT, allow_pickle=True)[k].item() for k in np.load(OUT).files} if os.path.exists(OUT) else {}
    ph = {"Baseline": sw.phonons(D2, "Baseline"), MOL: sw.phonons(D2, MOL)}
    for key in sorted(src.files, key=float):
        if key in db:
            continue
        V = float(key) / 1e3; U = src[key].item()["U"]; E = sw.window(V, 0.25); t = time.time(); row = {}
        for m in ("Baseline", MOL):
            r = sw.scba(V, U, E, ph[m], None, tol=1e-7, max_iter=4000)
            row[m] = dict(I=0.5 * (r.I_right - r.I_left), I_err=0.5 * abs(r.I_right + r.I_left),
                          iters=r.iters_used, conv=r.converged)
        dI = row[MOL]["I"] - row["Baseline"]["I"]
        db[key] = dict(V=V, Baseline=row["Baseline"], mol=row[MOL], dI=dI,
                       dI_err=float(np.hypot(row[MOL]["I_err"], row["Baseline"]["I_err"])), secs=time.time() - t)
        np.savez(OUT, **{k: np.array(v, dtype=object) for k, v in db.items()})
        print(f"  V={float(key):5.0f} mV  I_BL {row['Baseline']['I']*1e9:9.4f}  I_{MOL} {row[MOL]['I']*1e9:9.4f} nA  "
              f"dI {dI*1e9:+8.4f} +- {db[key]['dI_err']*1e9:.4f} ({100*dI/max(abs(row['Baseline']['I']),1e-30):+6.2f} %)"
              f"{'' if row['Baseline']['conv'] and row[MOL]['conv'] else '  SCBA NC'}  [{db[key]['secs']:.0f}s]", flush=True)
    print("[done]", flush=True)
