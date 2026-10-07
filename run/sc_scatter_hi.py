"""High-bias refinement of run/sc_scatter_sweep.py around the main resonance.

At 1152 mV with 30 % mixing the outer loop jumped between two states
(U_well ~270 and ~340-370 meV, current 15-185 nA) and never settled. Here:
10 % mixing, up to 40 outer iterations, 16 mV steps, starting from the
converged 1024 mV state; then the same range swept back down, so a genuine
bistability (both directions converge, to different states) can be told
apart from an iteration that will not settle.

Usage: python -u run/sc_scatter_hi.py [D2=0.015] [V_hi=1.6] [mix=0.1]
"""
from __future__ import annotations

import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np

import run.sc_scatter_sweep as sw

D2 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.015
V_HI = float(sys.argv[2]) if len(sys.argv) > 2 else 1.6
sw.OUTER_MIX = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1
sw.OUTER_MAX = 40
V_SEED, DV, LO_MARGIN = 1.024, 0.016, 0.25
BASE = os.path.join(os.path.dirname(sw.bc.OUT), "..", "2026-10-05", f"sc_scatter_sweep_D2{D2:g}.npz")


def run(direction, Vs, U, dn, ph):
    out = os.path.abspath(os.path.join(os.path.dirname(BASE),
                                       f"sc_scatter_hi_D2{D2:g}_{direction}.npz"))
    db = ({k: np.load(out, allow_pickle=True)[k].item() for k in np.load(out).files}
          if os.path.exists(out) else {})
    sigma = None
    for V in Vs:
        key = f"{V*1e3:.0f}"
        if key in db:
            U, dn, sigma = db[key]["U"], db[key]["dn"], None
            continue
        t = time.time()
        print(f"  [{direction}] V = {V*1e3:.0f} mV", flush=True)
        U, dn, sigma, hist, conv, I_bal, I_dev_bal, final = sw.solve_bias(
            float(V), U, dn, ph, LO_MARGIN, sigma)
        ratio = final["I_mode"] / max(hist[-1]["I_bal_win_mode"], 1e-30)
        db[key] = dict(V=float(V), U=U, dn=dn, hist=hist, converged=conv,
                       U_well=float(U[sw.bc.well_sites].mean()), I_scba_mode=final["I_mode"],
                       I_scba_err=final["I_err"], final_scba=final, I_bal_mode=I_bal,
                       I_dev_bal=I_dev_bal, I_dev_scba_est=I_dev_bal * ratio, D2=D2,
                       mix=sw.OUTER_MIX, secs=time.time() - t)
        np.savez(out, **{k: np.array(v, dtype=object) for k, v in db.items()})
        print(f"  [{direction}] V = {V*1e3:.0f} mV {'converged' if conv else 'NOT converged'} in "
              f"{len(hist)} outer: U_well {db[key]['U_well']*1e3:+.2f} meV  I_SCBA "
              f"{final['I_mode']*1e9:.3f} +- {final['I_err']*1e9:.3f} nA/mode  "
              f"[{db[key]['secs']/60:.1f} min]", flush=True)
    return U, dn


if __name__ == "__main__":
    seed = np.load(os.path.abspath(BASE), allow_pickle=True)[f"{V_SEED*1e3:.0f}"].item()
    assert seed["converged"], "seed bias must be converged"
    ph = sw.phonons(D2)
    up = np.round(np.arange(V_SEED + DV, V_HI + DV / 2, DV), 6)
    print(f"[setup] D2={D2:g}  mix={sw.OUTER_MIX}  outer max {sw.OUTER_MAX}  up {up[0]}-{up[-1]} V "
          f"step {DV}, seeded from converged {V_SEED*1e3:.0f} mV", flush=True)
    U, dn = run("up", up, seed["U"], seed["dn"], ph)
    run("down", up[::-1][1:], U, dn, ph)
    print("[done]", flush=True)
