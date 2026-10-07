"""Self-consistent I-V of the spacer device with phonon scattering, fully coupled.

The scattering density is re-solved inside every Poisson step
(sc_scatter_sweep.solve_bias_coupled), which replaced the outer-loop schemes:
those lagged (offset) or amplified (ratio) the correction near resonance, and
even where they "converged" they settled 7 meV / 17 % off the coupled
solution at 1184 mV. Up-sweep 0 -> V_hi, then down-sweep back to V_lo_down,
each bias warm-started from the previous one; resumable.

Usage: python -u run/sc_coupled_sweep.py [D2=0.015] [V_hi=1.6] [dV=0.032] [V_lo_down=0.64]
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
DV = float(sys.argv[3]) if len(sys.argv) > 3 else 0.032
V_LO_DOWN = float(sys.argv[4]) if len(sys.argv) > 4 else 0.64
LO_MARGIN = 0.25
OUTDIR = os.path.abspath(os.path.join(os.path.dirname(sw.bc.OUT), "..", "2026-10-07"))


def sweep(direction, Vs, U, sigma, ph):
    out = os.path.join(OUTDIR, f"sc_coupled_D2{D2:g}_{direction}.npz")
    db = ({k: np.load(out, allow_pickle=True)[k].item() for k in np.load(out).files}
          if os.path.exists(out) else {})
    for V in Vs:
        key = f"{V*1e3:.0f}"
        if key in db:
            U, sigma = db[key]["U"], None
            continue
        U, dn, sigma, info, final, I_bal, I_dev_bal = sw.solve_bias_coupled(
            float(V), U, ph, LO_MARGIN, sigma)
        db[key] = dict(V=float(V), U=U, dn=dn, info=info, final_scba=final,
                       converged=bool(info["poisson_conv"]),
                       U_well=float(U[sw.bc.well_sites].mean()), I_scba_mode=final["I_mode"],
                       I_scba_err=final["I_err"], I_bal_mode=I_bal, I_dev_bal=I_dev_bal, D2=D2)
        np.savez(out, **{k: np.array(v, dtype=object) for k, v in db.items()})
        print(f"  [{direction}] V = {V*1e3:5.0f} mV  {'conv' if info['poisson_conv'] else 'NOT CONV'} "
              f"({info['poisson_iters']} Poisson it)  U_well {db[key]['U_well']*1e3:+8.2f} meV  "
              f"I {final['I_mode']*1e9:9.3f} +- {final['I_err']*1e9:6.3f} nA/mode"
              f"{'' if final['conv'] else ' (final SCBA NC)'}  [{info['secs']/60:.1f} min]", flush=True)
    return U, sigma


if __name__ == "__main__":
    os.makedirs(OUTDIR, exist_ok=True)
    ph = sw.phonons(D2)
    up = np.round(np.arange(0.0, V_HI + DV / 2, DV), 6)
    down = np.round(np.arange(V_HI - DV, V_LO_DOWN - DV / 2, -DV), 6)
    U0 = np.load(sw.bc.OUT, allow_pickle=True)["up_0"].item()["U"]   # ballistic SC at 0 V
    print(f"[setup] coupled scattering-Poisson, D2={D2:g}, up 0-{V_HI} V, down to {V_LO_DOWN} V, "
          f"step {DV*1e3:.0f} mV", flush=True)
    U, sigma = sweep("up", up, U0, None, ph)
    sweep("down", down, U, sigma, ph)
    print("[done]", flush=True)
