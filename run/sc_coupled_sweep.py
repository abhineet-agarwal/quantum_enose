"""Self-consistent I-V of the spacer device with phonon scattering, fully coupled.

The scattering density is re-solved inside every Poisson step
(sc_scatter_sweep.solve_bias_coupled), which replaced the outer-loop schemes:
those lagged (offset) or amplified (ratio) the correction near resonance, and
even where they "converged" they settled 7 meV / 17 % off the coupled
solution at 1184 mV. Each bias is warm-started from the previous one; the
run is resumable. A --seed lets a finer sweep start from a converged state
of an earlier run, keeping the sweep history (up/down) that matters inside a
bistable window.

Examples:
  python -u run/sc_coupled_sweep.py --mol Baseline --v1 1.6 --down-to 0.64
  python -u run/sc_coupled_sweep.py --mol Mol_A --v1 1.6
  python -u run/sc_coupled_sweep.py --mol Mol_A --v0 1.04 --v1 1.376 --dv 0.016 \\
      --seed results/2026-10-07/sc_coupled_D20.015_Mol_A_up.npz:1024 --tag fine
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np

import run.sc_scatter_sweep as sw

LO_MARGIN = 0.25
OUTDIR = os.path.abspath(os.path.join(os.path.dirname(sw.bc.OUT), "..", "2026-10-07"))


def out_path(args, direction):
    mol = "" if args.mol == "Baseline" else f"_{args.mol}"    # Baseline keeps its original name
    tag = f"_{args.tag}" if args.tag else ""
    return os.path.join(OUTDIR, f"sc_coupled_D2{args.d2:g}{mol}{tag}_{direction}.npz")


def sweep(args, direction, Vs, U, sigma, ph):
    out = out_path(args, direction)
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
                       converged=bool(info["poisson_conv"]), molecule=args.mol,
                       U_well=float(U[sw.bc.well_sites].mean()), I_scba_mode=final["I_mode"],
                       I_scba_err=final["I_err"], I_bal_mode=I_bal, I_dev_bal=I_dev_bal, D2=args.d2)
        np.savez(out, **{k: np.array(v, dtype=object) for k, v in db.items()})
        print(f"  [{args.mol} {direction}] V = {V*1e3:5.0f} mV  "
              f"{'conv' if info['poisson_conv'] else 'NOT CONV'} ({info['poisson_iters']} Poisson it)  "
              f"U_well {db[key]['U_well']*1e3:+8.2f} meV  I {final['I_mode']*1e9:9.3f} +- "
              f"{final['I_err']*1e9:6.3f} nA/mode{'' if final['conv'] else ' (final SCBA NC)'}  "
              f"[{info['secs']/60:.1f} min]", flush=True)
    return U, sigma


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mol", default="Baseline")
    ap.add_argument("--d2", type=float, default=0.015)
    ap.add_argument("--v0", type=float, default=0.0)
    ap.add_argument("--v1", type=float, default=1.6)
    ap.add_argument("--dv", type=float, default=0.032)
    ap.add_argument("--down-to", type=float, default=None,
                    help="after the up-sweep, sweep back down to this bias")
    ap.add_argument("--seed", default=None,
                    help="file.npz:key -- start from that converged state instead of 0 V")
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)
    ph = sw.phonons(args.d2, args.mol)
    if args.seed:
        f, key = args.seed.rsplit(":", 1)
        s = np.load(f, allow_pickle=True)[key].item()
        assert s["converged"], f"seed {args.seed} is not converged"
        U0 = s["U"]
    else:
        U0 = np.load(sw.bc.OUT, allow_pickle=True)["up_0"].item()["U"]  # ballistic SC at 0 V
    up = np.round(np.arange(args.v0, args.v1 + args.dv / 2, args.dv), 6)
    print(f"[setup] coupled scattering-Poisson, {args.mol}, D2={args.d2:g}, up {up[0]}-{up[-1]} V "
          f"step {args.dv*1e3:.0f} mV" + (f", down to {args.down_to} V" if args.down_to else "")
          + (f", seeded from {args.seed}" if args.seed else ""), flush=True)
    U, sigma = sweep(args, "up", up, U0, None, ph)
    if args.down_to is not None:
        down = np.round(np.arange(args.v1 - args.dv, args.down_to - args.dv / 2, -args.dv), 6)
        sweep(args, "down", down, U, sigma, ph)
    print("[done]", flush=True)


if __name__ == "__main__":
    main()
