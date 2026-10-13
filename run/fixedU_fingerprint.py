"""Molecular current at specified fixed, archived Baseline potentials.

This checks transport only; it does not certify the supplied Poisson potential,
energy/bias refinement, or a device-level sensitivity. Failed pairs are saved
as diagnostics and cause a nonzero exit. Every invocation creates fresh output.

Usage:
  python run/fixedU_fingerprint.py Mol_A .015 --biases 576 1504 --dE .0005
  python run/fixedU_fingerprint.py Mol_A .015 --biases 576 --dE .0001 --tol 1e-9
"""
from __future__ import annotations
import argparse
import datetime as dt
import hashlib
import json
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import run.sc_scatter_sweep as sw
from core.validation import current_diagnostics


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("molecule", nargs="?", default="Mol_A")
    ap.add_argument("D2", nargs="?", type=float, default=.015)
    ap.add_argument("--source", type=Path)
    ap.add_argument("--output", type=Path)
    ap.add_argument("--biases", nargs="+", type=int, help="selected biases in mV")
    ap.add_argument("--dE", type=float, default=.0005, help="energy spacing in eV")
    ap.add_argument("--tol", type=float, default=1e-9)
    ap.add_argument("--max-iter", type=int, default=4000)
    ap.add_argument("--anderson-depth", type=int, default=0,
                    help="0: linear mixing (default); >=2: guarded Anderson history")
    ap.add_argument("--lo-margin", type=float, default=.25)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    if args.dE <= 0 or args.tol <= 0 or args.max_iter < 1 or args.lo_margin <= 0:
        ap.error("Require positive dE, tol, max-iter and lo-margin")
    source = args.source or Path("results/2026-10-07") / f"sc_coupled_D2{args.D2:g}_up.npz"
    output = args.output or Path("results") / dt.date.today().isoformat() / (
        "fixedU_fingerprint_" + dt.datetime.now().strftime("%H%M%S_%f") + ".json")
    output.parent.mkdir(parents=True, exist_ok=True)
    settings = vars(args).copy()
    settings.update(source=str(source), output=str(output), device=sw.bc.DEVICE,
                    Ef=sw.bc.EF, T=sw.bc.T_K,
                    current_rtol=.002, current_atol=1e-14, solver_validation_version=3)
    report = dict(status="running", scope="transport at fixed archived potentials only",
                  source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  settings=settings, currents=[], pairs=[], potentials={},
                  spectrum_grid_validated=False)
    # Reserve the output before any solve; never overwrite a prior run or resume failures.
    with output.open("x") as f:
        json.dump(report, f, indent=2)
    sw.DE_SCBA = args.dE
    ph = {m: sw.phonons(args.D2, m) for m in ("Baseline", args.molecule)}
    with np.load(source, allow_pickle=True) as src:
        keys = [str(v) for v in args.biases] if args.biases else sorted(src.files, key=float)
        for key in keys:
            state = src[key].item()
            U = state["U"]; V = float(key)/1000
            report["potentials"][key] = dict(U_eV=U.tolist(), source_info=state.get("info", {}))
            E = sw.window(V, args.lo_margin)
            pair = {}
            for mol in ("Baseline", args.molecule):
                start = time.time()
                r = sw.scba(V, U, E, ph[mol], tol=args.tol, max_iter=args.max_iter,
                            require_current=True, anderson_depth=args.anderson_depth)
                row = dict(V=V, dE=args.dE, molecule=mol, I_L=r.I_left, I_R=r.I_right,
                           I_mode=.5*(r.I_right-r.I_left), iters=r.iters_used,
                           converged=r.converged, physical_valid=r.physical_valid,
                           self_energy_residual=r.self_energy_residual,
                           mixer_fallbacks=r.mixer_fallbacks,
                           E_min=float(E[0]), E_max=float(E[-1]), seconds=time.time()-start,
                           **{k: v.item() for k, v in current_diagnostics(
                               r.I_left, r.I_right, r.converged).items()})
                report["currents"].append(row); pair[mol] = row
                output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
                if not args.quiet:
                    print(f"V={V:g} {mol} I={row['I_mode']*1e9:.9f} nA "
                          f"valid={row['current_valid']} imbalance={row['conservation_rel']:.3g}", flush=True)
            a, b = pair["Baseline"], pair[args.molecule]
            delta = b["I_mode"]-a["I_mode"]
            report["pairs"].append(dict(V=V, delta_I_A=delta,
                delta_I_percent=100*delta/abs(a["I_mode"]) if a["I_mode"] else None,
                numerical_valid=a["current_valid"] and b["current_valid"]))
    report["status"] = "complete"
    report["passed"] = all(p["numerical_valid"] for p in report["pairs"])
    output.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(f"{'PASS' if report['passed'] else 'FAIL'}: fixed-potential current checks; {output}")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
