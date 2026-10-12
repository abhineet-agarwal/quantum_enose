"""Check molecular delta-I and delta-d2I on independent energy/bias refinements.

Usage:
  python run/check_spectrum_convergence.py --centers 0.096 0.56 --molecules Mol_A

Uses the fixed-linear-potential, single-mode model. A pass applies only to the
tested stencils, molecule, window and tolerances; it does not validate a full
spectrum, transverse integration, or a device model. Exit 1 means the numerical
gate failed; the JSON records the measurements and each failed criterion.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np

from core.poisson import fermi_level_from_density
from core.scba_rank1_keldysh import run_rank1_keldysh_single_bias_lowmem
from core.validation import current_diagnostics, central_second_difference, refinement_agrees
from run.run_rank1_sweep import (
    build_stack, build_phonon_modes, gaussian_chi, emitter_barrier_center_nm,
    linear_bias_profile, _M0, MATERIALS,
)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--device", default="ZnO_MgZnO_symmetric")
    ap.add_argument("--molecules", nargs="+", default=["Mol_A"])
    ap.add_argument("--centers", nargs="+", type=float, default=[.096, .56])
    ap.add_argument("--energy-steps", nargs="+", type=float, default=[.002, .0005, .0001])
    ap.add_argument("--bias-steps", nargs="+", type=float, default=[.016, .008])
    ap.add_argument("--E-min", type=float, default=-.25)
    ap.add_argument("--E-max", type=float, default=.5)
    ap.add_argument("--T", type=float, default=300.)
    group = ap.add_mutually_exclusive_group()
    group.add_argument("--Ef", type=float, default=.02)
    group.add_argument("--donor-density", type=float, help="fixed ionized density, m^-3")
    ap.add_argument("--scba-tol", type=float, default=1e-7)
    ap.add_argument("--scba-max-iter", type=int, default=500)
    ap.add_argument("--anderson-depth", type=int, default=0,
                    help="SCBA mixer history; 0 = linear mixing. Anderson/DIIS returns "
                         "wrong-sign, unconverged currents at D^2 >= 0.01 and Mol_A has "
                         "D0^2 = 0.1; with it the 0.5 meV / 16 mV result here moved 6 %%")
    ap.add_argument("--rtol", type=float, default=.05)
    ap.add_argument("--delta-I-atol", type=float, default=1e-13, help="A")
    ap.add_argument("--delta-D-atol", type=float, default=1e-10, help="A/V^2")
    ap.add_argument("--output", type=Path, default=None)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    if any(len(v) < 2 or any(x <= 0 for x in v) or any(a <= b for a, b in zip(v, v[1:]))
           for v in [args.energy_steps, args.bias_steps]):
        ap.error("Energy and bias steps must each contain at least two decreasing positive values")
    if args.T <= 0 or args.E_max <= args.E_min or args.rtol <= 0:
        ap.error("Require T > 0, E-max > E-min and rtol > 0")
    output = args.output or Path("results") / ("validation_" + dt.datetime.now().strftime("%Y%m%d_%H%M%S_%f")) / "spectrum.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    # Reserve before doing expensive work; never replace a previous report.
    with output.open("x") as stream:
        json.dump({"status": "running"}, stream)
    H, UB, t0, Np, z, bounds = build_stack(args.device, .2e-9)
    chi = gaussian_chi(z, emitter_barrier_center_nm(bounds), .3)
    kt = .02585*args.T/300
    ef = args.Ef if args.donor_density is None else fermi_level_from_density(
        args.donor_density, kt, MATERIALS["ZnO"]["m_eff"]*_M0)
    biases = sorted({round(c+s*h, 12) for c in args.centers for h in args.bias_steps for s in (-1, 0, 1)})
    currents, rows, spectra = {}, [], []
    for de in args.energy_steps:
        E = np.arange(args.E_min, args.E_max + de/2, de)
        for mol in dict.fromkeys(["Baseline", *args.molecules]):
            D, hi, nb, ch, chl = build_phonon_modes(mol, de, kt, chi, 1, 1, UB=UB)
            for V in biases:
                r = run_rank1_keldysh_single_bias_lowmem(
                    V=V, E_grid=E, H_z=H, UB=UB, bias_profile=linear_bias_profile(V, Np, 1, 1),
                    t0=t0, Ef=ef, kT=kt, chi_diag=ch, D0_sq_per_mode=D,
                    hnu_idx_per_mode=hi, N_bose_per_mode=nb, chi_per_mode=chl,
                    max_iter=args.scba_max_iter, tol=args.scba_tol, mix=.3, method="rgf",
                    anderson_depth=args.anderson_depth)
                diag = current_diagnostics(r.I_left, r.I_right, r.converged)
                row = dict(dE=de, molecule=mol, V=V, I_L=r.I_left, I_R=r.I_right,
                           scba_converged=r.converged, iters=r.iters_used,
                           **{k: v.item() for k, v in diag.items()})
                rows.append(row)
                currents[de, mol, V] = r.I_right
                if not args.quiet:
                    print(f"dE={de*1e3:g} meV {mol} V={V:.4f}: {r.I_right*1e9:.6f} nA "
                          f"valid={row['current_valid']}", flush=True)
        for mol in args.molecules:
            for c in args.centers:
                for h in args.bias_steps:
                    delta = [currents[de, mol, round(c+s*h, 12)]
                             - currents[de, "Baseline", round(c+s*h, 12)] for s in (-1, 0, 1)]
                    spectra.append(dict(dE=de, molecule=mol, center=c, dV=h,
                                        delta_I_A=delta, delta_D_A_V2=float(central_second_difference(delta, h))))
    checks = []
    def compare(a, b, kind):
        current_ok = bool(np.all(refinement_agrees(a["delta_I_A"], b["delta_I_A"],
                                                   args.rtol, args.delta_I_atol))) if kind == "energy" else True
        spectrum_ok = bool(refinement_agrees(a["delta_D_A_V2"], b["delta_D_A_V2"],
                                             args.rtol, args.delta_D_atol))
        return dict(kind=kind, molecule=b["molecule"], center=b["center"],
                    coarse=a, fine=b, delta_I_passed=current_ok, delta_D_passed=spectrum_ok,
                    passed=current_ok and spectrum_ok)
    for mol in args.molecules:
        for c in args.centers:
            def find(de, dv):
                return next(x for x in spectra if x["molecule"] == mol and x["center"] == c
                            and x["dE"] == de and x["dV"] == dv)
            # Energy convergence at each tested voltage resolution.
            for h in args.bias_steps:
                checks.append(compare(find(args.energy_steps[-2], h), find(args.energy_steps[-1], h), "energy"))
            # Bias convergence at the finest tested energy resolution.
            checks.append(compare(find(args.energy_steps[-1], args.bias_steps[-2]),
                                  find(args.energy_steps[-1], args.bias_steps[-1]), "bias"))
    passed = all(r["current_valid"] for r in rows) and all(c["passed"] for c in checks)
    config = vars(args).copy(); config["output"] = str(output); config["Ef_used"] = ef
    report = dict(status="complete", passed=passed,
                  scope="selected fixed-potential stencils; single transverse mode per spin",
                  settings=config, currents=rows, spectra=spectra, checks=checks,
                  untested=["energy-window extent", "spatial grid", "transverse coupling/integration",
                            "other biases", "electrostatic model", "experimental detectability"])
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(f"{'PASS' if passed else 'FAIL'}: molecular-spectrum refinement checks; {output}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
