"""Does the molecule's electrostatic feedback change its fingerprint?

At selected biases, Baseline and the molecule are solved fully
self-consistently (coupled scattering-Poisson, paired from the same converged
Baseline potential, tight tolerances, plain predictor -- which matched Newton
to 0.03 meV at 960 mV) and the resulting dI is compared with the
fixed-potential dI of run/fixedU_fingerprint.py.

Usage: python -u run/feedback_check.py [Mol_A] [D2=0.015] [V_mV ...]
"""
from __future__ import annotations
import os, sys, time
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np
import run.sc_scatter_sweep as sw

MOL = sys.argv[1] if len(sys.argv) > 1 else "Mol_A"
D2 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.015
VS = [int(x) for x in sys.argv[3:]] or [320, 544, 576, 608, 960, 1472, 1504, 1536]
TIGHT = dict(scba_tol=1e-5, poisson_tol=1e-5, final_tol=1e-7, final_max_iter=4000)
D = os.path.abspath(os.path.join(os.path.dirname(sw.bc.OUT), "..", "2026-10-07"))
OUT = os.path.join(D, f"feedback_check_D2{D2:g}_{MOL}_validation_v2.npz")

if __name__ == "__main__":
    sw.bc.POISSON_MAX_ITER = 300
    base = np.load(os.path.join(D, f"sc_coupled_D2{D2:g}_up.npz"), allow_pickle=True)
    fixed = np.load(os.path.join(D, f"fixedU_fingerprint_D2{D2:g}_{MOL}.npz"), allow_pickle=True)
    db = {k: np.load(OUT, allow_pickle=True)[k].item() for k in np.load(OUT).files} if os.path.exists(OUT) else {}
    ph = {"Baseline": sw.phonons(D2, "Baseline"), MOL: sw.phonons(D2, MOL)}
    for Vm in VS:
        key = str(Vm)
        if key in db:
            continue
        V = Vm / 1e3; U0 = base[key].item()["U"]; t = time.time(); row = {}
        for m in ("Baseline", MOL):
            U, _, _, info, fin, _, _ = sw.solve_bias_coupled(V, U0, ph[m], 0.25, None, **TIGHT)
            row[m] = dict(U_well=float(U[sw.bc.well_sites].mean()), I=fin["I_mode"], I_err=fin["I_err"],
                          poisson_iters=info["poisson_iters"], poisson_conv=info["poisson_conv"],
                          dU_final=info["dU_final"], final_conv=fin["conv"],
                          numerical_valid=info["numerical_valid"], poisson_residual_eV=info["poisson_residual_eV"])
        dI = row[MOL]["I"] - row["Baseline"]["I"]; dIe = float(row[MOL]["I_err"] + row["Baseline"]["I_err"])
        f = fixed[key].item()
        db[key] = dict(V=V, rows=row, error_kind="sum_of_half_contact_defects_not_sigma", dI_sc=dI, dI_sc_err=dIe, dI_fixed=f["dI"], dI_fixed_err=f["dI_err"],
                       dU_well=row[MOL]["U_well"] - row["Baseline"]["U_well"], secs=time.time() - t)
        np.savez(OUT, **{k: np.array(v, dtype=object) for k, v in db.items()})
        conv = all(row[m]["numerical_valid"] for m in row)
        print(f"  V={Vm:5d} mV  dI self-consistent {dI*1e9:+8.4f} (contact half-defects sum {dIe*1e9:.4f}) nA ({100*dI/abs(row['Baseline']['I']):+6.2f} %)   "
              f"fixed-U {f['dI']*1e9:+8.4f} ({100*f['dI']/abs(f['Baseline']['I']):+6.2f} %)   "
              f"dU_well {db[key]['dU_well']*1e3:+6.2f} meV  {'conv' if conv else 'NOT CONV'} "
              f"(Poisson {row['Baseline']['poisson_iters']}/{row[MOL]['poisson_iters']} it)  [{db[key]['secs']/60:.1f} min]", flush=True)
    print("[done]", flush=True)
