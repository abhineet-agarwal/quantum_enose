"""Figures for the IEEE Sensors 2026 paper (Figs. 3a-c and 4).

Recovered from the inline script that made the 2026-09-14 versions (session
transcript), with two changes:

* currents are the transverse-integrated DEVICE current (`I_device`), not the
  single-mode `I_R` the earlier figures plotted;
* d2I/dV2 is the raw second difference np.diff(I, 2) / dV**2, which is what the
  paper states (Fig. 2: "diff(diff(I))/dV^2 -- no smoothing") and what the
  SISPAD figures used. The 2026-09-14 versions used np.gradient twice, i.e. a
  2*dV = 32 mV stencil that smooths the features the paper discusses.

Temperature data must come from a FIXED-DOPING sweep (E_F moving with T, see
--Ef in run_rank1_sweep.py); the 300 K file is shared with Fig. 3.

Usage:
  python run/make_ieee_sensors_figs.py --figs iv d2 deltaD temp --out DIR
"""
from __future__ import annotations

import argparse
import glob
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA_300 = "results/2026-09-26"
DATA_T = {10: "results/2026-09-30", 77: "results/2026-09-30",
          150: "results/2026-09-30", 300: "results/2026-09-26"}
MOLS = [("Baseline", "Baseline", "#444444", "-"),
        ("Mol_A", "Mol\\_A (100 meV)", "#1f77b4", "-"),
        ("Mol_B", "Mol\\_B (180 meV)", "#d62728", "-"),
        ("Mol_AB", "Mol\\_AB", "#2ca02c", "--")]
TEMPS = [(10, "#1f77b4"), (77, "#2ca02c"), (150, "#ff7f0e"), (300, "#d62728")]
LW = 1.5
FIELD = "I_device"
ISCALE, IUNIT = 1e3, "mA"

plt.rcParams.update({"font.size": 8, "axes.labelsize": 8.5, "legend.fontsize": 6.6,
                     "xtick.labelsize": 7.5, "ytick.labelsize": 7.5})


def load(data_dir, mol, T_K):
    files = sorted(glob.glob(os.path.join(
        data_dir, f"iets_ZnO_MgZnO_symmetric_{mol}_0-800mV_{T_K}K_rank1scba_*.npz")))
    if not files:
        raise SystemExit(f"no {mol} {T_K} K sweep in {data_dir}")
    a = np.load(files[-1], allow_pickle=True)
    if not bool(np.all(a["converged"])):
        raise SystemExit(f"{files[-1]}: not every bias point converged")
    return np.asarray(a["V"]), np.asarray(a[FIELD])


def d2(V, I):
    """Raw second difference, no smoothing; plotted at the interior points."""
    dV = V[1] - V[0]
    return V[1:-1], np.diff(I, 2) / dV ** 2


def fig_iv(out):
    f, a = plt.subplots(figsize=(4.0, 2.9))
    for m, lab, c, ls in MOLS:
        V, I = load(DATA_300, m, 300)
        a.plot(V * 1e3, I * ISCALE, ls, color=c, lw=LW, label=lab)
    a.set_xlabel("Bias Voltage (mV)"); a.set_ylabel(f"Current ({IUNIT})")
    a.legend(frameon=False); a.grid(alpha=.25); f.tight_layout(pad=.3)
    f.savefig(os.path.join(out, "fig4_IV.png"), dpi=300); plt.close(f)


def fig_d2(out):
    f, a = plt.subplots(figsize=(4.0, 2.9))
    for m, lab, c, ls in MOLS:
        Vc, dd = d2(*load(DATA_300, m, 300))
        a.plot(Vc * 1e3, dd, ls, color=c, lw=LW, label=lab)
    a.axhline(0, color="k", lw=.5)
    a.set_xlabel("Bias Voltage (mV)"); a.set_ylabel(r"$d^2I/dV^2$ (A/V$^2$)")
    a.legend(frameon=False); a.grid(alpha=.25); f.tight_layout(pad=.3)
    f.savefig(os.path.join(out, "fig5_d2IdV2.png"), dpi=300); plt.close(f)


def fig_deltaD(out):
    Vc, base = d2(*load(DATA_300, "Baseline", 300))
    f, a = plt.subplots(figsize=(7.0, 2.5))
    for m, lab, c, ls in MOLS[1:]:
        _, dd = d2(*load(DATA_300, m, 300))
        a.plot(Vc * 1e3, dd - base, ls, color=c, lw=LW, label=lab)
    a.axhline(0, color="k", lw=.5)
    a.set_xlabel("Bias Voltage (mV)"); a.set_ylabel(r"$\Delta D$ (A/V$^2$)")
    a.legend(frameon=False, ncol=3); a.grid(alpha=.25); f.tight_layout(pad=.3)
    f.savefig(os.path.join(out, "fig7_deltaD.png"), dpi=300); plt.close(f)


def fig_temp(out):
    f, ax = plt.subplots(1, 2, figsize=(7.0, 2.7))
    for T, c in TEMPS:
        V, I = load(DATA_T[T], "Mol_A", T)
        ax[0].plot(V * 1e3, I * ISCALE, "-", color=c, lw=LW, label=f"{T}K")
        Vc, dd = d2(V, I)
        ax[1].plot(Vc * 1e3, dd, "-", color=c, lw=LW, label=f"{T}K")
    ax[0].set_xlabel("Bias (mV)"); ax[0].set_ylabel(f"$I$ ({IUNIT})")
    ax[0].set_title("(a)", fontsize=9)
    ax[1].axhline(0, color="k", lw=.5)
    ax[1].set_xlabel("Bias (mV)"); ax[1].set_ylabel(r"$d^2I/dV^2$ (A/V$^2$)")
    ax[1].set_title("(b)", fontsize=9)
    # the caption refers to this inset (0-260 mV), invisible on the main axes
    ins = ax[1].inset_axes([0.10, 0.58, 0.40, 0.38])
    for T, c in TEMPS:
        Vc, dd = d2(*load(DATA_T[T], "Mol_A", T))
        m = Vc <= 0.26
        ins.plot(Vc[m] * 1e3, dd[m], "-", color=c, lw=1.0)
    ins.axhline(0, color="k", lw=.4); ins.tick_params(labelsize=5.5)
    ins.set_title(r"0–260 mV ($\times$ zoom)", fontsize=5.8)
    for a in ax:
        a.legend(frameon=False); a.grid(alpha=.25)
    f.tight_layout(pad=.3)
    f.savefig(os.path.join(out, "temp.png"), dpi=300); plt.close(f)


FIGS = {"iv": fig_iv, "d2": fig_d2, "deltaD": fig_deltaD, "temp": fig_temp}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--figs", nargs="+", choices=sorted(FIGS), default=sorted(FIGS))
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    for name in args.figs:
        FIGS[name](args.out)
        print(f"  wrote {name}")


if __name__ == "__main__":
    main()
