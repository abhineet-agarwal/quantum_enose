"""
Generate the SISPAD -9 paper data figures (Fig 3a/b/c + ΔI), reconstructed.

This is the recovered "correct" generator. The original was overwritten on
2026-04-28 by a rebroken rewrite (the surviving run/regen_sispad_figs.py), which
used the 201-pt grid, the analytic d2I field, a serif/black color scheme, and
renumbered the outputs to fig4_IV..fig8. This script reproduces the figures as
they actually appear in SISPAD_2026_Abhineet-9.pdf:

  fig1_IV.png       I-V at 300 K              -> paper Fig 3a (Overleaf: fig4_IV)
  fig2_d2IdV2.png   numerical d2I/dV2 at 300K -> paper Fig 3b (Overleaf: fig5_d2IdV2)
  fig3_deltaI.png   ΔI = I_mol - I_BL  (nA)
  fig4_deltaD.png   ΔD discrimination metric  -> paper Fig 3c (Overleaf: fig7_deltaD)
  fig5_temp.png     Mol_A 10/77/150/300 K     -> paper Fig 4 (Overleaf: temp.png)

Method — consistent across ALL panels (intentionally; the published paper mixed
numerical Fig 3b with analytic Fig 4b, which we do not replicate):
  * 51-point 0-800 mV grid (16 mV spacing) for every figure, including the
    temperature sweep (Mol_A at 10/77/150 K re-run at 51 pts to match Fig 3).
  * default matplotlib colors (C0..C3 = blue/orange/green/red), with titles.
  * d2I/dV2 by RAW second difference  d2 = np.diff(I, 2) / dV**2  (NO smoothing,
    plotted at the interior bias points V[1:-1]). Matches fig2's +125 / -170
    µA/V² scale at 300 K; np.gradient**2 (the broken path) gives ±55.
The npz `d2I` field (per-bias analytic) is NOT used — it underpredicts at 300 K
and was the source of the earlier inconsistency.

Usage:
  python run/make_paper_figs.py [--data-dir DIR] [--out-dir DIR]
"""
from __future__ import annotations

import argparse
import glob
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MOLS = ["Baseline", "Mol_A", "Mol_B", "Mol_AB"]
IV_LABEL = {
    "Baseline": "Baseline (bulk LO)",
    "Mol_A": "Mol A (100 meV)",
    "Mol_B": "Mol B (180 meV)",
    "Mol_AB": "Mol AB (both)",
}
D_LABEL = {"Mol_A": r"$\Delta$Mol_A", "Mol_B": r"$\Delta$Mol_B", "Mol_AB": r"$\Delta$Mol_AB"}


def load_51pt(data_dir, mol):
    """Load the clean 51-point (16 mV) 300 K sweep for `mol` (any run date)."""
    best = None
    for f in sorted(glob.glob(os.path.join(
            data_dir, f"iets_ZnO_MgZnO_symmetric_{mol}_0-800mV_300K_rank1scba_*.npz"))):
        a = np.load(f, allow_pickle=True)
        if a["V"].size <= 60:          # 51-pt files; skip the 201-pt grid
            best = (a["V"], a["I_R"])
    return best


def numerical_d2(V, I):
    """Raw second difference d2I/dV2 (no smoothing); returns (V_interior, d2)."""
    dV = V[1] - V[0]
    return V[1:-1], np.diff(I, 2) / dV ** 2


TEMPS = [10, 77, 150, 300]
TEMP_LABEL = {10: "10 K", 77: "77 K", 150: "150 K", 300: "300 K"}


def load_temp_51pt(data_dir, T_K):
    """Load the 51-pt Mol_A sweep at temperature T_K (any run date)."""
    best = None
    for f in sorted(glob.glob(os.path.join(
            data_dir, f"iets_ZnO_MgZnO_symmetric_Mol_A_0-800mV_{T_K}K_rank1scba_*.npz"))):
        a = np.load(f, allow_pickle=True)
        if a["V"].size <= 60:          # 51-pt grid, consistent with Fig 3
            best = (a["V"], a["I_R"])
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="results/sispad_scba_2026-04-14")
    ap.add_argument("--out-dir", default="results/sispad_scba_2026-04-14")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    data = {m: load_51pt(args.data_dir, m) for m in MOLS}
    missing = [m for m, d in data.items() if d is None]
    if missing:
        raise SystemExit(f"missing 51-pt data for: {missing} in {args.data_dir}")

    def save(fig, stem):
        for ext in ("png", "pdf"):
            p = os.path.join(args.out_dir, f"{stem}.{ext}")
            fig.savefig(p, dpi=150)
            print(f"  wrote {p}")
        plt.close(fig)

    # Fig 1 — I-V
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, m in enumerate(MOLS):
        V, I = data[m]
        ax.plot(V * 1e3, I * 1e9, color=f"C{i}", label=IV_LABEL[m])
    ax.set_title("I-V at 300 K (SCBA, Anderson mixing)")
    ax.set_xlabel("Bias (mV)"); ax.set_ylabel("Current (nA)")
    ax.legend(loc="upper left"); ax.grid(alpha=0.3)
    fig.tight_layout(); save(fig, "fig1_IV")

    # Fig 2 — numerical d2I/dV2
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, m in enumerate(MOLS):
        V, I = data[m]
        Vc, d2 = numerical_d2(V, I)
        ax.plot(Vc * 1e3, d2, color=f"C{i}", label=IV_LABEL[m])
    ax.set_title(r"Numerical $d^2I/dV^2$ at 300 K (SCBA)")
    ax.set_xlabel("Bias (mV)"); ax.set_ylabel(r"$d^2I/dV^2$ (A/V$^2$)")
    ax.legend(loc="upper left"); ax.grid(alpha=0.3)
    fig.tight_layout(); save(fig, "fig2_d2IdV2")

    # Fig 3 — ΔI (nA)
    Vb, Ib = data["Baseline"]
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, m in enumerate(["Mol_A", "Mol_B", "Mol_AB"]):
        V, I = data[m]
        ax.plot(V * 1e3, (I - Ib) * 1e9, color=f"C{i}", label=D_LABEL[m])
    ax.set_title("Current difference from Baseline (SCBA)")
    ax.set_xlabel("Bias (mV)"); ax.set_ylabel(r"$\Delta$I (nA)")
    ax.legend(loc="upper left"); ax.grid(alpha=0.3)
    fig.tight_layout(); save(fig, "fig3_deltaI")

    # Fig 4 — ΔD discrimination
    Vc, d2b = numerical_d2(Vb, Ib)
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, m in enumerate(["Mol_A", "Mol_B", "Mol_AB"]):
        V, I = data[m]
        _, d2 = numerical_d2(V, I)
        ax.plot(Vc * 1e3, d2 - d2b, color=f"C{i}", label=D_LABEL[m])
    ax.set_title(r"IETS discrimination $\Delta$D(V) = d$^2$I$_{mol}$/dV$^2$ $-$ d$^2$I$_{BL}$/dV$^2$")
    ax.set_xlabel("Bias (mV)"); ax.set_ylabel(r"$\Delta$D (A/V$^2$)")
    ax.legend(loc="upper left"); ax.grid(alpha=0.3)
    fig.tight_layout(); save(fig, "fig4_deltaD")

    # Fig 5 — temperature dependence (Mol_A), 2-panel like paper Fig 4.
    # Same 51-pt grid + raw np.diff(I,2)/dV2 as Fig 3 (consistent, not analytic).
    temp = {T: load_temp_51pt(args.data_dir, T) for T in TEMPS}
    miss_T = [T for T, d in temp.items() if d is None]
    if miss_T:
        print(f"[warn] fig5_temp skipped — missing 51-pt Mol_A data at T={miss_T}K "
              f"(run run_rank1_sweep.py --T {miss_T} --V-points 51 --molecules Mol_A)")
    else:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
        for i, T in enumerate(TEMPS):
            V, I = temp[T]
            ax1.plot(V * 1e3, I * 1e9, color=f"C{i}", label=TEMP_LABEL[T])
            Vc, d2 = numerical_d2(V, I)
            ax2.plot(Vc * 1e3, d2, color=f"C{i}", label=TEMP_LABEL[T])
        ax1.set_title("Mol A: I-V vs Temperature")
        ax1.set_xlabel("Bias (mV)"); ax1.set_ylabel("Current (nA)")
        ax1.legend(loc="upper left"); ax1.grid(alpha=0.3)
        ax1.text(0.02, 0.98, "(a)", transform=ax1.transAxes, va="top", fontweight="bold")
        ax2.set_title(r"Mol A: numerical $d^2I/dV^2$ vs Temperature")
        ax2.set_xlabel("Bias (mV)"); ax2.set_ylabel(r"$d^2I/dV^2$ (A/V$^2$)")
        ax2.legend(loc="upper left"); ax2.grid(alpha=0.3)
        ax2.text(0.02, 0.98, "(b)", transform=ax2.transAxes, va="top", fontweight="bold")
        fig.tight_layout(); save(fig, "fig5_temp")


if __name__ == "__main__":
    main()
