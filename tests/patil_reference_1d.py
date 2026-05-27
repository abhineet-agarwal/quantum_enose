"""
Patil 1D reference port.

Line-by-line faithful port of ``papers/rtd2modes_1d.m`` — Patil's MATLAB
reference for "The role of inelastic scattering in resonant tunneling
heterostructures". Parameters are copied verbatim from the MATLAB file;
see ``docs/PATIL_BENCHMARK.md`` for the parameter table and physics notes.

This module exposes two entry points:

* ``run_patil_reference(...)`` — drives the full SCBA self-consistency and
  returns a ``PatilResult`` dataclass with I(V), dI/dV, d²I/dV², and the
  coherent/incoherent current decomposition.
* A CLI (``python -m tests.patil_reference_1d``) that runs with Patil's
  exact defaults and saves the reference curves to
  ``tests/patil_reference_1d.npz``, plus an IETS/I-V plot.

Design rule: **do not "modernize" the MATLAB algorithm**. Where Patil made
a choice (only-diagonal phonon self-energy, Newton-mixer with it=0.5, etc.)
we match it exactly. Any physics corrections we want to apply go in our
multimode solver, not here — here we preserve the reference.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from typing import Optional

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


# ---------------------------------------------------------------------------
# Physical constants (match MATLAB header exactly).
# ---------------------------------------------------------------------------

HBAR = 1.06e-34            # J·s
Q = 1.6e-19                # C
M_E = 9.1e-31              # kg
M_STAR_GAAS = 0.067 * M_E  # kg, GaAs effective mass
IE = (Q * Q) / (2 * np.pi * HBAR)   # current prefactor — MATLAB convention


# ---------------------------------------------------------------------------
# Result container.
# ---------------------------------------------------------------------------


@dataclass
class PatilResult:
    V: np.ndarray          # (NV,) bias grid in V
    I_total: np.ndarray    # (NV,) total current in A, MATLAB II
    I_coherent: np.ndarray # (NV,) coherent part, MATLAB IIco
    I_incoherent: np.ndarray  # (NV,) incoherent part, MATLAB IInonco
    dIdV: np.ndarray       # (NV-1,) finite-difference dI/dV
    d2IdV2: np.ndarray     # (NV-2,) finite-difference d²I/dV²
    V_dIdV: np.ndarray     # corresponding V axis for dI/dV
    V_d2IdV2: np.ndarray
    # Parameters used (for downstream provenance).
    params: dict


# ---------------------------------------------------------------------------
# The main driver.
# ---------------------------------------------------------------------------


def run_patil_reference(
    NS: int = 1,
    NC: int = 40,
    ND: int = 1,
    Nb: int = 15,
    Vb: float = 0.6,
    t0: float = 5.2,
    Ef: float = 0.02,
    kT: float = 0.025,
    NV: int = 76,
    V_max: float = 0.3,
    dE: float = 0.005,
    E_min: float = -0.2,
    E_max: float = 0.8,
    D0: float = 0.0,
    Dnu: tuple[float, ...] = (0.1, 0.1),
    hnu_units: tuple[int, ...] = (18, 35),
    scba_tol: float = 1e-5,
    mix_it: float = 0.5,
    max_outer_iter: int = 200,
    verbose: bool = True,
) -> PatilResult:
    """Run Patil's MATLAB RTD 1D reference.

    Default arguments reproduce the ``%%% one peak in TM vs E %%%`` block of
    ``papers/rtd2modes_1d.m`` (lines 22-26) and the phonon settings
    (lines 44-48). The returned ``PatilResult`` contains the SCBA-converged
    current and IETS curves.
    """
    # --- Grid setup (lines 22-40) ---
    Np = NS + NC + ND
    UB = np.concatenate([
        np.zeros(NS),
        Vb * np.ones(Nb),
        np.zeros(NC - 2 * Nb),
        Vb * np.ones(Nb),
        np.zeros(ND),
    ])
    assert UB.shape == (Np,), f"UB shape {UB.shape} != Np {Np}"

    # Tight-binding Hamiltonian T[i,j] = 2t0 δ_ij − t0 (δ_{i,j±1}) + UB·δ_ij.
    T = 2.0 * t0 * np.eye(Np)
    T -= t0 * np.eye(Np, k=1)
    T -= t0 * np.eye(Np, k=-1)
    T += np.diag(UB)

    # --- Bias grid (line 42) ---
    VV = np.linspace(0.0, V_max, NV)
    dV = VV[1] - VV[0] if NV > 1 else 0.0

    # --- Energy grid (line 55) ---
    E = np.arange(E_min, E_max + 0.5 * dE, dE)
    NE = E.size

    # --- Phonon mode setup (lines 44-48) ---
    Dnu_arr = np.asarray(Dnu, dtype=float)
    hnu_units_arr = np.asarray(hnu_units, dtype=int)
    Nph = Dnu_arr.size
    if hnu_units_arr.size != Nph:
        raise ValueError("hnu_units must match Dnu in length")
    Nhnu = 1.0 / (np.exp(dE * hnu_units_arr / kT) - 1.0)  # Bose-Einstein
    hnu_eV = dE * hnu_units_arr                           # diagnostic

    if verbose:
        print(f"Patil 1D reference")
        print(f"  Np = {Np} (NS={NS}, NC={NC}, ND={ND}), Nb = {Nb}, Vb = {Vb} eV")
        print(f"  t0 = {t0} eV, Ef = {Ef} eV, kT = {kT} eV")
        print(f"  NV = {NV}, V ∈ [0, {V_max}]")
        print(f"  NE = {NE}, E ∈ [{E_min}, {E_max}], dE = {dE}")
        print(f"  Phonons: D² = {Dnu_arr.tolist()} eV², ℏω = {hnu_eV.tolist()} eV")
        print(f"  SCBA tol = {scba_tol}, mix it = {mix_it}")

    # --- Output buffers ---
    II = np.zeros(NV)
    II3 = np.zeros(NV)     # MATLAB II3 — total phonon-related trace
    II4 = np.zeros(NV)     # MATLAB II4 — second contact current
    IIco = np.zeros(NV)    # Coherent current
    IInonco = np.zeros(NV) # Incoherent current

    zplus = 1j * 1e-12

    # --- Main V-loop (line 50) ---
    for iV, V in enumerate(VV):
        mu1 = Ef + V / 2
        mu2 = Ef - V / 2
        # Applied potential profile (line 52): V/2 on left, linear drop, −V/2 on right
        U1 = V * np.concatenate([
            0.5 * np.ones(NS),
            np.linspace(0.5, -0.5, NC),
            -0.5 * np.ones(ND),
        ])

        # Fermi functions on the energy grid.
        f1 = 1.0 / (1.0 + np.exp((E - mu1) / kT))
        f2 = 1.0 / (1.0 + np.exp((E - mu2) / kT))

        # Allocate per-energy arrays — shape (Np, Np, NE), matching MATLAB.
        sigin1 = np.zeros((Np, Np, NE), dtype=complex)
        sigout1 = np.zeros((Np, Np, NE), dtype=complex)
        sigin2 = np.zeros((Np, Np, NE), dtype=complex)
        sigout2 = np.zeros((Np, Np, NE), dtype=complex)
        sigin = np.zeros((Np, Np, NE), dtype=complex)
        sigout = np.zeros((Np, Np, NE), dtype=complex)
        n = np.zeros((Np, Np, NE), dtype=complex)
        p = np.zeros((Np, Np, NE), dtype=complex)
        gam1 = np.zeros((Np, Np, NE), dtype=complex)
        gam2 = np.zeros((Np, Np, NE), dtype=complex)
        gamp = np.zeros((Np, Np, NE), dtype=complex)
        G = np.zeros((Np, Np, NE), dtype=complex)
        A = np.zeros((Np, Np, NE), dtype=complex)

        # --- SCBA outer loop (line 82) ---
        change = 1.0
        outer = 0
        while change > scba_tol and outer < max_outer_iter:
            outer += 1

            # --- Inner energy loop (line 83) ---
            for k in range(NE):
                # Left contact surface self-energy.
                ck = 1.0 - (E[k] + zplus - U1[0] - UB[0]) / (2 * t0)
                ka = np.arccos(ck)
                sig1_loc = -t0 * np.exp(1j * ka)
                sig1 = np.zeros((Np, Np), dtype=complex)
                sig1[0, 0] = sig1_loc
                gam1[:, :, k] = 1j * (sig1 - sig1.conj().T)

                # Right contact surface self-energy.
                ck = 1.0 - (E[k] + zplus - U1[-1] - UB[-1]) / (2 * t0)
                ka = np.arccos(ck)
                sig2_loc = -t0 * np.exp(1j * ka)
                sig2 = np.zeros((Np, Np), dtype=complex)
                sig2[-1, -1] = sig2_loc
                gam2[:, :, k] = 1j * (sig2 - sig2.conj().T)

                # Contact in/out scattering.
                sigin1[:, :, k] = f1[k] * gam1[:, :, k]
                sigin2[:, :, k] = f2[k] * gam2[:, :, k]
                sigout1[:, :, k] = (1 - f1[k]) * gam1[:, :, k]
                sigout2[:, :, k] = (1 - f2[k]) * gam2[:, :, k]

                # Total phonon broadening (Fix-4 causal combination, Patil
                # uses Σ_ph^R = (i/2) (Σ_in + Σ_out) — this is the original
                # MATLAB sign convention; matches when fed through the
                # Green's function below).
                gamp[:, :, k] = sigin[:, :, k] + sigout[:, :, k]

                # Retarded Green's function at energy E[k].
                M = ((E[k] + zplus) * np.eye(Np)
                     - T - np.diag(U1) - sig1 - sig2
                     + 0.5j * gamp[:, :, k])
                G[:, :, k] = np.linalg.inv(M)

                A[:, :, k] = 1j * (G[:, :, k] - G[:, :, k].conj().T)

                # Electron correlation (line 95): n = G [f1 Γ1 + f2 Γ2 + Σ_in] G†
                Sigma_in_total = (f1[k] * gam1[:, :, k]
                                  + f2[k] * gam2[:, :, k]
                                  + sigin[:, :, k])
                n[:, :, k] = np.real(
                    G[:, :, k] @ Sigma_in_total @ G[:, :, k].conj().T
                )
                p[:, :, k] = A[:, :, k] - n[:, :, k]

            # --- Phonon self-energy update (lines 98-116) ---
            # Diagonal-only local Einstein vibron (MATLAB ne(nn,nn,...) pattern).
            siginnew = D0 * n
            sigoutnew = D0 * p

            for iph in range(Nph):
                inu = int(hnu_units_arr[iph])
                if inu >= NE:
                    continue
                ne_arr = np.zeros_like(n)
                na_arr = np.zeros_like(n)
                pe_arr = np.zeros_like(p)
                pa_arr = np.zeros_like(p)

                # Diagonal-only energy shift (MATLAB: ne(nn,nn,1:NE-inu) = n(nn,nn,inu+1:NE))
                for nn in range(Np):
                    ne_arr[nn, nn, : NE - inu] = n[nn, nn, inu:NE]
                    na_arr[nn, nn, inu:NE] = n[nn, nn, : NE - inu]
                    pe_arr[nn, nn, : NE - inu] = p[nn, nn, inu:NE]
                    pa_arr[nn, nn, inu:NE] = p[nn, nn, : NE - inu]

                siginnew += (Nhnu[iph] + 1) * Dnu_arr[iph] * ne_arr \
                          + Nhnu[iph] * Dnu_arr[iph] * na_arr
                sigoutnew += Nhnu[iph] * Dnu_arr[iph] * pe_arr \
                           + (Nhnu[iph] + 1) * Dnu_arr[iph] * pa_arr

            # Convergence check (MATLAB line 117-118).
            change = (np.sum(np.abs(siginnew - sigin))
                      + np.sum(np.abs(sigoutnew - sigout)))

            # Mixer (line 119-120).
            sigin = (1 - mix_it) * sigin + mix_it * siginnew
            sigout = (1 - mix_it) * sigout + mix_it * sigoutnew

        # --- Current evaluation (lines 124-136) ---
        I1 = 0.0
        I2 = 0.0
        I3 = 0.0
        Ico = 0.0
        Inco = 0.0
        for k in range(NE):
            I1 += np.real(np.trace(
                sigout2[:, :, k] @ n[:, :, k]
                - sigin2[:, :, k] @ p[:, :, k]
            ))
            I2 += np.real(np.trace(
                sigout1[:, :, k] @ n[:, :, k]
                - sigin1[:, :, k] @ p[:, :, k]
            ))
            I3 += np.real(np.trace(
                sigout[:, :, k] @ n[:, :, k]
                - sigin[:, :, k] @ p[:, :, k]
            ))
            Ico += np.real(np.trace(
                sigin2[:, :, k] @ G[:, :, k] @ gam1[:, :, k] @ G[:, :, k].conj().T
                - gam2[:, :, k] @ G[:, :, k] @ sigin1[:, :, k] @ G[:, :, k].conj().T
            ))
            Inco += np.real(np.trace(
                sigin2[:, :, k] @ G[:, :, k] @ gamp[:, :, k] @ G[:, :, k].conj().T
                - gam2[:, :, k] @ G[:, :, k] @ sigin[:, :, k] @ G[:, :, k].conj().T
            ))

        II[iV] = I1 * dE * IE
        II3[iV] = I3 * dE * IE
        II4[iV] = I2 * dE * IE
        IIco[iV] = Ico * dE * IE
        IInonco[iV] = Inco * dE * IE

        if verbose:
            print(f"  V = {V:5.3f} V  scba_iter = {outer:3d}  I = {II[iV]:+.3e} A")

    # --- Derivatives (lines 148-149) ---
    if NV >= 2:
        dIdV = np.diff(II) / dV
        V_dIdV = VV[1:]
    else:
        dIdV = np.zeros(0)
        V_dIdV = np.zeros(0)
    if NV >= 3:
        d2IdV2 = np.diff(dIdV) / dV
        V_d2IdV2 = VV[2:]
    else:
        d2IdV2 = np.zeros(0)
        V_d2IdV2 = np.zeros(0)

    params = dict(
        NS=NS, NC=NC, ND=ND, Np=Np, Nb=Nb, Vb=Vb, t0=t0,
        Ef=Ef, kT=kT, NV=NV, V_max=V_max, dE=dE, E_min=E_min, E_max=E_max,
        D0=D0, Dnu=list(Dnu_arr), hnu_units=list(hnu_units_arr),
        hnu_eV=list(hnu_eV), Nhnu=list(Nhnu),
        scba_tol=scba_tol, mix_it=mix_it, max_outer_iter=max_outer_iter,
    )
    return PatilResult(
        V=VV, I_total=II, I_coherent=IIco, I_incoherent=IInonco,
        dIdV=dIdV, d2IdV2=d2IdV2, V_dIdV=V_dIdV, V_d2IdV2=V_d2IdV2,
        params=params,
    )


# ---------------------------------------------------------------------------
# CLI — run with Patil's MATLAB defaults, save npz, emit plot.
# ---------------------------------------------------------------------------


def _save_and_plot(result: PatilResult, out_npz: str,
                   out_png: Optional[str]) -> None:
    np.savez(out_npz,
             V=result.V, I_total=result.I_total,
             I_coherent=result.I_coherent, I_incoherent=result.I_incoherent,
             dIdV=result.dIdV, d2IdV2=result.d2IdV2,
             V_dIdV=result.V_dIdV, V_d2IdV2=result.V_d2IdV2,
             **{k: np.asarray(v) if isinstance(v, (list, tuple)) else v
                for k, v in result.params.items()})
    print(f"wrote {out_npz}")

    if out_png is None:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print(f"matplotlib not available — skipping {out_png}")
        return

    fig, axes = plt.subplots(3, 1, figsize=(7, 9), sharex=True)
    axes[0].plot(result.V, result.I_total * 1e6, label="I_total")
    axes[0].plot(result.V, result.I_coherent * 1e6, "--", label="I_coh")
    axes[0].plot(result.V, result.I_incoherent * 1e6, ":", label="I_incoh")
    axes[0].set_ylabel("I [µA]")
    axes[0].legend()
    axes[0].grid(True)
    axes[0].set_title("Patil 1D reference (rtd2modes_1d.m port)")

    axes[1].plot(result.V_dIdV, result.dIdV * 1e6)
    axes[1].set_ylabel("dI/dV [µS]")
    axes[1].grid(True)

    axes[2].plot(result.V_d2IdV2, result.d2IdV2 * 1e6)
    axes[2].set_ylabel("d²I/dV² [µS/V]")
    axes[2].set_xlabel("V [V]")
    axes[2].grid(True)
    # Mark expected IETS peak positions at V = ℏω/e (± sign for symmetric bias).
    for hw in result.params["hnu_eV"]:
        axes[2].axvline(hw, color="k", alpha=0.3, lw=0.8, ls=":")

    fig.tight_layout()
    fig.savefig(out_png, dpi=120)
    print(f"wrote {out_png}")


def main():
    out_dir = os.path.dirname(os.path.abspath(__file__))
    out_npz = os.path.join(out_dir, "patil_reference_1d.npz")
    out_png = os.path.join(out_dir, "patil_reference_1d.png")
    result = run_patil_reference(verbose=True)
    _save_and_plot(result, out_npz, out_png)


if __name__ == "__main__":
    main()
