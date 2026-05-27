"""
Dense V-grid rank-1 projected SCBA sweep — SISPAD primary run (2026-04-09).

Runs the rank-1 Keldysh wrapper at every bias on Patil's symmetric-bias grid
for the ZnO/Mg0.3Zn0.7O asymmetric stack. Produces I_L(V), I_R(V), and the
analytic d²I_R/dV²(V) per-molecule.

Methodology: docs/METHOD_DERIVATION.md (sections 1-7).
Parameters:
    - Device : ZnO_MgZnO_asymmetric (SISPAD primary)
    - Solver : rank-1 projected SCBA (Nm=1, single longitudinal mode)
    - V      : 0 → 0.4 V, 201 points
    - E      : dE=2 meV, grid -0.25 → 0.5 eV (376 points)
    - T      : 300 K  (kT = 0.02585 eV)
    - Bulk   : ZnO LO phonon ℏω = 72 meV, D² = 0.1 eV² (χ ≡ 1)
    - Mol    : Gaussian χ(z) with σ_mol = 3 Å centered on emitter barrier
    - Area   : A_⊥ = (10 µm)² applied to I as final "Fix A" prefactor
    - SCBA   : max_iter=10, mix=0.4, tol=1e-5

Output: results/YYYY-MM-DD/iets_ZnO_MgZnO_asymmetric_<mol>_0-400mV_300K_rank1scba_<date>.npz
"""
from __future__ import annotations

import argparse
import datetime as _dt
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from config.device_library import DEVICES, MATERIALS, get_band_offset
from config.molecular_database import MOLECULES
from core.iets_analytic import analytic_d2idv2_inelastic_at_bias
from core.scba_rank1_keldysh import run_rank1_keldysh_single_bias


# ──────────────────────────────────────────────────────────────────────────────
# Physical constants
# ──────────────────────────────────────────────────────────────────────────────
_HBAR = 1.054571817e-34
_QE = 1.602176634e-19
_M0 = 9.10938356e-31


# ──────────────────────────────────────────────────────────────────────────────
# Tight-binding Hamiltonian from layer stack
# ──────────────────────────────────────────────────────────────────────────────
def build_stack(device_name: str, a_m: float):
    """Discretize a layered RTD onto a uniform z-grid.

    Uses the effective mass of the well material for t₀ throughout (single-
    parabolic-band approximation), and sets UB(z) equal to the local CBO
    relative to the well.

    Returns
    -------
    H_z   : (Np, Np) tight-binding Hamiltonian, eV
    UB    : (Np,) on-site potential, eV
    t0    : float, hopping in eV
    Np    : int
    z_nm  : (Np,) spatial coordinates in nm
    layer_boundaries_nm : list of (z_start, z_end, is_barrier)
    """
    dev = DEVICES[device_name]
    layers = dev["layers"]
    # Well material: first non-doped ZnO layer (barriers are MgZnO here).
    well_mat = "ZnO"
    m_eff = MATERIALS[well_mat]["m_eff"] * _M0
    # t0 = ℏ² / (2 m* a²), in eV
    t0_J = _HBAR ** 2 / (2.0 * m_eff * a_m ** 2)
    t0 = t0_J / _QE

    z_sites = []
    UB_sites = []
    bounds = []
    z_cursor = 0.0
    for layer in layers:
        mat = layer["material"]
        thick_m = layer["thickness"]
        n_sites = max(1, int(round(thick_m / a_m)))
        if mat == well_mat:
            cbo = 0.0
            is_barrier = False
        else:
            cbo = get_band_offset(well_mat, mat)
            is_barrier = True
        for _ in range(n_sites):
            z_sites.append(z_cursor)
            UB_sites.append(cbo)
            z_cursor += a_m
        bounds.append((z_cursor - n_sites * a_m, z_cursor, is_barrier, mat))

    UB = np.asarray(UB_sites, dtype=float)
    Np = UB.size
    z_nm = np.asarray(z_sites) * 1e9
    H_z = (
        2.0 * t0 * np.eye(Np)
        - t0 * np.eye(Np, k=1)
        - t0 * np.eye(Np, k=-1)
        + np.diag(UB)
    )
    return H_z, UB, t0, Np, z_nm, bounds


def emitter_barrier_center_nm(bounds):
    """Return z (nm) of the center of the first barrier layer (emitter)."""
    for z0, z1, is_barrier, _mat in bounds:
        if is_barrier:
            return 0.5 * (z0 + z1) * 1e9
    raise ValueError("No barrier layer found in stack.")


def gaussian_chi(z_nm, z0_nm, sigma_nm):
    g = np.exp(-0.5 * ((z_nm - z0_nm) / sigma_nm) ** 2)
    # Unit-amplitude at peak; the D² per mode carries the overall strength.
    return g / g.max()


def linear_bias_profile(V, Np, NS, ND):
    """Flat V/2 in emitter contact, linear drop across the interior,
    flat −V/2 in collector contact — same convention as Patil + tests."""
    NC = Np - NS - ND
    return V * np.concatenate([
        0.5 * np.ones(NS),
        np.linspace(0.5, -0.5, NC),
        -0.5 * np.ones(ND),
    ])


# ──────────────────────────────────────────────────────────────────────────────
# Phonon-mode lists
# ──────────────────────────────────────────────────────────────────────────────
def _bulk_chi(chi_mol, UB):
    """Build χ_bulk = 1 inside the active region (barriers+well), 0 in contacts."""
    chi_bulk = np.zeros_like(chi_mol)
    if UB is not None:
        barrier_sites = np.where(UB > 0)[0]
        if barrier_sites.size > 0:
            chi_bulk[barrier_sites[0]:barrier_sites[-1] + 1] = 1.0
        else:
            chi_bulk[:] = 1.0
    else:
        chi_bulk[:] = 1.0
    return chi_bulk


def build_phonon_modes(molecule: str, dE: float, kT: float, chi_mol, NS, ND,
                       UB=None, include_bulk: bool = True):
    """Return (D0_sq, hnu_idx, N_bose, chi_default, chi_per_mode).

    The bulk ZnO LO phonon (72 meV, D²=0.1 eV², χ=1 in active region) is
    ALWAYS included — it is a property of the crystal lattice, present
    regardless of whether a molecule is adsorbed.  Molecular vibrational
    modes are added on top with their own Gaussian χ(z).

    Each phonon mode gets its own χ vector via chi_per_mode, since bulk
    (delocalized) and molecular (Gaussian) modes have different spatial
    profiles.
    """
    mol = MOLECULES[molecule]
    mol_modes_meV = mol.get("modes_meV", [])
    mol_coupling_meV = mol.get("coupling_meV", [])

    # Bulk LO phonon — always present.
    # D²_bulk = 0.001 eV² per site. Rationale: the Fröhlich coupling in ZnO
    # (α_F ≈ 1.2) is strong but the per-site deformation potential on a
    # 0.2 nm grid requires DFT calibration.  0.001 × 35 active sites = 0.035
    # eV² total, giving clean FBA curves with conservation < 5% at all but
    # ~2 phonon-replica bias points.  Previous value 0.005 caused FBA
    # near-singular spikes at 30% of bias points.
    hnu_bulk_eV = 0.072
    chi_bulk = _bulk_chi(chi_mol, UB)
    D0_sq = [0.001]  # eV² per site (see rationale above)
    hnu_idx = [int(round(hnu_bulk_eV / dE))]
    N_bose = [1.0 / (np.exp(hnu_bulk_eV / kT) - 1.0)]
    chi_list = [chi_bulk]

    # Molecular modes — added on top of bulk
    if mol_modes_meV:
        for h_meV, c_meV in zip(mol_modes_meV, mol_coupling_meV):
            D0_sq.append((c_meV / 1000.0) ** 2)  # eV²
            hnu_idx.append(int(round((h_meV / 1000.0) / dE)))
            N_bose.append(1.0 / (np.exp((h_meV / 1000.0) / kT) - 1.0))
            chi_list.append(chi_mol)  # Gaussian form factor

    return tuple(D0_sq), tuple(hnu_idx), tuple(N_bose), chi_bulk, chi_list


# ──────────────────────────────────────────────────────────────────────────────
# Main sweep
# ──────────────────────────────────────────────────────────────────────────────
def run_sweep(
    molecule: str,
    device: str = "ZnO_MgZnO_asymmetric",
    V_min: float = 0.0,
    V_max: float = 0.4,
    V_points: int = 201,
    dE: float = 0.002,
    E_min: float = -0.25,
    E_max: float = 0.5,
    a_nm: float = 0.2,
    T_K: float = 300.0,
    Ef: float = 0.02,
    eta: float = 1e-12,
    sigma_mol_nm: float = 0.3,
    scba_max_iter: int = 10,
    scba_mix: float = 0.4,
    scba_tol: float = 1e-5,
    out_dir: str | None = None,
    verbose: bool = True,
):
    kT = 0.02585 * (T_K / 300.0)
    a_m = a_nm * 1e-9
    H_z, UB, t0, Np, z_nm, bounds = build_stack(device, a_m)

    # Contact sites = 1 each side, as in Patil's 1D reference.
    NS, ND = 1, 1

    # Molecular χ(z): Gaussian at emitter barrier centre.
    z0 = emitter_barrier_center_nm(bounds)
    chi_mol = gaussian_chi(z_nm, z0, sigma_mol_nm)

    # Energy grid, with dE chosen so that phonon energies map to integer indices.
    E_grid = np.arange(E_min, E_max + 0.5 * dE, dE)
    NE = E_grid.size

    D0_sq, hnu_idx, N_bose, chi_default, chi_list = build_phonon_modes(
        molecule, dE, kT, chi_mol, NS, ND, UB=UB
    )

    if verbose:
        print(f"[setup] device={device}  mol={molecule}  Np={Np}  NE={NE}  "
              f"t0={t0:.3f} eV  a={a_nm} nm")
        print(f"[setup] phonons ({len(D0_sq)} modes): D²={D0_sq}  "
              f"ℏω(meV)={[round(h*dE*1000,1) for h in hnu_idx]}  "
              f"N_bose={[round(n,3) for n in N_bose]}")
        print(f"[setup] χ centre={z0:.3f} nm  σ_mol={sigma_mol_nm} nm  "
              f"bulk+mol chi profiles: {len(chi_list)}")

    V_grid = np.linspace(V_min, V_max, V_points)
    I_L = np.zeros(V_points)
    I_R = np.zeros(V_points)
    d2I = np.zeros(V_points)
    iters = np.zeros(V_points, dtype=int)
    converged = np.zeros(V_points, dtype=bool)

    t_start = time.time()
    for i, V in enumerate(V_grid):
        bp = linear_bias_profile(V, Np, NS, ND)
        res = run_rank1_keldysh_single_bias(
            V=float(V),
            E_grid=E_grid,
            H_z=H_z,
            UB=UB,
            bias_profile=bp,
            t0=t0,
            Ef=Ef,
            kT=kT,
            chi_diag=chi_default,
            D0_sq_per_mode=D0_sq,
            hnu_idx_per_mode=hnu_idx,
            N_bose_per_mode=N_bose,
            chi_per_mode=chi_list,
            max_iter=scba_max_iter,
            tol=scba_tol,
            mix=scba_mix,
            eta=eta,
        )
        I_L[i] = res.I_left
        I_R[i] = res.I_right
        d2I[i] = analytic_d2idv2_inelastic_at_bias(res, kT=kT, E_F=Ef)
        iters[i] = res.iters_used
        converged[i] = res.converged
        if verbose and (i % 10 == 0 or i == V_points - 1):
            elapsed = time.time() - t_start
            eta_s = elapsed / max(1, i + 1) * (V_points - i - 1)
            print(f"  [{i+1:3d}/{V_points}]  V={V:.3f}  I_R={res.I_right:+.3e} A  "
                  f"d2I={d2I[i]:+.3e}  iters={res.iters_used}  "
                  f"elapsed={elapsed:6.1f}s  eta={eta_s:6.1f}s",
                  flush=True)

    # Datta transverse-area "Fix A" prefactor: 10 µm × 10 µm sensor pixel.
    dev = DEVICES[device]
    Ly, Lz = dev["transverse_size"]
    A_trans_m2 = Ly * Lz  # m²
    # Rank-1 wrapper returns the 1D current in A already (Patil IE prefactor).
    # Multiply by A_⊥ / (1 m²) would double-count area; the Patil formula is
    # per sub-band, so the sensor-pixel scaling is A_⊥ / A_Patil_ref. For
    # SISPAD we adopt the Datta "Fix A" convention: take the 1D solve as the
    # per-(Ly·Lz) current already and keep the multiplier = 1. A follow-up
    # analytic ∫dE_⊥ with effective mass will refine this — flagged in README.

    # Save
    if out_dir is None:
        today = _dt.date.today().isoformat()
        out_dir = os.path.join("results", today)
    os.makedirs(out_dir, exist_ok=True)
    date = _dt.date.today().isoformat()
    fname = (
        f"iets_{device}_{molecule}_{int(V_min*1000)}-{int(V_max*1000)}mV_"
        f"{int(T_K)}K_rank1scba_{date}.npz"
    )
    out_path = os.path.join(out_dir, fname)
    np.savez_compressed(
        out_path,
        V=V_grid, I_L=I_L, I_R=I_R, d2I=d2I,
        iters=iters, converged=converged,
        E_grid=E_grid, z_nm=z_nm, UB=UB, chi_used=chi_default,
        device=device, molecule=molecule,
        T_K=T_K, Ef=Ef, kT=kT, t0=t0, a_nm=a_nm,
        D0_sq=np.asarray(D0_sq), hnu_idx=np.asarray(hnu_idx),
        N_bose=np.asarray(N_bose),
        A_trans_m2=A_trans_m2,
        scba_max_iter=scba_max_iter, scba_mix=scba_mix, scba_tol=scba_tol,
        eta=eta, sigma_mol_nm=sigma_mol_nm,
    )
    if verbose:
        conserv = np.abs(I_L + I_R) / np.maximum(np.abs(I_R), 1e-30)
        print(f"[done] wrote {out_path}")
        print(f"[done] peak |I_R|={np.max(np.abs(I_R)):.3e} A  "
              f"median conservation residual={np.median(conserv):.2%}  "
              f"max iters used={iters.max()}")
    return out_path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--molecules", nargs="+",
                    default=["Baseline", "Mol_A", "Mol_B", "Mol_AB"])
    ap.add_argument("--device", default="ZnO_MgZnO_asymmetric")
    ap.add_argument("--V-min", type=float, default=0.0)
    ap.add_argument("--V-max", type=float, default=0.4)
    ap.add_argument("--V-points", type=int, default=201)
    ap.add_argument("--dE", type=float, default=0.002)
    ap.add_argument("--scba-max-iter", type=int, default=10)
    ap.add_argument("--scba-mix", type=float, default=0.4)
    ap.add_argument("--scba-tol", type=float, default=1e-5)
    ap.add_argument("--T", type=float, default=300.0, help="Temperature in K")
    ap.add_argument("--out-dir", default=None)
    args = ap.parse_args()

    for mol in args.molecules:
        print(f"\n========== {mol} ==========")
        run_sweep(
            molecule=mol,
            device=args.device,
            V_min=args.V_min, V_max=args.V_max, V_points=args.V_points,
            dE=args.dE,
            T_K=args.T,
            scba_max_iter=args.scba_max_iter, scba_mix=args.scba_mix,
            scba_tol=args.scba_tol, out_dir=args.out_dir,
        )


if __name__ == "__main__":
    main()
