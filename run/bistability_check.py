"""Is the self-consistent Poisson-NEGF solution unique at each bias?

An RTD can be genuinely bistable near NDR (charge stored in the well feeds back
on the potential; Goldman et al. 1987), so path dependence is not automatically
a numerical fault. This solves each bias three ways -- cold start (linear
profile), warm-started on an upward sweep, warm-started on a downward sweep --
at a tight Poisson tolerance, and compares the converged potentials and the
ballistic device current. Solutions that agree => unique; solutions that differ
while BOTH are converged => bistability; anything unconverged => no verdict.

Ballistic throughout (density on a fine grid via the banded route, current
from the banded transmission), so nothing (NE, Np, Np) is ever formed: phonons
move the density by 0.003 % (core.poisson.ballistic_transverse_density).

Usage:
  python -u run/bistability_check.py [cold|up|down|all] [V0 V1 dV]
"""
from __future__ import annotations

import os
import sys
import time

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np

from config.device_library import DEVICES, MATERIALS
from core.poisson import ballistic_transmission
from core.poisson_negf import build_electrostatics, run_self_consistent_bias
from core.scba_rank1_keldysh import landauer_current_1mode, tsu_esaki_current
from run.run_rank1_sweep import _M0, build_stack

DEVICE = "ZnO_MgZnO_symmetric_spacer"
T_K, EF = 300.0, 0.0493          # E_F neutral for N_D = 1e25 m^-3 at 300 K
A_M = 0.2e-9
DE_DENSITY = 0.0002              # density needs <= 0.2 meV (van Hove edge)
DE_CURRENT = 0.0002
# The predictor contracts at ~0.98/iteration: a 5 meV step left the potential
# 356 meV from the fixed point at 640 mV. 1e-5 eV is reached in ~400 iterations
# cold and holds to 0.3 meV when continued to 800.
POISSON_TOL, POISSON_MAX_ITER = 1e-5, 1000
OUT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results",
                                   "2026-10-02", "bistability_spacer.npz"))

kT = 0.02585 * T_K / 300.0
m_eff = MATERIALS["ZnO"]["m_eff"] * _M0
H_z, UB, t0, Np, z_nm, bounds = build_stack(DEVICE, A_M)
eps_r, N_D, cmask = build_electrostatics(DEVICE, A_M)
E_density = np.arange(-0.6, 0.5, DE_DENSITY)
E_current = np.arange(-0.6, 0.5, DE_CURRENT)
E_transport = np.arange(-0.6, 0.5 + 0.001, 0.002)   # unused without the final SCBA
area = float(np.prod(DEVICES[DEVICE]["transverse_size"]))
well = np.where(UB > 0.5 * UB.max())[0]
gaps = np.where(np.diff(well) > 1)[0]
well_sites = np.arange(well[gaps[0]] + 1, well[gaps[0] + 1])


def solve(V, U_init):
    t = time.time()
    sc = run_self_consistent_bias(
        V=V, E_grid=E_transport, H_z=H_z, UB=UB, t0=t0, Ef=EF, kT=kT,
        chi_diag=np.zeros(Np), D0_sq_per_mode=[], hnu_idx_per_mode=[],
        N_bose_per_mode=[], chi_per_mode=None, eps_r=eps_r, N_D=N_D, a_m=A_M,
        density_prefactor=None, contact_mask=cmask, NS=1, ND=1,
        scba_max_iter=1, scba_mix=0.3, scba_tol=1e-4, eta=1e-12,
        poisson_max_iter=POISSON_MAX_ITER, poisson_tol=POISSON_TOL,
        kT_screen=0.002, bc_scheme="neumann", density_mode="physical",
        m_eff_kg=m_eff, density_E_grid=E_density, U_init=U_init,
        final_scba=False)
    T = ballistic_transmission(E_current, H_z, UB, sc.U, t0)
    mu_L, mu_R = EF + V / 2, EF - V / 2
    I_mode = landauer_current_1mode(E_current, T, mu_L, mu_R, kT)
    I_dev = tsu_esaki_current(E_current, T, mu_L, mu_R, kT, m_eff, area)
    return dict(U=sc.U, n=sc.n_e, it=sc.poisson_iters, conv=sc.poisson_converged,
                dU=sc.dU_final, I_mode=I_mode, I_dev=I_dev, secs=time.time() - t,
                U_well=float(sc.U[well_sites].mean()))


def load():
    if os.path.exists(OUT):
        d = np.load(OUT, allow_pickle=True)
        return {k: d[k].item() for k in d.files}
    return {}


def save(db):
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    np.savez(OUT, **{k: np.array(v, dtype=object) for k, v in db.items()})


def run(path, V_list):
    db = load()
    prev = None
    for V in V_list:
        key = f"{path}_{V*1e3:.0f}"
        if key in db:
            prev = db[key]["U"]
            continue
        r = solve(float(V), None if path == "cold" else prev)
        db[key] = r
        save(db)
        prev = r["U"]
        print(f"  {path:4s} V={V*1e3:4.0f} mV  conv={r['conv']!s:5s} it={r['it']:3d} "
              f"dU={r['dU']*1e3:6.3f} meV  U_well={r['U_well']*1e3:+8.2f} meV  "
              f"I_mode={r['I_mode']*1e9:8.3f} nA  I_dev={r['I_dev']*1e3:8.2f} mA  "
              f"[{r['secs']:.0f}s]", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    V0, V1, dV = (float(x) for x in sys.argv[2:5]) if len(sys.argv) > 4 else (0.448, 0.768, 0.032)
    Vs = np.round(np.arange(V0, V1 + dV / 2, dV), 6)
    print(f"[setup] {DEVICE} Np={Np} Ef={EF*1e3:.1f} meV  density dE={DE_DENSITY*1e3:g} meV "
          f"({E_density.size} pts)  tol={POISSON_TOL*1e3:g} meV  V={V0}-{V1} step {dV}", flush=True)
    if mode in ("up", "all"):
        run("up", Vs)
    if mode in ("down", "all"):
        run("down", Vs[::-1])
    if mode in ("cold", "all"):
        run("cold", Vs)
