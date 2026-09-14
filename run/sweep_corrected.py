"""Corrected self-consistent Poisson-NEGF sweep (resumable, memory-lean).

Differences from the earlier ad-hoc runs, all of them deliberate corrections:

  * device  : ZnO_MgZnO_symmetric_spacer -- includes the 10 nm lightly-doped
              spacers of Akkala Table 3.1. Omitting them put ionized donors
              (1e25 m^-3) inside the barrier-reflection region, where quantum
              reflection holds n at ~2.4e24, producing spurious net POSITIVE
              charge ("emitter depletion") and a backward peak shift.
  * Ef      : 49.3 meV -- bulk charge neutrality n(Ef)=N_D for N_D=1e25 m^-3
              (N_c = 3.718e24 m^-3, m*=0.28, 300 K). The legacy 20 meV supports
              only half the specified doping.
  * donors  : used AS SPECIFIED (anchor=0). The old anchoring rescaled N_D down
              to whatever density the (under-resolved) model produced; with the
              correct Ef and an adequate grid that fudge is unnecessary.

NOTE on the energy grid: dE=2 meV under-resolves the DENSITY by ~25% (band-edge
1/sqrt(E) van Hove singularity, which the density weight ln(1+e^((mu-E)/kT))
maximises). The CURRENT is unaffected (0.3%) because its weight f_L-f_R vanishes
there. Use dE <= 0.5 meV when the density itself matters.

Usage:
  python -u run/sweep_corrected.py <device> <T_K> <mode> [mol] [npts] [vmin] [vmax] [Ef] [anchor] [dE]

  mode 0 = linear drop over the whole device (legacy "no Poisson" reference)
  mode 1 = self-consistent Poisson
  mode 2 = flat-band contacts (space-charge-free reference)
"""
import sys, os, gc

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
           'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_v, '2')

import numpy as np
from config.device_library import MATERIALS
from core.iets_analytic import analytic_d2idv2_inelastic_at_bias
from core.scba_rank1_keldysh import run_rank1_keldysh_single_bias
from core.poisson_negf import (
    build_electrostatics, compute_density_prefactor, contact_dirichlet,
    deep_contact_sites, flat_contact_profile, run_self_consistent_bias,
)
from core.poisson import physical_transverse_density
from run.run_rank1_sweep import (
    build_stack, emitter_barrier_center_nm, gaussian_chi, build_phonon_modes,
    linear_bias_profile, _M0, _EMIN_BIAS_MARGIN,
)

device = sys.argv[1]
T_K    = float(sys.argv[2])
mode   = int(sys.argv[3])
mol    = sys.argv[4] if len(sys.argv) > 4 else 'Baseline'
npts   = int(sys.argv[5]) if len(sys.argv) > 5 else 51
V_min  = float(sys.argv[6]) if len(sys.argv) > 6 else 0.0
V_max  = float(sys.argv[7]) if len(sys.argv) > 7 else 0.8
Ef     = float(sys.argv[8]) if len(sys.argv) > 8 else 0.0493
anchor = int(sys.argv[9]) if len(sys.argv) > 9 else 0
dE     = float(sys.argv[10]) if len(sys.argv) > 10 else 0.002

pois, flat = (mode == 1), (mode == 2)
tag = {0: 'nopois', 1: 'pois', 2: 'flatband'}[mode]
a_nm, eta, sigma_mol_nm = 0.2, 1e-12, 0.3
scba_max_iter, scba_mix, scba_tol = 10, 0.4, 1e-5
kT = 0.02585 * (T_K / 300.0)
a_m = a_nm * 1e-9

outdir = os.path.join(os.path.dirname(__file__), '..', 'results', '2026-08-05')
outdir = os.path.abspath(outdir)
os.makedirs(outdir, exist_ok=True)
_wt = f'_{int(V_min*1000)}-{int(V_max*1000)}mV'
_et = f'_Ef{int(round(Ef*1000))}'
_at = '' if anchor else '_NDspec'
_dt = '' if abs(dE - 0.002) < 1e-12 else f'_dE{dE*1e3:g}'
outf = f'{outdir}/corr_{device}_{mol}_{tag}{_wt}{_et}{_at}{_dt}.npz'

E_min = min(-0.25, -(0.5 * 0.8 + _EMIN_BIAS_MARGIN)) if pois else -0.25
H_z, UB, t0, Np, z_nm, bounds = build_stack(device, a_m)
NS, ND = 1, 1
chi_mol = gaussian_chi(z_nm, emitter_barrier_center_nm(bounds), sigma_mol_nm)
E_grid = np.arange(E_min, 0.5 + 0.5 * dE, dE)
D0_sq, hnu_idx, N_bose, chi_default, chi_list = build_phonon_modes(
    mol, dE, kT, chi_mol, NS, ND, UB=UB)
m_eff = MATERIALS["ZnO"]["m_eff"] * _M0

eps_r = N_D = cmask = density_prefactor = None
if pois or flat:
    eps_r, N_D, cmask = build_electrostatics(device, a_m)
print(f"[setup] {device} Np={Np} NE={E_grid.size} dE={dE*1e3:g} meV "
      f"Ef={Ef*1e3:.1f} meV mode={tag}", flush=True)

if pois:
    density_prefactor = compute_density_prefactor(
        E_grid=E_grid, H_z=H_z, UB=UB, t0=t0, Ef=Ef, kT=kT,
        chi_diag=chi_default, D0_sq_per_mode=D0_sq, hnu_idx_per_mode=hnu_idx,
        N_bose_per_mode=N_bose, chi_per_mode=chi_list, N_D=N_D,
        contact_mask=cmask, NS=NS, ND=ND, scba_max_iter=scba_max_iter,
        scba_mix=scba_mix, scba_tol=scba_tol, eta=eta)
    r0 = run_rank1_keldysh_single_bias(
        V=0.0, E_grid=E_grid, H_z=H_z, UB=UB,
        bias_profile=linear_bias_profile(0.0, Np, NS, ND), t0=t0, Ef=Ef, kT=kT,
        chi_diag=chi_default, D0_sq_per_mode=D0_sq, hnu_idx_per_mode=hnu_idx,
        N_bose_per_mode=N_bose, chi_per_mode=chi_list,
        max_iter=scba_max_iter, tol=scba_tol, mix=scba_mix, eta=eta)
    n0 = physical_transverse_density(r0.G_R, r0.Gam_L, r0.Gam_R, E_grid,
                                     Ef, Ef, kT, m_eff, a_m)
    nmod = float(np.mean(n0[deep_contact_sites(cmask)]))
    NDc = float(N_D[cmask].max())
    if anchor:
        N_D = N_D * (nmod / NDc)
        print(f"[setup] donors ANCHORED to model density {nmod:.3e}", flush=True)
    else:
        print(f"[setup] donors AS SPECIFIED N_D={NDc:.3e}; model equilibrium "
              f"contact n={nmod:.3e} (ratio {nmod/NDc:.3f})", flush=True)
    r0 = n0 = None
    gc.collect()

V_grid = np.linspace(V_min, V_max, npts)
if os.path.exists(outf):
    d = np.load(outf, allow_pickle=True)
    if d['V'].size == npts and np.allclose(d['V'], V_grid):
        I_R, d2I, U, n_e, done = (d['I_R'].copy(), d['d2I'].copy(), d['U'].copy(),
                                  d['n_e'].copy(), d['done'].copy())
        print(f"[resume] {int(done.sum())}/{npts}", flush=True)
    else:
        I_R = np.zeros(npts); d2I = np.zeros(npts)
        U = np.zeros((npts, Np)); n_e = np.zeros((npts, Np))
        done = np.zeros(npts, bool)
    d.close()
else:
    I_R = np.zeros(npts); d2I = np.zeros(npts)
    U = np.zeros((npts, Np)); n_e = np.zeros((npts, Np))
    done = np.zeros(npts, bool)

bar = UB > 0.5 * UB.max(); _i = np.where(bar)[0]
_g = np.where(np.diff(_i) > 1)[0]; b1 = _i[:_g[0] + 1]; b2 = _i[_g[0] + 1:]
well = np.arange(b1[-1] + 1, b2[0])
pre = np.arange(max(0, b1[0] - 40), b1[0])

U_warm = None
for i, V in enumerate(V_grid):
    if done[i]:
        if pois:
            U_warm = U[i].copy()
        continue
    if pois:
        sc = run_self_consistent_bias(
            V=float(V), E_grid=E_grid, H_z=H_z, UB=UB, t0=t0, Ef=Ef, kT=kT,
            chi_diag=chi_default, D0_sq_per_mode=D0_sq, hnu_idx_per_mode=hnu_idx,
            N_bose_per_mode=N_bose, chi_per_mode=chi_list, eps_r=eps_r, N_D=N_D,
            a_m=a_m, density_prefactor=density_prefactor, contact_mask=cmask,
            NS=NS, ND=ND, scba_max_iter=scba_max_iter, scba_mix=scba_mix,
            scba_tol=scba_tol, eta=eta, poisson_max_iter=60, poisson_tol=5e-3,
            kT_screen=0.002, bc_scheme='neumann', density_mode='physical',
            m_eff_kg=m_eff, U_init=U_warm)
        res = sc.result; U[i] = sc.U; n_e[i] = sc.n_e; U_warm = sc.U.copy()
        extra = (f"  pois={sc.poisson_iters}{'' if sc.poisson_converged else '!'}"
                 f"  pre/N_D={sc.n_e[pre].mean()/max(N_D[pre].mean(),1e-30):.3g}")
    else:
        if flat:
            _, dv = contact_dirichlet(cmask, float(V))
            bp = flat_contact_profile(float(V), cmask, dv)
        else:
            bp = linear_bias_profile(V, Np, NS, ND)
        U[i] = bp
        res = run_rank1_keldysh_single_bias(
            V=float(V), E_grid=E_grid, H_z=H_z, UB=UB, bias_profile=bp, t0=t0,
            Ef=Ef, kT=kT, chi_diag=chi_default, D0_sq_per_mode=D0_sq,
            hnu_idx_per_mode=hnu_idx, N_bose_per_mode=N_bose,
            chi_per_mode=chi_list, max_iter=scba_max_iter, tol=scba_tol,
            mix=scba_mix, eta=eta)
        extra = ""
    I_R[i] = res.I_right
    d2I[i] = analytic_d2idv2_inelastic_at_bias(res, kT=kT, E_F=Ef)
    done[i] = True
    np.savez(outf, V=V_grid, I_R=I_R, d2I=d2I, U=U, n_e=n_e, done=done,
             T_K=T_K, device=device, mode=mode, Ef=Ef, anchor=anchor, dE=dE)
    print(f"  [{i+1:3d}/{npts}] V={V:.3f} I_R={res.I_right:+.4e}{extra}", flush=True)
    res = None
    if pois:
        sc = None
    gc.collect()

print(f"[done] {int(done.sum())}/{npts} -> {outf}", flush=True)
