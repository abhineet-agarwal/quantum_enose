"""The chunked low-memory SCBA must reproduce the fast solver exactly."""
import numpy as np

from core.scba_rank1_keldysh import (run_rank1_keldysh_single_bias_fast,
                                     run_rank1_keldysh_single_bias_lowmem)
from run.run_rank1_sweep import (build_phonon_modes, build_stack, emitter_barrier_center_nm,
                                 gaussian_chi, linear_bias_profile)


def test_lowmem_matches_fast():
    a_m, dE, kT = 0.5e-9, 0.004, 0.02585
    H, UB, t0, Np, z, b = build_stack("ZnO_MgZnO_symmetric", a_m)
    chi = gaussian_chi(z, emitter_barrier_center_nm(b), 0.3)
    D0, hi, Nb, chd, chl = build_phonon_modes("Mol_A", dE, kT, chi, 1, 1, UB=UB)
    E = np.arange(-0.1, 0.5, dE)
    kw = dict(V=0.4, E_grid=E, H_z=H, UB=UB, bias_profile=linear_bias_profile(0.4, Np, 1, 1),
              t0=t0, Ef=0.02, kT=kT, chi_diag=chd, D0_sq_per_mode=D0, hnu_idx_per_mode=hi,
              N_bose_per_mode=Nb, chi_per_mode=chl, max_iter=60, tol=1e-6, mix=0.3)
    a = run_rank1_keldysh_single_bias_fast(**kw)
    lo = run_rank1_keldysh_single_bias_lowmem(chunk=17, **kw)   # chunk not dividing NE
    assert lo.iters_used == a.iters_used and lo.converged == a.converged
    n_fast = np.real(np.einsum("eii->ei", a.G_lesser))
    assert np.abs(lo.n_diag - n_fast).max() <= 1e-12 * np.abs(n_fast).max()
    assert abs(lo.I_right / a.I_right - 1) < 1e-12
    assert abs(lo.I_left / a.I_left - 1) < 1e-12
