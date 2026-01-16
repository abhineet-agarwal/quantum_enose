"""
Test Hybrid Mode Selection vs All-Inelastic SCBA

Compares:
1. All-inelastic (all 9 modes get full SCBA)
2. Hybrid (4 modes SCBA, 5 modes coherent)

Expected result: Hybrid gives ~2x speedup with similar physics.
"""

import numpy as np
import time
import sys
sys.path.insert(0, '.')

from config.device_library import get_device, get_material, MATERIALS
from config.molecular_database import get_molecule
from core.hamiltonian import discretize_device, build_hamiltonian
from core.self_energy import contact_self_energy_matrix, local_projection_operator, get_molecule_location
from core.transverse_modes import TransverseModes
from core.scba_solver import scba_iteration_multimode, compute_current_multimode
from core.scba_solver_hybrid import scba_iteration_hybrid, compute_current_hybrid
from core.hybrid_modes import select_hybrid_modes, print_mode_selection

print("="*70)
print("HYBRID MODE SELECTION TEST")
print("="*70)

# ============================================================================
# SETUP
# ============================================================================

# Device
device_name = "GaAs_AlAs_symmetric"
device = get_device(device_name)
grid = discretize_device(device, grid_spacing=0.3e-9)
H, t = build_hamiltonian(grid)

print(f"\nDevice: {device_name}")
print(f"Grid: {grid['Np']} points")

# Molecule
molecule_name = "Benzene"
molecule = get_molecule(molecule_name)
print(f"Molecule: {molecule_name}")

# Transverse modes
m0 = 9.10938356e-31
mat = get_material(grid['material'][0])
m_trans = mat['m_eff'] * m0

transverse_modes = TransverseModes(
    Ly=1.0e-6, Lz=1.0e-6,
    n_max=3, m_max=3,
    m_trans=m_trans
)
transverse_modes.compute_modes()

print(f"Transverse modes: {transverse_modes.Nm}")

# Phonon modes
local_sites, neighbor_radius = get_molecule_location(grid, device)
phonon_modes = []

# Bulk phonon
phonon_modes.append({
    'energy': 0.036,
    'coupling': 0.010,
    'is_local': False
})

# Molecular modes (just use 2 for speed)
mol_energies = molecule['modes_meV']
mol_couplings = molecule['coupling_meV']

for i in range(min(2, len(mol_energies))):
    energy_meV = mol_energies[i]
    coupling_meV = mol_couplings[i]

    # Mode-dependent coupling
    D_base = coupling_meV / 1000.0  # Convert meV to eV
    D_nm = np.zeros(transverse_modes.Nm)
    for im in range(transverse_modes.Nm):
        n, m = transverse_modes.modes[im]
        y_mol = 1.0e-6 / 2
        z_mol = 1.0e-6 / 2
        psi_nm = transverse_modes.get_wavefunction(n, m, y_mol, z_mol)
        psi_11 = transverse_modes.get_wavefunction(1, 1, y_mol, z_mol)
        D_nm[im] = D_base * (psi_nm**2 / psi_11**2)

    phonon_modes.append({
        'energy': energy_meV / 1000.0,  # Convert meV to eV
        'coupling': D_nm,
        'is_local': True,
        'local_sites': local_sites,
        'neighbor_radius': neighbor_radius
    })

print(f"Phonon modes: {len(phonon_modes)}")

# Energy grid (coarse for speed)
E_array = np.linspace(-0.1, 0.8, 30)
print(f"Energy grid: {len(E_array)} points")

# Contact self-energies
def Sigma1_func(E):
    return contact_self_energy_matrix(E, grid, t, contact='left')

def Sigma2_func(E):
    return contact_self_energy_matrix(E, grid, t, contact='right')

# Bias
V = 0.2
mu1 = +V/2
mu2 = -V/2
temperature = 300

print(f"Bias: {V} V, T = {temperature} K")

# ============================================================================
# TEST 1: ALL-INELASTIC MODE
# ============================================================================

print("\n" + "="*70)
print("TEST 1: ALL-INELASTIC (all 9 modes with SCBA)")
print("="*70)

start_all = time.time()

result_all = scba_iteration_multimode(
    E_array, H, Sigma1_func, Sigma2_func,
    mu1, mu2, temperature,
    phonon_modes, grid,
    transverse_modes,
    max_iter=30,
    tol=1e-2,
    mix=0.3,
    verbose=True
)

time_all = time.time() - start_all

I_all, I_vs_E_all, I_vs_mode_all = compute_current_multimode(
    result_all, E_array, mu1, mu2, temperature
)

print(f"\nAll-inelastic results:")
print(f"  Time: {time_all:.1f} seconds")
print(f"  Converged: {result_all['converged']}")
print(f"  Iterations: {result_all['iterations']}")
print(f"  Current: {I_all*1e9:.6f} nA")

# ============================================================================
# TEST 2: HYBRID MODE
# ============================================================================

print("\n" + "="*70)
print("TEST 2: HYBRID (6 modes SCBA, 3 modes coherent)")
print("="*70)

# Select modes - use 6 to get better accuracy
n_inelastic = 6
inelastic_modes, coherent_modes, importance_scores = select_hybrid_modes(
    transverse_modes, phonon_modes, temperature, n_inelastic=n_inelastic
)

print_mode_selection(transverse_modes, inelastic_modes, coherent_modes, importance_scores)

start_hybrid = time.time()

result_hybrid = scba_iteration_hybrid(
    E_array, H, Sigma1_func, Sigma2_func,
    mu1, mu2, temperature,
    phonon_modes, grid,
    transverse_modes,
    inelastic_modes, coherent_modes,
    max_iter=30,
    tol=1e-2,
    mix=0.3,
    verbose=True
)

time_hybrid = time.time() - start_hybrid

I_hybrid, I_vs_E_hybrid, I_vs_mode_hybrid = compute_current_hybrid(
    result_hybrid, E_array, mu1, mu2, temperature
)

print(f"\nHybrid results:")
print(f"  Time: {time_hybrid:.1f} seconds")
print(f"  Converged: {result_hybrid['converged']}")
print(f"  Iterations: {result_hybrid['iterations']}")
print(f"  Current: {I_hybrid*1e9:.6f} nA")

# ============================================================================
# COMPARISON
# ============================================================================

print("\n" + "="*70)
print("COMPARISON")
print("="*70)

speedup = time_all / time_hybrid if time_hybrid > 0 else 0
current_diff = abs(I_all - I_hybrid) / abs(I_all) * 100 if I_all != 0 else 0

print(f"\n  Metric               All-Inelastic    Hybrid         Difference")
print(f"  {'-'*65}")
print(f"  Time (s)             {time_all:<16.1f} {time_hybrid:<14.1f} {speedup:.2f}x faster")
print(f"  Current (nA)         {I_all*1e9:<16.6f} {I_hybrid*1e9:<14.6f} {current_diff:.1f}%")
print(f"  Converged            {str(result_all['converged']):<16} {str(result_hybrid['converged']):<14}")
print(f"  Iterations           {result_all['iterations']:<16} {result_hybrid['iterations']:<14}")

# Per-mode current comparison
print(f"\n  Per-mode current (pA) at V={V}V:")
print(f"  {'Mode':<8} {'All-Inel':<14} {'Hybrid':<14} {'Diff %':<10}")
print(f"  {'-'*50}")

for im in range(transverse_modes.Nm):
    I_mode_all = np.trapz(I_vs_mode_all[:, im], E_array) * 1e12
    I_mode_hybrid = np.trapz(I_vs_mode_hybrid[:, im], E_array) * 1e12
    if abs(I_mode_all) > 1e-15:
        diff = abs(I_mode_all - I_mode_hybrid) / abs(I_mode_all) * 100
    else:
        diff = 0

    mode_type = "SCBA" if im in inelastic_modes else "Coh"
    print(f"  {im:<3} ({mode_type:<4}) {I_mode_all:<14.6f} {I_mode_hybrid:<14.6f} {diff:<10.2f}")

print("\n" + "="*70)
print("SUMMARY")
print("="*70)

if speedup > 1.5:
    print(f"\n  Hybrid achieved {speedup:.1f}x speedup!")
else:
    print(f"\n  Speedup was only {speedup:.1f}x (expected ~{transverse_modes.Nm/n_inelastic:.1f}x)")

if current_diff < 10:
    print(f"  Current difference is {current_diff:.1f}% (acceptable)")
else:
    print(f"  WARNING: Current difference is {current_diff:.1f}% (may be too large)")

print("\n" + "="*70)
print("TEST COMPLETE")
print("="*70)
