"""
Test: Ballistic Multimode Transport

Verifies quasi-3D mode summation without phonon scattering.
This is a critical validation test for the multimode implementation.

Tests:
1. Equilibrium current (I = 0 at V = 0) for multimode
2. Positive current for forward bias
3. Quasi-3D current > 1D current (more conduction channels)
4. Per-mode contributions sum to total
5. Mode energies affect transmission
"""

import numpy as np
import sys
sys.path.insert(0, '.')

from config.device_library import get_device
from core.hamiltonian import discretize_device, build_hamiltonian
from core.self_energy import contact_self_energy_matrix
from core.green_functions import broadening_function
from core.scba_solver import (
    scba_iteration, compute_current,
    scba_iteration_multimode, compute_current_multimode
)
from core.transverse_modes import TransverseModes

print("\n" + "=" * 70)
print("TEST: Ballistic Multimode Transport (Quasi-3D)")
print("=" * 70)

# ============================================================================
# SETUP: Device and Grid
# ============================================================================

print("\n[SETUP] Configuring GaAs/AlAs RTD...")
print("-" * 70)

device = get_device("GaAs_AlAs_symmetric")
grid = discretize_device(device, grid_spacing=0.15e-9)
H, t = build_hamiltonian(grid)
Np = grid['Np']

print(f"Device: {Np} grid points, {grid['x'][-1]*1e9:.1f} nm")

# ============================================================================
# SETUP: Transverse Modes
# ============================================================================

print("\n[SETUP] Creating transverse modes (3×3 = 9 modes)...")
print("-" * 70)

m0 = 9.10938356e-31  # kg
m_eff = 0.067 * m0   # GaAs

transverse_modes = TransverseModes(
    Ly=1.0e-6,   # 1 µm
    Lz=1.0e-6,   # 1 µm
    n_max=3,
    m_max=3,
    m_trans=m_eff
)
transverse_modes.compute_modes()

print(f"Transverse modes: {transverse_modes.Nm}")
print(f"Mode energies: {transverse_modes.energies[0]*1000:.3f} to {transverse_modes.energies[-1]*1000:.3f} meV")

# ============================================================================
# SETUP: Energy Grid and Contact Self-Energies
# ============================================================================

print("\n[SETUP] Energy grid and contacts...")
print("-" * 70)

# Coarse energy grid for fast testing
E_array = np.linspace(-0.2, 1.0, 50)  # 50 points
NE = len(E_array)

print(f"Energy grid: {NE} points from {E_array[0]:.2f} to {E_array[-1]:.2f} eV")

# Contact self-energies
def Sigma1(E):
    return contact_self_energy_matrix(E, grid, t, 'left')

def Sigma2(E):
    return contact_self_energy_matrix(E, grid, t, 'right')

# ============================================================================
# TEST 1: Equilibrium Current (V = 0)
# ============================================================================

print("\n" + "=" * 70)
print("TEST 1: Equilibrium Current (Multimode)")
print("=" * 70)

V = 0.0
mu1 = V / 2.0
mu2 = -V / 2.0
temperature = 300  # K

print(f"Bias: V = {V} V (mu1 = {mu1:.3f} eV, mu2 = {mu2:.3f} eV)")

# Run multimode SCBA (ballistic - no phonons)
result_3d = scba_iteration_multimode(
    E_array, H, Sigma1, Sigma2,
    mu1, mu2, temperature,
    phonon_modes=[],  # No phonon scattering
    grid=grid,
    transverse_modes=transverse_modes,
    max_iter=1,  # Should converge immediately for ballistic
    verbose=False
)

# Compute current
I_3d, I_vs_E_3d, I_vs_mode_3d = compute_current_multimode(
    result_3d, E_array, mu1, mu2, temperature
)

print(f"\nResult: I = {I_3d*1e12:.6f} pA")

if np.abs(I_3d) < 1e-15:  # Essentially zero
    print("✓ PASS: Current is zero at equilibrium")
else:
    print(f"✗ FAIL: Expected I = 0, got I = {I_3d*1e12:.3f} pA")

# ============================================================================
# TEST 2: Forward Bias Current (V = 0.1 V)
# ============================================================================

print("\n" + "=" * 70)
print("TEST 2: Forward Bias Current (Multimode)")
print("=" * 70)

V = 0.1
mu1 = V / 2.0
mu2 = -V / 2.0

print(f"Bias: V = {V} V (mu1 = {mu1:.3f} eV, mu2 = {mu2:.3f} eV)")

# Run multimode SCBA
result_3d = scba_iteration_multimode(
    E_array, H, Sigma1, Sigma2,
    mu1, mu2, temperature,
    phonon_modes=[],
    grid=grid,
    transverse_modes=transverse_modes,
    max_iter=1,
    verbose=False
)

# Compute current
I_3d, I_vs_E_3d, I_vs_mode_3d = compute_current_multimode(
    result_3d, E_array, mu1, mu2, temperature
)

print(f"\nResult: I = {I_3d*1e9:.3f} nA")

if I_3d > 0:
    print("✓ PASS: Current is positive for forward bias")
else:
    print(f"✗ FAIL: Expected I > 0, got I = {I_3d*1e9:.3f} nA")

# ============================================================================
# TEST 3: 1D vs Quasi-3D Comparison
# ============================================================================

print("\n" + "=" * 70)
print("TEST 3: 1D vs Quasi-3D Comparison")
print("=" * 70)

# Run 1D SCBA for comparison
result_1d = scba_iteration(
    E_array, H, Sigma1, Sigma2,
    mu1, mu2, temperature,
    phonon_modes=[],
    grid=grid,
    max_iter=1,
    verbose=False
)

# Pre-compute broadening arrays for 1D current
Gamma1_array = np.zeros((Np, Np, NE), dtype=float)
Gamma2_array = np.zeros((Np, Np, NE), dtype=float)

for iE, E in enumerate(E_array):
    Gamma1_array[:, :, iE] = broadening_function(Sigma1(E))
    Gamma2_array[:, :, iE] = broadening_function(Sigma2(E))

# Compute 1D current
I_1d, I_vs_E_1d = compute_current(
    result_1d, Gamma1_array, Gamma2_array, E_array,
    mu1, mu2, temperature
)

print(f"1D current:       I = {I_1d*1e9:.3f} nA")
print(f"Quasi-3D current: I = {I_3d*1e9:.3f} nA")
print(f"Ratio (3D/1D):    {I_3d/I_1d:.2f}x")

if I_3d > I_1d:
    print("✓ PASS: Quasi-3D current > 1D current (more channels)")
else:
    print(f"✗ FAIL: Expected quasi-3D > 1D, got {I_3d/I_1d:.2f}x")

# ============================================================================
# TEST 4: Per-Mode Contributions Sum to Total
# ============================================================================

print("\n" + "=" * 70)
print("TEST 4: Per-Mode Contribution Sum")
print("=" * 70)

# Check: sum of per-mode currents should equal total
I_sum_check = np.trapz(I_vs_E_3d, E_array) * (2 * 1.602176634e-19**2 / (2 * np.pi * 1.054571817e-34))

print(f"Total current:           I = {I_3d*1e9:.3f} nA")
print(f"Sum of I_vs_E:           I = {I_sum_check*1e9:.3f} nA")

# Check that I_vs_E sums correctly
if np.allclose(I_vs_E_3d, np.sum(I_vs_mode_3d, axis=1)):
    print("✓ PASS: Per-mode contributions sum correctly")
else:
    print("✗ FAIL: Per-mode sum does not match total")

# ============================================================================
# TEST 5: Mode Energy Analysis
# ============================================================================

print("\n" + "=" * 70)
print("TEST 5: Mode Energy Analysis")
print("=" * 70)

# Compute total current contribution from each mode
mode_currents = np.trapz(I_vs_mode_3d, E_array, axis=0) * (2 * 1.602176634e-19**2 / (2 * np.pi * 1.054571817e-34))

print("\nMode contributions:")
print("  Mode  (n,m)  Energy(meV)  Current(nA)  Fraction")
print("-" * 60)

for im in range(transverse_modes.Nm):
    n, m = transverse_modes.modes[im]
    E_mode = transverse_modes.energies[im] * 1000  # meV
    I_mode = mode_currents[im] * 1e9  # nA
    fraction = mode_currents[im] / I_3d if I_3d > 0 else 0

    print(f"  {im:4d}  ({n},{m})   {E_mode:6.3f}       {I_mode:6.3f}     {fraction:5.1%}")

# Check that lowest energy mode contributes most
dominant_mode = np.argmax(mode_currents)
if dominant_mode == 0:  # (1,1) mode should be dominant
    print("\n✓ PASS: Lowest energy mode (1,1) contributes most")
else:
    print(f"\n⚠ WARNING: Mode {dominant_mode} dominates instead of mode 0")

# ============================================================================
# SUMMARY
# ============================================================================

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print(f"\nKey Results:")
print(f"  • Equilibrium current: {I_3d*1e12:.3e} pA (should be ~0)")
print(f"  • Forward bias (V=0.1V): {I_3d*1e9:.2f} nA")
print(f"  • 1D current: {I_1d*1e9:.2f} nA")
print(f"  • Quasi-3D enhancement: {I_3d/I_1d:.2f}x")
print(f"  • Number of modes: {transverse_modes.Nm}")
print(f"  • Mode energy range: {transverse_modes.energies[0]*1000:.3f} to {transverse_modes.energies[-1]*1000:.3f} meV")

print("\n" + "=" * 70)
print("✓ All tests completed!")
print("=" * 70 + "\n")
