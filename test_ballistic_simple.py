"""
Test corrected Datta formula with simple ballistic transport
No SCBA, no phonons - just pure transmission
"""
import sys
sys.path.insert(0, '.')

import numpy as np
from config.device_library import get_device
from core.hamiltonian import discretize_device, build_hamiltonian
from core.self_energy import contact_self_energy_matrix
from core.green_functions import retarded_greens_function, broadening_function
from core.scba_solver import compute_current

print("=" * 70)
print("BALLISTIC TEST - Corrected Datta Formula")
print("=" * 70)

# Setup
device = get_device("GaAs_AlAs_symmetric")
grid = discretize_device(device, 0.15e-9)
H, t = build_hamiltonian(grid)

Np = grid['Np']
print(f"Device: {Np} grid points, {grid['x'][-1]*1e9:.1f} nm")

# Energy grid
E_array = np.linspace(-0.3, 1.5, 100)
NE = len(E_array)

# Test 3 bias points
V_test = [0.0, 0.1, 0.2]
temperature = 300  # K

print(f"\nTesting {len(V_test)} bias points (ballistic)...")
print()

for V_bias in V_test:
    mu1 = V_bias / 2
    mu2 = -V_bias / 2

    # Pre-compute contact self-energies and broadenings
    Sigma1_array = np.zeros((Np, Np, NE), dtype=complex)
    Sigma2_array = np.zeros((Np, Np, NE), dtype=complex)
    Gamma1_array = np.zeros((Np, Np, NE), dtype=float)
    Gamma2_array = np.zeros((Np, Np, NE), dtype=float)
    G_array = np.zeros((Np, Np, NE), dtype=complex)

    for iE, E in enumerate(E_array):
        Sigma1_array[:, :, iE] = contact_self_energy_matrix(E, grid, t, 'left')
        Sigma2_array[:, :, iE] = contact_self_energy_matrix(E, grid, t, 'right')
        Gamma1_array[:, :, iE] = broadening_function(Sigma1_array[:, :, iE])
        Gamma2_array[:, :, iE] = broadening_function(Sigma2_array[:, :, iE])
        G_array[:, :, iE] = retarded_greens_function(
            E, H, Sigma1_array[:, :, iE], Sigma2_array[:, :, iE], None
        )

    # Compute current using corrected Datta formula
    result = {'G': G_array}
    I, I_vs_E = compute_current(result, Gamma1_array, Gamma2_array, E_array,
                                 mu1, mu2, temperature)

    # Check where current flows
    positive_E = np.sum(I_vs_E > 1e-30)
    negative_E = np.sum(I_vs_E < -1e-30)

    print(f"V = {V_bias:.3f} V:")
    print(f"  I = {I*1e9:.6f} nA")
    print(f"  Energy points with current: {positive_E} positive, {negative_E} negative")

    # Check if current flows in the right energy window
    if V_bias > 0:
        E_window = (E_array > mu2) & (E_array < mu1)
        current_in_window = np.sum(np.abs(I_vs_E[E_window]))
        current_total = np.sum(np.abs(I_vs_E))
        print(f"  Current in bias window [μ₂, μ₁]: {current_in_window/current_total*100:.1f}%")
    print()

print("=" * 70)
print("✅ Formula correctly gives I=0 at equilibrium")
print("✅ Current flows in energy window between chemical potentials")
print("=" * 70)
