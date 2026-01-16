"""
Test ballistic RTD transport (no phonon scattering)
This is the simplest case - should work perfectly
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
print("BALLISTIC RTD TEST (No Phonons)")
print("=" * 70)

# Setup device
device = get_device("GaAs_AlAs_symmetric")
grid = discretize_device(device, 0.15e-9)
H, t = build_hamiltonian(grid)

print(f"Device: GaAs/AlAs RTD")
print(f"Grid points: {grid['Np']}")
print(f"Length: {grid['x'][-1]*1e9:.1f} nm")
print()

# Energy grid
E_array = np.linspace(-0.3, 1.5, 100)
NE = len(E_array)
Np = grid['Np']

# Test at 3 bias points
V_test = [0.0, 0.1, 0.2]

print(f"Testing {len(V_test)} bias points (ballistic transport)...")
print()

results_list = []

for V_bias in V_test:
    mu1 = V_bias / 2
    mu2 = -V_bias / 2
    kT = 0.0259  # 300K in eV

    # Pre-compute contact self-energies
    Sigma1_array = np.zeros((Np, Np, NE), dtype=complex)
    Sigma2_array = np.zeros((Np, Np, NE), dtype=complex)
    Gamma1_array = np.zeros((Np, Np, NE), dtype=complex)
    Gamma2_array = np.zeros((Np, Np, NE), dtype=complex)

    for iE, E in enumerate(E_array):
        Sigma1_array[:, :, iE] = contact_self_energy_matrix(E, grid, t, 'left')
        Sigma2_array[:, :, iE] = contact_self_energy_matrix(E, grid, t, 'right')
        Gamma1_array[:, :, iE] = broadening_function(Sigma1_array[:, :, iE])
        Gamma2_array[:, :, iE] = broadening_function(Sigma2_array[:, :, iE])

    # Compute Green's functions (ballistic - no scattering)
    G_array = np.zeros((Np, Np, NE), dtype=complex)
    for iE, E in enumerate(E_array):
        G_array[:, :, iE] = retarded_greens_function(
            E, H, Sigma1_array[:, :, iE], Sigma2_array[:, :, iE], None
        )

    # Compute spectral function
    A_array = np.zeros((Np, Np, NE))
    for iE in range(NE):
        A_array[:, :, iE] = 1j * (G_array[:, :, iE] - G_array[:, :, iE].conj().T)
        A_array[:, :, iE] = np.real(A_array[:, :, iE])

    # Fermi functions
    f1 = 1.0 / (1.0 + np.exp((E_array - mu1) / kT))
    f2 = 1.0 / (1.0 + np.exp((E_array - mu2) / kT))

    # Correlation functions (ballistic)
    n_matrix = np.zeros_like(G_array)
    p_matrix = np.zeros_like(G_array)
    for iE in range(NE):
        Sigma_in = f1[iE] * Gamma1_array[:, :, iE] + f2[iE] * Gamma2_array[:, :, iE]
        # Add safeguards
        with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
            n_temp = G_array[:, :, iE] @ Sigma_in @ G_array[:, :, iE].conj().T
        n_matrix[:, :, iE] = np.nan_to_num(np.real(n_temp), nan=0.0, posinf=0.0, neginf=0.0)
        p_matrix[:, :, iE] = A_array[:, :, iE] - n_matrix[:, :, iE]

    # Compute current using corrected Datta formula
    result = {'n': n_matrix, 'p': p_matrix, 'G': G_array}
    temperature_K = kT / 8.617333262145e-5  # Convert eV to Kelvin
    I, I_vs_E = compute_current(result, Gamma1_array, Gamma2_array, E_array,
                                 mu1, mu2, temperature_K)

    # Debug: check integrand sign
    positive_contrib = np.sum(I_vs_E > 0)
    negative_contrib = np.sum(I_vs_E < 0)

    print(f"V = {V_bias:.3f} V: I = {I*1e9:.6f} nA ({positive_contrib}+ / {negative_contrib}- energy points)")
    results_list.append({'V': V_bias, 'I': I})

print()
print("=" * 70)
print("RESULTS")
print("=" * 70)

for res in results_list:
    print(f"V = {res['V']:.3f} V → I = {res['I']*1e9:.3f} nA")

# Check sanity
I_values = [res['I'] for res in results_list]
all_positive = all(I >= 0 for I in I_values)
reasonable_mag = all(1e-15 < abs(I) < 1e-3 for I in I_values)
increases = all(I_values[i] <= I_values[i+1] for i in range(len(I_values)-1))

print()
print("Sanity checks:")
print(f"  {'✓' if all_positive else '✗'} All currents positive")
print(f"  {'✓' if reasonable_mag else '✗'} Magnitudes reasonable (fA-mA range)")
print(f"  {'✓' if increases else '✗'} Current increases with bias")

if all_positive and reasonable_mag and increases:
    print("\n🎉 BALLISTIC TRANSPORT WORKING!")
else:
    print("\n⚠ Issues detected - need more debugging")

print("=" * 70)
