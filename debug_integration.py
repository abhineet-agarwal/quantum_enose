"""
Debug: Check current integration step-by-step
"""
import sys
sys.path.insert(0, '.')

import numpy as np
from config.device_library import get_device
from core.hamiltonian import discretize_device, build_hamiltonian
from core.self_energy import contact_self_energy_matrix
from core.green_functions import retarded_greens_function, broadening_function

print("=" * 70)
print("DEBUG: Current Integration")
print("=" * 70)

# Setup
device = get_device("GaAs_AlAs_symmetric")
grid = discretize_device(device, 0.15e-9)
H, t = build_hamiltonian(grid)
Np = grid['Np']

# Coarse energy grid for debugging
E_array = np.linspace(0.0, 1.0, 10)
NE = len(E_array)
dE = E_array[1] - E_array[0]

# Bias
V_bias = 0.1
mu1 = V_bias / 2
mu2 = -V_bias / 2
kT = 0.0259  # 300K in eV

print(f"Device: {Np} points")
print(f"Energy grid: {NE} points from {E_array[0]:.2f} to {E_array[-1]:.2f} eV (dE={dE:.3f} eV)")
print(f"Bias: V = {V_bias} V (μ₁={mu1:.3f} eV, μ₂={mu2:.3f} eV)")
print()

# Compute for each energy
print("Energy    T(E)        f₁-f₂      Integrand")
print("-" * 50)

integrand_values = []
for iE, E in enumerate(E_array):
    # Compute transmission
    Sigma1 = contact_self_energy_matrix(E, grid, t, 'left')
    Sigma2 = contact_self_energy_matrix(E, grid, t, 'right')
    G = retarded_greens_function(E, H, Sigma1, Sigma2, None)
    Gamma1 = broadening_function(Sigma1)
    Gamma2 = broadening_function(Sigma2)

    with np.errstate(all='ignore'):
        A2 = G @ Gamma2 @ G.conj().T
    A2 = np.nan_to_num(A2, nan=0.0, posinf=0.0, neginf=0.0)

    with np.errstate(all='ignore'):
        T = np.real(np.trace(Gamma1 @ A2))
    T = np.nan_to_num(T, nan=0.0, posinf=0.0, neginf=0.0)

    # Fermi functions
    f1 = 1.0 / (1.0 + np.exp((E - mu1) / kT))
    f2 = 1.0 / (1.0 + np.exp((E - mu2) / kT))

    # Integrand
    integrand = T * (f1 - f2)
    integrand_values.append(integrand)

    print(f"{E:6.3f}    {T:.3e}   {f1-f2:+.3e}   {integrand:.3e}")

print()

# Integrate
I_integral = np.trapz(integrand_values, E_array)  # Units: dimensionless × eV = eV
print(f"∫ T(E)(f₁-f₂) dE = {I_integral:.6e} eV")
print()

# Convert to current
q = 1.602176634e-19  # C
h = 6.62607015e-34  # J·s
hbar = 1.054571817e-34  # J·s

I_Landauer = (2 * q**2 / h) * I_integral  # Landauer formula (with spin factor 2)
I_Datta_h = (q / h) * I_integral  # What I'm using
I_Datta_hbar = (q**2 / hbar) * I_integral  # What Datta paper might mean

print("Current conversion:")
print(f"  Landauer (2q²/h): I = {I_Landauer:.3e} A = {I_Landauer*1e9:.3f} nA")
print(f"  Using (q/h):      I = {I_Datta_h:.3e} A = {I_Datta_h*1e9:.3f} nA")
print(f"  Using (q²/ℏ):     I = {I_Datta_hbar:.3e} A = {I_Datta_hbar*1e9:.3f} nA")
print()

# Expected ballistic current
G_ballistic = 2 * q**2 / h  # Conductance quantum
I_expected = G_ballistic * V_bias * 4e-4  # × transmission ~ 4e-4
print(f"Expected for T~4e-4: I ~ {I_expected:.3e} A = {I_expected*1e9:.3f} nA")

print("=" * 70)
