"""
Debug: Check transmission coefficient values
"""
import sys
sys.path.insert(0, '.')

import numpy as np
from config.device_library import get_device
from core.hamiltonian import discretize_device, build_hamiltonian
from core.self_energy import contact_self_energy_matrix
from core.green_functions import retarded_greens_function, broadening_function

print("=" * 70)
print("DEBUG: Transmission Coefficient")
print("=" * 70)

# Setup
device = get_device("GaAs_AlAs_symmetric")
grid = discretize_device(device, 0.15e-9)
H, t = build_hamiltonian(grid)
Np = grid['Np']

# Single energy point
E = 0.5  # eV (above barrier)

# Compute Green's function
Sigma1 = contact_self_energy_matrix(E, grid, t, 'left')
Sigma2 = contact_self_energy_matrix(E, grid, t, 'right')
G = retarded_greens_function(E, H, Sigma1, Sigma2, None)

# Broadenings
Gamma1 = broadening_function(Sigma1)
Gamma2 = broadening_function(Sigma2)

# Compute transmission using Datta formula: T = Tr[Γ₁A₂] where A₂ = GΓ₂G†
with np.errstate(all='ignore'):
    A2 = G @ Gamma2 @ G.conj().T
A2 = np.nan_to_num(A2, nan=0.0, posinf=0.0, neginf=0.0)

with np.errstate(all='ignore'):
    transmission = np.trace(Gamma1 @ A2)
transmission = np.nan_to_num(transmission, nan=0.0, posinf=0.0, neginf=0.0)
T = np.real(transmission)

print(f"Energy: E = {E} eV")
print(f"Transmission: T(E) = {T:.6e}")
print()

# Check matrix norms
print("Matrix diagnostics:")
print(f"  ||Γ₁||_max = {np.abs(Gamma1).max():.6e}")
print(f"  ||Γ₂||_max = {np.abs(Gamma2).max():.6e}")
print(f"  ||G||_max = {np.abs(G).max():.6e}")
print(f"  ||A₂||_max = {np.abs(A2).max():.6e}")
print()

# For a 1D ballistic conductor, T should be ~ 1
# For RTD with barriers, T should be << 1
print("Expected transmission:")
print("  - 1D ballistic: T ~ 1")
print("  - RTD with barriers: T << 1 (exponentially suppressed)")
print()

if T > 1e10:
    print(f"⚠️ Transmission is HUGE ({T:.2e}) - numerical overflow!")
elif T > 10:
    print(f"⚠️ Transmission is too large ({T:.2e}) - should be ≤ 1 for single mode")
elif T > 1:
    print(f"⚠️ Transmission > 1 ({T:.2e}) - unphysical!")
elif T < 1e-10:
    print(f"⚠️ Transmission is very small ({T:.2e}) - tunneling heavily suppressed")
else:
    print(f"✓ Transmission is reasonable: T = {T:.6e}")

print("=" * 70)
