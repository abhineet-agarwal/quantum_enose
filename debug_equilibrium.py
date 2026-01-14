"""
Debug: Check if current is zero at equilibrium
At V=0 (mu1 = mu2), detailed balance requires I = 0
"""
import sys
sys.path.insert(0, '.')

import numpy as np
from config.device_library import get_device
from core.hamiltonian import discretize_device, build_hamiltonian
from core.self_energy import contact_self_energy_matrix
from core.green_functions import retarded_greens_function, broadening_function

print("Equilibrium Test: V=0, mu1=mu2")
print("=" * 70)

# Setup
device = get_device("GaAs_AlAs_symmetric")
grid = discretize_device(device, 0.15e-9)
H, t = build_hamiltonian(grid)

# Single energy point test
E = 0.5  # eV
mu = 0.0  # equilibrium
kT = 0.0259
f = 1.0 / (1.0 + np.exp((E - mu) / kT))

# Green's function
Sigma1 = contact_self_energy_matrix(E, grid, t, 'left')
Sigma2 = contact_self_energy_matrix(E, grid, t, 'right')
G = retarded_greens_function(E, H, Sigma1, Sigma2, None)

# Broadenings
Gamma1 = broadening_function(Sigma1)
Gamma2 = broadening_function(Sigma2)

# Spectral function
A = 1j * (G - G.conj().T)
A = np.real(A)

# At equilibrium: f1 = f2 = f
Sigma_in = f * (Gamma1 + Gamma2)

# Correlation functions
with np.errstate(all='ignore'):
    n = G @ Sigma_in @ G.conj().T
n = np.nan_to_num(np.real(n), nan=0.0, posinf=0.0, neginf=0.0)
p = A - n

# Current integrand terms
term1 = np.trace(Gamma1 @ n)
term2 = np.trace(Gamma2 @ p)
integrand = np.real(term1 - term2)

print(f"Energy: E = {E} eV")
print(f"Fermi function: f = {f:.6f}")
print()
print(f"Tr[Γ₁ n] = {term1.real:.6e}")
print(f"Tr[Γ₂ p] = {term2.real:.6e}")
print(f"Integrand = Tr[Γ₁ n - Γ₂ p] = {integrand:.6e}")
print()

# Expected: at equilibrium, this should be ~0
if abs(integrand) < 1e-6:
    print("✓ Integrand ≈ 0 at equilibrium (correct!)")
else:
    print(f"✗ Integrand = {integrand:.6e} ≠ 0 (BUG!)")

# Additional check: detailed balance
# At equilibrium: Tr[Γ₁ n] should equal Tr[Γ₂ p]
print()
print("Detailed Balance Check:")
print(f"  Tr[Γ₁ n] = {term1.real:.6e}")
print(f"  Tr[Γ₂ p] = {term2.real:.6e}")
print(f"  Ratio: {(term1/term2).real:.6f}")

if np.isclose(term1, term2, rtol=0.01):
    print("  ✓ Terms balance (within 1%)")
else:
    print("  ✗ Terms don't balance - formula error!")

print("=" * 70)
