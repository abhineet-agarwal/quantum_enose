"""
Debug Green's function sign issue
"""
import sys
sys.path.append('.')

import numpy as np
from config.device_library import get_device
from core.hamiltonian import discretize_device, build_hamiltonian
from core.green_functions import retarded_greens_function, spectral_function
from core.self_energy import contact_self_energy_matrix

print("=" * 70)
print("DEBUGGING GREEN'S FUNCTION SIGN")
print("=" * 70)

# Setup
device = get_device("In2O3_Al2O3_symmetric")
grid = discretize_device(device, 0.2e-9)
H, t = build_hamiltonian(grid)

# Compute Green's function
E = 0.1  # eV
Sigma_L = contact_self_energy_matrix(E, grid, t, 'left')
Sigma_R = contact_self_energy_matrix(E, grid, t, 'right')
G = retarded_greens_function(E, H, Sigma_L, Sigma_R, None)

print(f"\n[CHECK 1] Green's Function Properties")
print(f"  G[0,0] = {G[0, 0]}")
print(f"  Re[G[0,0]] = {G[0, 0].real:.6f}")
print(f"  Im[G[0,0]] = {G[0, 0].imag:.6f}")
print(f"  Tr[Re(G)] = {np.trace(G).real:.6f}")
print(f"  Tr[Im(G)] = {np.trace(G).imag:.6f}")

# Check if Im[G] is negative (as it should be for retarded GF)
im_G = G.imag
print(f"\n[CHECK 2] Imaginary Part of G")
print(f"  Mean Im[G] = {im_G.mean():.6f}")
print(f"  Min Im[G] = {im_G.min():.6f}")
print(f"  Max Im[G] = {im_G.max():.6f}")
print(f"  Are all diagonal Im[G] negative? {np.all(np.diag(G).imag < 0)}")

# Compute spectral function manually
print(f"\n[CHECK 3] Spectral Function Computation")
A_formula1 = 1j * (G - G.conj().T)
print(f"  Using A = i(G - G†):")
print(f"    Tr(A) = {np.trace(A_formula1).real:.6f}")

A_formula2 = -2 * G.imag
print(f"  Using A = -2 Im(G):")
print(f"    Tr(A) = {np.trace(A_formula2):.6f}")

# Check if G - G† is purely imaginary
diff = G - G.conj().T
print(f"\n[CHECK 4] G - G† Properties")
print(f"  Max |Re(G - G†)| = {np.abs(diff.real).max():.2e}")
print(f"  Max |Im(G - G†)| = {np.abs(diff.imag).max():.2e}")
print(f"  Is (G - G†) anti-Hermitian? {np.allclose(diff, -diff.conj().T)}")

# Check eigenvalues of A
A = spectral_function(G)
eigenvalues = np.linalg.eigvalsh(A)
print(f"\n[CHECK 5] Spectral Function A Properties")
print(f"  Tr(A) = {np.trace(A):.6f}")
print(f"  Min eigenvalue = {eigenvalues.min():.6f}")
print(f"  Max eigenvalue = {eigenvalues.max():.6f}")
print(f"  Negative eigenvalues: {np.sum(eigenvalues < -1e-10)}")
print(f"  Positive eigenvalues: {np.sum(eigenvalues > 1e-10)}")

print("\n" + "=" * 70)
if eigenvalues.min() < -1e-6:
    print("⚠ BUG CONFIRMED: Spectral function has negative eigenvalues!")
    print("   This violates the fundamental property A ≥ 0")
else:
    print("✓ Spectral function is positive semi-definite")
print("=" * 70)
