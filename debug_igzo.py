"""
Debug IGZO device setup - identify why we got zero current
"""

import numpy as np
import sys
sys.path.insert(0, '.')

from config.device_library import get_device, get_material, MATERIALS
from core.hamiltonian import discretize_device, build_hamiltonian

print("="*70)
print("DEBUGGING IGZO DEVICE")
print("="*70)

# Test 1: Check material properties
print("\n[1] Material Properties:")
print("-"*70)
for mat_name in ["IGZO", "Al2O3"]:
    mat = get_material(mat_name)
    print(f"  {mat_name}:")
    print(f"    m_eff = {mat['m_eff']:.3f} m₀")
    print(f"    Ec = {mat['Ec']:.3f} eV")

# Test 2: Get device and discretize
print("\n[2] Device Discretization:")
print("-"*70)
device = get_device("IGZO_Al2O3_symmetric")
print(f"Device: {device['description']}")
print(f"Layers:")
for i, layer in enumerate(device['layers']):
    print(f"  {i}: {layer['material']} - {layer['thickness']*1e9:.1f} nm")

grid = discretize_device(device, grid_spacing=0.3e-9)
print(f"\nGrid points: {grid['Np']}")
print(f"Grid spacing: {grid['a']*1e9:.2f} nm")

# Test 3: Check effective mass profile
print("\n[3] Effective Mass Profile:")
print("-"*70)
m0 = 9.10938356e-31
m_eff_m0 = grid['m_eff'] / m0
print(f"  Edge (site 0):  m_eff = {m_eff_m0[0]:.3f} m₀")
print(f"  Center:         m_eff = {m_eff_m0[grid['Np']//2]:.3f} m₀")
print(f"  Edge (site -1): m_eff = {m_eff_m0[-1]:.3f} m₀")
print(f"  Min m_eff: {m_eff_m0.min():.3f} m₀")
print(f"  Max m_eff: {m_eff_m0.max():.3f} m₀")

# Test 4: Build Hamiltonian and check hopping
print("\n[4] Hamiltonian and Hopping:")
print("-"*70)
H, t = build_hamiltonian(grid)
print(f"  Hopping t range: {t.min():.4f} to {t.max():.4f} eV")
print(f"  t[0] (edge): {t[0]:.4f} eV")
print(f"  Diagonal range: {np.diag(H).min():.4f} to {np.diag(H).max():.4f} eV")

# Test 5: Contact self-energies
print("\n[5] Contact Self-Energies:")
print("-"*70)
from core.self_energy import contact_self_energy_matrix, contact_self_energy_1d

E_test = 0.1
Sigma_L = contact_self_energy_matrix(E_test, grid, t, contact='left')
Sigma_R = contact_self_energy_matrix(E_test, grid, t, contact='right')

print(f"  At E = {E_test} eV:")
print(f"    Σ_L[0,0] = {Sigma_L[0,0]:.6f}")
print(f"    Σ_R[-1,-1] = {Sigma_R[-1,-1]:.6f}")

# Test at multiple energies
print(f"\n  Contact self-energy vs energy:")
for E in [0.0, 0.5, 1.0, 2.0, 3.0]:
    sigma = contact_self_energy_1d(E, grid['Ec'][0], t[0])
    print(f"    E = {E:.1f} eV: Σ = {sigma:.6f}")

# Test 6: Green's function
print("\n[6] Green's Function Test:")
print("-"*70)
from core.green_functions import retarded_greens_function, spectral_function, broadening_function

Sigma1 = contact_self_energy_matrix(E_test, grid, t, contact='left')
Sigma2 = contact_self_energy_matrix(E_test, grid, t, contact='right')

G = retarded_greens_function(E_test, H, Sigma1, Sigma2, SigmaS=None, eta=1e-4)
A = spectral_function(G)
Gamma1 = broadening_function(Sigma1)
Gamma2 = broadening_function(Sigma2)

print(f"  At E = {E_test} eV:")
print(f"    |G| max: {np.abs(G).max():.6e}")
print(f"    Trace(A): {np.trace(A):.6e}")
print(f"    Trace(Γ₁): {np.trace(Gamma1):.6e}")
print(f"    Trace(Γ₂): {np.trace(Gamma2):.6e}")

# Transmission
A2 = G @ Gamma2 @ G.conj().T
T = np.real(np.trace(Gamma1 @ A2))
print(f"    Transmission T: {T:.6e}")

# Test 7: Scan energies for transmission
print("\n[7] Transmission vs Energy:")
print("-"*70)
E_scan = np.linspace(-0.5, 4.0, 50)
T_scan = []

for E in E_scan:
    Sigma1 = contact_self_energy_matrix(E, grid, t, contact='left')
    Sigma2 = contact_self_energy_matrix(E, grid, t, contact='right')
    G = retarded_greens_function(E, H, Sigma1, Sigma2, SigmaS=None, eta=1e-4)
    Gamma1 = broadening_function(Sigma1)
    Gamma2 = broadening_function(Sigma2)
    A2 = G @ Gamma2 @ G.conj().T
    T = np.real(np.trace(Gamma1 @ A2))
    T_scan.append(T)

T_scan = np.array(T_scan)
print(f"  Energy range: {E_scan[0]:.2f} to {E_scan[-1]:.2f} eV")
print(f"  Transmission range: {T_scan.min():.6e} to {T_scan.max():.6e}")

# Find resonance
peak_idx = np.argmax(T_scan)
print(f"  Peak transmission at E = {E_scan[peak_idx]:.3f} eV, T = {T_scan[peak_idx]:.6e}")

print("\n" + "="*70)
print("DEBUG COMPLETE")
print("="*70)
