"""
Test core physics modules with minimal parameters
"""
import sys
sys.path.append('.')

import numpy as np
from config.device_library import get_device
from core.hamiltonian import discretize_device, build_hamiltonian
from core.green_functions import retarded_greens_function, spectral_function
from core.self_energy import contact_self_energy_matrix

print("=" * 70)
print("TESTING CORE PHYSICS MODULES")
print("=" * 70)

# Test 1: Device discretization
print("\n[TEST 1] Device Discretization")
print("-" * 70)

try:
    device = get_device("In2O3_Al2O3_symmetric")
    dx = 0.2e-9  # 0.2 nm grid spacing (coarse for testing)

    grid = discretize_device(device, dx)
    x_grid = grid['x']
    Np = grid['Np']

    print(f"✓ Grid points: {Np}")
    print(f"  - Grid spacing: {dx*1e9:.2f} nm")
    print(f"  - Total length: {x_grid[-1]*1e9:.2f} nm")

except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback
    traceback.print_exc()

# Test 2: Hamiltonian construction
print("\n[TEST 2] Hamiltonian Construction")
print("-" * 70)

try:
    H, t = build_hamiltonian(grid)

    print(f"✓ Hamiltonian shape: {H.shape}")
    print(f"  - Is Hermitian: {np.allclose(H, H.conj().T)}")
    print(f"  - Hopping parameter: {t}")
    print(f"  - Diagonal range: [{np.diag(H).real.min():.4f}, {np.diag(H).real.max():.4f}] eV")

except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback
    traceback.print_exc()

# Test 3: Contact self-energy
print("\n[TEST 3] Contact Self-Energy")
print("-" * 70)

try:
    # Test at a few energies
    E_test = np.array([0.0, 0.1, 0.2])  # eV

    Sigma_L_list = []
    Sigma_R_list = []

    for E in E_test:
        Sigma_L = contact_self_energy_matrix(E, grid, t, 'left')
        Sigma_R = contact_self_energy_matrix(E, grid, t, 'right')
        Sigma_L_list.append(Sigma_L)
        Sigma_R_list.append(Sigma_R)

    print(f"✓ Contact self-energies computed for {len(E_test)} energies")
    print(f"  - Sigma_L shape: {Sigma_L.shape}")
    print(f"  - Sigma_R shape: {Sigma_R.shape}")
    print(f"  - Sigma_L[0,0] at E=0.1: {Sigma_L_list[1][0, 0]}")
    print(f"  - Has imaginary part: {np.abs(Sigma_L_list[1].imag).max() > 1e-10}")

except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback
    traceback.print_exc()

# Test 4: Green's function (without phonons, just contacts)
print("\n[TEST 4] Green's Function (Ballistic)")
print("-" * 70)

try:
    # Compute Green's function at one energy
    E_single = 0.1  # eV
    Sigma_L_single = contact_self_energy_matrix(E_single, grid, t, 'left')
    Sigma_R_single = contact_self_energy_matrix(E_single, grid, t, 'right')
    Sigma_S = None  # No scattering

    G = retarded_greens_function(E_single, H, Sigma_L_single, Sigma_R_single, Sigma_S)

    print(f"✓ Green's function computed")
    print(f"  - G shape: {G.shape}")
    print(f"  - G[0,0]: {G[0, 0]}")
    print(f"  - Max |G|: {np.abs(G).max():.4e}")
    print(f"  - Is finite: {np.all(np.isfinite(G))}")

    # Test spectral function
    A = spectral_function(G)
    trace_A = np.trace(A).real

    print(f"✓ Spectral function computed")
    print(f"  - Tr(A): {trace_A:.4f}")
    print(f"  - Is positive: {np.all(A.real >= -1e-10)}")
    print(f"  - Is finite: {np.all(np.isfinite(A))}")

except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback
    traceback.print_exc()

# Test 5: Check for numerical issues
print("\n[TEST 5] Numerical Stability Check")
print("-" * 70)

try:
    # Test Green's function at multiple energies to check for instabilities
    E_range = np.linspace(0.0, 0.3, 10)
    issues = []

    for E in E_range:
        Sigma_L_test = contact_self_energy_matrix(E, grid, t, 'left')
        Sigma_R_test = contact_self_energy_matrix(E, grid, t, 'right')

        G_test = retarded_greens_function(E, H, Sigma_L_test, Sigma_R_test, None)

        if not np.all(np.isfinite(G_test)):
            issues.append(f"Non-finite values at E={E:.3f} eV")
        if np.abs(G_test).max() > 1e10:
            issues.append(f"Very large values at E={E:.3f} eV: {np.abs(G_test).max():.2e}")

    if issues:
        print(f"⚠ Found {len(issues)} potential issues:")
        for issue in issues:
            print(f"  - {issue}")
    else:
        print(f"✓ No numerical stability issues found")
        print(f"  - Tested {len(E_range)} energy points")
        print(f"  - All Green's functions are finite and reasonable")

except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 70)
print("CORE PHYSICS TESTS COMPLETE")
print("=" * 70)
