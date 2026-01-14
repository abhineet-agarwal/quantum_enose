"""
Minimal single molecule simulation with reduced parameters for testing
"""
import sys
sys.path.append('.')

import numpy as np
import warnings

# Show all warnings (don't convert to errors)
warnings.filterwarnings('always')

from config.molecular_database import get_molecule
from config.device_library import get_device
from core.hamiltonian import discretize_device, build_hamiltonian
from core.self_energy import contact_self_energy_matrix
from core.scba_solver import scba_iteration, compute_current
from run.run_single_molecule import build_phonon_modes, SimulationConfig

print("=" * 70)
print("MINIMAL SINGLE MOLECULE SIMULATION")
print("=" * 70)

# Setup with MINIMAL parameters
print("\n[SETUP] Creating minimal configuration")
device = get_device("In2O3_Al2O3_symmetric")
molecule = get_molecule("Benzene")

# Very coarse grid
grid_spacing = 0.3e-9  # 0.3 nm (coarse!)
grid = discretize_device(device, grid_spacing)
H, t = build_hamiltonian(grid)

print(f"  Device: In2O3_Al2O3_symmetric")
print(f"  Molecule: Benzene ({len(molecule['modes_meV'])} modes)")
print(f"  Grid points: {grid['Np']}")
print(f"  Grid spacing: {grid_spacing*1e9:.2f} nm")

# Minimal simulation config
config = SimulationConfig()
config.E_points = 30  # Very few energy points
config.V_points = 3   # Only 3 bias points
config.scba_max_iter = 5  # Few iterations
config.scba_tolerance = 1e-2  # Loose tolerance
config.verbose = True

# Build phonon modes
phonon_modes = build_phonon_modes(molecule, config, grid, device, "Benzene")
print(f"  Phonon modes: {len(phonon_modes)}")

# Contact self-energy functions
def Sigma1(E):
    return contact_self_energy_matrix(E, grid, t, 'left')

def Sigma2(E):
    return contact_self_energy_matrix(E, grid, t, 'right')

# Test at one bias point
print("\n[TEST] Running SCBA at V = 0.2 V")
print("-" * 70)

try:
    # Energy array
    E_min = -0.2
    E_max = 0.6
    E_array = np.linspace(E_min, E_max, config.E_points)

    # Bias point
    V_bias = 0.2
    mu1 = V_bias / 2
    mu2 = -V_bias / 2

    print(f"  Energy range: [{E_min:.2f}, {E_max:.2f}] eV ({config.E_points} points)")
    print(f"  Bias: V = {V_bias:.2f} V")
    print(f"  Chemical potentials: μ1 = {mu1:.2f}, μ2 = {mu2:.2f} eV")
    print(f"  Temperature: 300 K")
    print()

    # Run SCBA iteration
    result = scba_iteration(
        E_array, H, Sigma1, Sigma2,
        mu1=mu1, mu2=mu2, temperature=300,
        phonon_modes=phonon_modes, grid=grid,
        max_iter=config.scba_max_iter,
        tol=config.scba_tolerance,
        mix=0.3,
        verbose=True
    )

    print()
    conv_status = "✓ Converged" if result['converged'] else "⚠ Not fully converged"
    print(f"{conv_status} after {result['iterations']} iterations")
    print(f"  Final residual: {result['residual']:.2e}")
    print(f"  Tolerance: {config.scba_tolerance:.2e}")

    # Check for NaN/Inf in results
    G = result['G']
    has_nan = np.any(~np.isfinite(G))
    has_large = np.abs(G).max() > 1e10

    if has_nan:
        print(f"\n⚠ WARNING: Green's function contains NaN or Inf values")
    elif has_large:
        print(f"\n⚠ WARNING: Green's function has very large values: {np.abs(G).max():.2e}")
    else:
        print(f"\n✓ Green's function is numerically stable")
        print(f"  Max |G|: {np.abs(G).max():.2e}")

    print("\n" + "=" * 70)
    if not has_nan:
        print("SUCCESS: Minimal simulation completed!")
    else:
        print("PARTIAL SUCCESS: Simulation ran but has numerical issues")
    print("=" * 70)

except RuntimeWarning as e:
    print(f"\n✗ RUNTIME WARNING: {e}")
    import traceback
    traceback.print_exc()

except Exception as e:
    print(f"\n✗ FAILED: {e}")
    import traceback
    traceback.print_exc()
