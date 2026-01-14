"""
Quick Test - Validate Complete Pipeline

Fast test to verify all components working together.
"""

import numpy as np
import sys
import os

# Setup path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

print("\n" + "="*70)
print("QUANTUM E-NOSE PIPELINE - QUICK TEST")
print("="*70)

# ============================================================================
# Test 1: Configuration modules
# ============================================================================

print("\n[TEST 1/5] Configuration modules...")

try:
    from config.device_library import get_device, DEVICES
    from config.molecular_database import get_molecule, MOLECULES
    
    device = get_device("In2O3_Al2O3_symmetric")
    molecule = get_molecule("Benzene")
    
    device_name = "In2O3_Al2O3_symmetric"
    molecule_name = "Benzene"
    print(f"  ✓ Loaded device: {device_name}")
    print(f"  ✓ Loaded molecule: {molecule_name} ({len(molecule['modes_meV'])} modes)")
    print(f"  ✓ Available devices: {len(DEVICES)}")
    print(f"  ✓ Available molecules: {len(MOLECULES)}")
except Exception as e:
    print(f"  ✗ FAILED: {e}")
    sys.exit(1)

# ============================================================================
# Test 2: Core physics modules
# ============================================================================

print("\n[TEST 2/5] Core physics modules...")

try:
    from core.hamiltonian import discretize_device, build_hamiltonian
    from core.green_functions import retarded_greens_function, spectral_function
    from core.self_energy import contact_self_energy_matrix
    
    grid = discretize_device(device, grid_spacing=0.15e-9)
    H, t = build_hamiltonian(grid)
    
    print(f"  ✓ Grid: {grid['Np']} points")
    print(f"  ✓ Hamiltonian: {H.shape}")
    
    # Test Green's function
    E_test = 0.5
    S1 = contact_self_energy_matrix(E_test, grid, t, 'left')
    S2 = contact_self_energy_matrix(E_test, grid, t, 'right')
    G = retarded_greens_function(E_test, H, S1, S2)
    A = spectral_function(G)
    
    print(f"  ✓ Green's function computed")
    print(f"  ✓ Spectral function: Tr(A) = {np.trace(A).real:.3f}")
    
except Exception as e:
    print(f"  ✗ FAILED: {e}")
    sys.exit(1)

# ============================================================================
# Test 3: SCBA solver
# ============================================================================

print("\n[TEST 3/5] SCBA solver...")

try:
    from core.scba_solver import scba_iteration, compute_current, compute_iets
    from run.run_single_molecule import build_phonon_modes, SimulationConfig
    
    # Simple test with minimal parameters
    config = SimulationConfig()
    config.E_points = 50
    config.V_points = 5
    config.scba_max_iter = 10
    config.scba_tolerance = 1e-3
    config.verbose = False
    
    # Build phonon modes
    phonon_modes = build_phonon_modes(molecule, config, grid, device, "Benzene")
    
    # Contact self-energy functions
    def Sigma1(E):
        return contact_self_energy_matrix(E, grid, t, 'left')
    def Sigma2(E):
        return contact_self_energy_matrix(E, grid, t, 'right')
    
    # Energy grid
    E_array = np.linspace(-0.3, 1.0, config.E_points)
    
    # Run SCBA for one bias point
    result = scba_iteration(
        E_array, H, Sigma1, Sigma2,
        mu1=0.1, mu2=-0.1, temperature=300,
        phonon_modes=phonon_modes, grid=grid,
        max_iter=config.scba_max_iter,
        tol=config.scba_tolerance,
        mix=0.4, verbose=False
    )
    
    conv_str = "✓" if result['converged'] else "~"
    print(f"  {conv_str} SCBA: {result['iterations']} iterations")
    print(f"  ✓ Phonon modes: {len(phonon_modes)}")
    
except Exception as e:
    print(f"  ✗ FAILED: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# ============================================================================
# Test 4: Run scripts
# ============================================================================

print("\n[TEST 4/5] Run scripts...")

try:
    from run.run_single_molecule import SimulationConfig as SimConfig
    from run.run_batch import BatchConfig, select_molecules
    
    sim_config = SimConfig()
    batch_config = BatchConfig()
    
    test_molecules = select_molecules("Aromatic")
    
    print(f"  ✓ Single simulation config created")
    print(f"  ✓ Batch config created")
    print(f"  ✓ Aromatic class: {len(test_molecules)} molecules")
    
except Exception as e:
    print(f"  ✗ FAILED: {e}")
    sys.exit(1)

# ============================================================================
# Test 5: Analysis tools
# ============================================================================

print("\n[TEST 5/5] Analysis tools...")

try:
    from analysis.iets_analysis import (
        extract_iets_fingerprint, compute_similarity
    )
    
    # Create synthetic spectrum
    V_test = np.linspace(0, 0.5, 100)
    d2IdV2_test = np.exp(-((V_test - 0.2)/0.03)**2)
    
    # Extract fingerprint
    fp = extract_iets_fingerprint(V_test, d2IdV2_test)
    
    # Compute similarity with itself (should be 1.0)
    sim, dist = compute_similarity(fp, fp, method='euclidean')
    
    print(f"  ✓ Fingerprint: {fp['n_peaks']} peaks")
    print(f"  ✓ Self-similarity: {sim:.3f} (expected ~1.0)")
    
except Exception as e:
    print(f"  ✗ FAILED: {e}")
    sys.exit(1)

# ============================================================================
# Summary
# ============================================================================

print("\n" + "="*70)
print("ALL TESTS PASSED! ✓✓✓")
print("="*70)
print("\nYour quantum e-nose simulation framework is READY!")
print("\nNext steps:")
print("  1. Run single molecule: python run/run_single_molecule.py")
print("  2. Run batch (test): python run/run_batch.py test")
print("  3. Run full batch: python run/run_batch.py all")
print("\nResults will be saved to ./results/ and ./batch_results/")
print("="*70 + "\n")
