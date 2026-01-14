"""
Run IETS Simulation: Benzene with Quasi-3D Mode Summation

Tests the complete quasi-3D implementation with mode-dependent molecular coupling.
This is the first real molecular sensing simulation with correct spatial selectivity!

Configuration:
- Device: GaAs/AlAs symmetric RTD
- Molecule: Benzene (5 vibrational modes)
- Transverse modes: 3×3 = 9 modes
- Mode-dependent molecular coupling: ENABLED
"""

import numpy as np
import sys
import time
sys.path.insert(0, '.')

from run.run_single_molecule import SimulationConfig, run_iets_simulation

print("\n" + "=" * 70)
print("IETS SIMULATION: Benzene on GaAs/AlAs RTD (Quasi-3D)")
print("=" * 70)

# ============================================================================
# CONFIGURATION
# ============================================================================

config = SimulationConfig()

# Device
config.device_name = "GaAs_AlAs_symmetric"
config.grid_spacing = 0.3e-9  # Coarser grid for faster test (~70 points)

# Molecule
config.molecule_name = "Benzene"

# Phonons
config.bulk_phonon_energy = 0.036  # 36 meV (GaAs LO phonon)
config.bulk_phonon_coupling = 0.010  # 10 meV
config.molecular_coupling_scale = 1.0

# Energy grid (coarser for speed)
config.E_min = -0.2
config.E_max = 1.0
config.E_points = 50  # Reduced for speed

# Bias sweep (fewer points for initial test)
config.V_min = 0.0
config.V_max = 0.3  # Reduced from 0.5V
config.V_points = 7  # Just 7 points: 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3V

# Temperature
config.temperature = 300  # K

# SCBA parameters
config.scba_max_iter = 50  # Increased for better convergence
config.scba_tolerance = 1e-2  # Relaxed tolerance (0.01)
config.scba_mixing = 0.3

# Quasi-3D parameters
config.use_multimode = True   # ENABLE quasi-3D!
config.Ly = 1.0e-6           # 1 µm × 1 µm
config.Lz = 1.0e-6
config.n_max = 3             # 3×3 = 9 modes
config.m_max = 3

# Output
config.verbose = True

print("\nConfiguration:")
print(f"  Device: {config.device_name}")
print(f"  Molecule: {config.molecule_name}")
print(f"  Quasi-3D: {config.n_max}×{config.m_max} = {config.n_max*config.m_max} transverse modes")
print(f"  Energy grid: {config.E_points} points ({config.E_min:.1f} to {config.E_max:.1f} eV)")
print(f"  Bias sweep: {config.V_points} points ({config.V_min:.1f} to {config.V_max:.1f} V)")
print(f"  SCBA: max {config.scba_max_iter} iterations, tol={config.scba_tolerance:.0e}")

# ============================================================================
# RUN SIMULATION
# ============================================================================

print("\n" + "=" * 70)
print("RUNNING SIMULATION...")
print("=" * 70)

start_time = time.time()

try:
    results = run_iets_simulation(config)

    elapsed = time.time() - start_time

    print("\n" + "=" * 70)
    print("SIMULATION COMPLETE!")
    print("=" * 70)
    print(f"\nTotal time: {elapsed:.1f} seconds ({elapsed/60:.2f} minutes)")
    print(f"Time per bias point: {elapsed/config.V_points:.1f} seconds")

    # ========================================================================
    # ANALYZE RESULTS
    # ========================================================================

    print("\n" + "=" * 70)
    print("RESULTS ANALYSIS")
    print("=" * 70)

    V_array = results['V_array']
    I_array = results['I_array']
    dIdV = results['dIdV']
    d2IdV2 = results['d2IdV2']

    print(f"\nI-V Curve:")
    print(f"  Voltage range: {V_array[0]:.3f} to {V_array[-1]:.3f} V")
    print(f"  Current range: {I_array.min()*1e9:.3f} to {I_array.max()*1e9:.3f} nA")

    # Check equilibrium
    I_at_zero = I_array[0]
    print(f"\nEquilibrium check:")
    print(f"  I(V=0) = {I_at_zero*1e12:.6f} pA")
    if abs(I_at_zero) < 1e-15:
        print(f"  ✓ Perfect equilibrium (I = 0 at V = 0)")
    else:
        print(f"  ⚠ Non-zero current at equilibrium")

    # Show I-V data
    print(f"\n  V (V)     I (nA)    dI/dV (µS)   d²I/dV² (µS/V)")
    print(f"  " + "-" * 50)
    for i in range(len(V_array)):
        print(f"  {V_array[i]:.3f}    {I_array[i]*1e9:7.3f}   {dIdV[i]*1e6:8.3f}    {d2IdV2[i]*1e6:9.3f}")

    # IETS features
    print(f"\nIETS Spectrum (d²I/dV²):")
    print(f"  Range: {d2IdV2.min()*1e6:.3f} to {d2IdV2.max()*1e6:.3f} µS/V")

    # Find peaks (positive d²I/dV²)
    positive_features = d2IdV2 > 0
    n_peaks = np.sum(positive_features)
    print(f"  Positive features: {n_peaks}/{len(d2IdV2)} points")

    if n_peaks > 0:
        # Find voltage of maximum d²I/dV²
        peak_idx = np.argmax(d2IdV2)
        peak_V = V_array[peak_idx]
        peak_val = d2IdV2[peak_idx]
        print(f"  Maximum at: V = {peak_V:.3f} V, d²I/dV² = {peak_val*1e6:.3f} µS/V")

    # Convergence statistics
    print(f"\nSCBA Convergence:")
    scba_results = results['scba_results']
    converged_count = sum(1 for r in scba_results if r['converged'])
    print(f"  Converged bias points: {converged_count}/{len(scba_results)}")

    avg_iterations = np.mean([r['iterations'] for r in scba_results])
    print(f"  Average iterations: {avg_iterations:.1f}")

    # Check if any have mode information
    if 'transverse_modes' in scba_results[0]:
        modes = scba_results[0]['transverse_modes']
        print(f"\nTransverse Modes:")
        print(f"  Total: {modes.Nm}")
        print(f"  Energies: {modes.energies[0]*1000:.3f} to {modes.energies[-1]*1000:.3f} meV")

    # ========================================================================
    # SAVE RESULTS
    # ========================================================================

    print(f"\n" + "=" * 70)
    print("SAVING RESULTS")
    print("=" * 70)

    output_file = results.get('output_file', 'results/output.csv')
    print(f"\nResults saved to: {output_file}")

    print(f"\nResult dictionary keys:")
    for key in results.keys():
        if isinstance(results[key], np.ndarray):
            print(f"  {key}: array {results[key].shape}")
        else:
            print(f"  {key}: {type(results[key]).__name__}")

    # ========================================================================
    # PERFORMANCE ASSESSMENT
    # ========================================================================

    print(f"\n" + "=" * 70)
    print("PERFORMANCE ASSESSMENT")
    print("=" * 70)

    print(f"\nComputational Cost:")
    print(f"  Total time: {elapsed:.1f} seconds")
    print(f"  Time per bias point: {elapsed/config.V_points:.1f} seconds")
    print(f"  Estimated for 21 bias points: {elapsed/config.V_points*21:.1f} seconds ({elapsed/config.V_points*21/60:.1f} minutes)")

    if elapsed/config.V_points < 10:
        print(f"\n  ✓ Performance is GOOD (< 10 sec/point)")
        print(f"  Recommendation: Can run full simulations without hybrid optimization")
    elif elapsed/config.V_points < 30:
        print(f"\n  ⚠ Performance is ACCEPTABLE (10-30 sec/point)")
        print(f"  Recommendation: Hybrid mode selection would help but not critical")
    else:
        print(f"\n  ⚠ Performance is SLOW (> 30 sec/point)")
        print(f"  Recommendation: Implement hybrid mode selection (7-8x speedup)")

    # ========================================================================
    # PHYSICS CHECK
    # ========================================================================

    print(f"\n" + "=" * 70)
    print("PHYSICS VALIDATION")
    print("=" * 70)

    print(f"\nChecking physical correctness:")

    # 1. Equilibrium
    if abs(I_at_zero) < 1e-15:
        print(f"  ✓ Detailed balance: I = 0 at V = 0")
    else:
        print(f"  ✗ Equilibrium violated: I ≠ 0 at V = 0")

    # 2. Monotonic increase (approximately)
    if all(I_array[i] >= I_array[i-1] - 1e-12 for i in range(1, len(I_array))):
        print(f"  ✓ Current increases with bias")
    else:
        print(f"  ⚠ Current not monotonic (may have NDR)")

    # 3. Positive currents at forward bias
    if all(I_array[1:] >= 0):
        print(f"  ✓ Forward bias gives positive current")
    else:
        print(f"  ✗ Unexpected negative currents")

    # 4. Reasonable magnitude
    I_max_nA = I_array.max() * 1e9
    if 0.1 < I_max_nA < 1000:
        print(f"  ✓ Current magnitude reasonable ({I_max_nA:.1f} nA)")
    else:
        print(f"  ⚠ Current magnitude unusual ({I_max_nA:.1f} nA)")

    print(f"\n" + "=" * 70)
    print("✓ SIMULATION SUCCESSFUL!")
    print("=" * 70 + "\n")

except Exception as e:
    elapsed = time.time() - start_time
    print(f"\n" + "=" * 70)
    print("✗ SIMULATION FAILED")
    print("=" * 70)
    print(f"\nError after {elapsed:.1f} seconds:")
    print(f"  {type(e).__name__}: {e}")

    import traceback
    print(f"\nFull traceback:")
    traceback.print_exc()

    print(f"\n" + "=" * 70)
