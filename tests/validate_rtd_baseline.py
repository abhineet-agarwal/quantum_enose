"""
Validate RTD Baseline - GaAs/AlAs Double Barrier

This script validates the quantum transport simulation against known RTD physics.
Expected behavior:
- Resonant tunneling peak at V ~ 0.1-0.3 V
- Current on the order of nA-μA
- Possible negative differential resistance (NDR)
- Smooth I-V characteristic
"""
import sys
sys.path.insert(0, '.')

import numpy as np
from run.run_single_molecule import SimulationConfig, run_iets_simulation, save_results

print("=" * 70)
print("RTD BASELINE VALIDATION - GaAs/AlAs Double Barrier")
print("=" * 70)
print()
print("Expected RTD Physics:")
print("  ✓ Resonant tunneling peak at V ~ 0.1-0.3 V")
print("  ✓ Peak current ~ nA to μA range")
print("  ✓ Smooth I-V curve")
print("  ✓ Possible negative differential resistance (NDR)")
print()
print("=" * 70)
print()

# Create config for baseline validation
config = SimulationConfig()

# Use validated GaAs/AlAs device
config.device_name = "GaAs_AlAs_symmetric"
config.molecule_name = "Baseline"  # NO molecule - just device physics

# Reasonable parameters (not too slow, not too coarse)
config.grid_spacing = 0.15e-9  # 0.15 nm
config.E_points = 100  # Sufficient resolution
config.V_points = 21   # 0.025 V steps from 0 to 0.5 V
config.V_min = 0.0
config.V_max = 0.5

# GaAs bulk phonon (36 meV - well known)
config.bulk_phonon_energy = 0.036  # 36 meV (GaAs LO phonon)
config.bulk_phonon_coupling = 0.010  # Weak coupling

# SCBA settings
config.scba_max_iter = 20
config.scba_tolerance = 1e-3
config.verbose = True

print(f"Device: {config.device_name}")
print(f"Molecule: {config.molecule_name} (no molecule - pure device)")
print(f"Grid spacing: {config.grid_spacing*1e9:.2f} nm")
print(f"Energy points: {config.E_points}")
print(f"Bias points: {config.V_points} (0 to {config.V_max} V)")
print(f"Bulk phonon: {config.bulk_phonon_energy*1000:.1f} meV (GaAs LO)")
print()
print("Running simulation (estimated time: 5-10 minutes)...")
print("=" * 70)
print()

# Run simulation
results = run_iets_simulation(config)

# Save results
filepath = save_results(results)

# Analyze results
V_array = results['V_array']
I_array = results['I_array']
dIdV = results['dIdV']

# Find peak
peak_idx = np.argmax(I_array)
V_peak = V_array[peak_idx]
I_peak = I_array[peak_idx]

# Check for NDR (negative dI/dV after peak)
ndr_region = dIdV[peak_idx+1:] < 0
has_ndr = np.any(ndr_region)

# Find valley if NDR exists
if has_ndr:
    valley_idx = peak_idx + 1 + np.argmin(I_array[peak_idx+1:])
    I_valley = I_array[valley_idx]
    V_valley = V_array[valley_idx]
    pvr = I_peak / I_valley if I_valley > 0 else np.inf
else:
    I_valley = I_peak
    V_valley = V_peak
    pvr = 1.0

print()
print("=" * 70)
print("VALIDATION RESULTS")
print("=" * 70)
print()
print("📊 RTD Characteristics:")
print(f"  Peak voltage: V_peak = {V_peak:.3f} V")
print(f"  Peak current: I_peak = {I_peak*1e9:.3f} nA")
print(f"  Peak conductance: dI/dV_max = {dIdV.max():.3e} S")

if has_ndr:
    print(f"\n✓ Negative Differential Resistance (NDR) detected!")
    print(f"  Valley voltage: V_valley = {V_valley:.3f} V")
    print(f"  Valley current: I_valley = {I_valley*1e9:.3f} nA")
    print(f"  Peak-to-Valley Ratio (PVR) = {pvr:.2f}")
else:
    print(f"\n○ No NDR detected (may need higher bias or different structure)")

print(f"\n⏱ Simulation time: {results['elapsed_time']:.1f} seconds")
print(f"\n💾 Results saved to: {filepath}")

# Sanity checks
print()
print("=" * 70)
print("SANITY CHECKS")
print("=" * 70)

checks_passed = 0
checks_total = 0

# Check 1: Peak voltage in reasonable range
checks_total += 1
if 0.05 < V_peak < 0.6:
    print("✓ Peak voltage in reasonable range (0.05-0.6 V)")
    checks_passed += 1
else:
    print(f"⚠ Peak voltage ({V_peak:.3f} V) outside expected range")

# Check 2: Current magnitude reasonable
checks_total += 1
if 1e-12 < I_peak < 1e-3:
    print(f"✓ Peak current magnitude reasonable ({I_peak*1e9:.1f} nA)")
    checks_passed += 1
else:
    print(f"⚠ Peak current ({I_peak:.2e} A) outside typical range")

# Check 3: Positive current
checks_total += 1
if np.all(I_array >= 0):
    print("✓ Current is positive (no unphysical negative current)")
    checks_passed += 1
else:
    print("⚠ Negative current detected (check physics)")

# Check 4: Monotonic increase up to peak
checks_total += 1
monotonic = np.all(np.diff(I_array[:peak_idx+1]) >= -1e-15)
if monotonic:
    print("✓ Current increases monotonically to peak")
    checks_passed += 1
else:
    print("○ Current has oscillations before peak (numerical noise?)")

# Check 5: Conductance positive at low bias
checks_total += 1
if dIdV[0] > 0 and dIdV[1] > 0:
    print("✓ Positive differential conductance at low bias")
    checks_passed += 1
else:
    print("⚠ Unexpected conductance behavior at low bias")

print()
print("=" * 70)
if checks_passed == checks_total:
    print(f"🎉 ALL CHECKS PASSED ({checks_passed}/{checks_total}) - RTD physics validated!")
elif checks_passed >= checks_total - 1:
    print(f"✓ VALIDATION SUCCESSFUL ({checks_passed}/{checks_total}) - Minor issues only")
else:
    print(f"⚠ PARTIAL VALIDATION ({checks_passed}/{checks_total}) - Review results")
print("=" * 70)

# Expected vs Actual comparison
print()
print("Expected RTD Literature Values (typical double-barrier RTD):")
print("  - Peak voltage: 0.1-0.3 V")
print("  - Peak current density: 10-100 kA/cm² (depends on doping)")
print("  - PVR: 2-10 (for high-quality RTDs)")
print("  - NDR region: present in most symmetric RTDs")
print()
print(f"This Simulation:")
print(f"  - Peak voltage: {V_peak:.3f} V")
area_cm2 = (1e-6 * 1e-6) * 1e4  # 1 μm² in cm²
J_peak = I_peak / area_cm2 / 1000  # kA/cm²
print(f"  - Peak current density: {J_peak:.1f} kA/cm²")
print(f"  - PVR: {pvr:.2f}" if has_ndr else "  - PVR: N/A (no NDR)")
print(f"  - NDR: {'Yes' if has_ndr else 'No'}")
print()
print("=" * 70)
