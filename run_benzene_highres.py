"""
Run High-Resolution IETS: Benzene with Quasi-3D

41 bias points (0 to 0.4V, ΔV = 0.01V) to resolve individual vibrational modes.

Expected Benzene peaks:
- Mode 1: 49.5 meV  (~0.050 V)
- Mode 2: 79.0 meV  (~0.079 V)
- Mode 3: 134.4 meV (~0.134 V)
- Mode 4: 184.1 meV (~0.184 V)
- Mode 5: 395.4 meV (~0.395 V)
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import time
sys.path.insert(0, '.')

from run.run_single_molecule import SimulationConfig, run_iets_simulation

print("\n" + "=" * 70)
print("HIGH-RESOLUTION IETS: Benzene Quasi-3D")
print("=" * 70)

# ============================================================================
# CONFIGURATION
# ============================================================================

config = SimulationConfig()

# Device
config.device_name = "GaAs_AlAs_symmetric"
config.grid_spacing = 0.3e-9  # 0.3 nm (optimized)

# Molecule
config.molecule_name = "Benzene"

# Phonons
config.bulk_phonon_energy = 0.036  # 36 meV (GaAs LO phonon)
config.bulk_phonon_coupling = 0.010  # 10 meV
config.molecular_coupling_scale = 1.0

# Energy grid
config.E_min = -0.2
config.E_max = 1.0
config.E_points = 50  # Keep coarse for speed

# Bias sweep - HIGH RESOLUTION
config.V_min = 0.0
config.V_max = 0.4  # Extended to 0.4V to see mode 5 (395 meV)
config.V_points = 41  # ΔV = 0.01V - should resolve individual modes

# Temperature
config.temperature = 300  # K

# SCBA parameters
config.scba_max_iter = 50
config.scba_tolerance = 1e-2
config.scba_mixing = 0.3

# Quasi-3D parameters
config.use_multimode = True
config.Ly = 1.0e-6
config.Lz = 1.0e-6
config.n_max = 3
config.m_max = 3

# Output
config.verbose = True

print("\nConfiguration:")
print(f"  Device: {config.device_name}")
print(f"  Molecule: {config.molecule_name}")
print(f"  Quasi-3D: {config.n_max}×{config.m_max} = {config.n_max*config.m_max} transverse modes")
print(f"  Energy grid: {config.E_points} points ({config.E_min:.1f} to {config.E_max:.1f} eV)")
print(f"  Bias sweep: {config.V_points} points ({config.V_min:.1f} to {config.V_max:.1f} V)")
print(f"  Resolution: ΔV = {(config.V_max-config.V_min)/(config.V_points-1)*1000:.1f} mV")
print(f"  SCBA: max {config.scba_max_iter} iterations, tol={config.scba_tolerance:.0e}")

# Estimated time
est_time = 41 * 86.2  # 41 points × 86.2 sec/point
print(f"\nEstimated time: {est_time:.0f} seconds ({est_time/60:.1f} minutes)")

# ============================================================================
# RUN SIMULATION
# ============================================================================

print("\n" + "=" * 70)
print("RUNNING HIGH-RESOLUTION SIMULATION...")
print("=" * 70)
print("(This will take ~60 minutes)")

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
    # EXTRACT RESULTS
    # ========================================================================

    V_array = results['V_array']
    I_array = results['I_array']
    dIdV = results['dIdV']
    d2IdV2 = results['d2IdV2']

    # ========================================================================
    # ANALYSIS
    # ========================================================================

    print("\n" + "=" * 70)
    print("RESULTS ANALYSIS")
    print("=" * 70)

    print(f"\nI-V Curve:")
    print(f"  Voltage range: {V_array[0]:.3f} to {V_array[-1]:.3f} V")
    print(f"  Current range: {I_array.min()*1e9:.3f} to {I_array.max()*1e9:.3f} nA")
    print(f"  I(V=0) = {I_array[0]*1e12:.6f} pA")

    print(f"\nIETS Spectrum (d²I/dV²):")
    print(f"  Range: {d2IdV2.min()*1e6:.3f} to {d2IdV2.max()*1e6:.3f} µS/V")

    # Find all local maxima in d²I/dV²
    from scipy.signal import find_peaks
    peaks, properties = find_peaks(d2IdV2, prominence=1e-6)  # Adjust prominence as needed

    print(f"\n  Detected peaks: {len(peaks)}")
    if len(peaks) > 0:
        print(f"\n  Peak analysis:")
        print(f"  {'Peak':<6} {'Voltage (V)':<15} {'Energy (meV)':<15} {'d²I/dV² (µS/V)'}")
        print(f"  {'-'*60}")
        for i, peak_idx in enumerate(peaks):
            V_peak = V_array[peak_idx]
            E_peak = V_peak * 1000  # meV
            d2I_peak = d2IdV2[peak_idx] * 1e6
            print(f"  {i+1:<6} {V_peak:<15.3f} {E_peak:<15.1f} {d2I_peak:<15.3f}")

    # Compare with expected Benzene modes
    benzene_modes = [49.5, 79.0, 134.4, 184.1, 395.4]  # meV
    print(f"\n  Expected Benzene modes:")
    for i, E_mode in enumerate(benzene_modes):
        V_mode = E_mode / 1000.0
        print(f"    Mode {i+1}: {E_mode:6.1f} meV (~{V_mode:.3f} V)")

    # Convergence
    print(f"\nSCBA Convergence:")
    scba_results = results['scba_results']
    converged_count = sum(1 for r in scba_results if r['converged'])
    print(f"  Converged bias points: {converged_count}/{len(scba_results)}")
    avg_iterations = np.mean([r['iterations'] for r in scba_results])
    print(f"  Average iterations: {avg_iterations:.1f}")

    # ========================================================================
    # CREATE HIGH-RESOLUTION PLOTS
    # ========================================================================

    print("\n" + "=" * 70)
    print("CREATING HIGH-RESOLUTION PLOTS")
    print("=" * 70)

    fig = plt.figure(figsize=(16, 10))

    # Plot 1: I-V Curve
    ax1 = plt.subplot(2, 3, 1)
    ax1.plot(V_array, I_array*1e9, '-', linewidth=1.5, color='blue')
    ax1.set_xlabel('Bias Voltage (V)', fontsize=12)
    ax1.set_ylabel('Current (nA)', fontsize=12)
    ax1.set_title('I-V Curve (High-Res, 41 points)', fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.axhline(y=0, color='k', linestyle='--', alpha=0.3, linewidth=0.5)

    # Plot 2: dI/dV
    ax2 = plt.subplot(2, 3, 2)
    ax2.plot(V_array, dIdV*1e6, '-', linewidth=1.5, color='orange')
    ax2.set_xlabel('Bias Voltage (V)', fontsize=12)
    ax2.set_ylabel('dI/dV (µS)', fontsize=12)
    ax2.set_title('Differential Conductance', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3)

    # Plot 3: d²I/dV² (IETS) - Main result
    ax3 = plt.subplot(2, 3, 3)
    ax3.plot(V_array, d2IdV2*1e6, '-', linewidth=2, color='red')
    ax3.set_xlabel('Bias Voltage (V)', fontsize=12)
    ax3.set_ylabel('d²I/dV² (µS/V)', fontsize=12)
    ax3.set_title('IETS Spectrum (High Resolution)', fontsize=13, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.axhline(y=0, color='k', linestyle='--', alpha=0.3, linewidth=0.5)

    # Add vertical lines for expected Benzene modes
    colors_modes = ['green', 'purple', 'brown', 'pink', 'cyan']
    for i, E_mode in enumerate(benzene_modes):
        V_mode = E_mode / 1000.0
        if V_mode <= config.V_max:
            ax3.axvline(x=V_mode, color=colors_modes[i], linestyle='--',
                       alpha=0.5, linewidth=1.5,
                       label=f'Mode {i+1}: {E_mode:.1f} meV')
    ax3.legend(fontsize=9, loc='upper left')

    # Mark detected peaks
    if len(peaks) > 0:
        ax3.plot(V_array[peaks], d2IdV2[peaks]*1e6, 'rx', markersize=10,
                markeredgewidth=2, label='Detected peaks')

    # Plot 4: IETS zoomed (0-0.2V) to see low-energy modes
    ax4 = plt.subplot(2, 3, 4)
    mask_zoom = V_array <= 0.2
    ax4.plot(V_array[mask_zoom], d2IdV2[mask_zoom]*1e6, '-', linewidth=2, color='red')
    ax4.set_xlabel('Bias Voltage (V)', fontsize=12)
    ax4.set_ylabel('d²I/dV² (µS/V)', fontsize=12)
    ax4.set_title('IETS: 0-200 meV (Modes 1-4)', fontsize=13, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.axhline(y=0, color='k', linestyle='--', alpha=0.3, linewidth=0.5)

    # Add mode markers
    for i in range(4):  # First 4 modes
        V_mode = benzene_modes[i] / 1000.0
        ax4.axvline(x=V_mode, color=colors_modes[i], linestyle='--',
                   alpha=0.5, linewidth=1.5, label=f'{benzene_modes[i]:.1f} meV')
    ax4.legend(fontsize=9)

    # Plot 5: I-V log scale
    ax5 = plt.subplot(2, 3, 5)
    I_positive = np.maximum(I_array, 1e-15)
    ax5.semilogy(V_array, I_positive*1e9, '-', linewidth=1.5, color='blue')
    ax5.set_xlabel('Bias Voltage (V)', fontsize=12)
    ax5.set_ylabel('Current (nA, log scale)', fontsize=12)
    ax5.set_title('I-V (log scale)', fontsize=13, fontweight='bold')
    ax5.grid(True, alpha=0.3)

    # Plot 6: Normalized IETS (divide by conductance)
    ax6 = plt.subplot(2, 3, 6)
    # Normalized IETS: (d²I/dV²)/(dI/dV)
    dIdV_safe = np.maximum(dIdV, 1e-12)  # Avoid division by zero
    iets_normalized = d2IdV2 / dIdV_safe
    ax6.plot(V_array, iets_normalized, '-', linewidth=2, color='darkred')
    ax6.set_xlabel('Bias Voltage (V)', fontsize=12)
    ax6.set_ylabel('(d²I/dV²)/(dI/dV) [V⁻¹]', fontsize=12)
    ax6.set_title('Normalized IETS', fontsize=13, fontweight='bold')
    ax6.grid(True, alpha=0.3)
    ax6.axhline(y=0, color='k', linestyle='--', alpha=0.3, linewidth=0.5)

    # Add mode markers
    for i, E_mode in enumerate(benzene_modes):
        V_mode = E_mode / 1000.0
        if V_mode <= config.V_max:
            ax6.axvline(x=V_mode, color=colors_modes[i], linestyle='--',
                       alpha=0.3, linewidth=1)

    plt.suptitle('High-Resolution Benzene IETS (Quasi-3D, 9 modes)',
                 fontsize=15, fontweight='bold', y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.99])

    output_file = 'benzene_highres_results.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nPlots saved to: {output_file}")

    # Save data
    np.savez('benzene_highres_data.npz',
             V_array=V_array,
             I_array=I_array,
             dIdV=dIdV,
             d2IdV2=d2IdV2,
             peaks=peaks if len(peaks) > 0 else np.array([]),
             elapsed_time=elapsed)
    print(f"Data saved to: benzene_highres_data.npz")

    print("\n" + "=" * 70)
    print("✓ HIGH-RESOLUTION IETS COMPLETE!")
    print("=" * 70)

except Exception as e:
    elapsed = time.time() - start_time
    print(f"\n" + "=" * 70)
    print("✗ SIMULATION FAILED")
    print("=" * 70)
    print(f"\nError after {elapsed:.1f} seconds:")
    print(f"  {type(e).__name__}: {e}")

    import traceback
    traceback.print_exc()
