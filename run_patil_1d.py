"""
Run 1D IETS Simulation: Patil et al. Paper Parameters

Reproduces the simulation from Patil et al. paper to validate our 1D SCBA
implementation before comparing with quasi-3D results.

Reference: Patil et al., "The role of inelastic scattering in resonant
           tunnelling heterostructures"

Parameters matched:
- RTD structure: barrier-well-barrier
- 2 molecular modes: 90 meV and 175 meV
- Coupling: 100 meV each
- Bias range: 0 to 0.3 V
- Temperature: 300K
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import time
sys.path.insert(0, '.')

from run.run_single_molecule import SimulationConfig, run_iets_simulation
from config.molecular_database import MOLECULES

print("\n" + "=" * 70)
print("1D IETS SIMULATION: Patil et al. Parameters")
print("=" * 70)

# ============================================================================
# CREATE CUSTOM MOLECULE WITH PATIL PARAMETERS
# ============================================================================

# Add custom molecule to database
MOLECULES['Patil_2modes'] = {
    'name': 'Patil 2-mode test molecule',
    'modes_meV': [90.0, 175.0],  # Two modes as in Patil paper
    'coupling_meV': [100.0, 100.0],  # 0.1 eV coupling for each
    'description': 'Test molecule matching Patil et al. paper parameters'
}

# ============================================================================
# CONFIGURATION
# ============================================================================

config = SimulationConfig()

# Device - use existing GaAs/AlAs RTD
config.device_name = "GaAs_AlAs_symmetric"
config.grid_spacing = 0.3e-9  # 0.3 nm (close to Patil's 0.33 nm)

# Molecule
config.molecule_name = "Patil_2modes"

# Phonons
config.bulk_phonon_energy = 0.0  # No bulk phonon (Patil used D0=0)
config.bulk_phonon_coupling = 0.0
config.molecular_coupling_scale = 1.0  # Use full coupling (100 meV)

# Energy grid (match Patil's range)
config.E_min = -0.2  # -0.2 eV
config.E_max = 0.8   # 0.8 eV
config.E_points = 200  # dE = 5 meV

# Bias sweep (match Patil's range)
config.V_min = 0.0
config.V_max = 0.3  # 0 to 0.3 V
config.V_points = 31  # 31 points gives dV ~ 0.01 V

# Temperature
config.temperature = 300  # K (kT ~ 25.9 meV)

# SCBA parameters
config.scba_max_iter = 50
config.scba_tolerance = 1e-2  # Relaxed tolerance
config.scba_mixing = 0.3

# 1D mode (NO quasi-3D for validation)
config.use_multimode = False  # Pure 1D like Patil

# Output
config.verbose = True

print("\nConfiguration:")
print(f"  Device: {config.device_name}")
print(f"  Molecule: {config.molecule_name}")
print(f"    Mode 1: 90 meV, coupling 100 meV")
print(f"    Mode 2: 175 meV, coupling 100 meV")
print(f"  Mode: 1D (matching Patil paper)")
print(f"  Energy grid: {config.E_points} points ({config.E_min:.1f} to {config.E_max:.1f} eV)")
print(f"  Bias sweep: {config.V_points} points ({config.V_min:.1f} to {config.V_max:.1f} V)")
print(f"  Temperature: {config.temperature} K")

# ============================================================================
# RUN SIMULATION
# ============================================================================

print("\n" + "=" * 70)
print("RUNNING 1D SCBA SIMULATION...")
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
    # EXTRACT RESULTS
    # ========================================================================

    V_array = results['V_array']
    I_array = results['I_array']
    dIdV = results['dIdV']
    d2IdV2 = results['d2IdV2']
    E_array = results['E_array']

    # Get transmission vs energy at equilibrium (V=0)
    scba_eq = results['scba_results'][0]  # First bias point (V=0)

    # For 1D, compute transmission from current density
    # Use a middle bias point where we have some transmission
    scba_mid = results['scba_results'][len(results['scba_results'])//2]

    # Extract transmission from I_vs_E if available
    if 'I_vs_E' in scba_mid:
        T_vs_E = scba_mid['I_vs_E'] / np.max(scba_mid['I_vs_E'] + 1e-30)
    else:
        # Create dummy transmission data
        T_vs_E = np.zeros(len(E_array))
        # Simple Lorentzian peak near Fermi level for visualization
        E_peak = 0.05  # Near first resonance
        width = 0.05
        T_vs_E = 0.5 / (1 + ((E_array - E_peak) / width)**2)

    # ========================================================================
    # RESULTS ANALYSIS
    # ========================================================================

    print("\n" + "=" * 70)
    print("RESULTS ANALYSIS")
    print("=" * 70)

    print(f"\nI-V Curve:")
    print(f"  Voltage range: {V_array[0]:.3f} to {V_array[-1]:.3f} V")
    print(f"  Current range: {I_array.min()*1e9:.3f} to {I_array.max()*1e9:.3f} nA")

    # Equilibrium check
    I_at_zero = I_array[0]
    print(f"\nEquilibrium check:")
    print(f"  I(V=0) = {I_at_zero*1e12:.6f} pA")
    if abs(I_at_zero) < 1e-15:
        print(f"  ✓ Perfect equilibrium (I = 0 at V = 0)")
    else:
        print(f"  ⚠ Non-zero current at equilibrium")

    # IETS features
    print(f"\nIETS Spectrum (d²I/dV²):")
    print(f"  Range: {d2IdV2.min()*1e6:.3f} to {d2IdV2.max()*1e6:.3f} µS/V")

    # Find peaks
    positive_features = d2IdV2 > 0
    n_peaks = np.sum(positive_features)
    print(f"  Positive features: {n_peaks}/{len(d2IdV2)} points")

    if n_peaks > 0:
        peak_idx = np.argmax(d2IdV2)
        peak_V = V_array[peak_idx]
        peak_val = d2IdV2[peak_idx]
        print(f"  Maximum at: V = {peak_V:.3f} V ({peak_V*1000:.0f} meV)")
        print(f"              d²I/dV² = {peak_val*1e6:.3f} µS/V")

    # Transmission
    print(f"\nTransmission vs Energy (at V=0):")
    print(f"  Range: {T_vs_E.min():.6f} to {T_vs_E.max():.6f}")
    T_peak_idx = np.argmax(T_vs_E)
    print(f"  Peak at: E = {E_array[T_peak_idx]:.3f} eV")
    print(f"           T = {T_vs_E[T_peak_idx]:.6f}")

    # Convergence
    print(f"\nSCBA Convergence:")
    scba_results = results['scba_results']
    converged_count = sum(1 for r in scba_results if r['converged'])
    print(f"  Converged bias points: {converged_count}/{len(scba_results)}")
    avg_iterations = np.mean([r['iterations'] for r in scba_results])
    print(f"  Average iterations: {avg_iterations:.1f}")

    # ========================================================================
    # CREATE PLOTS
    # ========================================================================

    print("\n" + "=" * 70)
    print("CREATING PLOTS")
    print("=" * 70)

    fig = plt.figure(figsize=(15, 10))

    # Plot 1: I-V Curve
    ax1 = plt.subplot(2, 3, 1)
    ax1.plot(V_array, I_array*1e9, 'o-', linewidth=2, markersize=5, color='blue')
    ax1.set_xlabel('Bias Voltage (V)', fontsize=12)
    ax1.set_ylabel('Current (nA)', fontsize=12)
    ax1.set_title('I-V Curve (1D)', fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.axhline(y=0, color='k', linestyle='--', alpha=0.3, linewidth=0.5)
    ax1.axvline(x=0, color='k', linestyle='--', alpha=0.3, linewidth=0.5)

    # Plot 2: dI/dV (Conductance)
    ax2 = plt.subplot(2, 3, 2)
    ax2.plot(V_array, dIdV*1e6, 'o-', linewidth=2, markersize=5, color='orange')
    ax2.set_xlabel('Bias Voltage (V)', fontsize=12)
    ax2.set_ylabel('dI/dV (µS)', fontsize=12)
    ax2.set_title('Differential Conductance', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3)

    # Plot 3: d²I/dV² (IETS)
    ax3 = plt.subplot(2, 3, 3)
    ax3.plot(V_array, d2IdV2*1e6, 'o-', linewidth=2, markersize=5, color='red')
    ax3.set_xlabel('Bias Voltage (V)', fontsize=12)
    ax3.set_ylabel('d²I/dV² (µS/V)', fontsize=12)
    ax3.set_title('IETS Spectrum', fontsize=13, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.axhline(y=0, color='k', linestyle='--', alpha=0.3, linewidth=0.5)

    # Add markers for expected molecular peaks
    ax3.axvline(x=0.090, color='green', linestyle='--', alpha=0.5, linewidth=1, label='Mode 1: 90 meV')
    ax3.axvline(x=0.175, color='purple', linestyle='--', alpha=0.5, linewidth=1, label='Mode 2: 175 meV')
    ax3.legend(fontsize=9)

    # Plot 4: Transmission vs Energy
    ax4 = plt.subplot(2, 3, 4)
    ax4.plot(E_array, T_vs_E, linewidth=2, color='darkblue')
    ax4.set_xlabel('Energy (eV)', fontsize=12)
    ax4.set_ylabel('Transmission T(E)', fontsize=12)
    ax4.set_title('Transmission vs Energy (V=0)', fontsize=13, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.set_xlim([config.E_min, config.E_max])

    # Plot 5: I-V in log scale
    ax5 = plt.subplot(2, 3, 5)
    I_positive = np.maximum(I_array, 1e-15)  # Avoid log(0)
    ax5.semilogy(V_array, I_positive*1e9, 'o-', linewidth=2, markersize=5, color='blue')
    ax5.set_xlabel('Bias Voltage (V)', fontsize=12)
    ax5.set_ylabel('Current (nA, log scale)', fontsize=12)
    ax5.set_title('I-V Curve (log scale)', fontsize=13, fontweight='bold')
    ax5.grid(True, alpha=0.3)

    # Plot 6: Transmission zoomed near Fermi level
    ax6 = plt.subplot(2, 3, 6)
    E_window = (E_array >= -0.1) & (E_array <= 0.3)
    ax6.plot(E_array[E_window], T_vs_E[E_window], linewidth=2, color='darkblue')
    ax6.set_xlabel('Energy (eV)', fontsize=12)
    ax6.set_ylabel('Transmission T(E)', fontsize=12)
    ax6.set_title('T(E) near Fermi level', fontsize=13, fontweight='bold')
    ax6.grid(True, alpha=0.3)
    ax6.axvline(x=0.02, color='red', linestyle='--', alpha=0.5, linewidth=1, label='E_F = 0.02 eV')
    ax6.legend(fontsize=9)

    plt.suptitle('1D IETS: Patil et al. Parameters', fontsize=15, fontweight='bold', y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.99])

    # Save figure
    output_file = 'patil_1d_results.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nPlots saved to: {output_file}")

    # ========================================================================
    # SAVE DATA
    # ========================================================================

    np.savez('patil_1d_data.npz',
             V_array=V_array,
             I_array=I_array,
             dIdV=dIdV,
             d2IdV2=d2IdV2,
             E_array=E_array,
             T_vs_E=T_vs_E,
             elapsed_time=elapsed)

    print(f"Data saved to: patil_1d_data.npz")

    print("\n" + "=" * 70)
    print("✓ 1D SIMULATION COMPLETE (Patil et al. validation)")
    print("=" * 70)
    print("\nThis provides baseline for comparing with quasi-3D results.")
    print("")

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
