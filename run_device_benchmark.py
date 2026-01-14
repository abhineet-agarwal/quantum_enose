"""
Device Benchmark: Compare all device stacks

Tests multiple RTD structures with extended bias range to:
1. Find resonance/peak current
2. Measure on/off ratio
3. Identify best device for gas sensing

Uses 1D transport (no quasi-3D) for faster comparison.
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import sys
sys.path.insert(0, '.')

from config.device_library import get_device, get_material, list_all_devices
from core.hamiltonian import discretize_device, build_hamiltonian
from core.self_energy import contact_self_energy_matrix
from core.green_functions import retarded_greens_function, spectral_function, broadening_function

print("="*70)
print("DEVICE BENCHMARK")
print("="*70)

# ============================================================================
# DEVICE SELECTION
# ============================================================================

# Test these devices
devices_to_test = [
    "GaAs_AlAs_symmetric",
    "GaAs_AlAs_asymmetric",
    "GaAs_AlAs_thin",
    "GaAs_AlAs_wide_well",
    "InGaAs_InAlAs_symmetric",
]

print(f"\nDevices to benchmark: {len(devices_to_test)}")
for d in devices_to_test:
    print(f"  - {d}")

# ============================================================================
# BENCHMARK FUNCTION
# ============================================================================

def benchmark_device(device_name, V_array, E_array, temperature=300, verbose=True):
    """
    Benchmark a single device: compute I-V and transmission.
    Uses 1D ballistic transport (no phonons) for speed.
    """

    device = get_device(device_name)
    grid = discretize_device(device, grid_spacing=0.3e-9)
    H, t = build_hamiltonian(grid)

    Np = grid['Np']
    NE = len(E_array)
    NV = len(V_array)

    if verbose:
        print(f"\n  Device: {device_name}")
        print(f"  Grid: {Np} points")

    # Physical constants
    q = 1.602176634e-19
    h_planck = 6.62607015e-34
    kB_eV = 8.617333262145e-5
    kT = kB_eV * temperature

    # Pre-compute contact self-energies
    Sigma1_array = np.zeros((Np, Np, NE), dtype=complex)
    Sigma2_array = np.zeros((Np, Np, NE), dtype=complex)
    Gamma1_array = np.zeros((Np, Np, NE))
    Gamma2_array = np.zeros((Np, Np, NE))

    for iE, E in enumerate(E_array):
        Sigma1_array[:, :, iE] = contact_self_energy_matrix(E, grid, t, 'left')
        Sigma2_array[:, :, iE] = contact_self_energy_matrix(E, grid, t, 'right')
        Gamma1_array[:, :, iE] = broadening_function(Sigma1_array[:, :, iE])
        Gamma2_array[:, :, iE] = broadening_function(Sigma2_array[:, :, iE])

    # Compute transmission vs energy (at zero bias)
    T_vs_E = np.zeros(NE)
    for iE, E in enumerate(E_array):
        G = retarded_greens_function(E, H, Sigma1_array[:, :, iE],
                                     Sigma2_array[:, :, iE], SigmaS=None, eta=1e-4)
        with np.errstate(all='ignore'):
            A2 = G @ Gamma2_array[:, :, iE] @ G.conj().T
            T_vs_E[iE] = np.real(np.trace(Gamma1_array[:, :, iE] @ A2))
        T_vs_E[iE] = np.nan_to_num(T_vs_E[iE], nan=0, posinf=0, neginf=0)

    # Find resonance energy
    peak_idx = np.argmax(T_vs_E)
    E_resonance = E_array[peak_idx]
    T_peak = T_vs_E[peak_idx]

    if verbose:
        print(f"  Resonance: E = {E_resonance:.3f} eV, T = {T_peak:.6e}")

    # Compute I-V curve
    I_array = np.zeros(NV)

    for iV, V in enumerate(V_array):
        mu1 = +V / 2.0
        mu2 = -V / 2.0

        # Fermi functions
        f1 = 1.0 / (1.0 + np.exp((E_array - mu1) / kT))
        f2 = 1.0 / (1.0 + np.exp((E_array - mu2) / kT))

        # Current density: I = (2q²/h) ∫ T(E) [f1(E) - f2(E)] dE
        I_integrand = T_vs_E * (f1 - f2)
        I_array[iV] = (2 * q**2 / h_planck) * np.trapezoid(I_integrand, E_array)

    # Compute derivatives
    dIdV = np.gradient(I_array, V_array)
    d2IdV2 = np.gradient(dIdV, V_array)

    # Find peak current and voltage
    I_peak = np.max(I_array)
    V_peak = V_array[np.argmax(I_array)]
    I_min = np.min(I_array[V_array > V_peak]) if np.any(V_array > V_peak) else I_peak
    on_off_ratio = I_peak / I_min if I_min > 0 else 0

    if verbose:
        print(f"  I_peak = {I_peak*1e9:.3f} nA at V = {V_peak:.3f} V")
        print(f"  On/off ratio: {on_off_ratio:.2f}")

    return {
        'device_name': device_name,
        'V_array': V_array,
        'I_array': I_array,
        'dIdV': dIdV,
        'd2IdV2': d2IdV2,
        'E_array': E_array,
        'T_vs_E': T_vs_E,
        'E_resonance': E_resonance,
        'T_peak': T_peak,
        'I_peak': I_peak,
        'V_peak': V_peak,
        'on_off_ratio': on_off_ratio,
        'Np': Np
    }

# ============================================================================
# RUN BENCHMARKS
# ============================================================================

# Energy grid (wider range to find resonances)
E_array = np.linspace(-0.5, 2.0, 200)

# Bias sweep
V_array = np.linspace(0.0, 1.0, 21)

temperature = 300

print(f"\nParameters:")
print(f"  Energy grid: {len(E_array)} points ({E_array[0]:.1f} to {E_array[-1]:.1f} eV)")
print(f"  Bias sweep: {len(V_array)} points ({V_array[0]:.2f} to {V_array[-1]:.2f} V)")
print(f"  Temperature: {temperature} K")

results = {}
start_total = time.time()

for device_name in devices_to_test:
    print("\n" + "-"*70)
    try:
        result = benchmark_device(device_name, V_array, E_array, temperature, verbose=True)
        results[device_name] = result
    except Exception as e:
        print(f"  ERROR: {e}")
        results[device_name] = None

elapsed_total = time.time() - start_total
print(f"\n\nTotal benchmark time: {elapsed_total:.1f} seconds")

# ============================================================================
# COMPARISON TABLE
# ============================================================================

print("\n" + "="*70)
print("BENCHMARK RESULTS")
print("="*70)

print(f"\n{'Device':<30} {'E_res (eV)':<12} {'T_peak':<12} {'I_peak (nA)':<12} {'V_peak (V)':<10}")
print("-"*80)

for device_name, result in results.items():
    if result is not None:
        print(f"{device_name:<30} {result['E_resonance']:<12.3f} {result['T_peak']:<12.2e} "
              f"{result['I_peak']*1e9:<12.3f} {result['V_peak']:<10.3f}")
    else:
        print(f"{device_name:<30} {'ERROR':<12}")

# ============================================================================
# PLOTS
# ============================================================================

print("\n" + "="*70)
print("CREATING PLOTS")
print("="*70)

# Filter valid results
valid_results = {k: v for k, v in results.items() if v is not None}

if len(valid_results) > 0:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    colors = plt.cm.tab10(np.linspace(0, 1, len(valid_results)))

    # Plot 1: Transmission vs Energy
    ax1 = axes[0, 0]
    for i, (name, result) in enumerate(valid_results.items()):
        ax1.semilogy(result['E_array'], result['T_vs_E'], '-',
                     color=colors[i], linewidth=1.5, label=name[:20])
    ax1.set_xlabel('Energy (eV)', fontsize=12)
    ax1.set_ylabel('Transmission', fontsize=12)
    ax1.set_title('Transmission vs Energy', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim([1e-15, 1])

    # Plot 2: I-V Curves
    ax2 = axes[0, 1]
    for i, (name, result) in enumerate(valid_results.items()):
        ax2.plot(result['V_array'], result['I_array']*1e9, '-o',
                 color=colors[i], linewidth=1.5, markersize=4, label=name[:20])
    ax2.set_xlabel('Bias Voltage (V)', fontsize=12)
    ax2.set_ylabel('Current (nA)', fontsize=12)
    ax2.set_title('I-V Curves', fontsize=13, fontweight='bold')
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    # Plot 3: dI/dV
    ax3 = axes[1, 0]
    for i, (name, result) in enumerate(valid_results.items()):
        ax3.plot(result['V_array'], result['dIdV']*1e6, '-o',
                 color=colors[i], linewidth=1.5, markersize=4, label=name[:20])
    ax3.set_xlabel('Bias Voltage (V)', fontsize=12)
    ax3.set_ylabel('dI/dV (µS)', fontsize=12)
    ax3.set_title('Differential Conductance', fontsize=13, fontweight='bold')
    ax3.legend(fontsize=8)
    ax3.grid(True, alpha=0.3)

    # Plot 4: Bar chart of key metrics
    ax4 = axes[1, 1]
    device_names_short = [n[:15] for n in valid_results.keys()]
    I_peaks = [r['I_peak']*1e9 for r in valid_results.values()]

    bars = ax4.bar(device_names_short, I_peaks, color=colors[:len(valid_results)])
    ax4.set_ylabel('Peak Current (nA)', fontsize=12)
    ax4.set_title('Peak Current Comparison', fontsize=13, fontweight='bold')
    ax4.tick_params(axis='x', rotation=45)

    for bar, val in zip(bars, I_peaks):
        ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                 f'{val:.1f}', ha='center', fontsize=9)

    plt.suptitle('Device Stack Benchmark (1D Ballistic Transport)',
                 fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.98])

    output_file = 'device_benchmark_results.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nPlots saved to: {output_file}")

    # Save data
    np.savez('device_benchmark_data.npz',
             device_names=list(valid_results.keys()),
             E_array=E_array,
             V_array=V_array,
             **{f"{n.replace('-', '_').replace(' ', '_')}_T": r['T_vs_E'] for n, r in valid_results.items()},
             **{f"{n.replace('-', '_').replace(' ', '_')}_I": r['I_array'] for n, r in valid_results.items()})
    print(f"Data saved to: device_benchmark_data.npz")

print("\n" + "="*70)
print("DEVICE BENCHMARK COMPLETE")
print("="*70)

# Summary recommendations
print("\n" + "="*70)
print("RECOMMENDATIONS")
print("="*70)

if len(valid_results) > 0:
    # Find best device by current
    best_current = max(valid_results.items(), key=lambda x: x[1]['I_peak'])
    # Find device with lowest resonance energy
    best_resonance = min(valid_results.items(), key=lambda x: x[1]['E_resonance'])

    print(f"\n1. Highest current: {best_current[0]}")
    print(f"   I_peak = {best_current[1]['I_peak']*1e9:.2f} nA")

    print(f"\n2. Lowest resonance (best for low-bias): {best_resonance[0]}")
    print(f"   E_res = {best_resonance[1]['E_resonance']*1000:.1f} meV")

    print("\n3. For gas sensing applications:")
    print("   - Use thin barriers for high current (better SNR)")
    print("   - Use asymmetric design for position selectivity")
    print("   - Wide well gives lower resonance (low-bias operation)")

print("\n" + "="*70)
