"""
Oxide RTD Baseline Characterization

Run baseline simulations (no molecule) on all BEOL-compatible oxide RTD devices.
This characterizes the I-V curves, transmission, and identifies optimal operating points.

Key considerations for manufacturability:
- ALD minimum thickness: ~1 nm (but 1.5+ nm safer for uniformity)
- Barrier thickness: 1.5-2.5 nm
- Well thickness: 2-4 nm (determines resonance energy)
- Contact thickness: 8-15 nm
"""

import numpy as np
import sys
import os
import time
import argparse
import matplotlib.pyplot as plt

# Setup path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

from config.device_library import MATERIALS, DEVICES, create_custom_device, get_device, get_band_offset, BAND_OFFSETS
from core.hamiltonian import discretize_device, build_hamiltonian
from core.self_energy import contact_self_energy_matrix
from core.green_functions import (
    retarded_greens_function, spectral_function, broadening_function
)

# Physical constants
hbar = 1.054571817e-34  # J*s
m0 = 9.10938356e-31     # kg
q = 1.602176634e-19     # C
kB = 1.380649e-23       # J/K

# ============================================================================
# OPTIMIZED OXIDE DEVICE CONFIGURATIONS
# ============================================================================

def calculate_resonance_energy(well_width_nm, m_eff_rel):
    """Calculate ground state energy in quantum well"""
    L = well_width_nm * 1e-9
    m_eff = m_eff_rel * m0
    E1 = (np.pi**2 * hbar**2) / (2 * m_eff * L**2)
    return E1 / q  # Return in eV


def print_device_physics():
    """Print expected physics for oxide devices"""
    print("\n" + "="*70)
    print("OXIDE RTD PHYSICS ANALYSIS")
    print("="*70)

    # Well materials
    wells = ['In2O3', 'IGZO', 'ZnO', 'SnO2']
    barriers = ['Al2O3', 'HfO2']

    print("\nWell material properties:")
    print(f"  {'Material':<10} {'m*/m0':<8} {'Phonon (meV)':<15}")
    print(f"  {'-'*35}")
    for mat in wells:
        m = MATERIALS[mat]
        print(f"  {mat:<10} {m['m_eff']:<8.3f} {m['phonon_energy_meV']:<15.1f}")

    print("\nBarrier material properties:")
    print(f"  {'Material':<10} {'m*/m0':<8} {'Ec (eV)':<10} {'Notes':<30}")
    print(f"  {'-'*60}")
    for mat in barriers:
        m = MATERIALS[mat]
        print(f"  {mat:<10} {m['m_eff']:<8.3f} {m['Ec']:<10.2f} {m.get('notes', '')[:30]}")

    print("\nResonance energies for different well widths:")
    print(f"  {'Well (nm)':<12}", end="")
    for mat in wells:
        print(f"{mat:<10}", end="")
    print()
    print(f"  {'-'*52}")

    for L in [2.0, 2.5, 3.0, 3.5, 4.0, 5.0]:
        print(f"  {L:<12.1f}", end="")
        for mat in wells:
            E1 = calculate_resonance_energy(L, MATERIALS[mat]['m_eff'])
            print(f"{E1*1000:<10.1f}", end="")  # meV
        print()
    print("  (energies in meV)")

    print("\nExpected bias for resonance (V ≈ 2×E1):")
    print(f"  {'Well (nm)':<12}", end="")
    for mat in wells:
        print(f"{mat:<10}", end="")
    print()
    print(f"  {'-'*52}")

    for L in [2.0, 2.5, 3.0, 3.5, 4.0, 5.0]:
        print(f"  {L:<12.1f}", end="")
        for mat in wells:
            E1 = calculate_resonance_energy(L, MATERIALS[mat]['m_eff'])
            V_res = 2 * E1  # Approximate
            print(f"{V_res:<10.2f}", end="")
        print()
    print("  (voltage in V)")
    print("="*70 + "\n")


# Optimized oxide device configurations
OXIDE_DEVICES_OPTIMIZED = {
    # ========== In2O3/Al2O3 variants ==========
    "In2O3_Al2O3_thin_barrier": {
        "description": "In2O3/Al2O3 with thin barriers (1.5nm) for higher current",
        "layers": [
            {"material": "In2O3", "thickness": 10e-9, "doping": 1e25},
            {"material": "Al2O3", "thickness": 1.5e-9, "doping": 0},  # Thin
            {"material": "In2O3", "thickness": 3.0e-9, "doping": 0},  # 3nm well
            {"material": "Al2O3", "thickness": 1.5e-9, "doping": 0},
            {"material": "In2O3", "thickness": 10e-9, "doping": 1e25}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
    },

    "In2O3_Al2O3_wide_well": {
        "description": "In2O3/Al2O3 with wide well (4nm) for lower resonance",
        "layers": [
            {"material": "In2O3", "thickness": 10e-9, "doping": 1e25},
            {"material": "Al2O3", "thickness": 1.5e-9, "doping": 0},
            {"material": "In2O3", "thickness": 4.0e-9, "doping": 0},  # Wide well
            {"material": "Al2O3", "thickness": 1.5e-9, "doping": 0},
            {"material": "In2O3", "thickness": 10e-9, "doping": 1e25}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
    },

    # ========== In2O3/HfO2 (lower barrier) ==========
    "In2O3_HfO2_optimized": {
        "description": "In2O3/HfO2 optimized (lower barrier height)",
        "layers": [
            {"material": "In2O3", "thickness": 10e-9, "doping": 1e25},
            {"material": "HfO2", "thickness": 1.5e-9, "doping": 0},
            {"material": "In2O3", "thickness": 3.0e-9, "doping": 0},
            {"material": "HfO2", "thickness": 1.5e-9, "doping": 0},
            {"material": "In2O3", "thickness": 10e-9, "doping": 1e25}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
    },

    # ========== IGZO variants ==========
    "IGZO_Al2O3_thin_barrier": {
        "description": "IGZO/Al2O3 with thin barriers",
        "layers": [
            {"material": "IGZO", "thickness": 10e-9, "doping": 8e24},
            {"material": "Al2O3", "thickness": 1.5e-9, "doping": 0},
            {"material": "IGZO", "thickness": 3.0e-9, "doping": 0},
            {"material": "Al2O3", "thickness": 1.5e-9, "doping": 0},
            {"material": "IGZO", "thickness": 10e-9, "doping": 8e24}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
    },

    "IGZO_HfO2_optimized": {
        "description": "IGZO/HfO2 optimized (amorphous + lower barrier)",
        "layers": [
            {"material": "IGZO", "thickness": 10e-9, "doping": 8e24},
            {"material": "HfO2", "thickness": 1.5e-9, "doping": 0},
            {"material": "IGZO", "thickness": 3.5e-9, "doping": 0},
            {"material": "HfO2", "thickness": 1.5e-9, "doping": 0},
            {"material": "IGZO", "thickness": 10e-9, "doping": 8e24}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
    },

    # ========== ZnO variants ==========
    "ZnO_Al2O3_thin_barrier": {
        "description": "ZnO/Al2O3 with thin barriers",
        "layers": [
            {"material": "ZnO", "thickness": 10e-9, "doping": 1e25},
            {"material": "Al2O3", "thickness": 1.5e-9, "doping": 0},
            {"material": "ZnO", "thickness": 3.0e-9, "doping": 0},
            {"material": "Al2O3", "thickness": 1.5e-9, "doping": 0},
            {"material": "ZnO", "thickness": 10e-9, "doping": 1e25}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
    },

    # ========== SnO2 variants ==========
    "SnO2_Al2O3_thin_barrier": {
        "description": "SnO2/Al2O3 with thin barriers",
        "layers": [
            {"material": "SnO2", "thickness": 10e-9, "doping": 1e25},
            {"material": "Al2O3", "thickness": 1.5e-9, "doping": 0},
            {"material": "SnO2", "thickness": 3.0e-9, "doping": 0},
            {"material": "Al2O3", "thickness": 1.5e-9, "doping": 0},
            {"material": "SnO2", "thickness": 10e-9, "doping": 1e25}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
    },
}


# ============================================================================
# BASELINE SIMULATION (BALLISTIC - NO MOLECULE)
# ============================================================================

def run_ballistic_simulation(device, grid_spacing=0.12e-9, E_points=400,
                             V_min=0.0, V_max=5.0, V_points=51,
                             temperature=300, verbose=True):
    """
    Run ballistic (coherent) transport simulation with proper bias treatment.

    Includes bias-induced band bending (essential for NDR in high-barrier oxides).
    """

    # Discretize device
    grid = discretize_device(device, grid_spacing=grid_spacing)
    Np = grid['Np']
    x = grid['x']
    L = x[-1]

    # Base Hamiltonian (at V=0)
    H0, t = build_hamiltonian(grid)

    if verbose:
        print(f"  Grid: {Np} points, {L*1e9:.1f} nm total")

    # Energy range - need to cover resonance and barriers
    # For oxides with ~2-3 eV barrier, need wide range
    E_min = -0.5
    E_max = max(4.0, np.max(grid['Ec']) + 0.5)  # Cover barrier + margin
    E_array = np.linspace(E_min, E_max, E_points)
    dE = E_array[1] - E_array[0]

    # Thermal energy
    kT = kB * temperature / q  # in eV

    # Fermi function
    def fermi(E, mu):
        x_val = (E - mu) / kT
        x_val = np.clip(x_val, -500, 500)
        return 1.0 / (1.0 + np.exp(x_val))

    # Use a larger eta for numerical stability with high barriers
    eta = 1e-5

    # Compute transmission at V=0 (for analysis)
    T_vs_E = np.zeros(E_points)

    def Sigma1(E):
        return contact_self_energy_matrix(E, grid, t, 'left')

    def Sigma2(E):
        return contact_self_energy_matrix(E, grid, t, 'right')

    for iE, E in enumerate(E_array):
        try:
            S1 = Sigma1(E)
            S2 = Sigma2(E)
            Gamma1 = broadening_function(S1)
            Gamma2 = broadening_function(S2)

            H_eff = (E + 1j*eta) * np.eye(Np) - H0 - S1 - S2
            G = np.linalg.solve(H_eff, np.eye(Np))

            T = np.real(np.trace(Gamma1 @ G @ Gamma2 @ G.conj().T))
            T_vs_E[iE] = max(0, min(T, 10))
        except:
            T_vs_E[iE] = 0

    # Bias sweep with proper band bending
    V_array = np.linspace(V_min, V_max, V_points)
    I_array = np.zeros(V_points)

    for iV, V in enumerate(V_array):
        if verbose and iV % 10 == 0:
            print(f"  V = {V*1000:.0f} mV...")

        # LINEAR POTENTIAL DROP across device (bias-induced band bending)
        # This is ESSENTIAL for NDR in RTDs
        U_bias = -V * (x / L)

        # Build Hamiltonian with bias
        H = H0.copy()
        for i in range(Np):
            H[i, i] += U_bias[i]

        # Contact self-energies with bias shift
        # Left contact at reference potential
        # Right contact shifted down by V
        def Sigma1_V(E):
            return contact_self_energy_matrix(E, grid, t, 'left')

        def Sigma2_V(E):
            return contact_self_energy_matrix(E + V, grid, t, 'right')

        # Chemical potentials
        mu1 = 0.0   # Emitter (reference)
        mu2 = -V    # Collector (shifted by bias)

        # Current calculation via Landauer formula
        I = 0.0
        for iE, E in enumerate(E_array):
            try:
                S1 = Sigma1_V(E)
                S2 = Sigma2_V(E)
                Gamma1 = broadening_function(S1)
                Gamma2 = broadening_function(S2)

                H_eff = (E + 1j*eta) * np.eye(Np) - H - S1 - S2
                G = np.linalg.solve(H_eff, np.eye(Np))

                T_E = np.real(np.trace(Gamma1 @ G @ Gamma2 @ G.conj().T))
                T_E = max(0, min(T_E, 10))
            except:
                T_E = 0

            f1 = fermi(E, mu1)
            f2 = fermi(E, mu2)
            I += T_E * (f1 - f2) * dE

        # Convert to Amps (2e/h factor)
        I_array[iV] = (2 * q**2 / (2 * np.pi * hbar)) * I

    # Compute derivatives
    dIdV = np.gradient(I_array, V_array)
    d2IdV2 = np.gradient(dIdV, V_array)

    # Find peak and valley for NDR analysis
    I_peak_idx = np.argmax(I_array)
    I_peak = I_array[I_peak_idx]
    V_peak = V_array[I_peak_idx]

    has_ndr = np.any(dIdV < 0)
    if has_ndr and I_peak_idx < len(I_array) - 1:
        I_valley = np.min(I_array[I_peak_idx:])
        PVCR = I_peak / I_valley if I_valley > 0 else 1.0
        V_valley = V_array[I_peak_idx + np.argmin(I_array[I_peak_idx:])]
    else:
        I_valley = I_array[-1]
        PVCR = 1.0
        V_valley = V_array[-1]

    # Find resonance from transmission
    T_peaks = np.where(T_vs_E > 0.1 * np.max(T_vs_E))[0]
    if len(T_peaks) > 0:
        E_res = E_array[T_peaks[0]]
        T_res = T_vs_E[T_peaks[0]]
    else:
        E_res = 0
        T_res = 0

    return {
        'E_array': E_array,
        'T_vs_E': T_vs_E,
        'V_array': V_array,
        'I_array': I_array,
        'dIdV': dIdV,
        'd2IdV2': d2IdV2,
        'grid': grid,
        'I_peak': I_peak,
        'V_peak': V_peak,
        'I_valley': I_valley,
        'V_valley': V_valley,
        'has_ndr': has_ndr,
        'PVCR': PVCR,
        'E_res': E_res,
        'T_res': T_res
    }


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Run baseline simulations on oxide RTD devices',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--v-max', type=float, default=5.0,
                       help='Maximum bias voltage (V) - use 5V for oxide devices')
    parser.add_argument('--v-points', type=int, default=51,
                       help='Number of bias points')
    parser.add_argument('--e-points', type=int, default=300,
                       help='Number of energy points')
    parser.add_argument('--output-dir', '-o', type=str, default='results',
                       help='Output directory')
    parser.add_argument('--physics', action='store_true',
                       help='Print physics analysis and exit')
    parser.add_argument('--device', '-d', type=str, default=None,
                       help='Run specific device only')

    args = parser.parse_args()

    # Print physics analysis
    if args.physics:
        print_device_physics()
        return

    print("\n" + "="*70)
    print("OXIDE RTD BASELINE CHARACTERIZATION")
    print("="*70)
    print_device_physics()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Select devices to run
    if args.device:
        if args.device in OXIDE_DEVICES_OPTIMIZED:
            devices_to_run = {args.device: OXIDE_DEVICES_OPTIMIZED[args.device]}
        elif args.device in DEVICES:
            devices_to_run = {args.device: DEVICES[args.device]}
        else:
            print(f"Error: Device '{args.device}' not found")
            return
    else:
        devices_to_run = OXIDE_DEVICES_OPTIMIZED

    # Storage for results
    all_results = {}

    # Run simulations
    total_devices = len(devices_to_run)
    start_time = time.time()

    for idx, (name, device) in enumerate(devices_to_run.items()):
        print(f"\n[{idx+1}/{total_devices}] {name}")
        print(f"  {device['description']}")
        print("-" * 60)

        # Identify well and barrier materials
        well_mat = device['layers'][0]['material']
        barrier_mat = None
        for layer in device['layers']:
            mat_name = layer['material']
            if mat_name != well_mat:
                barrier_mat = mat_name
                break

        # Get and print band offset from literature
        if barrier_mat:
            try:
                cbo = get_band_offset(well_mat, barrier_mat)
                print(f"  Band offset: {well_mat}/{barrier_mat} CBO = {cbo:.2f} eV (from literature)")
            except (ValueError, KeyError):
                print(f"  Band offset: Using default MATERIALS values")

        # Print layer structure
        for i, layer in enumerate(device['layers']):
            thick_nm = layer['thickness'] * 1e9
            mat = layer['material']
            print(f"  Layer {i}: {mat:<8} {thick_nm:.1f} nm")

        # Run simulation
        sim_start = time.time()
        results = run_ballistic_simulation(
            device,
            E_points=args.e_points,
            V_max=args.v_max,
            V_points=args.v_points
        )
        sim_time = time.time() - sim_start

        # Store results
        all_results[name] = results

        print(f"\n  Results:")
        print(f"    Resonance: E = {results['E_res']*1000:.1f} meV, T = {results['T_res']:.3f}")
        print(f"    Peak: I = {results['I_peak']*1e6:.2f} μA at V = {results['V_peak']:.2f} V")
        print(f"    NDR: {'YES' if results['has_ndr'] else 'NO'}, PVCR = {results['PVCR']:.2f}")
        print(f"    Time: {sim_time:.1f} s")

    # Save all results
    total_time = time.time() - start_time

    print("\n" + "="*70)
    print(f"ALL SIMULATIONS COMPLETE in {total_time:.1f} s")
    print("="*70)

    # Save data
    save_dict = {'device_names': list(all_results.keys())}
    for name, res in all_results.items():
        key = name.replace('-', '_').replace(' ', '_')
        save_dict[f'{key}_E'] = res['E_array']
        save_dict[f'{key}_T'] = res['T_vs_E']
        save_dict[f'{key}_V'] = res['V_array']
        save_dict[f'{key}_I'] = res['I_array']
        save_dict[f'{key}_dIdV'] = res['dIdV']
        save_dict[f'{key}_d2IdV2'] = res['d2IdV2']

    npz_file = os.path.join(args.output_dir, 'oxide_baseline_data.npz')
    np.savez(npz_file, **save_dict)
    print(f"\nData saved to: {npz_file}")

    # Generate comparison plots
    generate_comparison_plots(all_results, args.output_dir)

    # Print summary table
    print_summary_table(all_results)


def generate_comparison_plots(all_results, output_dir):
    """Generate comparison plots for all devices"""

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    colors = plt.cm.tab10(np.linspace(0, 1, len(all_results)))

    # Panel (a): Transmission vs Energy (full range for oxide barriers)
    ax1 = axes[0, 0]
    for idx, (name, res) in enumerate(all_results.items()):
        short_name = name.replace('_thin_barrier', '').replace('_optimized', '').replace('_wide_well', '-wide')
        ax1.semilogy(res['E_array'], res['T_vs_E'] + 1e-30, '-', color=colors[idx],
                    linewidth=1.5, label=short_name[:20])
    ax1.set_xlabel('Energy (eV)')
    ax1.set_ylabel('Transmission')
    # Auto-scale x-axis based on data
    E_max_plot = max(res['E_array'][-1] for res in all_results.values())
    ax1.set_xlim([-0.5, min(E_max_plot, 4.0)])
    ax1.set_ylim([1e-12, 2])
    ax1.legend(loc='lower right', fontsize=7)
    ax1.set_title('(a) Transmission at V=0')
    ax1.grid(True, alpha=0.3)

    # Panel (b): I-V curves (in V for oxides, not mV)
    ax2 = axes[0, 1]
    for idx, (name, res) in enumerate(all_results.items()):
        short_name = name.replace('_thin_barrier', '').replace('_optimized', '').replace('_wide_well', '-wide')
        ax2.plot(res['V_array'], res['I_array']*1e6, '-', color=colors[idx],
                linewidth=1.5, label=short_name[:20])
    ax2.set_xlabel('Bias Voltage (V)')
    ax2.set_ylabel('Current (μA)')
    ax2.legend(loc='upper left', fontsize=7)
    ax2.set_title('(b) I-V Characteristics')
    ax2.grid(True, alpha=0.3)

    # Panel (c): dI/dV
    ax3 = axes[1, 0]
    for idx, (name, res) in enumerate(all_results.items()):
        short_name = name.replace('_thin_barrier', '').replace('_optimized', '').replace('_wide_well', '-wide')
        ax3.plot(res['V_array'], res['dIdV']*1e6, '-', color=colors[idx],
                linewidth=1.5, label=short_name[:20])
    ax3.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
    ax3.set_xlabel('Bias Voltage (V)')
    ax3.set_ylabel('dI/dV (μS)')
    ax3.legend(loc='upper right', fontsize=7)
    ax3.set_title('(c) Differential Conductance')
    ax3.grid(True, alpha=0.3)

    # Panel (d): PVCR comparison (bar chart)
    ax4 = axes[1, 1]
    names = list(all_results.keys())
    pvcrs = [res.get('PVCR', 1.0) for res in all_results.values()]
    short_names = [n.replace('_thin_barrier', '\n(thin)').replace('_optimized', '\n(opt)').replace('_wide_well', '\n(wide)')
                   for n in names]

    bar_colors = ['green' if p > 1.5 else 'orange' if p > 1.1 else 'red' for p in pvcrs]
    bars = ax4.bar(range(len(names)), pvcrs, color=bar_colors)
    ax4.set_xticks(range(len(names)))
    ax4.set_xticklabels(short_names, rotation=45, ha='right', fontsize=7)
    ax4.set_ylabel('Peak-to-Valley Current Ratio')
    ax4.set_title('(d) PVCR Comparison')
    ax4.axhline(y=1, color='k', linestyle='--', linewidth=0.5)

    # Add value labels
    for bar, val in zip(bars, pvcrs):
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2, height + 0.05,
                f'{val:.2f}', ha='center', va='bottom', fontsize=7)

    plt.tight_layout()

    # Save
    png_file = os.path.join(output_dir, 'oxide_baseline_comparison.png')
    pdf_file = os.path.join(output_dir, 'oxide_baseline_comparison.pdf')
    plt.savefig(png_file, dpi=300)
    plt.savefig(pdf_file)
    plt.close()

    print(f"Plots saved to: {png_file}")


def print_summary_table(all_results):
    """Print summary table of all results"""

    print("\n" + "="*100)
    print("SUMMARY TABLE")
    print("="*100)
    print(f"{'Device':<30} {'E_res(meV)':<12} {'I_peak(μA)':<12} {'V_peak(V)':<10} {'NDR':<6} {'PVCR':<8}")
    print("-"*100)

    for name, res in all_results.items():
        E_res = res.get('E_res', 0) * 1000  # Convert to meV
        I_peak = res.get('I_peak', np.max(res['I_array'])) * 1e6  # μA
        V_peak = res.get('V_peak', res['V_array'][np.argmax(res['I_array'])])
        has_ndr = 'Yes' if res.get('has_ndr', np.any(res['dIdV'] < 0)) else 'No'
        pvcr = res.get('PVCR', 1.0)

        print(f"{name:<30} {E_res:<12.1f} {I_peak:<12.3f} {V_peak:<10.2f} {has_ndr:<6} {pvcr:<8.2f}")

    print("="*100)


if __name__ == "__main__":
    main()
