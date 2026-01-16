"""
RTD Baseline Characterization with Proper Bias Treatment

Run baseline (no molecule) simulations on RTD devices with:
- Bias-induced band bending (essential for NDR)
- Proper chemical potential handling
- Transmission analysis at multiple bias points
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

from config.device_library import DEVICES, get_device, get_band_offset, BAND_OFFSETS
from core.hamiltonian import discretize_device, build_hamiltonian
from core.self_energy import contact_self_energy_matrix
from core.green_functions import broadening_function

# Physical constants
hbar = 1.054571817e-34
m0 = 9.10938356e-31
q = 1.602176634e-19
kB = 1.380649e-23


def run_rtd_baseline(device, grid_spacing=0.1e-9, E_points=400,
                     V_min=0.0, V_max=0.5, V_points=51,
                     temperature=300, verbose=True):
    """
    Run baseline RTD simulation with proper bias treatment.

    Includes bias-induced band bending for correct NDR behavior.
    """

    # Discretize device
    grid = discretize_device(device, grid_spacing=grid_spacing)
    Np = grid['Np']
    x = grid['x']
    L = x[-1]

    # Base Hamiltonian (at V=0)
    H0, t = build_hamiltonian(grid)

    if verbose:
        print(f"  Grid: {Np} points, {L*1e9:.1f} nm")

    # Energy grid - cover below and above expected resonance
    E_min = -0.2
    E_max = max(0.8, np.max(grid['Ec']) + 0.3)  # Cover barrier + margin
    E_array = np.linspace(E_min, E_max, E_points)
    dE = E_array[1] - E_array[0]

    # Thermal energy
    kT = kB * temperature / q

    def fermi(E, mu):
        x_val = (E - mu) / kT
        x_val = np.clip(x_val, -500, 500)
        return 1.0 / (1.0 + np.exp(x_val))

    # Bias sweep
    V_array = np.linspace(V_min, V_max, V_points)
    I_array = np.zeros(V_points)
    eta = 1e-5

    # Store transmission at V=0 for analysis
    T_at_V0 = np.zeros(E_points)

    for iV, V in enumerate(V_array):
        if verbose and iV % 10 == 0:
            print(f"  V = {V*1000:.0f} mV...")

        # Linear potential drop across device (bias-induced band bending)
        U_bias = -V * (x / L)

        # Build Hamiltonian with bias
        H = H0.copy()
        for i in range(Np):
            H[i, i] += U_bias[i]

        # Contact self-energies
        # Left contact at reference potential
        # Right contact shifted down by V
        def Sigma1(E):
            return contact_self_energy_matrix(E, grid, t, 'left')

        def Sigma2(E):
            return contact_self_energy_matrix(E + V, grid, t, 'right')

        # Chemical potentials
        mu1 = 0.0   # Emitter (reference)
        mu2 = -V    # Collector (shifted by bias)

        # Current calculation
        I = 0.0
        for iE, E in enumerate(E_array):
            S1 = Sigma1(E)
            S2 = Sigma2(E)
            Gamma1 = broadening_function(S1)
            Gamma2 = broadening_function(S2)

            H_eff = (E + 1j*eta) * np.eye(Np) - H - S1 - S2
            try:
                G = np.linalg.solve(H_eff, np.eye(Np))
                T_E = np.real(np.trace(Gamma1 @ G @ Gamma2 @ G.conj().T))
                T_E = max(0, min(T_E, 10))
            except:
                T_E = 0

            # Store transmission at V=0
            if iV == 0:
                T_at_V0[iE] = T_E

            f1 = fermi(E, mu1)
            f2 = fermi(E, mu2)
            I += T_E * (f1 - f2) * dE

        I_array[iV] = (2 * q**2 / (2 * np.pi * hbar)) * I

    # Compute derivatives
    dIdV = np.gradient(I_array, V_array)
    d2IdV2 = np.gradient(dIdV, V_array)

    # Find peak and valley
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
    T_peaks = np.where(T_at_V0 > 0.1 * np.max(T_at_V0))[0]
    if len(T_peaks) > 0:
        E_res = E_array[T_peaks[0]]
        T_res = T_at_V0[T_peaks[0]]
    else:
        E_res = 0
        T_res = 0

    return {
        'E_array': E_array,
        'T_vs_E': T_at_V0,
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


def main():
    parser = argparse.ArgumentParser(
        description='Run baseline RTD characterization with proper bias treatment',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--devices', '-d', type=str, default='GaAs',
                       help='Device filter: "GaAs", "oxide", "all", or specific device name')
    parser.add_argument('--v-max', type=float, default=0.5,
                       help='Maximum bias voltage (V)')
    parser.add_argument('--v-points', type=int, default=51,
                       help='Number of bias points')
    parser.add_argument('--e-points', type=int, default=400,
                       help='Number of energy points')
    parser.add_argument('--output-dir', '-o', type=str, default='results',
                       help='Output directory')

    args = parser.parse_args()

    # Select devices
    if args.devices == 'GaAs':
        device_names = [n for n in DEVICES.keys() if 'GaAs' in n or 'InGaAs' in n]
    elif args.devices == 'oxide':
        device_names = [n for n in DEVICES.keys() if any(x in n for x in ['In2O3', 'IGZO', 'ZnO', 'SnO2'])]
    elif args.devices == 'all':
        device_names = list(DEVICES.keys())
    elif args.devices in DEVICES:
        device_names = [args.devices]
    else:
        print(f"Unknown device selection: {args.devices}")
        print(f"Available: {list(DEVICES.keys())}")
        return

    print("\n" + "="*70)
    print("RTD BASELINE CHARACTERIZATION")
    print("="*70)
    print(f"Devices: {len(device_names)}")
    print(f"Bias range: 0 to {args.v_max*1000:.0f} mV ({args.v_points} points)")

    os.makedirs(args.output_dir, exist_ok=True)

    results = {}
    start_time = time.time()

    for idx, name in enumerate(device_names):
        print(f"\n[{idx+1}/{len(device_names)}] {name}")
        print("-" * 60)

        device = get_device(name)

        # Identify well and barrier materials
        well_mat = device['layers'][0]['material']
        barrier_mat = None
        for layer in device['layers']:
            mat_name = layer['material']
            if mat_name != well_mat:
                barrier_mat = mat_name
                break

        # Get and print band offset
        if barrier_mat:
            try:
                cbo = get_band_offset(well_mat, barrier_mat)
                print(f"  Band offset: {well_mat}/{barrier_mat} CBO = {cbo:.2f} eV (from literature)")
            except (ValueError, KeyError):
                print(f"  Band offset: Using default MATERIALS values")

        # Print structure
        for i, layer in enumerate(device['layers']):
            print(f"  Layer {i}: {layer['material']:<8} {layer['thickness']*1e9:.1f} nm")

        sim_start = time.time()
        res = run_rtd_baseline(
            device,
            E_points=args.e_points,
            V_max=args.v_max,
            V_points=args.v_points
        )
        sim_time = time.time() - sim_start

        results[name] = res

        print(f"\n  Results:")
        print(f"    Resonance: E = {res['E_res']*1000:.1f} meV, T = {res['T_res']:.3f}")
        print(f"    Peak: I = {res['I_peak']*1e6:.2f} uA at V = {res['V_peak']*1000:.0f} mV")
        print(f"    NDR: {'YES' if res['has_ndr'] else 'NO'}, PVCR = {res['PVCR']:.2f}")
        print(f"    Time: {sim_time:.1f} s")

    total_time = time.time() - start_time

    # Summary
    print("\n" + "="*80)
    print("SUMMARY TABLE")
    print("="*80)
    print(f"{'Device':<30} {'E_res(meV)':<12} {'I_peak(uA)':<12} {'V_peak(mV)':<12} {'NDR':<6} {'PVCR':<8}")
    print("-"*80)
    for name, res in results.items():
        print(f"{name:<30} {res['E_res']*1000:<12.1f} {res['I_peak']*1e6:<12.2f} {res['V_peak']*1000:<12.0f} {'Yes' if res['has_ndr'] else 'No':<6} {res['PVCR']:<8.2f}")
    print("="*80)
    print(f"Total time: {total_time:.1f} s")

    # Save data
    save_dict = {'device_names': list(results.keys())}
    for name, res in results.items():
        key = name.replace('-', '_').replace(' ', '_')
        save_dict[f'{key}_V'] = res['V_array']
        save_dict[f'{key}_I'] = res['I_array']
        save_dict[f'{key}_E'] = res['E_array']
        save_dict[f'{key}_T'] = res['T_vs_E']

    npz_file = os.path.join(args.output_dir, f'rtd_baseline_{args.devices}.npz')
    np.savez(npz_file, **save_dict)
    print(f"\nData saved to: {npz_file}")

    # Generate plots
    generate_plots(results, args.output_dir, args.devices)


def generate_plots(results, output_dir, prefix):
    """Generate comparison plots"""

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    colors = plt.cm.tab10(np.linspace(0, 1, len(results)))

    # I-V curves
    ax1 = axes[0, 0]
    for idx, (name, res) in enumerate(results.items()):
        short = name.replace('GaAs_AlAs_', '').replace('InGaAs_InAlAs_', 'InGaAs_')[:15]
        ax1.plot(res['V_array']*1000, res['I_array']*1e6, '-', color=colors[idx], lw=1.5, label=short)
    ax1.set_xlabel('Bias Voltage (mV)')
    ax1.set_ylabel('Current (μA)')
    ax1.legend(fontsize=7)
    ax1.set_title('(a) I-V Characteristics')
    ax1.grid(True, alpha=0.3)

    # Transmission
    ax2 = axes[0, 1]
    for idx, (name, res) in enumerate(results.items()):
        short = name.replace('GaAs_AlAs_', '').replace('InGaAs_InAlAs_', 'InGaAs_')[:15]
        ax2.semilogy(res['E_array']*1000, res['T_vs_E'], '-', color=colors[idx], lw=1.5, label=short)
    ax2.set_xlabel('Energy (meV)')
    ax2.set_ylabel('Transmission')
    ax2.set_xlim([-100, 500])
    ax2.set_ylim([1e-8, 2])
    ax2.legend(fontsize=7)
    ax2.set_title('(b) Transmission at V=0')
    ax2.grid(True, alpha=0.3)

    # dI/dV
    ax3 = axes[1, 0]
    for idx, (name, res) in enumerate(results.items()):
        short = name.replace('GaAs_AlAs_', '').replace('InGaAs_InAlAs_', 'InGaAs_')[:15]
        ax3.plot(res['V_array']*1000, res['dIdV']*1e3, '-', color=colors[idx], lw=1.5, label=short)
    ax3.axhline(y=0, color='k', ls='-', lw=0.5)
    ax3.set_xlabel('Bias Voltage (mV)')
    ax3.set_ylabel('dI/dV (mS)')
    ax3.legend(fontsize=7)
    ax3.set_title('(c) Differential Conductance')
    ax3.grid(True, alpha=0.3)

    # Bar chart of PVCR
    ax4 = axes[1, 1]
    names = list(results.keys())
    pvcrs = [res['PVCR'] for res in results.values()]
    short_names = [n.replace('GaAs_AlAs_', '').replace('InGaAs_InAlAs_', 'InGaAs_')[:12] for n in names]
    bar_colors = ['green' if p > 1.5 else 'red' for p in pvcrs]
    bars = ax4.bar(range(len(names)), pvcrs, color=bar_colors)
    ax4.set_xticks(range(len(names)))
    ax4.set_xticklabels(short_names, rotation=30, ha='right', fontsize=8)
    ax4.set_ylabel('Peak-to-Valley Current Ratio')
    ax4.set_title('(d) PVCR Comparison')
    ax4.axhline(y=1, color='k', ls='--', lw=0.5)
    for bar, val in zip(bars, pvcrs):
        ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                f'{val:.1f}', ha='center', va='bottom', fontsize=8)

    plt.tight_layout()
    png_file = os.path.join(output_dir, f'rtd_baseline_{prefix}.png')
    pdf_file = os.path.join(output_dir, f'rtd_baseline_{prefix}.pdf')
    plt.savefig(png_file, dpi=300)
    plt.savefig(pdf_file)
    plt.close()
    print(f"Plots saved to: {png_file}")


if __name__ == "__main__":
    main()
