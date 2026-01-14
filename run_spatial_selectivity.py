"""
Spatial Selectivity Study: IETS response vs molecule position

Tests how the molecular signature depends on where the molecule is placed:
1. Emitter barrier (highest sensitivity)
2. Quantum well center
3. Collector barrier (lower sensitivity)

Expected: Different positions give different IETS amplitudes,
demonstrating spatial sensitivity for targeted gas sensing.
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import sys
sys.path.insert(0, '.')

from config.device_library import get_device, get_material
from config.molecular_database import get_molecule
from core.hamiltonian import discretize_device, build_hamiltonian
from core.self_energy import contact_self_energy_matrix, local_projection_operator
from core.transverse_modes import TransverseModes
from core.scba_solver import scba_iteration_multimode, compute_current_multimode

print("="*70)
print("SPATIAL SELECTIVITY STUDY")
print("="*70)

# ============================================================================
# DEVICE SETUP
# ============================================================================

device_name = "GaAs_AlAs_symmetric"
device = get_device(device_name)
grid = discretize_device(device, grid_spacing=0.3e-9)
H, t = build_hamiltonian(grid)

print(f"\nDevice: {device_name}")
print(f"Grid: {grid['Np']} points, a = {grid['a']*1e9:.2f} nm")

# Show layer boundaries
print("\nLayer structure:")
for i, (start, end) in enumerate(grid['layer_boundaries']):
    mat = grid['material'][start]
    thickness = (end - start + 1) * grid['a']
    print(f"  Layer {i}: sites {start:2d}-{end:2d}, {mat:6s}, {thickness*1e9:.1f} nm")

# ============================================================================
# MOLECULE POSITIONS
# ============================================================================

# Get layer boundaries
layers = grid['layer_boundaries']
# Layer 0: Emitter contact
# Layer 1: Emitter barrier (AlAs)
# Layer 2: Quantum well (GaAs)
# Layer 3: Collector barrier (AlAs)
# Layer 4: Collector contact

positions = [
    {
        'name': 'Emitter Barrier',
        'layer_idx': 1,
        'description': 'On emitter side barrier - highest sensitivity expected'
    },
    {
        'name': 'Quantum Well',
        'layer_idx': 2,
        'description': 'In quantum well center - moderate sensitivity'
    },
    {
        'name': 'Collector Barrier',
        'layer_idx': 3,
        'description': 'On collector side barrier - lower sensitivity expected'
    }
]

# Determine molecule sites for each position
for pos in positions:
    start, end = layers[pos['layer_idx']]
    pos['local_sites'] = [(start + end) // 2]  # Center of layer
    pos['neighbor_radius'] = 2
    print(f"\n{pos['name']}:")
    print(f"  Site: {pos['local_sites'][0]}, Layer: {pos['layer_idx']}")
    print(f"  {pos['description']}")

# ============================================================================
# SIMULATION PARAMETERS
# ============================================================================

molecule_name = "Benzene"
molecule = get_molecule(molecule_name)
print(f"\nMolecule: {molecule_name}")
print(f"  Modes: {molecule['modes_meV'][:3]} meV (first 3)")

# Transverse modes (use smaller grid for speed)
m0 = 9.10938356e-31
mat = get_material(grid['material'][0])
m_trans = mat['m_eff'] * m0

transverse_modes = TransverseModes(
    Ly=1.0e-6, Lz=1.0e-6,
    n_max=3, m_max=3,
    m_trans=m_trans
)
transverse_modes.compute_modes()
print(f"Transverse modes: {transverse_modes.Nm}")

# Energy grid
E_array = np.linspace(-0.1, 0.8, 30)
print(f"Energy grid: {len(E_array)} points")

# Bias sweep (5 points for quick test)
V_array = np.linspace(0.0, 0.2, 5)
print(f"Bias sweep: {len(V_array)} points ({V_array[0]:.2f} to {V_array[-1]:.2f} V)")

temperature = 300

# Contact self-energies
def Sigma1_func(E):
    return contact_self_energy_matrix(E, grid, t, contact='left')

def Sigma2_func(E):
    return contact_self_energy_matrix(E, grid, t, contact='right')

# ============================================================================
# RUN SIMULATIONS FOR EACH POSITION
# ============================================================================

results_by_position = {}

for pos in positions:
    print("\n" + "="*70)
    print(f"SIMULATING: {pos['name']}")
    print("="*70)

    # Build phonon modes for this position
    phonon_modes = []

    # Bulk phonon (position-independent)
    phonon_modes.append({
        'energy': 0.036,
        'coupling': 0.010,
        'is_local': False
    })

    # Molecular modes (position-dependent)
    local_sites = pos['local_sites']
    neighbor_radius = pos['neighbor_radius']

    mol_energies = molecule['modes_meV']
    mol_couplings = molecule['coupling_meV']

    for i in range(min(2, len(mol_energies))):  # Use 2 modes for speed
        energy_meV = mol_energies[i]
        coupling_meV = mol_couplings[i]

        D_base = coupling_meV / 1000.0  # meV to eV
        D_nm = np.zeros(transverse_modes.Nm)
        for im in range(transverse_modes.Nm):
            n, m = transverse_modes.modes[im]
            y_mol = 1.0e-6 / 2
            z_mol = 1.0e-6 / 2
            psi_nm = transverse_modes.get_wavefunction(n, m, y_mol, z_mol)
            psi_11 = transverse_modes.get_wavefunction(1, 1, y_mol, z_mol)
            D_nm[im] = D_base * (psi_nm**2 / psi_11**2)

        phonon_modes.append({
            'energy': energy_meV / 1000.0,
            'coupling': D_nm,
            'is_local': True,
            'local_sites': local_sites,
            'neighbor_radius': neighbor_radius
        })

    # Run bias sweep
    I_array = np.zeros(len(V_array))
    scba_results = []

    start_time = time.time()

    for iV, V in enumerate(V_array):
        mu1 = +V / 2.0
        mu2 = -V / 2.0

        print(f"  V = {V:.3f} V ... ", end='', flush=True)

        result = scba_iteration_multimode(
            E_array, H, Sigma1_func, Sigma2_func,
            mu1, mu2, temperature,
            phonon_modes, grid,
            transverse_modes,
            max_iter=20,
            tol=1e-2,
            mix=0.3,
            verbose=False
        )

        I, I_vs_E, I_vs_mode = compute_current_multimode(
            result, E_array, mu1, mu2, temperature
        )

        I_array[iV] = I
        scba_results.append({
            'converged': result['converged'],
            'iterations': result['iterations']
        })

        print(f"I = {I*1e9:.4f} nA, converged = {result['converged']}")

    elapsed = time.time() - start_time

    # Compute derivatives
    dIdV = np.gradient(I_array, V_array)
    d2IdV2 = np.gradient(dIdV, V_array)

    results_by_position[pos['name']] = {
        'V_array': V_array,
        'I_array': I_array,
        'dIdV': dIdV,
        'd2IdV2': d2IdV2,
        'scba_results': scba_results,
        'elapsed': elapsed,
        'position': pos
    }

    print(f"\n  Completed in {elapsed:.1f} seconds")
    print(f"  Current at V=0.2V: {I_array[-1]*1e9:.4f} nA")

# ============================================================================
# ANALYSIS
# ============================================================================

print("\n" + "="*70)
print("COMPARISON")
print("="*70)

print(f"\n{'Position':<20} {'I @ 0.2V (nA)':<15} {'Relative':<12} {'Time (s)':<10}")
print("-"*60)

I_ref = results_by_position['Emitter Barrier']['I_array'][-1]
for pos_name, data in results_by_position.items():
    I_final = data['I_array'][-1]
    relative = I_final / I_ref if I_ref != 0 else 0
    print(f"{pos_name:<20} {I_final*1e9:<15.4f} {relative:<12.3f} {data['elapsed']:<10.1f}")

# ============================================================================
# PLOT RESULTS
# ============================================================================

print("\n" + "="*70)
print("CREATING PLOTS")
print("="*70)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

colors = {'Emitter Barrier': 'blue', 'Quantum Well': 'green', 'Collector Barrier': 'red'}

# Plot 1: I-V curves
ax1 = axes[0, 0]
for pos_name, data in results_by_position.items():
    ax1.plot(data['V_array'], data['I_array']*1e9, 'o-',
             label=pos_name, color=colors[pos_name], linewidth=2, markersize=6)
ax1.set_xlabel('Bias Voltage (V)', fontsize=12)
ax1.set_ylabel('Current (nA)', fontsize=12)
ax1.set_title('I-V Curves: Position Dependence', fontsize=13, fontweight='bold')
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: dI/dV
ax2 = axes[0, 1]
for pos_name, data in results_by_position.items():
    ax2.plot(data['V_array'], data['dIdV']*1e6, 'o-',
             label=pos_name, color=colors[pos_name], linewidth=2, markersize=6)
ax2.set_xlabel('Bias Voltage (V)', fontsize=12)
ax2.set_ylabel('dI/dV (µS)', fontsize=12)
ax2.set_title('Differential Conductance', fontsize=13, fontweight='bold')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Plot 3: d²I/dV² (IETS)
ax3 = axes[1, 0]
for pos_name, data in results_by_position.items():
    ax3.plot(data['V_array'], data['d2IdV2']*1e6, 'o-',
             label=pos_name, color=colors[pos_name], linewidth=2, markersize=6)
ax3.set_xlabel('Bias Voltage (V)', fontsize=12)
ax3.set_ylabel('d²I/dV² (µS/V)', fontsize=12)
ax3.set_title('IETS Signal: Position Dependence', fontsize=13, fontweight='bold')
ax3.legend()
ax3.grid(True, alpha=0.3)
ax3.axhline(y=0, color='k', linestyle='--', alpha=0.3)

# Plot 4: Relative sensitivity bar chart
ax4 = axes[1, 1]
pos_names = list(results_by_position.keys())
I_values = [results_by_position[n]['I_array'][-1]*1e9 for n in pos_names]
bar_colors = [colors[n] for n in pos_names]

bars = ax4.bar(pos_names, I_values, color=bar_colors, edgecolor='black')
ax4.set_ylabel('Current at 0.2V (nA)', fontsize=12)
ax4.set_title('Position Sensitivity Comparison', fontsize=13, fontweight='bold')
ax4.tick_params(axis='x', rotation=15)

# Add relative values on bars
for bar, val in zip(bars, I_values):
    ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
             f'{val:.2f}', ha='center', fontsize=10)

plt.suptitle(f'Spatial Selectivity Study: {molecule_name} on {device_name}',
             fontsize=14, fontweight='bold', y=0.995)
plt.tight_layout(rect=[0, 0, 1, 0.98])

output_file = 'spatial_selectivity_results.png'
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"\nPlots saved to: {output_file}")

# Save data
np.savez('spatial_selectivity_data.npz',
         **{f"{n.replace(' ', '_')}_V": results_by_position[n]['V_array'] for n in pos_names},
         **{f"{n.replace(' ', '_')}_I": results_by_position[n]['I_array'] for n in pos_names},
         **{f"{n.replace(' ', '_')}_dIdV": results_by_position[n]['dIdV'] for n in pos_names},
         **{f"{n.replace(' ', '_')}_d2IdV2": results_by_position[n]['d2IdV2'] for n in pos_names})
print(f"Data saved to: spatial_selectivity_data.npz")

print("\n" + "="*70)
print("SPATIAL SELECTIVITY STUDY COMPLETE")
print("="*70)
