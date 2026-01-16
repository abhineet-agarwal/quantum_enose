"""
Single Molecule IETS Simulation

Complete pipeline for simulating one molecule on an RTD device.
Computes full I-V curve and IETS spectrum.
"""

import numpy as np
import sys
import os
import time

# Setup path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

from config.device_library import get_device
from config.molecular_database import get_molecule
from core.hamiltonian import discretize_device, build_hamiltonian
from core.self_energy import (
    contact_self_energy_matrix, 
    get_molecule_location
)
from core.green_functions import broadening_function
from core.scba_solver import (
    scba_iteration, compute_current, compute_iets,
    scba_iteration_multimode, compute_current_multimode
)

# ============================================================================
# SIMULATION CONFIGURATION
# ============================================================================

class SimulationConfig:
    """Configuration for IETS simulation"""
    
    def __init__(self):
        # Device
        self.device_name = "In2O3_Al2O3_symmetric"
        self.grid_spacing = 0.12e-9  # 0.12 nm
        
        # Molecule
        self.molecule_name = "Benzene"
        
        # Material phonon (bulk)
        self.bulk_phonon_energy = 0.070  # 70 meV (In2O3)
        self.bulk_phonon_coupling = 0.010  # 10 meV
        
        # Molecular phonon coupling
        self.molecular_coupling_scale = 1.0  # Multiplier for molecule couplings
        
        # Energy grid
        self.E_min = -0.3  # eV
        self.E_max = 1.5   # eV
        self.E_points = 200
        
        # Bias sweep
        self.V_min = 0.0   # V
        self.V_max = 0.5   # V
        self.V_points = 26  # 0.02 V steps
        
        # Temperature
        self.temperature = 300  # K
        
        # SCBA parameters
        self.scba_max_iter = 50
        self.scba_tolerance = 1e-4
        self.scba_mixing = 0.3

        # Transverse mode parameters (Quasi-3D)
        self.use_multimode = False      # Enable quasi-3D mode summation
        self.Ly = 1.0e-6               # Transverse width y (m) - 1 µm default
        self.Lz = 1.0e-6               # Transverse width z (m) - 1 µm default
        self.n_max = 3                 # Max mode index in y direction
        self.m_max = 3                 # Max mode index in z direction
        self.m_trans = None            # Effective mass for transverse motion (defaults to device m*)

        # Hybrid mode selection (for speedup)
        self.use_hybrid = False        # Enable hybrid mode selection (7-8x speedup)
        self.n_inelastic = 4           # Number of modes for full SCBA (rest are coherent)

        # Output
        self.output_dir = "results"
        self.verbose = True

# ============================================================================
# PHONON MODE BUILDER
# ============================================================================

def build_phonon_modes(molecule, config, grid, device_config, molecule_name="unknown",
                       transverse_modes=None):
    """
    Build phonon mode list from molecule data

    Parameters:
    -----------
    molecule : dict
        Molecule data from database
    config : SimulationConfig
        Simulation configuration
    grid : dict
        Device grid
    device_config : dict
        Device configuration
    molecule_name : str
        Molecule name (for labeling)
    transverse_modes : TransverseModes or None
        If provided, compute mode-dependent molecular coupling

    Returns:
    --------
    phonon_modes : list of dict
        List of phonon modes for SCBA
    """

    phonon_modes = []

    # Add bulk phonon (global, material property)
    # Bulk phonons are mode-independent
    if config.bulk_phonon_energy > 0:
        phonon_modes.append({
            'energy': config.bulk_phonon_energy,
            'coupling': config.bulk_phonon_coupling,
            'is_local': False,
            'name': 'Bulk phonon'
        })

    # Add molecular vibrations (local)
    if molecule['modes_meV']:
        # Get molecule location (longitudinal)
        mol_sites, mol_radius = get_molecule_location(grid, device_config)

        # Molecule transverse position
        # For now, assume molecule at center of cross-section
        # TODO: Make this configurable via config.y_mol, config.z_mol
        if transverse_modes is not None:
            y_mol = transverse_modes.Ly / 2.0  # Center
            z_mol = transverse_modes.Lz / 2.0  # Center

        for i, (energy_meV, coupling_meV) in enumerate(
            zip(molecule['modes_meV'], molecule['coupling_meV'])
        ):
            # Base coupling strength
            D_base = (coupling_meV / 1000.0) * config.molecular_coupling_scale  # eV

            # Mode-dependent coupling for quasi-3D
            if transverse_modes is not None:
                # Compute |ψ_nm(y_mol, z_mol)|² for all modes
                psi_sq = np.zeros(transverse_modes.Nm)

                for im in range(transverse_modes.Nm):
                    n, m = transverse_modes.modes[im]
                    psi = transverse_modes.get_wavefunction(n, m, y_mol, z_mol)
                    psi_sq[im] = psi**2

                # Normalize to get relative weights [0, 1]
                # This preserves physics while keeping magnitudes reasonable
                max_psi_sq = np.max(psi_sq)
                if max_psi_sq > 0:
                    weights = psi_sq / max_psi_sq
                else:
                    weights = np.zeros(transverse_modes.Nm)

                # Apply normalized weights: D_nm = D_base × weight
                # Modes with nodes → weight=0 → no coupling
                # Mode with antinode → weight=1 → full coupling
                D_nm = D_base * weights

                coupling = D_nm  # Array[Nm]
            else:
                # 1D: scalar coupling
                coupling = D_base

            phonon_modes.append({
                'energy': energy_meV / 1000.0,  # Convert meV to eV
                'coupling': coupling,  # Scalar (1D) or Array[Nm] (quasi-3D)
                'is_local': True,
                'local_sites': mol_sites,
                'neighbor_radius': mol_radius,
                'name': f"{molecule_name}_mode_{i+1}"
            })

    return phonon_modes

# ============================================================================
# TRANSVERSE MODE SETUP (Quasi-3D)
# ============================================================================

def setup_transverse_modes(config, device):
    """
    Create TransverseModes object from configuration

    Parameters:
    -----------
    config : SimulationConfig
        Contains Ly, Lz, n_max, m_max, m_trans
    device : dict
        Device specification with layers

    Returns:
    --------
    transverse_modes : TransverseModes object
        Configured and computed mode object
    """
    from core.transverse_modes import TransverseModes
    from config.device_library import MATERIALS

    # Default to device effective mass if not specified
    m_trans = config.m_trans
    if m_trans is None:
        # Extract effective mass from first well/channel layer
        # Look for the material with lowest band edge (well material)
        m0 = 9.10938356e-31  # kg

        # Try to find well material (typically the first layer with Ec = 0 or lowest)
        if 'layers' in device:
            well_material = None
            min_Ec = float('inf')

            for layer in device['layers']:
                mat_name = layer['material']
                if mat_name in MATERIALS:
                    mat = MATERIALS[mat_name]
                    if mat['Ec'] < min_Ec:
                        min_Ec = mat['Ec']
                        well_material = mat_name

            if well_material:
                m_eff_rel = MATERIALS[well_material]['m_eff']
            else:
                # Fallback: use first layer
                m_eff_rel = MATERIALS[device['layers'][0]['material']]['m_eff']
        else:
            # No layers info, use GaAs default
            m_eff_rel = 0.067

        m_trans = m_eff_rel * m0

    # Create and compute modes
    modes = TransverseModes(
        config.Ly, config.Lz,
        config.n_max, config.m_max,
        m_trans
    )
    modes.compute_modes()

    return modes

# ============================================================================
# MAIN SIMULATION FUNCTION
# ============================================================================

def run_iets_simulation(config):
    """
    Run complete IETS simulation for one molecule
    
    Parameters:
    -----------
    config : SimulationConfig
        Simulation configuration
        
    Returns:
    --------
    results : dict
        Complete simulation results
    """
    
    if config.verbose:
        print("\n" + "="*70)
        print(f"IETS SIMULATION: {config.molecule_name} on {config.device_name}")
        print("="*70)
    
    start_time = time.time()
    
    # ========================================================================
    # SETUP
    # ========================================================================
    
    if config.verbose:
        print("\n[1/5] Setting up device and molecule...")
    
    # Load device
    device = get_device(config.device_name)
    grid = discretize_device(device, grid_spacing=config.grid_spacing)
    H, t = build_hamiltonian(grid)

    # Load molecule
    molecule = get_molecule(config.molecule_name)

    # Setup transverse modes FIRST (needed for mode-dependent coupling)
    if config.use_multimode:
        transverse_modes = setup_transverse_modes(config, device)
        if config.verbose:
            print(f"  Transverse modes: {transverse_modes.Nm} total ({config.n_max}×{config.m_max})")
            print(f"    Mode energies: {transverse_modes.energies[0]*1000:.3f} to {transverse_modes.energies[-1]*1000:.3f} meV")
            print(f"    Dimensions: Ly={config.Ly*1e6:.2f} µm, Lz={config.Lz*1e6:.2f} µm")
    else:
        transverse_modes = None
        if config.verbose:
            print(f"  Mode: 1D transport (no transverse modes)")

    # Build phonon modes (with mode-dependent molecular coupling if quasi-3D)
    phonon_modes = build_phonon_modes(molecule, config, grid, device, config.molecule_name,
                                      transverse_modes=transverse_modes)

    if config.verbose:
        print(f"  Device: {grid['Np']} grid points, {grid['x'][-1]*1e9:.1f} nm")
        print(f"  Molecule: {config.molecule_name} ({len(molecule['modes_meV'])} vibrational modes)")
        print(f"  Total phonon modes: {len(phonon_modes)} (1 bulk + {len(molecule['modes_meV'])} molecular)")

        # Print mode-dependent coupling info for quasi-3D
        if config.use_multimode and len(molecule['modes_meV']) > 0:
            print(f"  Mode-dependent molecular coupling:")
            # Show coupling for first molecular mode
            mol_mode = phonon_modes[1] if len(phonon_modes) > 1 else None
            if mol_mode and hasattr(mol_mode['coupling'], '__len__'):
                D_array = mol_mode['coupling']
                print(f"    First vibrational mode: D_min={np.min(D_array)*1000:.3f} meV, D_max={np.max(D_array)*1000:.3f} meV")
                # Count how many modes have significant coupling (>1% of max)
                n_coupled = np.sum(D_array > 0.01 * np.max(D_array))
                print(f"    Modes with >1% coupling: {n_coupled}/{transverse_modes.Nm}")

    # ========================================================================
    # ENERGY GRID
    # ========================================================================
    
    if config.verbose:
        print("\n[2/5] Creating energy grid...")
    
    E_array = np.linspace(config.E_min, config.E_max, config.E_points)
    dE = E_array[1] - E_array[0]
    
    if config.verbose:
        print(f"  Energy range: {config.E_min:.2f} to {config.E_max:.2f} eV")
        print(f"  Points: {config.E_points} (dE = {dE*1000:.1f} meV)")
    
    # ========================================================================
    # CONTACT SELF-ENERGIES
    # ========================================================================
    
    if config.verbose:
        print("\n[3/5] Computing contact self-energies...")
    
    def Sigma1(E):
        return contact_self_energy_matrix(E, grid, t, 'left')
    
    def Sigma2(E):
        return contact_self_energy_matrix(E, grid, t, 'right')
    
    # Pre-compute Gamma arrays for current calculation
    Np = grid['Np']
    Gamma1_array = np.zeros((Np, Np, config.E_points))
    Gamma2_array = np.zeros((Np, Np, config.E_points))
    
    for iE, E in enumerate(E_array):
        Gamma1_array[:, :, iE] = broadening_function(Sigma1(E))
        Gamma2_array[:, :, iE] = broadening_function(Sigma2(E))
    
    if config.verbose:
        print(f"  Contact broadenings computed for {config.E_points} energies")
    
    # ========================================================================
    # BIAS SWEEP (I-V CURVE)
    # ========================================================================
    
    if config.verbose:
        print("\n[4/5] Running bias sweep...")
        print(f"  Voltage: {config.V_min:.2f} to {config.V_max:.2f} V ({config.V_points} points)")
    
    V_array = np.linspace(config.V_min, config.V_max, config.V_points)
    I_array = np.zeros(config.V_points)
    scba_results = []
    
    for iV, V in enumerate(V_array):
        if config.verbose:
            print(f"\n  --- Bias {iV+1}/{config.V_points}: V = {V:.3f} V ---")
        
        # Chemical potentials (symmetric bias)
        mu1 = +V / 2.0
        mu2 = -V / 2.0
        
        # Run SCBA (choose between 1D, quasi-3D, or hybrid)
        if config.use_multimode:
            # Quasi-3D with optional hybrid mode selection
            if config.use_hybrid:
                # Hybrid: select important modes for SCBA, rest are coherent
                from core.hybrid_modes import select_hybrid_modes, print_mode_selection
                from core.scba_solver_hybrid import scba_iteration_hybrid, compute_current_hybrid

                # Select modes (only once at first bias point)
                if iV == 0:
                    inelastic_modes, coherent_modes, importance_scores = select_hybrid_modes(
                        transverse_modes, phonon_modes, config.temperature,
                        n_inelastic=config.n_inelastic
                    )
                    if config.verbose:
                        print_mode_selection(transverse_modes, inelastic_modes,
                                           coherent_modes, importance_scores)

                result = scba_iteration_hybrid(
                    E_array, H, Sigma1, Sigma2,
                    mu1, mu2, config.temperature,
                    phonon_modes, grid,
                    transverse_modes,
                    inelastic_modes, coherent_modes,
                    max_iter=config.scba_max_iter,
                    tol=config.scba_tolerance,
                    mix=config.scba_mixing,
                    verbose=config.verbose
                )

                # Compute current
                I, I_vs_E, I_vs_mode = compute_current_hybrid(
                    result, E_array, mu1, mu2, config.temperature
                )
            else:
                # All-inelastic: standard multimode SCBA
                result = scba_iteration_multimode(
                    E_array, H, Sigma1, Sigma2,
                    mu1, mu2, config.temperature,
                    phonon_modes, grid,
                    transverse_modes,
                    max_iter=config.scba_max_iter,
                    tol=config.scba_tolerance,
                    mix=config.scba_mixing,
                    verbose=config.verbose
                )

                # Compute current with mode summation
                I, I_vs_E, I_vs_mode = compute_current_multimode(
                    result, E_array, mu1, mu2, config.temperature
                )
        else:
            # 1D: standard SCBA
            result = scba_iteration(
                E_array, H, Sigma1, Sigma2,
                mu1, mu2, config.temperature,
                phonon_modes, grid,
                max_iter=config.scba_max_iter,
                tol=config.scba_tolerance,
                mix=config.scba_mixing,
                verbose=config.verbose
            )

            # Compute current (1D)
            I, I_vs_E = compute_current(result, Gamma1_array, Gamma2_array, E_array,
                                         mu1, mu2, config.temperature)

        I_array[iV] = I

        # Store result
        scba_results.append(result)

        if config.verbose:
            conv_str = "✓" if result['converged'] else "⚠"
            print(f"  {conv_str} I = {I*1e12:.3f} pA (iterations: {result['iterations']})")
    
    # ========================================================================
    # IETS CALCULATION
    # ========================================================================
    
    if config.verbose:
        print("\n[5/5] Computing IETS spectrum...")
    
    dIdV, d2IdV2 = compute_iets(V_array, I_array)
    
    if config.verbose:
        print(f"  dI/dV range: {dIdV.min():.3e} to {dIdV.max():.3e} S")
        print(f"  d²I/dV² range: {d2IdV2.min():.3e} to {d2IdV2.max():.3e} S/V")
        n_peaks = np.sum(d2IdV2 > 0)
        print(f"  Positive features: {n_peaks}/{len(d2IdV2)}")
    
    # ========================================================================
    # FINALIZE
    # ========================================================================
    
    elapsed = time.time() - start_time
    
    if config.verbose:
        print("\n" + "="*70)
        print(f"SIMULATION COMPLETE in {elapsed:.1f} seconds")
        print("="*70 + "\n")
    
    return {
        'config': config,
        'device': device,
        'molecule': molecule,
        'grid': grid,
        'phonon_modes': phonon_modes,
        'E_array': E_array,
        'V_array': V_array,
        'I_array': I_array,
        'dIdV': dIdV,
        'd2IdV2': d2IdV2,
        'scba_results': scba_results,
        'elapsed_time': elapsed
    }

# ============================================================================
# SAVE RESULTS
# ============================================================================

def save_results(results, filename=None):
    """Save simulation results to file"""
    
    import csv
    
    config = results['config']
    
    # Create output directory
    os.makedirs(config.output_dir, exist_ok=True)
    
    # Generate filename
    if filename is None:
        filename = f"{config.device_name}_{config.molecule_name}_IETS.csv"
    
    filepath = os.path.join(config.output_dir, filename)
    
    # Write CSV
    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Header
        writer.writerow(['# IETS Simulation Results'])
        writer.writerow(['# Device:', config.device_name])
        writer.writerow(['# Molecule:', config.molecule_name])
        writer.writerow(['# Temperature:', f'{config.temperature} K'])
        writer.writerow(['#'])
        writer.writerow(['V (V)', 'I (A)', 'dI/dV (S)', 'd2I/dV2 (S/V)'])
        
        # Data
        for V, I, g, iets in zip(
            results['V_array'], results['I_array'], 
            results['dIdV'], results['d2IdV2']
        ):
            writer.writerow([f'{V:.6e}', f'{I:.6e}', f'{g:.6e}', f'{iets:.6e}'])
    
    print(f"Results saved to: {filepath}")
    return filepath

# ============================================================================
# CLI ARGUMENT PARSING
# ============================================================================

def parse_args():
    """Parse command-line arguments"""
    import argparse

    parser = argparse.ArgumentParser(
        description='Run IETS simulation for a single molecule on an RTD device',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required arguments
    parser.add_argument('--device', '-d', type=str, default='GaAs_AlAs_symmetric',
                       help='Device name from device_library')
    parser.add_argument('--molecule', '-m', type=str, default='Benzene',
                       help='Molecule name from molecular_database')

    # Bias sweep
    parser.add_argument('--v-min', type=float, default=0.0,
                       help='Minimum bias voltage (V)')
    parser.add_argument('--v-max', type=float, default=0.5,
                       help='Maximum bias voltage (V)')
    parser.add_argument('--v-points', type=int, default=26,
                       help='Number of bias points')

    # Energy grid
    parser.add_argument('--e-min', type=float, default=-0.3,
                       help='Minimum energy (eV)')
    parser.add_argument('--e-max', type=float, default=1.5,
                       help='Maximum energy (eV)')
    parser.add_argument('--e-points', type=int, default=200,
                       help='Number of energy points')

    # Physics parameters
    parser.add_argument('--temperature', '-T', type=float, default=300,
                       help='Temperature (K)')
    parser.add_argument('--bulk-phonon', type=float, default=0.036,
                       help='Bulk phonon energy (eV)')
    parser.add_argument('--bulk-coupling', type=float, default=0.010,
                       help='Bulk phonon coupling (eV)')
    parser.add_argument('--mol-coupling-scale', type=float, default=1.0,
                       help='Molecular coupling scale factor')

    # SCBA parameters
    parser.add_argument('--scba-max-iter', type=int, default=50,
                       help='Maximum SCBA iterations')
    parser.add_argument('--scba-tol', type=float, default=1e-4,
                       help='SCBA convergence tolerance')
    parser.add_argument('--scba-mix', type=float, default=0.3,
                       help='SCBA mixing parameter')

    # Quasi-3D options
    parser.add_argument('--multimode', action='store_true',
                       help='Enable quasi-3D multimode transport')
    parser.add_argument('--Ly', type=float, default=1.0,
                       help='Transverse width y (um)')
    parser.add_argument('--Lz', type=float, default=1.0,
                       help='Transverse width z (um)')
    parser.add_argument('--n-max', type=int, default=3,
                       help='Max mode index in y')
    parser.add_argument('--m-max', type=int, default=3,
                       help='Max mode index in z')

    # Hybrid mode selection
    parser.add_argument('--hybrid', action='store_true',
                       help='Enable hybrid mode selection (requires --multimode)')
    parser.add_argument('--n-inelastic', type=int, default=4,
                       help='Number of inelastic modes for hybrid')

    # Output options
    parser.add_argument('--output-dir', '-o', type=str, default='results',
                       help='Output directory')
    parser.add_argument('--output-file', type=str, default=None,
                       help='Output filename (auto-generated if not specified)')
    parser.add_argument('--save-npz', action='store_true',
                       help='Save full results as .npz file')
    parser.add_argument('--quiet', '-q', action='store_true',
                       help='Quiet mode (minimal output)')

    # Utility
    parser.add_argument('--list-devices', action='store_true',
                       help='List available devices and exit')
    parser.add_argument('--list-molecules', action='store_true',
                       help='List available molecules and exit')

    return parser.parse_args()


def config_from_args(args):
    """Create SimulationConfig from parsed arguments"""
    config = SimulationConfig()

    # Device and molecule
    config.device_name = args.device
    config.molecule_name = args.molecule

    # Bias sweep
    config.V_min = args.v_min
    config.V_max = args.v_max
    config.V_points = args.v_points

    # Energy grid
    config.E_min = args.e_min
    config.E_max = args.e_max
    config.E_points = args.e_points

    # Physics
    config.temperature = args.temperature
    config.bulk_phonon_energy = args.bulk_phonon
    config.bulk_phonon_coupling = args.bulk_coupling
    config.molecular_coupling_scale = args.mol_coupling_scale

    # SCBA
    config.scba_max_iter = args.scba_max_iter
    config.scba_tolerance = args.scba_tol
    config.scba_mixing = args.scba_mix

    # Quasi-3D
    config.use_multimode = args.multimode
    config.Ly = args.Ly * 1e-6  # um to m
    config.Lz = args.Lz * 1e-6
    config.n_max = args.n_max
    config.m_max = args.m_max

    # Hybrid
    config.use_hybrid = args.hybrid
    config.n_inelastic = args.n_inelastic

    # Output
    config.output_dir = args.output_dir
    config.verbose = not args.quiet

    return config


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    args = parse_args()

    # Handle utility options
    if args.list_devices:
        from config.device_library import DEVICES
        print("\nAvailable devices:")
        for name in DEVICES.keys():
            print(f"  - {name}")
        sys.exit(0)

    if args.list_molecules:
        from config.molecular_database import MOLECULES, PERCEPTUAL_CLASSES
        print("\nAvailable molecules by class:")
        for cls, mols in PERCEPTUAL_CLASSES.items():
            print(f"\n  {cls}:")
            for mol in mols:
                print(f"    - {mol}")
        sys.exit(0)

    # Create configuration from args
    config = config_from_args(args)

    # Validate hybrid requires multimode
    if config.use_hybrid and not config.use_multimode:
        print("Warning: --hybrid requires --multimode, enabling multimode")
        config.use_multimode = True

    # Run simulation
    results = run_iets_simulation(config)

    # Save results
    save_results(results, filename=args.output_file)

    # Save .npz if requested
    if args.save_npz:
        npz_file = os.path.join(config.output_dir,
                                f"{config.device_name}_{config.molecule_name}_data.npz")
        np.savez(npz_file,
                 V_array=results['V_array'],
                 I_array=results['I_array'],
                 dIdV=results['dIdV'],
                 d2IdV2=results['d2IdV2'],
                 E_array=results['E_array'])
        print(f"Data saved to: {npz_file}")

    # Print summary
    print("\nSUMMARY:")
    print(f"  Device: {config.device_name}")
    print(f"  Molecule: {config.molecule_name}")
    print(f"  Mode: {'Quasi-3D' if config.use_multimode else '1D'}" +
          (f" (hybrid, {config.n_inelastic} inelastic)" if config.use_hybrid else ""))
    print(f"  Peak current: {results['I_array'].max()*1e12:.3f} pA")
    print(f"  Max conductance: {results['dIdV'].max():.3e} S")
    print(f"  IETS features: {np.sum(results['d2IdV2'] > 0)}")
    print(f"  Time: {results['elapsed_time']:.1f} s")
