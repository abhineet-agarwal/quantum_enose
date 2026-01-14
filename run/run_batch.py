"""
Batch Molecular Screening

Screen all molecules (or selected classes) and generate IETS database.
"""

import numpy as np
import sys
import os
import time
import csv
from datetime import datetime

# Setup path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

from config.molecular_database import (
    MOLECULES, PERCEPTUAL_CLASSES, get_molecule, get_class_molecules
)
from config.device_library import get_device, DEVICES
from run.run_single_molecule import SimulationConfig, run_iets_simulation

# ============================================================================
# BATCH CONFIGURATION
# ============================================================================

class BatchConfig:
    """Configuration for batch screening"""
    
    def __init__(self):
        # Device selection
        self.device_names = ["In2O3_Al2O3_symmetric"]  # Can add more
        
        # Molecule selection
        self.molecule_selection = "all"  # "all", "class_name", or list of names
        
        # Simulation parameters (copied from SimulationConfig)
        self.grid_spacing = 0.12e-9
        self.bulk_phonon_energy = 0.070  # 70 meV (In2O3)
        self.bulk_phonon_coupling = 0.010
        self.molecular_coupling_scale = 1.0
        
        # Reduced resolution for batch (faster)
        self.E_points = 150
        self.V_points = 21  # 0.025 V steps for 0-0.5 V
        
        # SCBA parameters (more tolerant for batch)
        self.scba_max_iter = 30
        self.scba_tolerance = 5e-4
        self.scba_mixing = 0.4
        
        # Output
        self.output_dir = "batch_results"
        self.save_individual = True
        self.save_summary = True
        self.verbose = False  # Quiet mode for batch

# ============================================================================
# MOLECULE SELECTION
# ============================================================================

def select_molecules(selection):
    """
    Select molecules based on criteria
    
    Parameters:
    -----------
    selection : str or list
        "all", class name, or list of molecule names
        
    Returns:
    --------
    molecule_names : list
        List of molecule names to simulate
    """
    
    if isinstance(selection, list):
        return selection
    
    if selection == "all":
        return list(MOLECULES.keys())
    
    if selection in PERCEPTUAL_CLASSES:
        return get_class_molecules(selection)
    
    raise ValueError(f"Unknown selection: {selection}")

# ============================================================================
# BATCH RUNNER
# ============================================================================

def run_batch_screening(batch_config):
    """
    Run batch screening of molecules
    
    Parameters:
    -----------
    batch_config : BatchConfig
        Batch configuration
        
    Returns:
    --------
    batch_results : dict
        Results for all simulations
    """
    
    print("\n" + "="*70)
    print("BATCH MOLECULAR SCREENING")
    print("="*70)
    
    start_time = time.time()
    
    # Select molecules
    molecule_names = select_molecules(batch_config.molecule_selection)
    
    print(f"\nDevices: {batch_config.device_names}")
    print(f"Molecules: {len(molecule_names)}")
    print(f"Total simulations: {len(batch_config.device_names) * len(molecule_names)}")
    
    # Create output directory
    os.makedirs(batch_config.output_dir, exist_ok=True)
    
    # Storage
    batch_results = {
        'config': batch_config,
        'simulations': [],
        'summary': {}
    }
    
    # Run simulations
    total_sims = len(batch_config.device_names) * len(molecule_names)
    sim_count = 0
    
    for device_name in batch_config.device_names:
        for molecule_name in molecule_names:
            sim_count += 1
            
            print(f"\n[{sim_count}/{total_sims}] {device_name} + {molecule_name}")
            print("-" * 70)
            
            # Create simulation config
            sim_config = SimulationConfig()
            sim_config.device_name = device_name
            sim_config.molecule_name = molecule_name
            sim_config.grid_spacing = batch_config.grid_spacing
            sim_config.bulk_phonon_energy = batch_config.bulk_phonon_energy
            sim_config.bulk_phonon_coupling = batch_config.bulk_phonon_coupling
            sim_config.molecular_coupling_scale = batch_config.molecular_coupling_scale
            sim_config.E_points = batch_config.E_points
            sim_config.V_points = batch_config.V_points
            sim_config.scba_max_iter = batch_config.scba_max_iter
            sim_config.scba_tolerance = batch_config.scba_tolerance
            sim_config.scba_mixing = batch_config.scba_mixing
            sim_config.output_dir = batch_config.output_dir
            sim_config.verbose = batch_config.verbose
            
            # Run simulation
            try:
                results = run_iets_simulation(sim_config)
                
                # Extract summary metrics
                summary = extract_summary_metrics(results)
                
                # Store
                batch_results['simulations'].append({
                    'device': device_name,
                    'molecule': molecule_name,
                    'results': results,
                    'summary': summary
                })
                
                # Save individual results
                if batch_config.save_individual:
                    filename = f"{device_name}_{molecule_name}.csv"
                    save_iets_csv(results, batch_config.output_dir, filename)
                
                # Print summary
                print(f"  ✓ I_peak = {summary['I_peak']*1e12:.2f} pA")
                print(f"  ✓ IETS features: {summary['n_iets_peaks']}")
                print(f"  ✓ Time: {results['elapsed_time']:.1f} s")
                
            except Exception as e:
                print(f"  ✗ FAILED: {e}")
                batch_results['simulations'].append({
                    'device': device_name,
                    'molecule': molecule_name,
                    'error': str(e)
                })
    
    # Generate summary
    if batch_config.save_summary:
        save_batch_summary(batch_results, batch_config.output_dir)
    
    # Total time
    elapsed = time.time() - start_time
    
    print("\n" + "="*70)
    print(f"BATCH COMPLETE in {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")
    print(f"Success: {sum(1 for s in batch_results['simulations'] if 'error' not in s)}/{total_sims}")
    print("="*70 + "\n")
    
    return batch_results

# ============================================================================
# SUMMARY METRICS
# ============================================================================

def extract_summary_metrics(results):
    """Extract key metrics from simulation results"""
    
    V = results['V_array']
    I = results['I_array']
    dIdV = results['dIdV']
    d2IdV2 = results['d2IdV2']
    
    # Find peaks
    I_peak = I.max()
    I_peak_V = V[I.argmax()]
    
    # Conductance
    G_peak = dIdV.max()
    G_peak_V = V[dIdV.argmax()]
    
    # IETS features (count positive peaks)
    n_iets_peaks = np.sum(d2IdV2 > 0.1 * d2IdV2.max())
    
    # IETS peak positions
    iets_peak_indices = np.where(d2IdV2 > 0.5 * d2IdV2.max())[0]
    iets_peak_voltages = V[iets_peak_indices] if len(iets_peak_indices) > 0 else []
    
    # Convergence
    n_converged = sum(1 for r in results['scba_results'] if r['converged'])
    convergence_rate = n_converged / len(results['scba_results'])
    
    return {
        'I_peak': I_peak,
        'I_peak_V': I_peak_V,
        'G_peak': G_peak,
        'G_peak_V': G_peak_V,
        'n_iets_peaks': n_iets_peaks,
        'iets_peak_voltages': iets_peak_voltages,
        'convergence_rate': convergence_rate
    }

# ============================================================================
# SAVE FUNCTIONS
# ============================================================================

def save_iets_csv(results, output_dir, filename):
    """Save IETS data to CSV"""
    
    filepath = os.path.join(output_dir, filename)
    
    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Header
        writer.writerow(['V (V)', 'I (A)', 'dI/dV (S)', 'd2I/dV2 (S/V)'])
        
        # Data
        for V, I, g, iets in zip(
            results['V_array'], results['I_array'],
            results['dIdV'], results['d2IdV2']
        ):
            writer.writerow([f'{V:.6e}', f'{I:.6e}', f'{g:.6e}', f'{iets:.6e}'])

def save_batch_summary(batch_results, output_dir):
    """Save batch summary to CSV"""
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"batch_summary_{timestamp}.csv"
    filepath = os.path.join(output_dir, filename)
    
    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Header
        writer.writerow([
            'Device', 'Molecule', 'Perceptual_Class',
            'I_peak (A)', 'I_peak_V (V)',
            'G_peak (S)', 'G_peak_V (V)',
            'N_IETS_peaks', 'Convergence_rate',
            'Status'
        ])
        
        # Data
        for sim in batch_results['simulations']:
            if 'error' in sim:
                writer.writerow([
                    sim['device'], sim['molecule'], '',
                    '', '', '', '', '', '', 'FAILED'
                ])
            else:
                mol = get_molecule(sim['molecule'])
                s = sim['summary']
                writer.writerow([
                    sim['device'], sim['molecule'], mol['perceptual_class'],
                    f"{s['I_peak']:.6e}", f"{s['I_peak_V']:.6e}",
                    f"{s['G_peak']:.6e}", f"{s['G_peak_V']:.6e}",
                    s['n_iets_peaks'], f"{s['convergence_rate']:.2f}",
                    'SUCCESS'
                ])
    
    print(f"\nBatch summary saved to: {filepath}")
    return filepath

# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    # Create configuration
    config = BatchConfig()
    
    # Command-line options
    if len(sys.argv) > 1:
        selection = sys.argv[1]
        if selection in PERCEPTUAL_CLASSES:
            config.molecule_selection = selection
            print(f"Running class: {selection}")
        elif selection == "test":
            # Quick test with 3 molecules
            config.molecule_selection = ["Benzene", "Toluene", "Naphthalene"]
            config.V_points = 11
            print("Running TEST mode (3 molecules)")
        else:
            config.molecule_selection = selection
    
    # Run batch
    batch_results = run_batch_screening(config)
    
    # Print final summary
    print("\nFINAL SUMMARY:")
    successful = [s for s in batch_results['simulations'] if 'error' not in s]
    print(f"  Total simulations: {len(batch_results['simulations'])}")
    print(f"  Successful: {len(successful)}")
    print(f"  Failed: {len(batch_results['simulations']) - len(successful)}")
    
    if successful:
        I_peaks = [s['summary']['I_peak'] for s in successful]
        print(f"  Current range: {min(I_peaks)*1e12:.2f} to {max(I_peaks)*1e12:.2f} pA")
