"""
Fast single molecule test with reduced parameters
"""
import sys
sys.path.insert(0, '.')

from run.run_single_molecule import SimulationConfig, run_iets_simulation, save_results

# Create config with reduced parameters for speed
config = SimulationConfig()

# Faster settings
config.device_name = "In2O3_Al2O3_symmetric"
config.molecule_name = "Benzene"
config.grid_spacing = 0.2e-9  # Coarser grid (faster)
config.E_points = 80  # Fewer energy points
config.V_points = 11  # Fewer bias points (0.05 V steps)
config.V_max = 0.5
config.scba_max_iter = 15  # More iterations for convergence
config.scba_tolerance = 5e-3  # Looser tolerance
config.verbose = True

print("=" * 70)
print("FAST IETS SIMULATION TEST")
print("=" * 70)
print(f"Device: {config.device_name}")
print(f"Molecule: {config.molecule_name}")
print(f"Grid spacing: {config.grid_spacing*1e9:.2f} nm")
print(f"Energy points: {config.E_points}")
print(f"Bias points: {config.V_points}")
print(f"Expected time: ~5-10 minutes")
print("=" * 70)
print()

# Run simulation
results = run_iets_simulation(config)

# Save results
filepath = save_results(results)

# Print summary
print("\n" + "=" * 70)
print("SIMULATION COMPLETE!")
print("=" * 70)
print(f"Device: {config.device_name}")
print(f"Molecule: {config.molecule_name}")
print(f"Peak current: {results['I_array'].max()*1e12:.3f} pA")
print(f"Max conductance: {results['dIdV'].max():.3e} S")
print(f"Time: {results['elapsed_time']:.1f} seconds")
print(f"\nResults saved to: {filepath}")
print("=" * 70)
