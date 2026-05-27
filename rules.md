# Quantum E-Nose Simulation Code Rules

## Core Principle: One Code, Many Configurations

All simulation scripts must follow the principle of **single codebase with reconfigurable parameters**. This ensures:
- Easy navigation and understanding
- Consistent methodology across simulations
- Simple parameter sweeps and comparisons
- Reproducible results

---

## Rule 1: Single Script Per Simulation Type

| Simulation Type | Script | Purpose |
|-----------------|--------|---------|
| Baseline (no molecule) | `run/run_baseline.py` | Ballistic RTD transport |
| IETS (with molecules) | `run/run_iets.py` | Molecular fingerprinting |
| Temperature study | `run/run_temperature.py` | Temperature dependence |

**DO NOT** create multiple versions (e.g., `run_iets_v2.py`, `run_iets_v3.py`). Instead, use command-line arguments or config files.

---

## Rule 2: Command-Line Configuration

All configurable parameters must be exposed via `argparse`:

```python
parser.add_argument('--device', '-d', type=str, required=True,
                    help='Device name from device_library.py')
parser.add_argument('--well-material', type=str, default=None,
                    help='Override well material')
parser.add_argument('--barrier-material', type=str, default=None,
                    help='Override barrier material')
parser.add_argument('--well-thickness', type=float, default=None,
                    help='Override well thickness (nm)')
parser.add_argument('--barrier-thickness', type=float, default=None,
                    help='Override barrier thickness (nm)')
parser.add_argument('--v-min', type=float, default=0.0,
                    help='Minimum bias voltage (V)')
parser.add_argument('--v-max', type=float, default=0.5,
                    help='Maximum bias voltage (V)')
parser.add_argument('--output-prefix', '-o', type=str, default=None,
                    help='Output file prefix (default: auto-generated)')
```

---

## Rule 3: Output File Naming and Folder Convention

All results are saved in a **dated subfolder** under `results/`:

```
results/YYYY-MM-DD/
    README.md                        ← auto-generated + human notes
    iets_{device}_{molecule}_{vrange}_{T}K_{solver}_{date}.npz
    iets_{device}_{molecule}_{vrange}_{T}K_{solver}_{date}.csv
    iets_{device}_{molecule}_{vrange}_{T}K_{solver}_{date}.png
    iets_{device}_{molecule}_{vrange}_{T}K_{solver}_{date}.pdf
```

The `solver` tag **must** encode the approximation level:
- `fba` — First Born Approximation (max_iter=1)
- `scba{N}` — Self-Consistent Born Approximation with N iterations

Examples:
```
results/2026-02-27/iets_In2O3_Ga2O3_kappa_symmetric_Baseline_0-400mV_300K_fba_2026-02-27.npz
results/2026-02-27/iets_In2O3_Ga2O3_kappa_symmetric_Mol_A_0-400mV_300K_fba_2026-02-27.npz
results/2026-03-01/iets_In2O3_Ga2O3_kappa_symmetric_Benzene_0-1500mV_300K_scba5_2026-03-01.npz
```

**NEVER overwrite existing result files.** Because the date is in the filename, re-running on the same day with the same settings would collide — in that case check for existence and raise an error or add a suffix.

Each `README.md` is auto-generated with a parameter table per run. **Edit it manually** to add observations, parameter changes, or reasons for the run.

---

## Rule 4: Device Library as Single Source of Truth

All material parameters and device configurations live in `config/device_library.py`.

To add a new device:
```python
DEVICES["MyNewDevice"] = {
    "description": "Description here",
    "layers": [...],
    "transverse_size": (1e-6, 1e-6),
    "molecule_location": "emitter_barrier"
}
```

**DO NOT** hardcode material parameters in simulation scripts.

---

## Rule 5: Molecular Database as Single Source

All molecular parameters live in `config/molecular_database.py`.

To add a new molecule:
```python
MOLECULES["NewMolecule"] = {
    "name": "New Molecule",
    "modes_meV": [...],
    "coupling_meV": [...],
    "description": "..."
}
```

---

## Rule 6: Consistent Plot Format

All plots must include:
1. Clear axis labels with units
2. Legend with device/molecule names
3. Grid lines (alpha=0.3)
4. Title indicating simulation type and key parameters
5. Consistent color scheme across related plots

Standard figure size: `(14, 10)` for 4-panel plots, `(12, 8)` for 2-panel.

---

## Rule 7: Data Saving

Every simulation must save (in the dated subfolder):
1. **Plot** (`.png`, 300 dpi + `.pdf`)
2. **Raw data** (`.npz`) with all computed arrays
3. **CSV** (`.csv`) — human-readable IV and IETS columns
4. **README.md** — auto-generated parameter table, edit to add notes

Data structure in `.npz`:
```python
np.savez(output_file,
    V_array=V_array,
    I_array=I_array,
    dIdV=dIdV,
    d2IdV2=d2IdV2,
    device_name=device_name,
    parameters=parameter_dict
)
```

---

## Rule 8: No Code Duplication

If functionality is needed in multiple scripts, extract it to a module:
- `core/transport.py` - Transport calculations
- `core/plotting.py` - Plotting utilities
- `core/analysis.py` - Analysis functions

---

## Rule 9: Verbose Mode

All scripts must support `--verbose` / `--quiet` flags:
```python
parser.add_argument('--verbose', '-v', action='store_true')
parser.add_argument('--quiet', '-q', action='store_true')
```

---

## Rule 10: Example Usage in Docstring

Every script must have usage examples in its docstring:

```python
"""
Baseline RTD Simulation

Usage:
    # Run with predefined device
    python run/run_baseline.py --device In2O3_Ga2O3_kappa_symmetric

    # Override thickness
    python run/run_baseline.py --device In2O3_Ga2O3_kappa_symmetric --barrier-thickness 1.5

    # Custom output name
    python run/run_baseline.py --device ZnO_MgZnO_2nm --output-prefix my_simulation
"""
```

---

## Summary

| Aspect | Rule |
|--------|------|
| Scripts | One per simulation type |
| Parameters | Command-line configurable |
| Materials | In device_library.py only |
| Molecules | In molecular_database.py only |
| Output folder | `results/YYYY-MM-DD/` per run date |
| Output naming | `{type}_{device}_{mol}_{vrange}_{T}K_{solver}_{date}.{ext}` |
| Plots | Consistent format, labeled, with legend |
| Data | Always save .npz + .csv + README alongside plots |

---

*Last updated: February 2026*
