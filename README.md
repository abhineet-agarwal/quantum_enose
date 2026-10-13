# Quantum Electronic Nose — 1D ZnO/MgZnO NEGF-SCBA

This repository studies vibrationally assisted transport in a ZnO/MgZnO RTD.
The implemented production model is **one transverse mode, one spin, with
local longitudinal phonon self-energies**. The separate rank-1 Dyson algebra
is tested, but a multimode Keldysh scattering solver is not implemented.

**Validation status (2026-10-13):** archived SISPAD currents reproduce within the historical gate's tolerance,
but the molecular second derivative fails energy/bias-grid refinement. The
69 nA peak is per mode per spin, not the current of a 10 µm × 10 µm pixel.
The newer coupled Poisson work uses an approximate windowed single-mode-ratio
charge correction. Neither path is a validated quantitative planar sensor model.
See [VALIDATION.md](VALIDATION.md) for measured checks and remaining work.

## Main paths

- `config/`: layer stacks, material parameters and synthetic molecular modes.
- `core/scba_rank1_keldysh.py`: 1D SCBA, dense/diagonal/RGF implementations,
  contact currents and independent-channel transverse integration.
- `core/poisson.py`, `core/poisson_negf.py`: electrostatics, charge assembly,
  fixed-doping contact filling and Poisson residual checks.
- `core/validation.py`: numerical diagnostics, spectral refinement and plot gates.
- `run/run_rank1_sweep.py`: fixed-potential or exploratory Poisson sweeps.
- `run/sc_coupled_sweep.py`: spacer-device Poisson with the approximate scattering charge.
- `run/check_spectrum_convergence.py`: compare molecular delta-I and delta-d2I
  under independent energy and voltage refinement.
- `results/sispad_scba_2026-04-14/`: immutable historical data, not validation targets.

`I_L` and `I_R` are per transverse mode per spin. `I_device`, when explicitly
computed with `--transverse`, sums **independent** transverse channels and
includes area and spin. It still omits localized intermode scattering.
`I_dev_elastic` integrates only the coherent transmission.

## Tests and historical reproduction

```bash
source venv/bin/activate
# A fresh environment also needs the test dependencies:
# python -m pip install -r requirements-dev.txt
python -m pytest tests/ -q
python run/verify_reproduction.py
python run/make_paper_figs.py --historical --out-dir results/historical_figures_new
```

The reproduction gate checks saved per-mode currents, not physical correctness.
`bash run/reproduce_sispad.sh` deliberately retains the old parameters, including
fixed Ef at all temperatures, for historical comparison. It explicitly permits
failed numerical checks and labels its figures. Existing result files and
figures are protected from replacement.

## Numerical validation

```bash
python run/check_spectrum_convergence.py --centers 0.096 0.56 --molecules Mol_A
```

The JSON report records raw currents, SCBA status, current conservation, delta-I,
delta-d2I and energy/bias refinement checks. Exit 1 means the gate failed; it is
not a reason to smooth the spectrum or choose a coarser grid. A pass covers only
the selected stencils and parameters. Energy-window, spatial-grid and model
validation remain separate requirements.

New sweeps require solver convergence and current conservation. Poisson also
requires converged inner charge solves and an unscreened equation defect below
tolerance; a small damped step alone no longer passes. Failing sweeps save
explicit diagnostic flags and exit nonzero. `--allow-unconverged` is available
for diagnostics or historical reproduction. Figure generators require spectral
validation or an explicit historical/diagnostic option and visible labelling.

SCBA now defaults to linear mixing in every implementation. Optional
`--anderson-depth 8` uses guarded extrapolation: negative/nonfinite rates are
rejected, and returned rates, populations and the self-energy residual are
checked before convergence is accepted. Final transport solves also continue
until the contact currents balance. The fixed-potential spacer check uses a
tighter 1e-9 self-energy tolerance and fresh output with explicit acceptance:

```bash
python run/fixedU_fingerprint.py Mol_A .015 --biases 576 1504 --dE .0005
```

Its archived input potentials are held fixed; a passing transport check does
not certify those potentials or the derivative. Spectrum checks also support
`--energy-offset .5` to test grid alignment; passing one alignment does not
establish independence from alignment.

For a temperature series at fixed fully ionized contact density:

```bash
for T in 10 77 150 300; do
  python run/run_rank1_sweep.py --device ZnO_MgZnO_symmetric --molecules Mol_A \
    --T "$T" --donor-density 5.03636495047155e24 --V-max 0.8 --V-points 51 \
    --out-dir results/fixed_doping_new
done
```

This sets Ef(T) from the parabolic-band Fermi integral. It does not account for
donor freeze-out. These are exploratory sweeps; choose and validate energy and
voltage grids before interpreting the derivative.

## Self-consistent work

```bash
python -u run/sc_coupled_sweep.py --mol Baseline --v1 1.6 --down-to 0.64
```

New coupled runs use fresh output directories. An explicit `--out-dir` can resume
a run certified by the current diagnostic schema; old/unconverged cached points
cannot be silently reused. A `--seed file.npz:key` is an initial potential only.

The charge correction is an approximation, not full transverse-integrated
SCBA. Contact imbalance is a numerical diagnostic, **not a statistical standard
deviation**, experimental noise level, or bound on model/grid error.
Older Poisson peak-shift claims were withdrawn; see
[docs/POISSON_STATUS.md](docs/POISSON_STATUS.md).
