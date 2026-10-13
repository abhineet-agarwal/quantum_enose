# Numerical audit and corrections — 2026-10-12

The implemented model is a single transverse channel per spin with local
longitudinal phonon SCBA. A full localized-molecule multimode Keldysh solver
has **not** been implemented or validated by this change. Independent-channel
transverse integration is a different approximation. The submitted PDFs and
archived NPZ/figure files are preserved; the working manuscript source has been
corrected to match the implemented model and its evidence.

## What changed

1. The Poisson loop and final output share one charge-assembly routine, including
   offset/callback scattering corrections and Thomas–Fermi terminal charge.
2. Poisson convergence requires both a small predictor step and an unscreened
   equation defect, plus convergence of any nested charge solve. The final
   returned density and potential are checked again. The defect is not a bound
   on nonlinear, grid or modelling error.
3. Sweep acceptance also checks contact-current conservation and transverse-node
   solves. Failed results save diagnostics and exit nonzero. Missing device
   currents are NaN, not a zero-current placeholder. Existing results are
   protected from replacement; coupled runs use fresh output locations and do
   not silently resume old-schema or failed cached points.
4. `--donor-density` obtains Ef(T) from the parabolic-band Fermi integral. Both
   fixed-Ef and fixed-doping policies are recorded in output metadata. The
   historical reproduction script deliberately keeps fixed Ef and says so.
5. A spectrum gate tests molecular delta-I and delta-d2I under independent energy
   and voltage refinement. Passing a solver test never sets spectral validation.
6. Figure generators require numerical/spectral certification or explicit,
   visibly labelled historical/diagnostic plotting. Uncomputed `I_device`
   placeholders are rejected even in diagnostic mode.
7. Documentation and the working manuscript withdraw unsupported 3D, whole-pixel
   current, conservation, electrostatic-shift, temperature and sigma claims.
   Contact imbalance is a diagnostic, not a statistical uncertainty.

## Executed checks

- `venv/bin/python -m pytest tests/ -q`: **89 passed, 6 subtests passed**.
  Regression tests exercise returned charge with both boundary schemes,
  falsely small damped steps, unfinished inner solves, fixed-doping filling,
  output preservation, failed-sweep diagnostics and uncomputed plot fields.
- `venv/bin/python run/verify_reproduction.py`: **PASS**, exact archived I_R
  at Baseline 0.32/0.56 V and Mol_A 0.56 V. This is historical reproduction only.
- Fresh approximate coupled Baseline solve at 576 mV: **numerical checks pass**.
  139 Poisson iterations; unscreened defect 9.559e-6 eV; predictor step
  9.402e-9 eV; final contact imbalance 0.0304%. Current is 0.976068 nA per mode
  per spin, versus the earlier 1.006464 nA. The well potential changes by
  +0.1603 meV. [Raw report](validation/poisson_576_2026-10-12.json).
- Molecular spectrum refinement at 560 mV: **FAIL** despite every individual
  SCBA/current check passing. [Raw report](validation/spectrum_2026-10-12.json).

Mol_A minus Baseline second derivative, in µA/V² per mode per spin:

| Energy spacing | 16 mV voltage stencil | 8 mV voltage stencil |
|---|---:|---:|
| 2 meV | 1.422141 | 16.597249 |
| 0.5 meV | 1.216411 | 5.001512 |
| 0.1 meV | 0.624855 | 5.431076 |

The 0.1 meV result is **not** established as converged. The large variation is
why the gate fails. The check uses SCBA tolerance 1e-7, maximum 500 iterations,
fixed Ef=20 meV, T=300 K, the 2/3/2 nm original stack, a=0.2 nm, energy window
[-0.25, 0.5] eV, and the original bulk/molecular couplings. The acceptance
threshold is 5% plus explicit absolute floors of 1e-13 A for delta-I and
1e-10 A/V² for delta-d2I. Checks compare the finest two energy grids at each
voltage resolution, and the finest two voltage stencils on the finest energy
grid. The scope is only the selected stencils, not a whole spectrum.

Reproduce the report in a fresh location:

```bash
venv/bin/python run/check_spectrum_convergence.py --centers 0.56 \
  --output results/new_validation/spectrum.json
```

Exit 1 and a completed JSON with `passed: false` are expected on these grids.

## Remaining scientific work

- Converge the molecular derivative across the bias range in energy, voltage,
  spatial spacing and energy-window extent; narrow resonances require more than
  visual smoothness. The present report does not establish a passing grid.
- Validate the approximate Poisson scattering charge against an actual transverse
  integration of scattering-dependent density. The 576 mV check establishes a
  fixed point of the approximate charge closure only.
- Implement and validate localized-molecule transverse-mode mixing if retaining
  the proposed rank-1 3D method. Simply summing independent channel currents does
  not provide that calculation or its molecule/area normalization.
- Run fixed-doping temperature sweeps on validated grids. Ef=20 meV at every T
  varies the implied density from about 5.04e24 m^-3 at 300 K to 1.91e24 at 10 K.
- Revisit bulk coupling calibration, omitted real phonon self-energy, material
  interfaces, molecular coverage, and experimental noise before sensitivity
  claims. None is estimated by a contact-current imbalance.

The updated manuscript is a working source revision. Its archival figures are
labelled as unvalidated; a new submission-ready PDF has not been generated.

## Refined spectrum study at 560 mV — 2026-10-12

`run/check_spectrum_convergence.py --centers 0.56 --molecules Mol_A --energy-steps
0.002 0.0005 0.0002 0.0001 0.00005 0.000025 --bias-steps 0.016 0.008 0.004 0.002
--scba-max-iter 4000` with **linear SCBA mixing** (`--anderson-depth 0`, now the
default: Mol_A has D0^2 = 0.1, where Anderson/DIIS is unreliable; with it the
0.5 meV / 16 mV value above read 1.216 instead of 1.149). Every one of the 108
individual solves passed its numerical checks; the refinement gate **FAILS**.
[Raw report](validation/spectrum_560_refined_2026-10-12.json).

Mol_A minus Baseline, per mode per spin:

| Energy spacing | dD, 16 mV | dD, 8 mV | dD, 4 mV | dD, 2 mV | dI at 560 mV |
|---|---:|---:|---:|---:|---:|
| 2 meV | 1.42 | 16.60 | 151.6 | 273.5 | 0.4233 nA |
| 0.5 meV | 1.15 | 4.73 | 40.65 | 71.10 | 0.4561 nA |
| 0.2 meV | 0.56 | 9.22 | 31.17 | 36.90 | 0.5364 nA |
| 0.1 meV | 0.62 | 5.41 | 17.02 | 19.79 | 0.5168 nA |
| 0.05 meV | 0.61 | 3.11 | 7.76 | 8.22 | 0.5401 nA |
| 0.025 meV | 0.57 | 2.22 | 4.78 | 3.44 | 0.5612 nA |

(dD in µA/V².) Within this single-mode model:

- **dI is stable to about ±5 % for spacing <= 0.2 meV (~0.54 nA)**; the published
  2 meV grid gives 0.42 nA.
- **dD with the published 16 mV stencil is stable to about ±5 % for spacing
  <= 0.1 meV (~0.59 µA/V²)**; the published 2 meV value, 1.42 µA/V², is 2.4x larger.
- **dD on 8/4/2 mV stencils is not converged even at 0.025 meV**: it keeps falling
  as the grid is refined, consistent with the 134 µeV-wide E1 resonance being
  sampled unevenly (I(V) is jagged on a few-mV scale at 0.5 meV). A local
  curvature would need spacing well below the E1 linewidth.

"Stable to ±5 %" is the spread of the last refinements, not an error bound.

## Energy-grid convergence of the spacer-device fingerprint — 2026-10-13

Fixed coupled-Baseline potential, Mol_A + bulk LO at D² = 0.015 on every
undoped site, linear SCBA mixing, tolerance 1e-7:

| spacing | Mol_A current at 576 mV |
|---|---:|
| 4 meV | 1.29319711 nA |
| 2 meV | 1.30377999 nA |
| 1 meV | 1.30388337 nA |
| 0.5 meV | 1.30388641 nA |
| 0.2, 0.1 meV | 1.3039 nA (dI/I = +32.36 % on all three of 0.5/0.2/0.1 meV) |

At 1504 mV dI/I = +41.14 % on all three of 0.5/0.2/0.1 meV. Unlike the
published device (weak bulk coupling confined to barriers + well, E1 134 µeV
wide), ZnO-level scattering on every undoped site broadens the narrow states to
many meV, so the current is converged in energy spacing by ~1 meV. This covers
the energy grid at fixed potential only; the single-mode model, the approximate
scattering-charge closure and bias-stencil (d²I) convergence are not covered.
