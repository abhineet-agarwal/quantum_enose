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
