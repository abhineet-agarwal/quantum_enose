# Implemented transport model and validation

The executable model is a longitudinal effective-mass tight-binding RTD with
local diagonal phonon SCBA in one transverse channel. At every bias the code
solves the retarded Green function and the contact-plus-phonon Keldysh equations,
then recomputes observables after the final self-energy update. Bulk LO and
synthetic molecular modes have separate longitudinal envelopes. The production
SCBA uses only the imaginary broadening part of the retarded phonon self-energy;
this is a model approximation, not a general causal self-energy reconstruction.

## Scope and units

- `I_L/I_R`: A per spin per transverse mode, prefactor q²/h with energy in eV.
  These currents have opposite signs; conservation requires I_L + I_R = 0.
- `I_device`: independent-channel transverse quadrature including spin and
  pixel area; NaN when not calculated. This is not the localized-molecule
  multimode scattering solution.
- `I_dev_elastic`: transverse integral of coherent transmission only.
- `core/scba_solver_rank1.py`: tested algebraic projection helpers. The actual
  multimode Keldysh reconstruction remains unimplemented.
- `d2I` in sweep NPZ files: a frozen-kernel analytic diagnostic, not the full
  voltage derivative. Historical figures use raw numerical second differences.

## Poisson charge and stopping criteria

The returned `n_e` is assembled by the same function as the iterative charge:
quantum density, optional offset/callback correction, then Thomas–Fermi
substitution on deep terminal sites. A `DensityCorrection` callback carries
both density and inner convergence status; unfinished inner solves cannot
certify the Poisson result.

Convergence requires a small predictor step and a small unscreened defect
`max(abs(solve_poisson(n(U)) - U))`, using the same boundary conditions and
interior clamps. This residual is independent of predictor damping. It is an
equation diagnostic in eV, not a nonlinear solution-error estimate. The final
charge and residual are checked again at the returned potential.

The coupled spacer model uses `n_window_3D * (n_SCBA_1mode/n_ballistic_1mode - 1)`.
It is self-consistent within this **approximate** windowed charge closure. Full
transverse scattering charge and finite-window convergence remain to be tested.
The older `physical` mode without a correction is a contact-injection density
approximation when fed a phonon-dressed retarded function; it omits phonon
in-scattering and must not be called exact SCBA charge.

## Numerical evidence

See [VALIDATION.md](VALIDATION.md). Historical reproduction is kept separate
from physical and grid validation. Passing unit tests verifies selected algebra
and implementation properties, not the whole scientific claim.

New sweep defaults use 500 SCBA iterations, tolerance 1e-7; Poisson defaults use
1000 iterations and 1e-5 eV step/unscreened-residual tolerances. These defaults
are safeguards, not a promise of convergence. Current conservation uses an
absolute 1e-14 A floor plus 0.2% relative tolerance.

The spectrum gate compares molecular current differences and second derivatives
on the finest two energy grids, and derivatives on the finest two voltage
stencils. It records all intermediate grids. It deliberately does not infer
grid convergence from a solver flag or from visually smooth curves.

## Temperature and output provenance

`--donor-density` solves Ef(T) at fixed fully ionized parabolic-band density.
Without it, `--Ef` stays fixed and the implied density changes with temperature.
NPZ metadata identifies this choice, current units, model scope, numerical
acceptance and the absence of spectral certification. Results are never silently
overwritten. Historical data and submitted PDFs remain unchanged.
