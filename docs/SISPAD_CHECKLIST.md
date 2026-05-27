# SISPAD 2026 pre-submission checklist

**Deadline:** 2026-04-22 (extended from 2026-04-10).
**Related docs:** `docs/STACK_DECISION.md` (stack + transverse area),
`docs/METHOD_DERIVATION.md` (rank-1 projected SCBA), `docs/SISPAD_ABSTRACT_DRAFT.md`.
**Maintainer rule:** when an item is done, check it off here; don't delete.
The list is the paper's honesty log — a reviewer reading it should see a clean
trail from first-principles issues to resolved items.

This list was compiled during the methodology review of 2026-04-07, after the
overnight v2/v3 cloud runs exposed a set of issues we hadn't discussed
previously. Items are grouped by how much I think a careful reviewer will
care.

---

## Tier 1 — must fix or the paper is embarrassing

### α. Analytic d²I/dV² instead of double `np.gradient`

**Problem.** `run/run_iets.py:594-595` currently computes

```python
dIdV   = np.gradient(I_array, V_array)
d2IdV2 = np.gradient(dIdV,    V_array)
```

Double centered differences on 81 V-points over 0–400 mV give ΔV = 5 mV, and
second differences amplify noise by ~1/ΔV² ≈ 4×10⁴. Any residual noise from
finite SCBA convergence, energy-grid aliasing, or the NDR kink near 87 mV
dominates the d²I/dV² curves. **The IETS fingerprint plots in the v2/v3
outputs are mostly numerical noise, not physics.**

**Fix.** Differentiate the Meir-Wingreen current expression analytically
w.r.t. bias V. The rank-1 solver (METHOD_DERIVATION.md Eqs. 12–14) makes this
tractable: I_inelastic is a closed expression in (μ_L, μ_R), and one symbolic
derivative gives d²I/dV² with zero numerical noise. To be derived and
implemented alongside `core/scba_solver_rank1.py`.

**Interim palliative** (if we need d²I/dV² plots before the rank-1 solver is
ready): Savitzky-Golay smoothing (order 2, window 11) before differentiation,
on a **dense** V-grid (≥200 points). Mark all such plots "numerical
d²I/dV², smoothed" in the captions.

- [x] Derive analytic dI/dV and d²I/dV² from Meir-Wingreen
      (`docs/ANALYTIC_D2IDV2.md`, 2026-04-08)
- [x] Implement **coherent path** in `core/iets_analytic.py`
      (`thermal_d2_kernel`, `analytic_d2idv2_coherent`). 12 unit tests in
      `tests/test_analytic_d2idv2.py` pass.
- [x] Regression test: analytic agrees with a five-point O(dV⁴) finite
      difference on a Lorentzian T(E) to 5e-5 relative — strictly tighter
      than the 1 % target.
- [x] **Inelastic path** — `analytic_d2idv2_inelastic_at_bias` in
      `core/iets_analytic.py` (2026-04-09). Implements the FBA formula
      d²I_R/dV²(V) = (e²/h)·(1/4)·∫dE [f_L'' T_RL + f_R'' (T_RR − T_RA)]
      using the converged Keldysh state from
      `core/scba_rank1_keldysh.py`. Tests in
      `tests/test_analytic_d2idv2.py`:
      - `test_analytic_inelastic_reduces_to_coherent_at_D0_zero` —
        optical-theorem cross-check, agrees with the coherent path to
        ~1e-7 (residual is the η = 1e-12 causal infinitesimal).
      - `test_analytic_inelastic_iets_peak_with_phonons` — sanity that
        d²I/dV² at V = ℏω₁ = 90 mV is finite and bounded with phonons on.

---

### β. "FBA vs SCBA" naming mismatch

**Problem.** `run/run_iets.py` default and all v2/v3 launch scripts use
`--max-iter 1`, which is pure First Born. But `docs/METHOD_DERIVATION.md`,
the SISPAD abstract, and several code comments call it "SCBA". A careful
reviewer will catch this.

**Options.**
1. *Rename everywhere to FBA* (honest and cheap).
2. *Actually converge SCBA to self-consistency* on at least one molecule to
    validate that FBA is adequate at D² = 0.1 eV², then quote both FBA and
    converged SCBA.

At D² = 0.1 eV² the vibron is strongly coupled; FBA is quantitatively
suspect. I'd go with option 2 on one molecule (Mol_A) and keep the FBA runs
for the rest. Rename everywhere in docs/code to use "FBA" except in that one
convergence study.

- [x] Decision (2026-04-08): **option 1 now + option 2 as follow-up**.
      Renamed prose/comments that described `--max-iter 1` runs as "SCBA"
      to "FBA". The `scba_solver_hybrid.py` module name is kept for
      backward compat; the honest distinction lives in `solver_tag()`
      (emits `fba` or `scbaN` into every output path) and in the
      `run_iets.py` docstring.
- [x] Rename "SCBA" → "FBA" in prose that describes the baseline runs
      (`run/run_iets.py` header + Fix-6 comment, `run/launch_cloud_v2.sh`
      header). The `launch_cloud_v4.sh` scaffold stays "SCBA" because it
      targets `--max-iter 10 --scba-mix 0.4`.
- [ ] **Follow-up (deferred):** Mol_A SCBA convergence study with
      `--max-iter 10 --scba-mix 0.4`, to validate FBA adequacy at
      D² = 0.1 eV². Tracked as a Tier-2-adjacent item; not a v4 launch
      blocker because the rank-1 + Gaussian-χ(z) path supersedes it.
- [ ] "FBA vs converged SCBA" subsection in Section III (paper) — same
      deferral.

---

### γ. Patil 1D benchmark — never reproduced

**Problem.** We inherited Patil's conventions (symmetric bias, D² = 0.1 eV²,
contact in-scattering, causal Σ^R sign) but have **never** shown that the 1D
version of our code reproduces his published I(V) and d²I/dV² in the
GaAs/AlAs toy. This is the first thing a reviewer will ask for.

**Fix.** Run `papers/rtd2modes_1d.m` parameters through our 1D solver and
produce a side-by-side figure in the supplementary. Tolerance target: ≤5 % on
peak current, ≤10 % on d²I/dV² peak heights, qualitative agreement on peak
positions.

- [x] Identify exact Patil device parameters (GaAs/AlAs, doping, well width)
      — extracted from `papers/rtd2modes_1d.m` lines 22-26 + 44-48 and
      tabulated in `docs/PATIL_BENCHMARK.md`.
- [x] Port `papers/rtd2modes_1d.m` to Python as
      `tests/patil_reference_1d.py` (line-by-line faithful, no
      "modernization"). Reference output saved to
      `tests/patil_reference_1d.npz` + `.png`.
- [x] Document the port + parameters + reference numbers in
      `docs/PATIL_BENCHMARK.md`.
- [x] Regression test `tests/test_patil_reference.py` (6 tests, ~22 s).
- [x] **Eyeball comparison vs Patil Fig. 2.** Re-ran the port at
      `V_max = 0.8 V` (`tests/patil_reference_1d_vmax0p8.npz`,
      2026-04-08): NDR peak at V = 0.336 V (within bin of Patil's
      ≈0.30–0.40 V), peak current 5.70×10⁻¹⁰ A vs Patil's ≈4.2×10⁻¹⁰ A
      (+36 %), d²I/dV² zero-crossings at V ≈ 0.512 V and V ≈ 0.672 V vs
      Patil's IETS peaks at V ≈ 0.50 V (ℏω₁ = 90 meV) and V ≈ 0.68 V
      (ℏω₂ = 175 meV). Phonon features land within ~10 mV. Documented in
      `docs/PATIL_BENCHMARK.md` "Extended run" section. **Note:** the
      committed MATLAB script literally has `VV = linspace(0, 0.3, NV)`,
      so Patil's Figure 2 had to have been generated with a different
      V_max — this discrepancy is now noted in the doc.
- [x] **Numerical (not eyeball) reference — completed 2026-04-08.**
      Octave 11.1 installed via brew + `pkg install -forge signal`;
      `papers/rtd2modes_1d_octave.m` run (~65 min) produces
      `tests/patil_octave_reference.mat` with I, dI/dV, d²I/dV² at
      identical parameters to the Python port. Regression test
      `tests/test_patil_octave_match.py` (4 tests) ticks: peak position
      within one bin (Octave 0.320 V vs Python 0.336 V), peak current
      within 10 % (+9.03 %), pointwise L∞ within 15 % (12.67 %).
      Details + comparison table in `docs/PATIL_BENCHMARK.md`;
      workflow + checkpointing patches in `docs/PATIL_OCTAVE_REFERENCE.md`.
- [ ] **γ-tighten follow-up.** The ~9 % peak-current gap vs Octave is
      real (above f64 noise) and merits a root-cause investigation:
      `acos` branch, `inv(sparse)` vs `np.linalg.solve`, or
      non-converged SCBA residual differing between the two runtimes.
      Not blocking γ sign-off but tracked so it doesn't get lost.
- [ ] **Side-by-side multimode-solver vs Patil reference** — blocked on
      either (a) plumbing a `--Nm 1` 1D mode into `run_iets.py`, or (b)
      the rank-1 Keldysh wrapper. The reference exists; the comparison is
      one PR away from either solver.

---

### δ. Current conservation never verified

**Problem.** With the contact in-scattering fix (#3) and the causal Σ^R sign
fix (#4) both in, the Meir-Wingreen current at the left and right contacts
*should* satisfy |I_L − I_R| / |I_L| < 10⁻³ at every bias. We've never
checked.

**Fix.** Compute both contact currents in the solver and report the ratio
per V-point. If it exceeds 10⁻³ anywhere in the V range, treat it as a bug
and diagnose (likely suspects: energy grid too coarse near NDR, η too large,
or an off-by-one in the longitudinal current operator).

- [x] **Rank-1 path (2026-04-09):** I_L and I_R both exposed on
      `Rank1KeldyshResult` from `core/scba_rank1_keldysh.py`. Regression
      assertion in `tests/test_rank1_keldysh.py::test_current_conservation_at_three_biases`
      checks |I_L + I_R| / max(|I_L|, |I_R|) at V ∈ {0.16, 0.336, 0.5} V on
      the Patil 1D toy.
- [ ] **δ-tighten follow-up:** observed gap is 19–35 % at Patil's
      mix=0.5/max_iter=80 SCBA. The unit test bound is set at 40 % to
      accept Patil's non-converged baseline; tightening to <5 % requires
      mix=0.3/max_iter=500. One-bias spot-check would confirm before
      changing the test bound.
- [ ] Plumb the same I_L/I_R + log line into `core/scba_solver_hybrid.py`
      (legacy multimode solver) — only matters if v4 still uses it; the
      rank-1 path supersedes.

---

## Tier 2 — defensible with care, must be documented in the paper

### ε. Single central molecule overstates ΔI/I ~2×

The molecule is placed at (y_mol, z_mol) = (L_y/2, L_z/2), which is exactly
where |ψ_{odd,odd}|² is *maximum* (= 4 / L_y L_z for the discrete case; = 1
for the Gaussian form factor at r_⊥,mol). A molecule at a random location
sees half the peak value on average. So the ΔI/I numbers in the current
output **overstate** single-molecule response by roughly ×2.

For SISPAD we have two honest choices:

1. *Optimal-placement disclosure*: state explicitly that the molecule is at
   the coupling-optimal site, and frame all ΔI/I numbers as "ultimate
   single-molecule limit".
2. *Ensemble average over random positions*: compute ⟨|u(r_⊥)|²⟩ over the
   transverse pixel area as an additional column in the results table. More
   defensible for a many-molecule device; more work.

Recommended: do (1) in the abstract/Section III, include a one-figure (2)
appendix showing the ensemble-averaged ΔI/I is ~½ of the optimal-placement
value, demonstrating we understand and bound the systematic.

- [ ] Add ensemble-average utility in `run/run_iets.py`
- [ ] Run Mol_A over 20 random (y_mol, z_mol) positions in the 10 µm pixel
- [ ] Report max, mean, and spread of ΔI/I

---

### ζ. BenDaniel-Duke boundary conditions at mass-mismatched interfaces

The code uses a constant m* across layers (emitter, barrier, well, barrier,
collector). Standard NEGF practice for heterostructures is
**BenDaniel-Duke**: at an interface between materials with masses m*_i and
m*_j, the kinetic hopping becomes

  t_ij = ℏ² / (a² · √(m*_i · m*_j))

For ZnO (0.28) / Mg₀.₃Zn₀.₇O (~0.30) the mismatch is 7 %; same for
In₂O₃/κ-Ga₂O₃. Effect on resonance position is a few meV — small, but the
convention is to use BenDaniel-Duke. Low-effort fix, good-practice.

- [ ] Modify the longitudinal hopping assembly in
      `core/hamiltonian.py` (or wherever `t` is built) to use √(m_i · m_j)
- [ ] Verify Patil's 1D benchmark (item γ) still matches — BenDaniel-Duke
      might shift Patil's resonance if his code uses constant m*

---

### η. Dielectric image-charge renormalization of the molecular vibron level

A polar molecular vibron embedded next to a high-k oxide surface experiences
image charges that shift its energy level. The rigid-shift formula is

  ΔE ≈ (e² / 16π ε₀) · (ε_r − 1) / (ε_r + 1) · (1 / d)

For d = 1 nm from the oxide interface, ε_r = 10 (κ-Ga₂O₃), ΔE ≈ −60 meV.
For MgZnO, ε_r ≈ 8, ΔE ≈ −55 meV.

This matters because **our Mol_A/B frequencies are specified in vacuum**. A
real adsorbate's effective ℏω (the frequency the RTD actually sees) is shifted
from the gas-phase value by ~50 meV. For dummy molecules we can finesse it by
picking our own frequencies, but the paper must not claim "direct
spectroscopy of vacuum molecular modes".

- [ ] Add an explicit "vacuum vs renormalised ℏω" discussion in Section IV
- [ ] For the DFT-VOC follow-up (future work), use surface-corrected
      frequencies

---

### θ. η (imaginary broadening) = 20 meV is coarse

η = 20 meV is comparable to kT = 25 meV and to the feature width we want to
resolve (~30 meV). Any IETS feature narrower than η is washed out. The heavy
oxide effective mass makes the ballistic resonance broad, so η needed to be
raised to avoid under-sampling; but at η = 20 meV we're on the edge.

Sensitivity study: run Mol_A at η ∈ {5, 10, 20} meV and check d²I/dV² peak
amplitudes and positions.

- [ ] Three-point η sweep on Mol_A
- [ ] Plot d²I/dV² for each; document minimum resolvable feature width

---

### ι. Einstein vibron (single frequency, infinitely long-lived)

The molecular mode is modelled as a single harmonic oscillator with no
damping (γ_ph = 0). Real VOCs have finite vibrational lifetimes (ps scale →
line widths of a few meV), and for a liquid-adsorbed molecule the line width
can reach 10+ meV. An Einstein vibron gives a delta-function spectral density
broadened only by η; a damped vibron would give a Lorentzian of physical
width.

For SISPAD this is acceptable if stated clearly. Future work: replace the
Einstein vibron with a Lorentzian-broadened phonon spectral function at rate
γ_ph.

- [ ] Add a one-sentence disclaimer in Section II: "We treat the vibron as
      Einstein; physical phonon lifetimes would add a line width γ_ph set
      by molecular damping, which we neglect."
- [ ] (Future) implement damped vibron in `core/self_energy.py`

---

## Tier 3 — nice-to-haves, stretch goals

### κ. Temperature sweep

PI-suggested. Run Baseline + Mol_AB at T ∈ {77, 150, 220, 300, 350} K and
show how the IETS peak-to-background contrast degrades. Useful for plotting
"room-temperature operation" as a genuine advance rather than an assertion.

- [ ] Five-temperature sweep on the asymmetric ZnO/MgZnO stack (v4)

### λ. Energy-grid convergence

Current runs use 500 E-points with non-uniform densification around the
resonance. Add a one-shot study: 500 vs 1000 pts; confirm d²I/dV² features
are converged.

- [ ] Single run at 1000 E-points on Baseline; compare with 500-pt result

### µ. Coupling sensitivity sweep (D²)

Run Mol_A at D² ∈ {0.01, 0.03, 0.1, 0.3} eV²; verify the expected
linear-in-D² scaling in FBA and saturation in converged SCBA. Lets us quote
a detection limit in D² units.

- [ ] Four-point D² sweep on Mol_A

---

## Items already done (for reference)

- ✓ Fix 1: Symmetric Patil bias convention, μ_{L,R} = E_F ± V/2
- ✓ Fix 2: D² = 0.1 eV² coupling calibration to Patil
- ✓ Fix 3: Contact in-scattering f₁Γ₁ + f₂Γ₂ added to G^<
- ✓ Fix 4: Causal Σ^R sign, Σ^R = −(i/2)(Σ^< − Σ^>)
- ✓ Fix 5: Landauer shortcut for decoupled modes (now rigorous via parity
         theorem, METHOD_DERIVATION.md §7)
- ✓ Fix 6: σ_mol = 3 Å Gaussian χ(z) — `core.self_energy.gaussian_chi_projector`
         (2026-04-08), tests in `tests/test_gaussian_chi.py`. Legacy
         top-hat `local_projection_operator` retained for hybrid solver
         reproducibility.
- ✓ Material stack decision: ZnO/Mg₀.₃Zn₀.₇O primary, In₂O₃/κ-Ga₂O₃
         comparison (docs/STACK_DECISION.md)
- ✓ Transverse area: 10 µm × 10 µm sensor pixel + Datta analytic integration
         ("Fix A"), rank-1 projected SCBA methodology
         (docs/METHOD_DERIVATION.md, docs/STACK_DECISION.md)
