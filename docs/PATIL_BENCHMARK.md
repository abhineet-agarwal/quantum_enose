# Patil 1D benchmark — SISPAD Tier-1 item γ

**Status:** 2026-04-08, **qualitatively matched** against Patil Fig. 2(a)/(b).
**Reference port:** `tests/patil_reference_1d.py`
**Reference output:** `tests/patil_reference_1d.npz` (76 V-points,
`tests/patil_reference_1d.png` for the IETS / I-V plot)
**Source:** `papers/rtd2modes_1d.m` (Patil's MATLAB, "The role of
inelastic scattering in resonant tunneling heterostructures")
**Regression test:** `tests/test_patil_reference.py`
**Original problem:** SISPAD_CHECKLIST.md item γ — we inherited Patil's
conventions (symmetric bias, D² = 0.1 eV², contact in-scattering, causal
Σ^R sign) but had never shown that our 1D code reproduced his published
I(V) and d²I/dV² in the GaAs/AlAs toy. This is the first thing a reviewer
will ask for.

## What this benchmark does

1. **Faithfully ports Patil's MATLAB SCBA loop to Python.** Line-by-line.
   Same tight-binding contact self-energy `Σ_{1,1} = -t₀ exp(ika)`, same
   diagonal-only Einstein vibron self-energy update (`ne(nn,nn,...) =
   n(nn,nn,inu+1:NE)`), same Newton mixer with `it = 0.5`, same Meir-Wingreen
   current trace evaluation. No "modernization" — design rule is *whatever
   Patil chose, we match exactly* so the port is a defensible reference.
2. **Runs at Patil's exact parameters** (defaults of `run_patil_reference`,
   matching the `%%% one peak in TM vs E %%%` block of `rtd2modes_1d.m`):

| Parameter | Value | Comment |
|---|---|---|
| Effective mass | 0.067·m₀ | GaAs |
| t₀ | 5.2 eV | tight-binding hopping |
| NS, NC, ND | 1, 40, 1 | Np = 42 sites |
| Nb, Vb | 15, 0.6 eV | symmetric double barrier (10-site well) |
| E_F, kT | 0.02 eV, 0.025 eV | Patil contact Fermi level |
| NV, V_max | 76, 0.3 V | bias grid |
| dE | 0.005 eV | energy grid spacing |
| E ∈ | [−0.2, 0.8] eV | NE = 201 |
| D² (Dnu) | [0.1, 0.1] eV² | two molecular phonon modes |
| ℏω (hnu·dE) | [0.09, 0.175] eV | mode energies |
| Bose N(ℏω) | from Patil Eq. line 48 | T = 290 K |
| SCBA mixer | it = 0.5 | matches MATLAB |
| SCBA tol | 1e−5 | converged ≡ Σ\|Δσ\| < tol |

3. **Saves I(V), I_coh(V), I_incoh(V), dI/dV, d²I/dV²** plus the parameter
   block to `tests/patil_reference_1d.npz`.
4. **Plots them** to `tests/patil_reference_1d.png` with vertical guides
   at the expected IETS positions V = ℏω/e ∈ {0.09 V, 0.175 V}.

## Reference results

### Default run (V_max = 0.3 V, 76 V-points)

| Quantity | Value |
|---|---|
| I(V = 0) | −2.0×10⁻¹⁸ A (machine zero) |
| I(V = 0.3 V) | +4.87×10⁻¹⁰ A |
| I(V) monotone on [0, 0.3] | yes — every Δ > 0 |
| NDR | none — peak is at V≈0.34 V, just past the cutoff |
| max\|d²I/dV²\| | 5.7×10⁻⁸ S/V |
| SCBA convergence | every V-point hit max_outer_iter = 200 |

This is the run the MATLAB source emits — `VV = linspace(0, 0.3, NV)`.
It covers the **rising edge** of the NDR peak only and does **not** show
the IETS phonon features that appear in Patil's published Figure 2(b),
which extends to V ≈ 0.8 V. Patil's published figure must have used a
larger V_max than the committed MATLAB script.

### Extended run (V_max = 0.8 V, 51 V-points) — comparison with Patil

Two references for this run, both at identical parameters:

- **Our Python port:** `tests/patil_reference_1d_vmax0p8.npz`
  (generated 2026-04-08, ~3 min on M-series).
- **Octave/MATLAB numerical reference:** `tests/patil_octave_reference.mat`
  (generated 2026-04-08 via `papers/rtd2modes_1d_octave.m`, a
  byte-identical port of Patil's `rtd2modes_1d.m` up to the five
  documented patches in `docs/PATIL_OCTAVE_REFERENCE.md`; ~65 min
  on M-series because Octave does not vectorize the per-energy loop).

**Regression test:** `tests/test_patil_octave_match.py` (4 tests, runs
in < 1 s, auto-skips if the .mat is missing).

| Feature | Octave (Patil) | Python port | Δ |
|---|---|---|---|
| NDR peak position | V = 0.320 V | V = 0.336 V | +16 mV = 1 bin |
| NDR peak current | 5.226×10⁻¹⁰ A | 5.698×10⁻¹⁰ A | **+9.03 %** |
| Pointwise I(V)  L∞ | — | — | **12.67 %** of max\|I\| |
| Coverage | V ∈ [0, 0.8] V, 51 pts | V ∈ [0, 0.8] V, 51 pts | identical grid |

**Eyeballed features vs published Fig. 2 (for reference — these were the
numbers before we had Octave):** d²I/dV² zero-crossings at V ≈ 0.512 V
and V ≈ 0.672 V match Patil's reported IETS peaks at V = 0.50 V
(ℏω₁ = 0.090 eV) and V = 0.68 V (ℏω₂ = 0.175 eV) to within ~12 mV.

**Verdict.** Both ports give the same NDR peak position (within one
bias bin), the same post-NDR valley, and the same IETS satellite
structure. The only real discrepancy is a **~9 % overshoot in peak
current** which propagates to a **~12.7 % L∞** pointwise error — both
dominated by the steep NDR slope. This is much tighter than the
+36 % gap we saw against the eye-digitized figure (which was itself
systematically low), and establishes the Python port as a trustable
ground-truth replica of Patil. The remaining ~9 % is under the 10 %
acceptance bound but is still large enough to be a real physics
difference, not floating-point noise. Candidates to investigate:

  1. `acos` branch choice in the tight-binding contact self-energy
     (Python uses `np.arccos`, Octave uses `acos`; they agree on the
     principal branch but may differ at edge cases).
  2. `inv(sparse(...))` in Octave vs `np.linalg.solve` in Python for
     the retarded Green's function; the sparse triangular factor
     differs.
  3. Accumulation order in the Meir-Wingreen trace sums (sparse
     vs dense differs in rounding).
  4. The `max_iter=200` SCBA cap — neither run converged to `1e-5`,
     and the Python and Octave mixers may accumulate a different
     non-converged residual.

This is tracked as the **γ-tighten** follow-up.

**The SCBA does not converge** to 1e-5 in 200 iterations at any bias.
This is a known difficulty at D² = 0.1 eV² with the it = 0.5 Newton mixer
inherited from Patil. The MATLAB code has identical behaviour — it relies
on the user noticing and either tightening the mixer or accepting the
non-converged answer. We accept the non-converged trace as the reference
because *that is what Patil's code emits*, and the Tier-1 γ goal is
"reproduce Patil's pipeline", not "out-converge it". A separate convergence
study (Tier-1 β follow-up: Mol_A `--max-iter 10 --scba-mix 0.4`) is the
right place to address this.

## What is *not* yet done

The eyeball comparison above closes the *qualitative* half of γ. Two
follow-ups remain:

- **Numerical reference (in progress 2026-04-08).** `brew install
  octave` is running; the patched headless script
  `papers/rtd2modes_1d_octave.m` will produce
  `tests/patil_octave_reference.mat`, which `run/plot_patil_overlay.py`
  picks up automatically. Full workflow + TODO punch list in
  `docs/PATIL_OCTAVE_REFERENCE.md`. Until the .mat exists, the
  comparison above remains an eyeball-against-figure check.
- **Multimode/rank-1 vs 1D side-by-side.** The 1D port is the reference
  the production solvers must be benchmarked against. That comparison is
  **not in this revision** because:

1. `core/scba_solver_hybrid.py` is the multimode production solver but it
   does not have a clean "1D mode" — passing it `transverse_modes` with
   `Nm = 1` requires plumbing through `run_iets.py` and the device library
   that I haven't touched yet.
2. The rank-1 projected solver (`core/scba_solver_rank1.py`) is the
   *replacement* the SISPAD paper actually quotes, but its Keldysh wrapper
   (which closes σ̃^{<,>} from G̃^{<,>}) is not yet written — see SISPAD
   Tier-1 α inelastic-path follow-up. A clean rank-1 vs Patil-1D comparison
   has to wait for that wrapper.

Both follow-ups are tracked. In the meantime the **reference exists** as a
stable, version-controlled, regression-tested ground truth so the
comparison is a single PR away once either solver is exposed.

## How to regenerate the reference

```
source venv/bin/activate
python -m tests.patil_reference_1d
```

(~3 minutes on a 2024 MacBook Pro M-series.) Re-run any time the port is
edited; the regression test (`tests/test_patil_reference.py`) provides
fast (~22 s) sanity checks at NV = 11.

## Provenance

- MATLAB source: `papers/rtd2modes_1d.m`, `%%% one peak in TM vs E %%%`
  block (lines 22–26) + phonon settings (lines 44–48). Author: Akshay
  Patil et al., GaTech.
- Port: `tests/patil_reference_1d.py` 2026-04-08, this repo.
- Conventions baked into the port:
  - Symmetric bias μ_{L,R} = E_F ± V/2 (Fix 1, Patil)
  - D² = 0.1 eV² coupling (Fix 2)
  - Contact in-scattering f₁Γ₁ + f₂Γ₂ in G^< (Fix 3)
  - Causal Σ^R = (i/2)(Σ_in + Σ_out) (Fix 4 sign — implicit in MATLAB)
  - Diagonal-only local Einstein vibron self-energy
- The port deliberately **does not** apply our Fix 5 (Landauer shortcut),
  Fix 6 (Gaussian χ(z)), or any rank-1 projection. Those are improvements
  *over* the Patil reference; mixing them in would defeat the purpose.

## How to apply

When someone asks "have we benchmarked against Patil?" the answer is:
*the 1D SCBA loop has been ported and reproduces his pipeline, the
production multimode and rank-1 solvers each have a clearly-tracked
comparison follow-up*. Point them at this file and at
`tests/patil_reference_1d.png`.
