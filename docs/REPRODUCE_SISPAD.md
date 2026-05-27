# Reproducing the SISPAD 2026 abstract

**Authoritative source:** `SISPAD_2026_Abhineet-9.pdf` (the `-9` revision; supersedes
`-5`). Title: *"NEGF Study of Room-Temperature IETS in a BEOL-Compatible
ZnO/Mg₀.₃Zn₀.₇O RTD for a Quantum Biomimetic Electronic Nose"*,
Agarwal & Ganguly, IIT Bombay.

This document pins exactly how to regenerate the abstract's results. Verified
**2026-05-27**: the current working tree reproduces the published numbers
bit-for-bit (see "Verification status" below).

---

## What the paper reports (ground truth)

| Quantity | Value |
|---|---|
| Device | ZnO / Mg₀.₃Zn₀.₇O **symmetric** RTD, 2 nm barriers / 3 nm well / 10 nm n⁺ contacts |
| Material | ΔEc = 0.47 eV, m\* = 0.28 m₀, ℏω_LO ≈ 72 meV, N_D = 10²⁵ m⁻³ |
| Pixel area | 10 µm × 10 µm = 10⁻⁷ cm² (Datta "Fix A" prefactor) |
| Bound state | E₁ = 149 meV |
| NDR | V_res ≈ 560 mV, Baseline peak |I_R| **≈ 69 nA**, PVR ≈ 6 (300 K) |
| Temp sweep (Mol_A) | peak rises 69 nA (300 K) → ~110 nA (10 K) |
| Analytes | Mol_A (ℏω=100 meV), Mol_B (180 meV), Mol_AB (both), Baseline (bulk 72 meV only) |
| Couplings | bulk D²_bulk = 0.001 eV²/site (χ≡1 in active region); molecular D₀² = 0.1 eV² (Gaussian χ(z) on emitter barrier, σ_mol = 3 Å) |

**Figures** (`-9`): Fig 3a = I-V, Fig 3b = d²I/dV², Fig 3c = ΔD discrimination,
Fig 4 = Mol_A temperature dependence.

---

## Method & exact parameters

Rank-1 projected NEGF-SCBA, single longitudinal mode (Nm=1).
Solver: `core/scba_rank1_keldysh.py` (Anderson/DIIS mixing, depth M=8, built in).
Driver: `run/run_rank1_sweep.py`.

| Parameter | Value | Note |
|---|---|---|
| device | `ZnO_MgZnO_symmetric` | driver default is *asymmetric* — must override |
| V grid | 0 → 0.8 V | driver default is 0–0.4 V — must override |
| V points | 201 (production figures) | abstract caption says 51; same physics, coarser |
| grid spacing | a = 0.2 nm → Np = 135 | |
| energy grid | dE = 2 meV, [−0.25, 0.5] eV → NE = 376 | |
| T | 300 K (+ 10/77/150 K for Fig 4) | |
| E_F | 20 meV | |
| SCBA mixing | mix (β) = **0.3** | driver default 0.4 — must override |
| SCBA max_iter | **100** | driver default 10 — must override |
| SCBA tol | **1e-4** | driver default 1e-5 — must override |
| convergence | 10–60 iters typical, current conservation < 0.2% | |

The driver's *defaults do not match the abstract* — always pass the overrides
(encoded in `run/reproduce_sispad.sh`).

---

## How to reproduce

**Fast gate (~25 s)** — confirm the current tree still produces the published
numbers without a full sweep:

```bash
python run/verify_reproduction.py
# re-runs Baseline/Mol_A at V=0.32/0.56 and checks |I_R| against golden npz
# PASS => 68.972 nA @ 560 mV etc., 0.00% diff
```

**Figures only (seconds)** — regenerate the paper data figures from the golden
data already in the repo:

```bash
python run/make_paper_figs.py
# reads the 51-pt results/sispad_scba_2026-04-14/*.npz, writes
#   fig1_IV, fig2_d2IdV2, fig3_deltaI, fig4_deltaD  (.png + .pdf)
```

**Full from-scratch sweep (~2 h at 51 pts)** — regenerate the .npz *and* the
data figures into a fresh dir:

```bash
bash run/reproduce_sispad.sh            # 51 points (paper grid), results/reproduce_<date>/
bash run/reproduce_sispad.sh 201        # denser grid (jittery raw d2I; ~14 h)
```

---

## Data & artifact locations

- **Golden npz:**
  `results/sispad_scba_2026-04-14/iets_ZnO_MgZnO_symmetric_*_0-800mV_*_rank1scba_*.npz`
  - `*_2026-04-14.npz` = **51-pt grid (the paper grid — use this)**;
    `*_2026-04-15.npz` = 201-pt grid (denser, but raw d2I is jitter-dominated).
  - `make_paper_figs.py` loads the 51-pt files (≤60 V-points).
- **Correct paper figures (Apr 14):** `results/sispad_scba_2026-04-14/`
  `fig1_IV`, `fig2_d2IdV2`, `fig3_deltaI`, `fig4_deltaD` (+ schematics
  `fig2_band_diagram`, `fig3_flowchart`, `fig1_table_stacks`). The Apr-14 PNG
  originals are preserved in `_correct_apr14_backup/`.
- **Overleaf rename mapping** (the paper `.tex` was edited on Overleaf):
  local `fig1_IV` → `fig4_IV`, `fig2_d2IdV2` → `fig5_d2IdV2`,
  `fig4_deltaD` → `fig7_deltaD`; band diagram/flowchart kept their names.
- **Do NOT use** `results/sispad_final_2026-04-14/*` (archived) — older 401-pt
  **FBA** run, not the SCBA production data.

---

## Provenance note: d²I/dV² and the lost generator

The correct paper figures (`fig1_IV`/`fig2_d2IdV2`/`fig3_deltaI`/`fig4_deltaD`,
Apr 14) were made by a generator that was **overwritten on 2026-04-28** by a
rebroken rewrite. That rewrite (the file that had been `run/regen_sispad_figs.py`,
now in `archive/run/`) introduced four regressions, all wrong:
1. used the **201-pt** grid instead of 51-pt (jittery I-V),
2. plotted the **analytic** `d2I` npz field (`core/iets_analytic`), smooth and
   ~±6 µA/V² — instead of numerical double-differentiation,
3. a serif/black color scheme instead of default matplotlib colors,
4. renumbered the outputs to `fig4_IV..fig8_temp_dependence`.

The **recovered** generator is `run/make_paper_figs.py`, verified to reproduce
the Apr-14 figures exactly. The d²I/dV² method is **raw second difference**
`d2 = np.diff(I, 2) / dV²` (NO smoothing), matching the paper's Fig 2 flowchart
and the +125 / −170 µA/V² scale in `fig2_d2IdV2` (the analytic field gives ±6;
`np.gradient²` gives ±55 — both wrong). The separate analytic-d²I accuracy issue
(under-captures dT/dV; `archive/docs/EXPLORATION_LOG.md`, 2026-04-12) is moot for
the figures, which never used it.

**Open item:** the paper's Fig 4 (temperature, `temp.png`) is not yet
regenerated — the surviving candidates disagree (Apr-14 `sispad_temp_sweep.png`
is 3-panel at ±100 µA/V²; the Apr-28 `fig8` is 2-panel at ±800), so the canonical
temp figure/method needs to be re-confirmed before adding it to `make_paper_figs.py`.

---

## Verification status (2026-05-27)

- `run/verify_reproduction.py` → PASS, 0.00% rel diff at V∈{0.32, 0.56} for
  Baseline and Mol_A (68.972 / 69.396 nA @ 560 mV).
- `run/make_paper_figs.py` → fig1_IV, fig2_d2IdV2, fig3_deltaI, fig4_deltaD
  reproduce the `SISPAD_2026_Abhineet-9.pdf` Fig 3 figures (51-pt, numerical
  d2I). The broken `regen_sispad_figs.py` was moved to `archive/run/`.
- Rank-1 pipeline is self-contained: it imports only `core/scba_rank1_keldysh`,
  `core/scba_solver_rank1`, `core/iets_analytic`, and the two config modules.
  The legacy `core/self_energy.py` (kept only because `tests/test_gaussian_chi.py`
  uses `gaussian_chi_projector`) and the archived `core/scba_solver_hybrid.py`
  are **not** on the reproduction path.
- `pytest tests/` → 47 pass, **4 pre-existing failures** in
  `tests/test_rank1_keldysh.py` (the Patil GaAs 1D benchmark, Tier-1 γ/δ — an
  open tightening item in `SISPAD_CHECKLIST.md`, not the ZnO production case and
  not used by any `-9` figure). The rank-1 Keldysh wrapper currently gives
  I_R≈1.7e-11 A at the Patil NDR peak vs the 5.7e-10 A reference; this regressed
  before the 2026-05-27 archive reorg and is unrelated to it.

## Repo layout note (2026-05-27 archive reorg)

Non-production code, scripts, old results, and superseded artifacts were moved
into `archive/` (preserving subtree structure: `archive/{core,run,tests,docs,
results,root_artifacts}` plus whole dirs `analysis/ debug/ scripts/ logs/
figures_publication/ data/ quantum_enose_poisson/`). Nothing was deleted; any
file can be restored by moving it back. Exploratory docs referenced elsewhere
(e.g. `EXPLORATION_LOG.md`) now live under `archive/docs/`.
