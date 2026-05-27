# Quantum Electronic Nose — ZnO/MgZnO RTD-IETS (NEGF-SCBA)

NEGF study of **room-temperature Inelastic Electron Tunneling Spectroscopy
(IETS)** in a BEOL-compatible **ZnO/Mg₀.₃Zn₀.₇O resonant tunneling diode (RTD)**,
operated as a quantum biomimetic electronic nose. Agarwal & Ganguly, IIT Bombay —
SISPAD 2026 (`SISPAD_2026_Abhineet-9.pdf`).

The transport problem is solved with a **rank-1 projected NEGF-SCBA**: a localized
molecular vibron whose transverse-mode coupling matrix is rank-1 (to O(10⁻⁶) for
σ_mol ≈ 3 Å) collapses the 3D phonon-dressed NEGF problem to a single projected 1D
SCBA loop. Derivation: `docs/METHOD_DERIVATION.md`.

---

## What's here (production pipeline)

```
quantum_enose/
├── core/
│   ├── scba_rank1_keldysh.py   # rank-1 Keldysh SCBA driver (Anderson/DIIS mixing)
│   ├── scba_solver_rank1.py    # rank-1 projected Dyson algebra
│   ├── iets_analytic.py        # analytic d²I/dV² (stored in npz; NOT used for figs)
│   └── self_energy.py          # gaussian_chi_projector (used by a test)
├── config/
│   ├── device_library.py       # ZnO_MgZnO_symmetric + other stacks
│   └── molecular_database.py   # Mol_A/B/AB (+ VOC database)
├── run/
│   ├── run_rank1_sweep.py      # bias/temperature sweep driver
│   ├── make_paper_figs.py      # the SISPAD paper figures (fig1_IV … fig5_temp)
│   ├── verify_reproduction.py  # fast (~25 s) reproduction gate
│   └── reproduce_sispad.sh     # full from-scratch sweep + figures
├── tests/                      # 6 production tests + Patil reference data
├── docs/                       # current method/stack/benchmark docs (see below)
├── results/sispad_scba_2026-04-14/   # the SISPAD paper results (figs + 51-pt npz)
├── papers/                     # references, LaTeX sources, poster
└── archive/                    # legacy code/scripts/old results (gitignored, local only)
```

## Device & key numbers (SISPAD `-9`)

- Symmetric ZnO / Mg₀.₃Zn₀.₇O RTD: 2 nm barriers / 3 nm well / 10 nm n⁺ contacts.
- ΔEc = 0.47 eV, m\* = 0.28 m₀, ℏω_LO ≈ 72 meV, 10 µm × 10 µm pixel.
- NDR at **V_res ≈ 560 mV**, Baseline peak **≈ 69 nA** (PVR ≈ 6) at 300 K.
- Analytes: Mol_A (ℏω=100 meV), Mol_B (180 meV), Mol_AB (both), Baseline (bulk LO only).
- Sweep: 0–800 mV, **51 points (16 mV)**, Anderson/DIIS mixing (M=8, β=0.3),
  current conservation < 0.2%. d²I/dV² = raw `np.diff(I,2)/dV²` (no smoothing).

## Reproduce

```bash
source venv/bin/activate

# Fast gate (~25 s): re-run key bias points, check vs golden data
python run/verify_reproduction.py

# Regenerate the paper figures from the golden npz (seconds)
python run/make_paper_figs.py
#   -> results/sispad_scba_2026-04-14/{fig1_IV,fig2_d2IdV2,fig3_deltaI,fig4_deltaD,fig5_temp}

# Full from-scratch sweep + figures (~2 h at 51 pts)
bash run/reproduce_sispad.sh
```

Full recipe, parameters, data locations, and the Overleaf figure-name mapping:
**`docs/REPRODUCE_SISPAD.md`**.

## Documentation

- `docs/REPRODUCE_SISPAD.md` — authoritative reproduction recipe (start here).
- `docs/METHOD_DERIVATION.md` — rank-1 projected SCBA derivation + parity theorem.
- `docs/STACK_DECISION.md` — why ZnO/Mg₀.₃Zn₀.₇O (vs In₂O₃/κ-Ga₂O₃).
- `docs/PATIL_BENCHMARK.md`, `docs/PATIL_OCTAVE_REFERENCE.md` — 1D GaAs benchmark
  (validated by `tests/test_patil_reference.py`, `tests/test_patil_octave_match.py`).
- `rules.md` — code conventions.

Legacy framework docs and the pre-submission draft/checklist live under
`archive/docs/`; the older multimode/hybrid solver and exploratory scripts are
under `archive/` (kept locally, not versioned).

## Tests

```bash
python -m pytest tests/ -q
```
47 pass. Four `tests/test_rank1_keldysh.py` cases (the Patil 1D benchmark γ/δ
tightening items) are known-failing and unrelated to the ZnO production results —
see `docs/REPRODUCE_SISPAD.md`.

## References

- Patil, Saha, Ganguly, *Sci. Rep.* 8:128 (2018) — RTD-IETS concept.
- Pandey et al., *Sci. Rep.* 11:11389 (2021) — biomimetic e-nose / VOC modes.
- Datta, *Quantum Transport: Atom to Transistor*, Cambridge — NEGF.

IIT Bombay, Dept. of Electrical Engineering.
