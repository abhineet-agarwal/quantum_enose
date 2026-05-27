# Technical Documentation — Rank-1 Projected NEGF-SCBA Pipeline

Current production codebase behind `SISPAD_2026_Abhineet-9.pdf`. For the full
mathematical derivation see `docs/METHOD_DERIVATION.md`; for the exact
reproduction recipe see `docs/REPRODUCE_SISPAD.md`. The pre-2026-04 multimode/
hybrid framework this file used to describe is preserved at
`archive/docs/TECHNICAL_DOCUMENTATION_legacy.md`.

---

## 1. Physical model

A double-barrier ZnO/Mg₀.₃Zn₀.₇O RTD biased into resonance acts as a narrow
energy filter (≪ kT), so phonon-assisted (inelastic) channels show up as
resolvable features in d²I/dV² even at 300 K. Two phonon families enter:

- **Bulk LO phonon** (ℏω = 72 meV), uniform in the transverse plane, present
  always (D²_bulk = 0.001 eV²/site, χ ≡ 1 over the active region).
- **Molecular vibron(s)** localized at the emitter barrier with an isotropic
  Gaussian form factor (σ_mol = 3 Å), coupling D₀² = 0.1 eV².

## 2. Rank-1 projection (the key idea)

In the transverse-mode basis the electron-vibron coupling factorizes as

    M = D₀ · χ(z) · |u⟩⟨u| ,   u_nm = ψ_nm(r⊥,mol)

i.e. **rank-1 in mode space** to O((σ_mol·k⊥)²) ≈ 10⁻⁶. The SCBA self-energy
inherits the rank-1 structure, so the full 3D phonon-dressed Dyson equation
collapses to a single projected **1D** Dyson equation for G̃ᴿ = ⟨u|Gᴿ|u⟩. The
transverse area enters only as a Datta prefactor (A = 10 µm × 10 µm). At a
symmetric pixel centre, parity exactly selects the (odd,odd) modes — the
coherent/inelastic current split is a theorem, not a heuristic
(`docs/METHOD_DERIVATION.md`).

## 3. Module map

| File | Role |
|------|------|
| `core/scba_rank1_keldysh.py` | Single-bias Keldysh-consistent SCBA: bare contact self-energies (Patil tight-binding leads), Gᴿ/G^</> on the energy grid, FBA/SCBA phonon update, **Anderson/DIIS mixing** (history M=8, β=mix), Meir-Wingreen current at *both* contacts (I_L, I_R for conservation). Entry point: `run_rank1_keldysh_single_bias(...)`. |
| `core/scba_solver_rank1.py` | Projected Dyson algebra (`build_projected_bare_g`, `projected_dyson_step`) used by the Keldysh driver. |
| `core/iets_analytic.py` | Closed-form analytic d²I/dV² (`analytic_d2idv2_inelastic_at_bias`). Stored in the npz for reference; **the published figures do not use it** (they use raw numerical `np.diff²`). |
| `core/self_energy.py` | `gaussian_chi_projector` (Gaussian χ(z)); retained for `tests/test_gaussian_chi.py`. Not on the figure path. |
| `config/device_library.py` | `ZnO_MgZnO_symmetric` (paper device) + comparison stacks. |
| `config/molecular_database.py` | `Mol_A`/`Mol_B`/`Mol_AB` (100/180 meV, coupling 316.2 meV → D₀²=0.1 eV²) + VOC database. |
| `run/run_rank1_sweep.py` | Builds the tight-binding stack, χ profiles and phonon mode lists; loops `run_rank1_keldysh_single_bias` over bias; writes `iets_<device>_<mol>_<Vrange>_<T>K_rank1scba_<date>.npz`. |
| `run/make_paper_figs.py` | Generates the paper figures from the 51-pt npz. |
| `run/verify_reproduction.py` | Fast gate: re-runs a few bias points, checks I_R vs golden. |
| `run/reproduce_sispad.sh` | Full sweep (4 species @300K + Mol_A @10/77/150K) + figures. |

## 4. SCBA loop (per bias)

```
seed Σ_ph = 0  (FBA)
repeat (max_iter=100, tol=1e-4):
    Gᴿ = [(E+iη)I − H − U_bias − Σ_L − Σ_R + (i/2)(Σ^<_ph+Σ^>_ph)]⁻¹
    G^<,> = Gᴿ (Σ^<,>_contacts + Σ^<,>_ph) Gᴬ           # Keldysh
    Σ^<,>_ph ← Einstein-vibron update from G^<,>          # per mode, energy-shifted
    mix Σ_ph via Anderson/DIIS (depth 8, β=0.3)
FINAL PASS: recompute Gᴿ, G^<,> with converged Σ_ph
I_L, I_R = Meir-Wingreen trace at each contact
```

Symmetric bias μ_{L,R} = E_F ± V/2, E_F = 20 meV. Causal Σᴿ_ph =
−(i/2)(Σ^< − Σ^>). Contact in-scattering f₁Γ₁ + f₂Γ₂ included.

## 5. Parameters (paper run)

| | value |
|---|---|
| grid | a = 0.2 nm → Np = 135 |
| energy | dE = 2 meV, [−0.25, 0.5] eV → NE = 376 |
| bias | 0–800 mV, 51 points (16 mV) |
| T | 300 K (+ 10/77/150 K for the temperature figure) |
| mixing | Anderson/DIIS, M = 8, β = 0.3, max_iter = 100, tol = 1e-4 |
| couplings | D²_bulk = 0.001 eV²/site; molecular D₀² = 0.1 eV² |
| conservation | |I_L + I_R| / max|I| < 0.2 % |

## 6. d²I/dV² convention

The published figures use **raw second difference** `d2 = np.diff(I,2)/dV²` on the
51-pt grid (no smoothing), plotted at the interior bias points. On this grid that
gives the +125 / −170 µA/V² scale of `fig2_d2IdV2`. The analytic field
(`iets_analytic.py`, ±6 µA/V² at 300 K) and `np.gradient²` (±55) do **not** match
and are not used. `fig5_temp` (temperature) uses the same numerical method for
consistency — this differs from the published Fig 4b (which used the analytic
field); update the paper figure to match.

## 7. Known limitations / open items

- Patil 1D GaAs benchmark (Tier-1 γ/δ): `tests/test_rank1_keldysh.py` has 4
  failing cases — the rank-1 wrapper gives I_R ≈ 1.7e-11 A at the Patil NDR peak
  vs a 5.7e-10 A reference. Not on the ZnO production path.
- First Born vs converged SCBA adequacy at D²=0.1 eV², self-consistent Poisson,
  BenDaniel-Duke interfaces, vibron lifetime broadening, and DFT-derived VOC
  spectra are future work (were tracked in the now-archived SISPAD checklist).
