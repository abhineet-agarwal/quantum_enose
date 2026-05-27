# SISPAD 2026 — primary device stack & transverse-area decision

**Date:** 2026-04-07
**Status:** Locked for the SISPAD submission (deadline 2026-04-22, extended
from 2026-04-10).
**Related docs:** `docs/METHOD_DERIVATION.md` (rank-1 projected SCBA),
`docs/REPRODUCE_SISPAD.md`. (The earlier abstract draft + checklist are now under
`archive/docs/`.)

> **Note (2026-05-27):** the decision below — ZnO/Mg₀.₃Zn₀.₇O symmetric stack,
> 10 µm pixel as a Datta prefactor — is current and implemented in
> `config/device_library.py` (`ZnO_MgZnO_symmetric`). Mentions of `run_iets.py`
> and `launch_cloud_v*.sh` are historical; those files are now under `archive/`.
**Supersedes:** an earlier (now-deleted) `TRANSVERSE_SIZE_DECISION.md` that was
written post-compaction around the *wrong* stack. The correct decision is
documented here from scratch.

---

## 1. What we're locking in

1. **Primary material stack:** **ZnO / Mg₀.₃Zn₀.₇O asymmetric RTD**
   (3 nm emitter barrier / 3 nm ZnO well / 1.5 nm collector barrier, n⁺ ZnO
   contacts). Device library entry: `ZnO_MgZnO_asymmetric`.
2. **Comparison stack:** **In₂O₃ / κ-Ga₂O₃ asymmetric RTD** (same geometry),
   used for Section III/IV side-by-side plots showing the numerical results are
   insensitive to the choice of oxide pair at this CBO. In-flight v3 cloud run
   on this stack stays at the legacy **1 µm × 1 µm** transverse box and is the
   comparison data.
3. **Primary-stack transverse area:** **A = 10 µm × 10 µm = 10⁻⁷ cm²**
   (sensor-pixel area entering Datta's prefactor — see §3).
4. **Transverse mode handling:** **"Fix A"** — replace the brute-force discrete
   mode sum with Datta's **analytic** transverse integration. The 10 × 10 µm²
   number is an *A prefactor* in the current expression, **not** a hard-wall
   quantization box that the discrete modes need to resolve.
5. **Methodology:** the **rank-1 projected SCBA** derived in
   `docs/METHOD_DERIVATION.md` replaces the current "4 inelastic + 5 coherent"
   hybrid solver. New file: `core/scba_solver_rank1.py`; new CLI flag:
   `--solver rank1`.

---

## 2. Why ZnO/Mg₀.₃Zn₀.₇O over In₂O₃/κ-Ga₂O₃

The κ-Ga₂O₃ stack and the MgZnO stack have essentially identical physics
numbers at this CBO:

| Criterion                    | In₂O₃ / κ-Ga₂O₃    | **ZnO / Mg₀.₃Zn₀.₇O** |
|------------------------------|---------------------|------------------------|
| CBO (eV)                     | 0.45 (1 XPS source) | **0.47 (multi-source, tunable as 1.57·x_Mg)** |
| Phase stability of barrier   | metastable κ-phase  | **stable wurtzite for x_Mg ≤ 0.3** |
| ALD process maturity         | emerging (2020+)    | **decades-old**        |
| Deposition T                 | < 400 °C            | **< 350 °C**           |
| *Experimentally demonstrated RTD with NDR* | **no**  | **yes** (Tampo et al., IEEE NANO 2019) |
| Effective mass (well)        | 0.30 m₀ (In₂O₃)     | 0.28 m₀ (ZnO)          |
| LO phonon (well)             | ~70 meV             | ~72 meV                |
| Tunable CBO                  | no                  | **yes** (via x_Mg)     |

So MgZnO wins on every axis a reviewer cares about (process maturity,
experimental grounding, tunability, phase stability), and matches on all the
physics axes that affect the simulation results. Swapping the primary stack is
a cosmetic change in the solver (CBO 0.45 → 0.47, m* 0.30 → 0.28, ℏω 70 →
72 meV) but a **substantive** change in the paper's credibility.

### Reviewer flags the swap addresses
1. *κ-phase metastability* — κ-Ga₂O₃ is orthorhombic and not thermodynamically
   stable; narrow growth window; long-term device-stack stability unknown.
   MgZnO is stable wurtzite at x_Mg ≤ 0.3.
2. *Single-source CBO* — the 0.45 eV number for In₂O₃/κ-Ga₂O₃ comes from one
   XPS paper (ACS Appl. Electron. Mater. 2020, 10.1021/acsaelm.0c00947) with
   no cross-check. MgZnO CBO is measured by multiple groups across the last
   decade and parameterised as 1.57·x_Mg.
3. *No demonstrated RTD* in the κ-Ga₂O₃ system (the ACS paper showed a 2DEG,
   not a double-barrier RTD). ZnO/Mg₀.₃₃Zn₀.₆₇O has a published working RTD
   with NDR.
4. *Native defect landscape* — β-Ga₂O₃ has Fe and V_Ga deep levels; κ-phase
   defect chemistry is open. ZnO defect physics is mature.

### What we keep from the In₂O₃/κ-Ga₂O₃ work
- The in-flight v3 cloud run (asymmetric In₂O₃/κ-Ga₂O₃, 81 V-pts, 1 µm) is
  **not killed**. It becomes Tier-2 comparison data for Section III of the
  paper: "the result is insensitive to the particular oxide pair at this CBO".
- The In₂O₃ device entries stay in `config/device_library.py`, tagged as
  "comparison stack" in their notes.

---

## 3. Transverse area: 10 µm² is a sensor-pixel, not a hard-wall box

### The old approach (wrong for large areas)

`run/run_iets.py` currently treats L_y × L_z as a hard-wall quantization box
and sums over a finite grid of discrete transverse modes (n_max = m_max = 3).
At L = 1 µm the transverse subband spacing is

  Δε = π²ℏ²/(2m*L²) ≈ 4 µeV   (m* = 0.30 m₀, L = 1 µm)

which is already ~6000× smaller than kT = 25 meV, so the discrete modes
approximate a continuum. At L = 10 µm the spacing drops another 100× (≈ 40
neV). Trying to resolve that continuum with 3×3 modes is nonsense — the
discrete mode basis is wildly inadequate at 10 µm, and brute-force raising
n_max to a sensible convergence is both wasteful and beside the point.

### Fix A — Datta's analytic transverse integration

The transverse integral over the continuum is analytic for a separable
Hamiltonian with parabolic subbands. The current through each longitudinal
resonance is (Datta, *Quantum Transport*, 2005, Eq. 9.3.14):

  I = (2e / h) · A · (m* kT / 2π ℏ²)
       · ∫ dE_z T(E_z) · { log(1 + exp((E_F − E_z + V/2) / kT))
                         − log(1 + exp((E_F − E_z − V/2) / kT)) }

where `A = L_y × L_z` enters only as a multiplicative prefactor in the current
density. There are **no discrete transverse modes in this expression**. The
transverse degree of freedom has been integrated out analytically; A sets the
device area and nothing else.

- The **absolute current** is proportional to A, so 10 µm² gives a nA-range
  peak current directly comparable to published oxide RTD I–V traces.
- The **elastic I(V) shape** is independent of A (up to the trivial scaling).
- The **inelastic fingerprint** (d²I/dV²) depends on the local molecular
  coupling density, which in the rank-1 derivation is controlled by σ_mol, not
  by L. So the A prefactor scales the overall signal without smearing the
  peaks.
- This means the 10 µm² choice is a **one-number** decision affecting an
  overall prefactor, not a knob whose change demands re-convergence of any
  discrete mode basis. That's why it's cheap, and it's why we were wrong to
  frame the earlier "1 µm → 10 µm" decision as changing the transverse box.

### Interaction with the rank-1 molecular coupling

The molecular form factor g(r) in `METHOD_DERIVATION.md` Eq. (1) is an
isotropic Gaussian with σ_mol = 3 Å. At the point molecule limit it sets a
rank-1 coupling operator in the transverse continuum rather than in the
discrete subband basis. The projected Green's function G̃^R = ⟨u|G^R|u⟩ uses

  u(r_⊥) = g(r_⊥ − r_⊥,mol) / ‖g‖₂

which is a fixed continuum function **independent of A**. The rank-1 derivation
carries over verbatim from the hard-wall case: the only thing that changes is
that the sum Σ_nm |u_nm|² becomes an integral ∫d²r_⊥ |u(r_⊥)|² = 1 (normalised
to unity, by construction). So at the level of G̃^R and σ̃^R, the 10 µm² area
**drops out entirely** — A re-enters only in the Datta prefactor for the
final current (elastic and inelastic alike).

This is the clean justification: sensor-pixel area is a prefactor, rank-1
projected SCBA is a molecular property, and the two are decoupled.

---

## 4. Why 10 µm² specifically (and not 100 nm² or 1 mm²)

- Typical e-nose **sensor pixel** in a printed/ALD gas-sensor array is 5–50 µm
  on a side. 10 µm is the midpoint and is what published oxide RTD mesas
  report.
- 100 nm × 100 nm would fit a CMOS-scale logic device but breaks the "pixel
  in a sensor array" narrative that the paper is being built around.
- 1 mm² would push the peak current into the µA range — achievable, but the
  molecular signal is drowned by the elastic background.
- The choice is a **pixel area**, not a device-physics parameter. It will be
  stated that way explicitly in Section II.

---

## 5. What this decision does and doesn't change

**Changes:**
- `config/device_library.py`: `ZnO_MgZnO_symmetric` updated to 10 µm notes;
  new `ZnO_MgZnO_asymmetric` entry added at 10 µm, Patil-style 3 nm / 1.5 nm
  barriers. In₂O₃ entries reverted to 1 µm and tagged "comparison stack".
- `run/launch_cloud_v3.sh` reverted (the in-flight In₂O₃ run stays at 1 µm).
- **New `run/launch_cloud_v4.sh`** targeting `ZnO_MgZnO_asymmetric`, 10 µm,
  and `--solver rank1` once that solver is implemented.
- **New `docs/SISPAD_CHECKLIST.md`** with the Tier-1 + Tier-2 open items.
- SISPAD abstract rewritten to lead with ZnO/MgZnO, rank-1 theorem, and Fix A.
- `docs/METHOD_DERIVATION.md` — unchanged; already describes the rank-1
  derivation. This doc just references it.

**Does not change:**
- The rank-1 derivation in `METHOD_DERIVATION.md` (still valid verbatim —
  transverse subband structure doesn't matter to the rank-1 argument).
- The six "fixes" from the overnight session (#1 bias convention, #2 D²
  calibration, #3 contact in-scattering, #4 causal Σ^R sign, #5 Landauer
  shortcut, #6 longitudinal vibron footprint) — all stay in, though #6 will
  be replaced by a proper σ_mol = 3 Å Gaussian in the rank-1 solver.
- The v3 in-flight In₂O₃ study — it runs to completion untouched and becomes
  Tier-2 comparison data.

---

## 6. Open items (see `docs/SISPAD_CHECKLIST.md` for full list)

The transverse-area change is cheap. The *methodology* changes that back it up
are not, and none are done yet:

1. **Implement `core/scba_solver_rank1.py`** per METHOD_DERIVATION.md §3–5.
2. **Unit test** `tests/test_rank1_vs_full.py` — 2-mode toy, verify rank-1
   agrees with brute-force coupled SCBA to machine precision.
3. **Replace Fix 6** with a σ_mol = 3 Å Gaussian χ(z). Recalibrate D₀ if signal
   drops.
4. **Analytic d²I/dV²** from differentiating Meir-Wingreen — replaces the
   double `np.gradient` that's currently dominating the IETS plots with
   numerical noise. (Tier-1 item α.)
5. **Patil 1D benchmark** — reproduce `papers/rtd2modes_1d.m` I(V) and
   d²I/dV² in our 1D solver. (Tier-1 item γ.)
6. **Current conservation** check. (Tier-1 item δ.)
7. **FBA vs SCBA** — rename or actually converge SCBA on one molecule.
   (Tier-1 item β.)

Only after these land does `launch_cloud_v4.sh` get activated.
