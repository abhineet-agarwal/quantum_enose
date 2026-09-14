# Suggested revisions — *IEEE Sensors 2026-10* draft

Source: `IEEE_Sensors_2026-10.pdf` (Agarwal & Ganguly, "Design of a Room-Temperature
Quantum Biomimetic Electronic Nose Using a CMOS BEOL Compatible Resonant
Tunneling Diode").

Items are ordered by severity. Each states the current text, what we now know,
and suggested replacement wording.

---

## A. Must fix — factual errors

### A1. Fig. 1 caption: "bound state E₁ = 149 meV" is not a resonance of this structure

**Current:** *"(i) V = 0: bound state E₁ = 149 meV; orange Gaussian χ(z) marks
molecular coupling. (ii) At resonance (V ≈ 560 mV): E₁ aligns with μ_L."*
The band-profile panel also annotates `E₁ = 149 meV`.

**Finding.** Ballistic transmission of the actual `ZnO_MgZnO_symmetric` stack
(Fisher–Lee `T = Tr[Γ_L G^R Γ_R G^A]`, 20 µeV grid) has **two** resonances:

| level | energy | FWHM | T_max |
|---|---|---|---|
| n = 1 | **78.74 meV** | 134 µeV | 0.995 |
| n = 2 | **295.82 meV** | broad | 1.000 |

149 meV is the *infinite-well* estimate `2t₀[1−cos(πa/L)]`, not an eigenvalue of
the finite 0.47 eV structure. An independent analytic finite-well solution gives
80 meV for the ground state, confirming n = 1 = 78.7 meV.

**This matters because the paper's own resonance condition is inconsistent with
it.** With `E₁ = 149 meV` and `E_F = 20 meV`, `V_res = 2(E₁−E_F) = 258 mV`,
which contradicts the reported 560 mV. The 560 mV NDR is the **n = 2**
resonance: `2(295.8 − 20) = 552 mV ≈ 560 mV` ✓.

**Suggested:** annotate the *two* levels, and state that transport at the NDR is
through the n = 2 resonance:

> (i) V = 0: quasi-bound states at E₁ = 79 meV and E₂ = 296 meV. The ground
> state is only 134 µeV wide and carries negligible current; the NDR at
> V ≈ 560 mV corresponds to alignment of μ_L with **E₂**
> (V_res = 2(E₂ − E_F) = 552 mV).

### A2. Peak-to-valley ratio is 12, not 6

**Current:** *"Baseline peak current is ∼69 nA with a peak-to-valley ratio (PVR)
of 6"* (§III-A), and *"Baseline peak ∼69 nA (PVR ≈6)"* (Fig. 3 caption).

**Finding.** From the archival run
(`results/sispad_scba_2026-04-14/...Baseline_0-800mV_300K...npz`):
peak 68.97 nA at 560 mV, valley 5.76 nA at 640 mV → **PVR = 11.98**.
The peak current (69 nA) is correct; only the PVR is wrong.

**Suggested:** replace "PVR of 6" with "PVR of ≈12" in both places, or state the
bias window used if a different valley definition was intended.

---

## B. Must fix — claims contradicted by the self-consistent work

### B1. §II "Self-consistent electrostatics" paragraph is wrong for *this* device

**Current:**
> *"In quantum transport studies of double-barrier RTDs, the dominant effect of
> adding a self-consistent Hartree coupling is found to be to renormalize the
> resonance bias position by tens of meV through charge accumulation in the well
> and the formation of a triangular emitter accumulation layer, while
> maintaining the qualitative I–V shape and the NDR feature [22]–[24]."*

**Findings.** We have since implemented self-consistent Poisson–NEGF and the
statement does not describe the simulated structure:

1. **No emitter accumulation layer forms in this device.** The simulated stack
   places n⁺ contacts (N_D = 10²⁵ m⁻³) *directly against the barrier*. Quantum
   reflection holds the electron density there at ≈ 2.35 × 10²⁴ m⁻³, so the
   uncompensated donor charge is **net positive (+7.65 × 10²⁴ m⁻³) —
   depletion**, not accumulation.
2. **The shift is not "tens of meV".** Against a space-charge-free flat-band
   reference the self-consistent peak moves **−96 mV** (528 → 432 mV);
   against the linear-drop reference, −128 mV. It is also *backward*, opposite
   to the direction the cited accumulation picture predicts.
3. **The cited references describe structures with spacers.** Akkala [24]
   Table 3.1 places **10 nm lightly-doped (10¹⁵ cm⁻³) spacers** between each
   doped lead and its barrier, and states their role explicitly: to keep ionized
   impurities away from the coherent region (§2.2), and that the emitter *notch*
   hosting the accumulation layer forms **in the emitter spacer** (§3.7).
   Adding those spacers to our stack flips the pre-barrier space charge to
   **−2.66 × 10²⁴ m⁻³ (accumulation)**, restoring the literature behaviour, with
   the electron density essentially unchanged (2.35 → 2.66 × 10²⁴, 13 %).

**The sign of the space charge is set by the donor placement, not by the
electrons.** So the paragraph is defensible only for a spacer-containing stack.

**Suggested — two options.**

*Option 1 (preferred): add the spacers to the device.* This makes the device
faithful to the cited literature, removes the contradiction, and is what a
fabricated structure would use anyway. Requires regenerating Figs. 3–4.

*Option 2: restate honestly for the present geometry:*

> A self-consistent Hartree treatment renormalizes the resonance bias position.
> For structures with lightly-doped spacers, emitter accumulation screens the
> barrier field and displaces the resonance to higher bias [22]–[24]. The stack
> simulated here has no spacer: ionized donors sit adjacent to the barrier,
> where quantum reflection suppresses the carrier density, giving a net positive
> space charge that concentrates the field and displaces the resonance to
> *lower* bias by ≈100 mV. Because IETS discrimination relies on the *relative*
> spacing of phonon-induced peaks, this bias renormalization shifts the
> measurement window without altering the molecular fingerprint structure.

The final sentence of the existing paragraph (relative spacing is preserved)
remains valid and should be kept in either option.

### B2. E_F = 20 meV is inconsistent with the quoted doping

**Current:** §II states contacts with `N_D = 10²⁵ m⁻³` and semi-infinite leads
with `E_F = 20 meV`.

**Finding.** Bulk charge neutrality `n(E_F) = N_D` for ZnO
(`N_c = 2[m*kT/2πħ²]^{3/2} = 3.718 × 10²⁴ m⁻³`, m* = 0.28, 300 K, full
ionization) requires **E_F = 49.3 meV**. E_F = 20 meV supports only
`5.04 × 10²⁴ m⁻³`, i.e. **half** the stated doping. Verified in-model: on a
converged grid the transverse-integrated NEGF density in the deep contact at
E_F = 49.3 meV is 9.98 × 10²⁴ m⁻³ = 0.998 N_D.

**Impact on the present results: none.** The calculations are not
self-consistent, so N_D never enters the transport; E_F is purely a contact
occupation parameter. Ballistic check on a resolved grid: the dominant peak
stays at 560 mV for both E_F values (peak current rises 110 → 162 nA).

**Suggested:** either quote `E_F = 49.3 meV` (consistent), or drop `N_D` from
§II and describe E_F as the contact occupation parameter. Leaving both numbers
as-is invites a reviewer to notice the factor-of-two mismatch.

---

## C. Should address — unresolved feature and an overstated claim

### C1. The n = 1 resonance is not resolved by the energy grid

The ground resonance (78.74 meV, **FWHM 134 µeV**) is ~15× narrower than the
2 meV energy grid and is therefore sampled erratically. It appears in the
archival data as the small **1.46 nA shoulder near 144 mV** visible in
Fig. 3(a), which is a numerical remnant rather than a converged feature; the
converged ballistic value is a 5–8 nA peak near 117–150 mV.

**Transport conclusions are unaffected:** on a 10 µeV grid the dominant n = 2
peak reproduces to **0.3 %** (109.3 → 109.62 nA), because a 134 µeV resonance
contributes negligible `∫T dE` even at unity transmission.

**Suggested:** add one sentence to §II (Convergence) —

> The n = 1 resonance at 79 meV has a 134 µeV linewidth and is not resolved on
> the 2 meV transport grid; it contributes < 5 % of the peak current, and the
> n = 2 transport peak is converged to 0.3 % against a 10 µeV reference grid.

— and either remove or explicitly label the 144 mV shoulder in Fig. 3(a).

### C2. "The 100 mV IETS peak ... remains clearly resolved at 300 K" is overstated

**Current (§III-C):** *"Critically, the 100 mV IETS peak in d²I/dV² narrows at
low T but remains clearly resolved at 300 K. This is the central practical claim
of the device."*

**Finding.** In the archival Mol_A data the largest |d²I/dV²| anywhere in the
60–160 mV window is **1.95 µA/V² at 300 K**, against a global maximum of
**58.4 µA/V²** — i.e. **3.3 % of full scale** (10 K: 2.99 vs 177.9, 1.7 %). On
the axes of Fig. 4(b) it is not visible.

Since this is described as the paper's central practical claim, it needs support
that the current figure does not provide.

**Suggested:** add a zoomed inset to Fig. 4(b) over 0–250 mV with its own
ordinate scale, and restate quantitatively, e.g.

> The molecular feature near ħω/e = 100 mV persists to 300 K with amplitude
> ≈ 2 µA/V², ~3 % of the resonance-derived background; it is resolved in the
> discrimination metric ΔD (Fig. 3c) where the elastic background cancels.

Note the discrimination metric ΔD, not the raw d²I/dV², is where the molecular
signature is actually separable — that is the stronger framing and is already
supported by Fig. 3(c).

---

## D. Worth raising — device feasibility

### D1. The analyte cannot reach the barrier in a planar vertical stack

χ(z) is centred on the emitter barrier, which in the simulated stack sits
**≈11 nm below the top surface** (and 46 nm in the longer Poisson geometry),
buried under the top contact. A physisorbed odorant cannot occupy that position.
The abstract's framing ("molecules adsorbed in or near a tunnel junction") is
standard for planar MIM IETS but not realizable for a buried RTD barrier.

**Suggested:** add one or two sentences to §IV (Conclusion / future work)
acknowledging the access problem and naming a route. A **nanopore array** is the
natural one: a pore etched through the stack exposes the barrier at the pore
wall, placing the molecule at exactly the z where χ(z) peaks, so the present
coupling model carries over unchanged.

Quantitatively, the signal is set by the fraction of current passing within
coupling range of a pore wall (≈3 nm), for 5 nm pores on a 10 × 10 µm mesa:

| pitch | pores | modulated fraction of I | pore area fraction |
|---|---|---|---|
| 100 nm | 10⁴ | 0.5 % | 0.002 |
| 50 nm | 4 × 10⁴ | 1.9 % | 0.008 |
| **20 nm** | 2.5 × 10⁵ | **12 %** | 0.05 |
| 10 nm | 10⁶ | 47 % | 0.20 |

A single pore modulates ~5 × 10⁻⁷ of the device current and is unmeasurable; a
~20 nm-pitch array modulates ~12 % while removing only 5 % of the mesa area.
This also means the transverse area is no longer a pure prefactor: it splits
into modulated and unmodulated fractions, and the unmodulated part dilutes ΔD.

---

## E. Optional — figure presentation

* Figs. 3 and 4 have been regenerated as **line curves without point markers**
  (`results/paper-figures/fig3_NO_poisson.png`,
  `fig4_temperature_NO_poisson.png`), which reads better at column width.
* Fig. 3(c): ΔD is dominated by structure above 600 mV where the SCBA valley is
  poorly converged. Consider restricting the abscissa to 0–600 mV, or marking
  the > 600 mV region as unconverged.
* Fig. 3(a): consider a log ordinate so the low-bias region and the 144 mV
  shoulder are legible alongside the 69 nA peak.

---

## Summary table

| # | Item | Severity | Affects results? |
|---|---|---|---|
| A1 | E₁ = 149 meV wrong; NDR is the n = 2 state at 296 meV | **high** | no — labelling/interpretation |
| A2 | PVR is 12, not 6 | **high** | no — quoted number |
| B1 | Emitter accumulation claim false without spacers; shift is −96 mV backward | **high** | yes — §II text, and device choice |
| B2 | E_F = 20 meV inconsistent with N_D = 10²⁵ m⁻³ (needs 49.3 meV) | medium | no — non-self-consistent |
| C1 | n = 1 resonance unresolved (134 µeV vs 2 meV grid) | medium | no — 0.3 % on the dominant peak |
| C2 | "100 mV peak clearly resolved" is 3 % of full scale | medium | claim needs restating |
| D1 | Analyte cannot reach a buried barrier; nanopore array route | medium | concept-level |
| E | Figure presentation | low | no |

**None of A1, A2, B2, C1 change any computed I–V.** The substantive scientific
change is **B1**, which is a statement about physics the paper asserts but the
simulated geometry does not exhibit.
