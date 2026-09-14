# Self-Consistent Electrostatics and the Resonance Peak Position

**Device:** ZnO/Mg₀.₃Zn₀.₇O double-barrier RTD, `ZnO_MgZnO_symmetric_long`
(30 nm n⁺ ZnO / 2 nm MgZnO / 3 nm ZnO / 2 nm MgZnO / 30 nm n⁺ ZnO),
CBO = 0.47 eV, m\* = 0.28 m₀, T = 300 K.
**Method:** rank-1 projected SCBA NEGF, optionally coupled self-consistently to a
1-D Poisson solver (Newton–Raphson, Dirichlet emitter / Neumann collector,
transverse-integrated physical density).

---

> **RESOLVED — see §0.** The backward shift reported below is an artifact of a
> missing device layer (the lightly-doped spacers of Akkala Table 3.1), not a
> physical disagreement with the emitter-accumulation picture. Sections 1–5
> describe the no-spacer structure and remain valid *for that structure*; the
> physical interpretation in §4 is superseded.

## 0. Root cause: the missing spacer layers

The reference structure (Akkala, Table 3.1) places a **10 nm lightly-doped
(10¹⁵ cm⁻³) spacer between each doped lead and its barrier**. The thesis states
their purpose explicitly: to keep ionized impurities away from the coherent
well region (§2.2), and to host the **emitter notch** — the triangular well in
which the accumulation layer's quasi-bound states form, captured by including
the emitter spacer in the NEGF quantum region (§3.7).

`ZnO_MgZnO_symmetric_long` omits them: donors sit directly against the barrier.
Evaluating the ballistic density at V = 0.432 V, E_F = 49.3 meV, on a converged
grid (dE = 0.2 meV), for the 8 nm immediately in front of the emitter barrier:

| structure | n (m⁻³) | N_D (m⁻³) | n/N_D | net charge N_D − n | state |
|---|---|---|---|---|---|
| no spacer (as simulated) | 2.346 × 10²⁴ | 1.00 × 10²⁵ | 0.235 | **+7.65 × 10²⁴** | depletion |
| with spacer (as in Table 3.1) | 2.660 × 10²⁴ | 1.00 × 10²¹ | 2660 | **−2.66 × 10²⁴** | accumulation |

**The electron density is essentially unchanged (2.35 vs 2.66 × 10²⁴, 13 %).**
What reverses the sign of the space charge is the *donor* density. Quantum
reflection from the barrier holds n near the barrier at ≈ 2.4 × 10²⁴ m⁻³; when
10²⁵ m⁻³ of donors are placed there, the uncompensated remainder is a large
*positive* charge, which reads as "emitter depletion." In the reference
structure the spacer carries essentially no donors, so the same electrons
constitute *negative* charge — the conventional accumulation layer.

The backward peak shift documented below therefore follows from placing ionized
donors in the barrier-reflection region, not from any physical difference in
emitter behaviour. Production runs should use `ZnO_MgZnO_symmetric_spacer`.

---

## 1. Summary

Turning on self-consistent electrostatics moves the resonant-tunnelling current
peak **to lower bias**. Against a space-charge-free flat-band reference the shift
is **−96 mV** (528 → 432 mV); against the naive linear-drop reference it is
−128 mV (560 → 432 mV). The shift arises because the self-consistent potential
places **more of the applied bias between the emitter and the quantum well** than
either fixed-profile model assumes, which lowers the resonance level E₁ and
brings it into alignment with the emitter chemical potential at a smaller applied
bias.

This is *opposite in sign* to the emitter-accumulation picture usually quoted for
RTDs (Akkala, §3.5), in which space charge screens the field and delays
resonance to *higher* bias. Section 4 shows the two pictures are not in conflict:
they describe different operating regimes, and this device is driven deep into
the regime where the accumulation layer cannot survive.

---

## 1a. Choice of doping and Fermi level

**Doping: N_D = 1 × 10¹⁹ cm⁻³ in both 30 nm leads, undoped elsewhere.**

The reference structure (Akkala, Table 3.1) is a GaAs/AlGaAs DBRTD with
10¹⁸ cm⁻³ leads. That number cannot be transferred to ZnO directly: the
conduction-band effective mass differs by a factor of four
(m\* = 0.28 m₀ for ZnO vs 0.067 m₀ for GaAs), so the effective density of states
`N_c = 2[m* kT/(2πħ²)]^{3/2}` is **8.5× larger** in ZnO. What must be preserved
is the *degeneracy* of the emitter — its ability to act as a reservoir — not the
absolute doping:

| system | m\*/m₀ | N_c (m⁻³) | n (m⁻³) | n/N_c | E_F |
|---|---|---|---|---|---|
| Akkala GaAs lead, 10¹⁸ cm⁻³ | 0.067 | 4.35 × 10²³ | 1 × 10²⁴ | 2.30 | +41.9 meV |
| **this work, ZnO lead, 10¹⁹ cm⁻³** | 0.280 | 3.72 × 10²⁴ | 1 × 10²⁵ | **2.69** | **+49.3 meV** |
| ZnO at 10¹⁸ cm⁻³ (literal transfer) | 0.280 | 3.72 × 10²⁴ | 1 × 10²⁴ | 0.27 | −31.5 meV |

10¹⁹ cm⁻³ in ZnO reproduces the reference degeneracy (n/N_c = 2.69 vs 2.30;
E_F = +49.3 vs +41.9 meV). A literal transfer of 10¹⁸ cm⁻³ would place E_F
31.5 meV *below* the band edge, giving a non-degenerate emitter that cannot
supply the resonant channel.

**Known deviation from the reference structure.** Table 3.1 includes 10 nm
lightly-doped (10¹⁵ cm⁻³) **spacers** between each doped lead and its barrier;
the present structure has none, and donors sit directly against the barrier.
Spacers separate the ionized-donor region from the tunnel barrier and are known
to control the emitter accumulation layer, so this is a material structural
difference and a candidate contributor to the emitter charge-state behaviour
reported in §3.2. The barrier/well thicknesses also differ (2/3/2 nm here vs
5/5/5 nm in Table 3.1).

**Fermi level: E_F = 49.3 meV above E_c**, fixed by bulk charge neutrality
`n(E_F) = N_D` with a Fermi–Dirac integral and full donor ionization,

    n(E_F) = N_c (2/√π) F_{1/2}(E_F/kT),   N_c = 3.718 × 10²⁴ m⁻³

This is confirmed *within the model*: on a converged energy grid
(E_min = −0.02 eV, dE = 0.05 meV) the transverse-integrated NEGF density in the
deep contact at E_F = 49.3 meV is 9.98 × 10²⁴ m⁻³, i.e. **0.998 × N_D**.

Earlier runs used E_F = 20 meV, which supports only half the specified doping;
see §6.2. An apparent discrepancy suggesting E_F ≈ 63 meV was traced to the
energy-grid convergence error of §6.4 and is not physical.

---

## 2. The resonance condition

Resonant current flows while the quasi-bound level lies inside the *occupied*
band of the emitter. In the symmetric-bias convention the emitter band edge sits
at `V/2` and its chemical potential at `μ_L = E_F + V/2`, so the level
`E₁ + U_well(V)` conducts for

    V/2  <  E₁ + U_well  <  E_F + V/2                                      (1)

equivalently, for a bias profile antisymmetric about the device centre
(`U_well = 0`, true for both the linear-drop and flat-band references),

    2(E₁ − E_F)  <  V  <  2E₁                                              (2)

Ballistic transmission at V = 0 gives **E₁ = 296 meV** for the n = 2 resonance,
so the conducting window is 552–592 mV for `E_F = 20 meV`. Note that Eq. (2)
brackets the window; the current *peak* lies inside it and is not in general at
either edge. Because the window is only 40 meV wide at this E_F, its lower edge
(552 mV) happens to lie close to the observed peak (560 mV) — a coincidence of
the narrow window, not an identity. Raising E_F to 49.3 meV widens the window to
493–592 mV and the peak remains at 560 mV (see §6.2), confirming that the lower
edge is an onset, not the peak.

Self-consistency breaks the antisymmetry. With `U_well ≠ 0` the whole window in
Eq. (2) is displaced rigidly:

    2(E₁ + U_well − E_F)  <  V  <  2(E₁ + U_well)                          (3)

so a well pulled *down* by |U_well| shifts the entire conducting window — and
with it the peak — *down* by 2|U_well|.

---

## 3. What the self-consistent solution actually does

### 3.1 Redistribution of the bias drop

Partitioning the applied bias between the emitter contact, the active region
(barriers + well) and the collector contact, at V = 432 mV:

| region | linear reference | self-consistent |
|---|---|---|
| emitter contact (30 nm) | 194 meV | 193 meV |
| **active region (7 nm)** | **44 meV** | **214 meV** |
| collector contact (30 nm) | 194 meV | 12 meV |

The linear model drops ~45 % of the bias across each *doped* contact, which a
degenerate semiconductor cannot sustain — it screens. The self-consistent
solution correctly confines the drop to the undoped active region. The
consequence for Eq. (3) is a well displaced below the antisymmetric value by

    ΔU_well = −71 meV        (at V = 432 mV, Neumann collector)

Predicted shift `2 × 71 = 142 mV`, measured **128 mV** relative to the linear
reference. The residual is accounted for by the flat-band reference itself
sitting 32 mV below the linear one (528 vs 560 mV).

### 3.2 Emitter charge state

The pre-barrier emitter charge state, measured as `n/N_D` averaged over the 8 nm
immediately in front of the emitter barrier:

| V (mV) | n/N_D | state |
|---|---|---|
| 0 | 1.105 | accumulation |
| 16 | 1.146 | accumulation |
| **32** | **0.925** | **crossover** |
| 80 | 0.581 | depletion |
| 352 | 0.516 | depletion |
| 432 | 0.399 | depletion (at the current peak) |

At equilibrium the emitter **accumulates**, exactly as the conventional picture
requires. The accumulation layer is weak (10–15 % above N_D) and collapses at
**V/2 ≈ 16 meV**, i.e. when the bias per contact becomes comparable to
`E_F = 20 meV` and to `kT = 25.9 meV`. Beyond that the emitter depletes
monotonically. The resulting positive space charge concentrates the field in
front of the emitter barrier, which is the microscopic origin of the
redistribution in §3.1.

---

## 4. Reconciliation with the accumulation picture

The disagreement is one of **operating regime**, traceable to the device's
resonance structure.

Ballistic T(E) at zero bias resolves two resonances:

| level | energy | linewidth (FWHM) | T_max | peak current |
|---|---|---|---|---|
| n = 1 | 78.74 meV | **134 µeV** | 0.9999 | 5.2 nA |
| n = 2 | 296 meV | broad | 0.990 | **109.6 nA** |

The ground state is confined far below the 0.47 eV barriers, so its tunnelling
rate — and hence its linewidth — is minute. Although it transmits perfectly on
resonance (T = 0.9999, confirming the structure is symmetric and the level is
genuine), a 134 µeV resonance contributes negligible `∫T dE`. Transport is
therefore dominated by n = 2, and the device operates at **V_peak ≈ 560 mV**,
where

    V/2 ≈ 280 meV ≈ 14 × E_F

A conventional RTD operates on its *ground* resonance at V/2 of a few × E_F,
comfortably within the regime where a degenerate emitter sustains an accumulation
layer. This device, by contrast, is driven an order of magnitude past that point,
where — as §3.2 shows directly — the accumulation layer no longer exists.

**The accumulation and depletion pictures are therefore both correct, in their
respective regimes.** The present device does accumulate at low bias and would
exhibit the conventional forward shift there; it depletes at its actual operating
point and exhibits a backward shift.

---

## 5. Numerical validation

* **Transmission.** A tridiagonal (banded) evaluation of
  `T = Γ_L Γ_R |G^R_{1,N}|²` reproduces the dense-matrix solver to
  7 × 10⁻¹⁴ relative error in T, and the ballistic current to 5 significant
  figures (1.040366 × 10⁻⁸ A, both).
* **Energy grid.** The 2 meV production grid does not resolve the 134 µeV n = 1
  resonance, but reproduces the dominant n = 2 peak to **0.3 %**
  (109.3 vs 109.62 nA on a 10 µeV grid). Transport conclusions are unaffected;
  the omitted n = 1 feature is a 5 nA shoulder near 150 mV.
* **Poisson sign conventions** verified by four controlled tests (electron sign,
  positive-charge sign, absolute units against an analytic parabola, and the
  U_well → E₁ coupling, found to be exactly 1:1).
* **Charge neutrality** holds in the deep contacts (`n/N_D = 1.005–1.037`).
* **Scattering** is irrelevant to the effect: switching the LO-phonon coupling
  off entirely changes `n/N_D` from 0.434 to 0.437 and ΔU_well from −89.8 to
  −86.8 meV.

---

## 6. Open items

These must be settled before the quantitative shift is quoted:

1. **Convergence-path dependence.** Self-consistent solutions started from
   different initial potentials reach different endpoints at the same bias
   (9.3 vs 90.7 nA at 432 mV; 21.1 vs 66.0 nA at 352 mV). The Poisson
   convergence threshold is 5 meV while the two endpoints differ by ~19 meV in
   U_well, so this is plausibly a tolerance artifact rather than genuine
   bistability. **Peak positions and shift magnitudes are provisional until
   this is resolved.** The qualitative mechanism of §3 is not affected — it
   rests on a smooth trend across 40+ bias points.
2. **Fermi level vs doping.** Doping is uniform per layer: leads
   `N_D = 1 × 10²⁵ m⁻³ (10¹⁹ cm⁻³)` over 30 nm each, barriers and well undoped,
   with no spacer between the doped lead and the barrier. The runs use
   `E_F = 20 meV`, but bulk charge neutrality
   (`n = N_c (2/√π) F_½(E_F/kT)`, `N_c = 3.718 × 10²⁴ m⁻³`, full ionization)
   requires **E_F = 49.3 meV**. The physical-density anchoring masks the
   inconsistency by rescaling the donor density down (to 3.795 × 10²⁴ m⁻³ at
   E_F = 20 meV — a factor 0.38).

   *Caveat:* the model's own contact density is consistently ≈ 0.77 × the
   analytic bulk value (3.795 vs 5.036 × 10²⁴ at 20 meV; 7.79 vs 10.0 × 10²⁴ at
   49.3 meV), so the *model-consistent* neutrality point is **E_F ≈ 63 meV**,
   not 49.3 meV. The origin of the 23 % shortfall (transverse-integration
   approximation, finite device length, or tight-binding vs parabolic DOS) has
   not been established and should be checked before E_F is changed in
   production.

   Impact on transport, quantified ballistically on the resolved grid:

   | E_F | dominant peak | peak current | onset 2(E₁−E_F) |
   |---|---|---|---|
   | 20.0 meV | 560 mV | 109.6 nA | 552 mV |
   | 49.3 meV | **560 mV** | **162.1 nA** | 493 mV |

   The **peak position is unchanged**; only the current rises (+48 %), because
   the wider occupied emitter band admits more carriers. A single-bias
   self-consistent check at the corrected E_F still shows emitter depletion
   (n/N_D = 0.450) and a *larger* well displacement (−115 meV), so the sign of
   the shift is not expected to change — if anything it grows.
3. **Boundary conditions.** The collector uses a zero-field (Neumann) condition,
   which suppresses the collector drop to 12 meV. A symmetric Dirichlet
   comparison should be reported alongside.

4. **Energy-grid convergence of the density (most serious).** The production
   grid (dE = 2 meV) under-resolves the density integral by ~25 %. Deep-contact
   density at E_F = 20 meV, against the analytic bulk value 5.036 × 10²⁴ m⁻³:

   | dE | n_deep | % of bulk |
   |---|---|---|
   | 2.0 meV (production) | 3.795 × 10²⁴ | 75.4 % |
   | 1.0 meV | 4.173 × 10²⁴ | 82.9 % |
   | 0.5 meV | 4.655 × 10²⁴ | 92.4 % |
   | 0.2 meV | 4.931 × 10²⁴ | 97.9 % |
   | 0.05 meV | 5.025 × 10²⁴ | 99.8 % |

   E_min is irrelevant (−0.60 and −0.10 eV give identical results at fixed dE),
   so this is resolution at the band edge, where the 1-D DOS carries a 1/√E van
   Hove singularity.

   **This affects the density but not the current.** The current weight
   `f_L − f_R` vanishes at the band edge and is finite only inside the transport
   window, so the current is converged to 0.3 % on the coarse grid (§5). The
   density weight `ln(1+e^{(μ−E)/kT})` is *maximal* at the band edge, directly
   on the singularity.

   Donor anchoring rescales N_D to the same under-resolved contact density, so
   deep-contact neutrality is preserved by construction; but the error is not
   uniform in z (it originates in a contact band-edge feature, while the active
   region has resonant structure), so reported `n/N_D` ratios carry an unknown
   bias. **Production Poisson runs should evaluate the density on dE ≤ 0.2 meV,
   or a grid refined near the band edge.**

---

## 7. Figures

| file | content |
|---|---|
| `thick_pois_vs_nopois_fig3.png` | I–V, dI/dV, d²I/dV² with and without Poisson |
| `thick_band_density_peak.png` | potential and carrier/donor profiles at the peak |
| `emitter_charge_vs_bias.png` | emitter charge state vs bias (accumulation → depletion) |
| `thick_nopois_baseline_4panel.png` | I–V, conductance, numerical and analytic d²I/dV² |
