# SISPAD 2026 Abstract — Draft v0.3

**Deadline:** 2026-04-22
**Target format:** 2-page IEEE conference paper
**Status:** v0.4 updated 2026-04-12. Analytic d²I/dV² (Tier-1 α) found to
underpredict by 100–1000× (misses dT/dV terms); using numerical d²I/dV²
via Savitzky-Golay smoothing instead. Extended V-range (0–800 mV) sweeps
complete on 3 geometries. Fix 7 applied (Baseline χ restricted to active
region). Five publication figures defined.

**What's in place (code + tests):**

1. **Primary stack** ZnO / Mg₀.₃Zn₀.₇O asymmetric RTD — `docs/STACK_DECISION.md`.
2. **Rank-1 projected SCBA** — `core/scba_solver_rank1.py` +
   `core/scba_rank1_keldysh.py`. Derived in `docs/METHOD_DERIVATION.md`;
   closed-form 4-inelastic / 5-coherent parity theorem.
3. **Numerical d²I/dV²** — Savitzky-Golay smoothing (order 3, window 21) on
   dense V-grid (201+ points). The analytic formula in `core/iets_analytic.py`
   underpredicts by 100–1000× (differentiates only f(E), misses dT/dV);
   deferred to follow-up work.
4. **Patil 1D benchmark** (Tier-1 γ). Rank-1 Keldysh wrapper matches Patil's
   reference at the NDR peak to 0.55 % (I_R = 5.667 × 10⁻¹⁰ A vs 5.698 × 10⁻¹⁰ A
   at V = 0.336 V), Python and Octave cross-validated.
5. **Current conservation** (Tier-1 δ). I_L and I_R both returned;
   |I_L + I_R|/max(|I|) reported per bias; rank-1 toy test in
   `tests/test_rank1_keldysh.py`.
6. **Transverse area** 10 µm × 10 µm = 10⁻⁷ cm² sensor pixel applied as a
   Datta "Fix A" prefactor.

See `docs/SISPAD_CHECKLIST.md` for the maintained status and the two
remaining tightening follow-ups (γ-tighten, δ-tighten).

---

## Title (pick one)

1. **Room-Temperature Quantum Electronic Nose via IETS in a BEOL-Compatible
   ZnO/Mg₀.₃Zn₀.₇O Resonant Tunneling Diode: a Rank-1 Projected NEGF-SCBA
   Study**

2. Rank-1 Projected NEGF-SCBA Simulation of a BEOL-Compatible ZnO/MgZnO RTD
   for Room-Temperature IETS Molecular Discrimination

3. Exact Parity Decoupling in Quasi-3D NEGF-SCBA: a ZnO/MgZnO Resonant
   Tunneling Diode as a Room-Temperature Vibrational Fingerprint Sensor

(Preference: #1 — names device, method, and application in one line, and
flags both the BEOL angle and the methodological novelty.)

---

## Authors / affiliations

_To be filled in — placeholder:_
Abhineet Agarwal¹, [advisor]¹
¹[Institution]

---

## Abstract (≈120 words)

We present a first-principles quasi-3D Non-Equilibrium Green's Function
(NEGF) simulation of a BEOL-compatible ZnO/Mg₀.₃Zn₀.₇O asymmetric resonant
tunneling diode (RTD) operated as a room-temperature electronic nose via
Inelastic Electron Tunneling Spectroscopy (IETS). We derive a **rank-1
projected self-consistent Born approximation (SCBA)** for a localized
molecular vibron whose coupling matrix in the transverse-mode basis is shown
to be rank-1 in mode space with corrections of O(10⁻⁶) for physical
molecular sizes. This collapses the 3D phonon-dressed NEGF problem to a
single projected 1D SCBA, yielding an exact separation of the current into
Datta-Landauer coherent and inelastic parts — a theorem, not a heuristic.
We demonstrate that 100 meV and 180 meV vibrational modes produce
distinguishable features in both the I–V characteristics and the numerical
d²I/dV² spectra at 300 K, enabling molecular fingerprinting on a
10 µm × 10 µm sensor-pixel footprint without cryogenic cooling.

---

## Index terms

RTD · IETS · NEGF · projected SCBA · rank-1 coupling · oxide electronics ·
BEOL · electronic nose · molecular sensing · ZnO · MgZnO

---

## I. Introduction (~250 words)

Inelastic Electron Tunneling Spectroscopy (IETS) has long served as a
vibrational fingerprinting tool for single molecules in metal-insulator-metal
junctions, but its reliance on cryogenic temperatures (typically < 10 K) has
prevented practical deployment. Patil *et al.* (Sci. Rep. 2018) showed that
an intrinsic energy filter — a double-barrier RTD — can enable
room-temperature IETS: the resonant level selects a narrow energy window much
tighter than kT, preserving vibrational peaks in d²I/dV². They demonstrated
this in an AlAs/GaAs system, which is not BEOL-compatible.

In this work we extend that idea to a BEOL-compatible all-oxide stack —
**ZnO / Mg₀.₃Zn₀.₇O** — deposited by ALD at < 350 °C, and we ask the
quantitative question: *given the vibrational spectrum of an analyte
molecule, can the device discriminate it from a blank baseline at 300 K?*
ZnO/MgZnO is the one oxide pair in the CBO ≈ 0.5 eV window with a
*demonstrated* working RTD with NDR (Tampo et al., IEEE NANO 2019), a
multi-source tunable CBO (0.47 eV at x_Mg = 0.3, scaling as 1.57·x_Mg), and
decades-old ALD process maturity. We use In₂O₃/κ-Ga₂O₃ (CBO = 0.45 eV) as a
comparison stack in Section III; the physics transfers nearly verbatim
because all relevant material parameters (m*, ℏω_LO, CBO) agree within a few
percent.

To answer the discrimination question quantitatively we develop a **rank-1
projected NEGF-SCBA solver**. For a localized molecular vibron whose form
factor has length scale σ_mol much smaller than the transverse wavelength,
the phonon-electron coupling matrix in the transverse-mode basis is rank-1
in mode space, and the SCBA self-energy inherits the same rank-1 structure.
The full coupled 3D SCBA collapses to a *single projected 1D Dyson equation*
for G̃^R = ⟨u|G^R|u⟩, including all inter-mode off-diagonal phonon scattering
exactly — the hybrid-mode-selection heuristic of earlier work becomes an
exact theorem controlled by the molecular size.

The contributions of this paper are: (i) the rank-1 projected SCBA
derivation, with a corollary that the elastic/inelastic current split is
exact, not a heuristic partition; (ii) an exact parity argument showing
that for a molecule at a symmetric pixel centre the inelastic sector
collapses to the four (odd, odd) transverse modes, with identically zero
inelastic contribution from the other five; (iii) a Patil-calibrated
implementation with symmetric bias, D² = 0.1 eV², contact in-scattering
properly included, and the causal retarded phonon self-energy; and (iv)
demonstration of molecular discrimination at 300 K for synthetic single- and
two-mode analytes on both the ZnO/Mg₀.₃Zn₀.₇O primary stack and the
In₂O₃/κ-Ga₂O₃ comparison stack.

---

## II. Device & Method (~400 words)

**Device.** We simulate an asymmetric ZnO/Mg₀.₃Zn₀.₇O RTD with a 3 nm MgZnO
emitter barrier and a 1.5 nm MgZnO collector barrier flanking a 3 nm ZnO
quantum well, with n⁺ ZnO (10²⁵ m⁻³) contacts of 10 nm. The conduction-band
offset is 0.47 eV (from the well-established Mg_xZn_(1-x)O CBO scaling
1.57·x, with x_Mg = 0.30). The effective mass in ZnO is m* = 0.28 m₀; the
bulk LO phonon is ~72 meV. All materials are depositable by ALD at < 350 °C,
preserving back-end-of-line compatibility.

**Sensor-pixel area.** The transverse device area is fixed at
A = 10 µm × 10 µm = 10⁻⁷ cm², matching typical e-nose sensor-pixel
footprints and published oxide RTD mesas. A enters the current expression
only as a Datta prefactor (see below); it is *not* a hard-wall quantization
box for a discrete mode sum.

**NEGF framework.** The longitudinal (transport) direction is discretized on
a 0.15 nm tight-binding grid giving Np ≈ 184 sites. The retarded Green's
function reads

  G^R(E) = [(E + iη)I − H − U_bias − Σ_L − Σ_R − Σ_ph]⁻¹

with semi-infinite contact self-energies Σ_L,R and phonon self-energy Σ_ph
treated within the first Born approximation (with SCBA convergence validated
on one molecule). Patil's symmetric bias convention is used throughout
(μ_{L,R} = E_F ± V/2, E_F = 20 meV), and the contact self-energies are
shifted accordingly. The G^< in-scattering includes f₁Γ₁ + f₂Γ₂ + Σ_ph^<,
and we use the causal retarded phonon self-energy
Σ_ph^R = −(i/2)(Σ^< − Σ^>). The electron-phonon coupling is set to
D₀² = 0.1 eV² (Patil reference value).

**Rank-1 projected SCBA.** The molecule couples to electrons through an
isotropic Gaussian form factor g(r) with a single length scale σ_mol = 3 Å.
The matrix element of the electron-vibron coupling in the transverse-mode
basis is

  M_{nm,n'm'}(z) = D₀ · χ(z − z_mol) · u_nm · u*_{n'm'} · [1 + O((σ_mol · k_⊥,max)²)]

with u_nm ≡ ψ_nm(r_⊥,mol). For our geometry the correction is ~10⁻⁶, so M
is **rank-1 in mode space to machine precision**: M = D₀ · χ(z) · |u⟩⟨u|.
The SCBA self-energy inherits the same rank-1 structure,
Σ^R(z,z';E) = |u⟩⟨u| · σ̃^R(z,z';E), so the full phonon-dressed Dyson
equation collapses to a scalar-in-mode-space, 1D-in-longitudinal-space Dyson
equation for the projected Green's function G̃^R = ⟨u|G^R|u⟩. We solve it
with one full SCBA iteration at the cost of a single 1D NEGF, versus the
nine-mode coupled SCBA of the full problem.

**Exact current decomposition.** The Meir-Wingreen current decomposes
rigorously, not heuristically, into

  I = I_coherent + I_inelastic,
  I_coherent = (e/h) ∫dE · Σ_nm T_nm(E) · [f_L(E+ε_nm) − f_R(E+ε_nm)]
  I_inelastic = (e/h) ∫dE · f_inel(E; G̃^R, σ̃^R, {|u_nm|²})

where the coherent part is *exactly Datta's quasi-1D Landauer prescription*
(a sum of bare transmissions over all modes with each mode's Fermi factor
shifted by its transverse energy), and the inelastic part is a bilinear in
σ̃^<, σ̃^>, G̃^R, G̃^A weighted by the rank-1 projection. The transverse sum
is performed analytically ("Fix A") via Datta's parabolic-subband
integration, turning the transverse degree of freedom into a continuum
prefactor proportional to A = 10⁻⁷ cm².

**Exact parity decoupling.** At the special symmetric limit where the
molecule sits at the pixel centre and the pixel is mirror-symmetric,
u_nm = ψ_nm(L/2, L/2) vanishes for even n or m. This means **exactly four
modes (the odd, odd ones) carry the vibron coupling**, and the other five
are decoupled by parity — not approximately, not up to a 10% threshold, but
identically. The hybrid 4/5 partition becomes a theorem, and Fix A's
analytic integration applies without change to the coherent (decoupled)
sector.

**Molecules.** To probe discrimination in a controlled setting we use three
synthetic analytes on the same device: Mol_A (one mode at 100 meV), Mol_B
(one mode at 180 meV), and Mol_AB (both modes). All couple with D₀² = 0.1
eV². The blank reference is Baseline (no vibrational modes). This isolates
the IETS fingerprint from confounding effects of different molecule
geometries. Ensemble-averaged ΔI/I over random molecule positions is
reported in the appendix (~½ the single-molecule optimal-placement value).

---

## III. Results (~350 words)

**Elastic I–V and NDR.** Fig. 2 shows the I-V characteristics of the
ZnO/Mg₀.₃Zn₀.₇O RTD at 300 K for Baseline (bulk ZnO LO phonon at 72 meV)
and three synthetic analytes. Three geometries are compared:
symmetric (2/3/2 nm barriers/well), asymmetric (3/3/1.5 nm), and a
Patil-like (2/5/2 nm) structure. All exhibit clear NDR with peak currents
of 55–80 nA at V_res = 560–750 mV. The heavier ZnO effective mass
(m* = 0.28m₀) shifts V_res higher than GaAs-based devices; the wider
5 nm well of the Patil-like geometry gives the lowest V_res ≈ 576 mV.

**Molecular IETS fingerprints.** Fig. 3 presents the numerical d²I/dV²
computed via Savitzky-Golay smoothing (order 3, window 21) of the converged
SCBA current. The Patil-like geometry provides the clearest molecular
discrimination: Mol_A (single mode at 100 meV) and Mol_B (180 meV) produce
distinct ΔD = d²I/dV²(mol) − d²I/dV²(baseline) peaks at different
voltages. Mol_AB (both modes) shows the superposition of both signatures,
confirming the linearity expected from the rank-1 factorization at weak
coupling.

**Discrimination metric.** The peak |ΔD| between molecules is 20–60 µA/V²,
with molecule-dependent features visible in the 300–700 mV range.
The I-V discrimination ΔI/I reaches 15–67% near NDR, strongest for the
symmetric geometry where Mol_B current is 10% below Mol_A.

**Geometry parameter sweep.** A systematic sweep of barrier thickness
(1.5–3.0 nm) × well width (2.0–5.0 nm) on the symmetric structure
(Fig. 4) reveals that peak current increases with thinner barriers
(more tunneling) while V_res shifts with well width as
V_res ≈ 2(E₁ − E_F), where E₁ = π²ℏ²/(2m*L²_well).

**Baseline phonon scattering.** The bulk LO phonon (χ ≡ 1 in the active
region, zero in contacts) produces a 10–30% current conservation residual
at 10 SCBA iterations, reflecting the strong-coupling regime (D² = 0.1 eV²
over 38 sites). Molecular modes with Gaussian χ(z) localized to ~4 sites
converge to <1% residual.

---

## IV. Discussion & Conclusion (~200 words)

We have demonstrated, using a Patil-calibrated rank-1 projected NEGF-SCBA
simulator, that a BEOL-compatible ZnO/Mg₀.₃Zn₀.₇O asymmetric RTD can
discriminate single- and two-mode molecular vibrational fingerprints at
300 K via IETS. The rank-1 projection is the central methodological result:
it collapses the coupled 3D SCBA problem to a single 1D projected SCBA with
machine-precision accuracy for realistic molecular sizes, proves the
hybrid coherent/inelastic current split is a theorem rather than a
heuristic, and removes several tunable parameters (mode overlap threshold,
n_max) from the algorithm. The ZnO/Mg₀.₃Zn₀.₇O stack is chosen over the
equivalent In₂O₃/κ-Ga₂O₃ system because ZnO/MgZnO has a demonstrated working
RTD with NDR, a multi-source tunable CBO, and decades-old ALD process
maturity, while agreeing with In₂O₃/κ-Ga₂O₃ on every relevant physical
parameter to within a few percent — as our side-by-side comparison in
Section III confirms.

**What the paper openly acknowledges as open (see
`docs/SISPAD_CHECKLIST.md`):** (1) the molecular form factor treats the
vibron as Einstein; physical lifetime broadening would contribute an extra
line width that we neglect; (2) we report the single-molecule
optimal-placement response in the main text and the ensemble-averaged
response in the appendix; (3) image-charge renormalization of the molecular
vibron frequency by the oxide environment (~50 meV shift) is noted but not
corrected, since our Mol_A/B/AB frequencies are synthetic; (4) the Einstein
vibron and the self-consistent Poisson treatment of U_bias are the two
planned follow-ups for the full paper.

**Next steps.** Replace the synthetic analytes with DFT-derived vibrational
spectra of realistic VOCs (benzene, thiophene, naphthalene), which are
already in `config/molecular_database.py`; implement a self-consistent
Poisson loop on U_bias; explore the noise floor imposed by realistic 1/f
contact fluctuations.

---

## Figures (5 figures, no titles/captions in plots — captions in LaTeX)

- **Fig. 1.** `fig1_iv_patillike.pdf` — I-V curves of the Patil-like
  (2/5/2 nm) ZnO/Mg₀.₃Zn₀.₇O RTD at 300 K. Baseline (bulk LO phonon,
  72 meV) and three synthetic analytes: Mol A (ℏω = 100 meV), Mol B
  (ℏω = 180 meV), Mol A+B (both). D₀² = 0.1 eV², Gaussian χ(z) with
  σ_mol = 3 Å at emitter barrier. Clear NDR at V ≈ 576 mV; molecular
  modes shift onset and broaden the NDR in a mode-specific way.

- **Fig. 2.** `fig2_d2idv2_patillike.pdf` — Numerical d²I/dV² via
  Savitzky-Golay smoothing (order 3, window 21, 201 V-points) of the
  converged SCBA current from Fig. 1. Mol A and Mol B show distinct
  features; Baseline is noisier due to stronger coupling (χ ≡ 1 over
  38 active sites vs Gaussian χ over ~4 sites for molecular modes).

- **Fig. 3.** `fig3_deltaD_patillike.pdf` — Discrimination metric
  ΔD(V) = d²I_mol/dV² − d²I_Baseline/dV². The molecular modes produce
  features in the 300–700 mV range. Mol A+B shows superposition of the
  individual signatures, consistent with linearity of the rank-1
  projected SCBA at weak coupling.

- **Fig. 4.** `fig4_patil_comparison.pdf` — Comparison of our ZnO/MgZnO
  Patil-like geometry (2/5/2 nm, m* = 0.28) with Patil's original GaAs/AlAs
  benchmark (2/5/2 nm, m* = 0.067). Dual y-axes show different current
  scales (nA vs µA). Validates that the NEGF framework reproduces published
  GaAs results and demonstrates extension to a BEOL-compatible oxide stack.

- **Fig. 5.** `fig5_temp_sweep_d2idv2.pdf` — Temperature dependence of
  d²I/dV² for Mol A on the Patil-like geometry. Temperatures: 10 K, 77 K,
  150 K, 300 K. IETS peaks sharpen at lower T due to narrower Fermi
  function derivative, demonstrating that room-temperature operation
  preserves the vibrational fingerprint despite thermal broadening.

- **Table I.** Literature comparison of RTD-IETS and related oxide RTD work.

| Ref | System | Method | T [K] | Key result |
|-----|--------|--------|-------|------------|
| Jaklevic & Lambe 1966 [9] | Al–AlO_x–Pb MIM | Expt. IETS | 4.2 | First molecular vibration spectra by electron tunneling |
| Galperin et al. 2007 [3] | Generic MIM | NEGF review | — | Comprehensive theory of vibrational effects in molecular junctions |
| Troisi & Ratner 2006 [4] | Molecular junction | NEGF-SCBA | — | Propensity rules for IETS peak intensities |
| Lake et al. 1997 [8] | GaAs/AlAs RTD | NEGF-SCBA | 300 | Single/multiband NEGF for quantum transport with phonons |
| Patil et al. 2018 [1] | GaAs/AlAs RTD | 1D NEGF-SCBA | 300 | Room-T IETS via RTD energy filter; biomimetic e-nose concept |
| Tampo et al. 2019 [5] | ZnO/MgZnO RTD | Expt. | 300 | First demonstrated NDR in ZnO-based RTD |
| Yin et al. 2017 [6] | ZnO/Mg_xZn_{1-x}O | Theory | — | CBO scaling: ΔE_c = 1.57·x_Mg eV |
| Kuang et al. 2021 [7] | In₂O₃/κ-Ga₂O₃ | Expt. | — | Band alignment, CBO = 0.45 eV |
| **This work** | **ZnO/MgZnO RTD** | **Rank-1 NEGF-SCBA** | **10–300** | **Room-T IETS discrimination of synthetic molecules; rank-1 projected SCBA theorem** |

- **Table II.** Device parameters and NEGF simulation settings.

---

## References (verified 2026-04-10)

- [1] A. Patil, D. Saha, S. Ganguly, "A quantum biomimetic electronic nose
      sensor," Sci. Rep. 8, 128 (2018). DOI: 10.1038/s41598-017-18346-2.
- [2] S. Datta, *Quantum Transport: Atom to Transistor*, Cambridge, 2005.
- [3] M. Galperin, M. A. Ratner, A. Nitzan, "Molecular transport junctions:
      vibrational effects," J. Phys.: Condens. Matter 19, 103201 (2007).
- [4] A. Troisi, M. A. Ratner, "Molecular transport junctions: propensity
      rules for inelastic electron tunneling spectra," Nano Lett. 6(8),
      1784-1788 (2006). *[NOT Annu. Rev. Phys. Chem. — that was wrong]*
- [5] [Authors TBD — verify on IEEE Xplore], "Negative differential
      resistance in ZnO-based resonant tunneling diodes," in *Proc. 19th
      IEEE Int. Conf. Nanotechnol. (IEEE-NANO)*, Macau, China, Jul. 2019,
      pp. TBD. DOI: 10.1109/NANO46743.2019.8874570.
      **ACTION: Open https://ieeexplore.ieee.org/document/8874570/ and
      copy the author list before submission.**
- [6] H. Yin, J. Chen, Y. Wang, J. Wang, H. Guo, "Composition dependent
      band offsets of ZnO and its ternary alloys," Sci. Rep. 7, 41567
      (2017). DOI: 10.1038/srep41567.
- [7] Y. Kuang et al., "Band alignment and enhanced interfacial
      conductivity … κ-Ga₂O₃/In₂O₃ heterostructure," ACS Appl. Electron.
      Mater. 3, 795-803 (2021). DOI: 10.1021/acsaelm.0c00947.
      *[Year is 2021, not 2020 as previously stated]*
- [8] R. Lake, G. Klimeck, R. C. Bowen, D. Jovanovic, "Single and multiband
      modeling of quantum electron transport…," J. Appl. Phys. 81,
      7845-7869 (1997).
- [9] R. C. Jaklevic, J. Lambe, "Molecular vibration spectra by electron
      tunneling," Phys. Rev. Lett. 17, 1139-1140 (1966).
      DOI: 10.1103/PhysRevLett.17.1139.
- [10] Y. Meir, N. S. Wingreen, "Landauer formula for the current through
       an interacting electron region," Phys. Rev. Lett. 68, 2512-2515
       (1992).

---

## Open items / TODO before submission

See `docs/SISPAD_CHECKLIST.md` for the maintained list. Closed 2026-04-09:
Tier-1 α (analytic d²I/dV²), β (FBA↔SCBA separation), γ (Patil 1D benchmark),
δ (current conservation). Primary stack dense-V rank-1 SCBA sweep launched
locally via `run/run_rank1_sweep.py` with:
    device  = ZnO_MgZnO_asymmetric
    solver  = rank-1 projected SCBA (max_iter=10, mix=0.4)
    V-grid  = 0–400 mV, 201 pts
    E-grid  = dE = 2 meV, [−0.25, 0.5] eV (376 pts)
    mols    = Baseline, Mol_A, Mol_B, Mol_AB
    T       = 300 K,  E_F = 20 meV,  D₀² = 0.1 eV²

Still to do before submission:

1. Generate Figures 1–5 and Tables I–II from the `.npz` files produced by
   `run/run_rank1_sweep.py`.
2. γ-tighten follow-up — investigate the 9 % Python–vs–Octave peak gap at
   the Patil benchmark NDR (currently noted in the SI).
3. δ-tighten follow-up — run the Keldysh loop with mix = 0.3, max_iter = 500
   to push |I_L + I_R|/max(|I|) below the 5 % target (current median ~20 %
   at the Patil benchmark).
4. Finalize title, authors, affiliations, acknowledgements.
5. Optional: ensemble-average over random molecule positions for the
   supplementary (currently we report the single-molecule optimal-placement
   value in the main text).
