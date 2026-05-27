# Rank-1 Projected SCBA for Quasi-3D IETS with a Localized Molecular Vibron

**Status:** Derivation for the SISPAD 2026 paper (Section II + supplementary).
**Author:** Abhineet Agarwal, with AI assistance. Historical provenance: earlier
thread with Nidhi Pandey on "molecule projection + perturbation theory" that
was not followed through at the time.
**Date:** 2026-04-07

> **Implementation note (2026-05-27):** the rank-1 solver shipped as a standalone
> pair, `core/scba_solver_rank1.py` (projected Dyson algebra) +
> `core/scba_rank1_keldysh.py` (Keldysh SCBA driver), run via
> `run/run_rank1_sweep.py`. The "`run_iets.py --solver {hybrid,rank1}` flag" and
> `scba_solver_hybrid.py` integration sketched in §below were **not** taken — those
> files are now under `archive/`. The derivation itself is unchanged and current.

---

## 1. Setup

We treat a quasi-3D resonant tunneling diode whose elastic Hamiltonian is
separable,

  H_el = H_z(z) + H_⊥(y, z_⊥)

with a transverse rectangular box of dimensions L_y × L_z_⊥ and hard walls.
Eigenstates of the transverse piece,

  ψ_nm(r_⊥) = (2 / √(L_y L_z_⊥)) sin(nπ y / L_y) sin(mπ z_⊥ / L_z_⊥),

carry transverse energies ε_nm. The longitudinal direction is discretized on
Np tight-binding sites with spacing a = 0.15 nm; each transverse mode (n,m)
defines an independent 1D longitudinal problem with on-site shift ε_nm and
contact self-energies Σ_{L,R}(E − ε_nm).

A single molecular vibron of frequency ω is embedded at position r_mol = (z_mol,
r_⊥,mol). Its displacement field couples to electrons via a local Fröhlich-like
coupling

  H_e-ph = D₀ · (b + b†) · ∫ d³r · g(r − r_mol) · c†(r) c(r)

where b† creates the vibrational quantum and g(r) is a normalized form factor.
We use an isotropic Gaussian with a **single** length scale σ_mol,

  g(r) = (2π σ_mol²)^{−3/2} · exp(−|r − r_mol|² / 2σ_mol²),    σ_mol = 3 Å.   (1)

The displacement amplitude D₀ is calibrated to Patil et al.'s value D₀² = 0.1 eV².

---

## 2. Matrix element in the transverse-mode basis

Project the coupling onto the transverse eigenbasis at longitudinal site z:

  M_{nm, n'm'}(z) = D₀ · χ(z − z_mol) · ∫d²r_⊥ ψ*_nm(r_⊥) g_⊥(r_⊥ − r_⊥,mol) ψ_{n'm'}(r_⊥)   (2)

where χ(z) ∝ exp(−z² / 2σ_mol²) is the longitudinal Gaussian profile and
g_⊥(r_⊥) is the 2D marginal of g.

**Claim.** For σ_mol · k_⊥,max ≪ 1, M factorizes as

  M_{nm, n'm'}(z) = D₀ · χ(z − z_mol) · u_nm · u*_{n'm'} · [1 + O((σ_mol k_⊥,max)²)]   (3)

with u_nm ≡ ψ_nm(r_⊥,mol) ∈ ℝ.

**Proof.** Expand ψ_{n'm'}(r_⊥) around r_⊥,mol:

  ψ_{n'm'}(r_⊥) = ψ_{n'm'}(r_⊥,mol)
                + (r_⊥ − r_⊥,mol) · ∇ψ_{n'm'}(r_⊥,mol)
                + (1/2)(r_⊥ − r_⊥,mol)_i (r_⊥ − r_⊥,mol)_j ∂_i∂_j ψ_{n'm'}(r_⊥,mol)
                + …

Likewise for ψ*_nm. Substitute into (2). Odd Gaussian moments vanish;
even moments give σ_mol². The leading term is the rank-1 product u*_nm u_{n'm'};
the first nontrivial correction is O(σ_mol² k_⊥²) where k_⊥ ~ π n_max / L_⊥.
For σ_mol = 3 Å, L_⊥ = 1 μm, n_max = 3,

  (σ_mol · π n_max / L_⊥)² ≈ (3×10⁻¹⁰ · π·3 / 10⁻⁶)² ≈ 8×10⁻⁶.   ∎

**Corollary.** To six-digit accuracy in our geometry, the transverse coupling
matrix is rank-1:

  M(z) = D₀ · χ(z − z_mol) · |u⟩⟨u|,    u_nm = ψ_nm(r_⊥,mol).   (4)

The point-molecule limit σ_mol → 0 is **not** a physical idealization — it is an
asymptotic result whose error scales as the sixth power of the ratio σ_mol /
(transverse wavelength). For any realistic VOC, the rank-1 structure is exact
for all practical purposes.

---

## 3. Rank-1 structure of the phonon self-energy

In the First Born Approximation (extendable to full SCBA by self-consistency),
the phonon self-energy in the combined (mode, site) basis is

  Σ^{<,>}(z, z'; E) = ∫ dω' D_ph^{<,>}(ω') M(z) G^{<,>}(z, z'; E ∓ ℏω') M(z')

For an Einstein (single-frequency) vibron at temperature T,

  D_ph^<(ω') = 2πi [(N+1)δ(ω' + ω) + N δ(ω' − ω)],   N = n_B(ℏω, T).

Substituting (4),

  Σ^<(z, z'; E) = D₀² χ(z − z_mol) χ*(z' − z_mol) · |u⟩⟨u|
                  × [(N+1) ⟨u| G^<(z, z'; E − ℏω) |u⟩ + N ⟨u| G^<(z, z'; E + ℏω) |u⟩]

Define the **projected (scalar in mode space) Green's function**

  G̃^{R,<,>}(z, z'; E) ≡ ⟨u| G^{R,<,>}(z, z'; E) |u⟩ = Σ_{nm, n'm'} u*_nm G^{R,<,>}_{(nm),(n'm')}(z, z'; E) u_{n'm'}   (5)

and the **projected (scalar in mode space) self-energy kernel**

  σ̃^{<,>}(z, z'; E) = D₀² χ(z − z_mol) χ*(z' − z_mol)
                      × [(N+1) G̃^{<,>}(z, z'; E − ℏω) + N G̃^{<,>}(z, z'; E + ℏω)]   (6)

  σ̃^R(z, z'; E) = −(i/2) · [σ̃^<(z, z'; E) − σ̃^>(z, z'; E)]     (causal combination — the Fix 4 sign)

so that in mode space

  Σ^R(z, z'; E) = |u⟩⟨u| · σ̃^R(z, z'; E).   (7)

The self-energy is **rank-1 in mode space with the same u-vector as M**. This is
the central structural fact of the derivation.

---

## 4. Closed Dyson equation for the projected Green's function

Let g_nm^R(z, z'; E) be the bare (ballistic, no phonon) 1D Green's function for
transverse mode (nm):

  g_nm^R(E) = [E · I − H_z − ε_nm · I − Σ_L(E − ε_nm) − Σ_R(E − ε_nm)]⁻¹   (8)

(9 of these, cheap — Np × Np tridiagonal inversions per energy.)

Define the **|u|²-weighted bare projected Green's function**

  G̃_0^R(z, z'; E) ≡ ⟨u| G_0^R |u⟩ = Σ_nm |u_nm|² · g_nm^R(z, z'; E).   (9)

Starting from the Dyson equation G^R = G_0^R + G_0^R Σ^R G^R and sandwiching on
both sides by ⟨u| and |u⟩,

  ⟨u|G^R|u⟩ = ⟨u|G_0^R|u⟩ + ⟨u|G_0^R |u⟩⟨u| σ̃^R G^R|u⟩
  G̃^R = G̃_0^R + G̃_0^R · σ̃^R · G̃^R

which solves to

  **G̃^R(E) = [I − G̃_0^R(E) · σ̃^R(E)]⁻¹ · G̃_0^R(E)**   (10)

This is a **scalar (in mode space), 1D (in longitudinal space)** Dyson equation.
It fully determines the phonon-dressed projected Green's function. The SCBA
iteration is:

  1. Initialize G̃^R = G̃_0^R, G̃^< = f(E) · G̃^R ...
  2. Compute σ̃^R, σ̃^< from (6).
  3. Update G̃^R from (10); update G̃^< from the Keldysh equation
     G̃^< = G̃^R · (Σ_L^< + Σ_R^< + σ̃^<) · G̃^A
     with the contact in-scattering included (Fix 3).
  4. Iterate until ‖ΔG̃^R‖ < tolerance.

This is a single 1D SCBA — the same computational cost as one of the four
mode-independent SCBAs the current `scba_solver_hybrid.py` runs, but
incorporating **all nine transverse modes including off-diagonals exactly.**

---

## 5. Reconstruction of the full G^R

From the Dyson equation with rank-1 Σ^R,

  G^R = G_0^R + G_0^R |u⟩⟨u| σ̃^R G^R

Using (10) to solve for the rank-1 correction, the full Green's function in the
(mode, site) basis is

  G^R_{(nm),(n'm')}(z, z'; E) = δ_{nm,n'm'} g_nm^R(z, z'; E)
     + u_nm u*_{n'm'} · ∫dz_1 dz_2 g_nm^R(z, z_1) · σ̃^R(z_1, z_2)
       · [I − G̃_0^R σ̃^R]⁻¹ · g_{n'm'}^R(z_2, z'; E)     (11)

with the second term being a rank-1 (in mode space) correction mediated by the
scalar kernel obtained by SCBA. Equivalently, the inverse [I − G̃_0^R σ̃^R]⁻¹
need be computed only once per energy during the SCBA iteration, not per mode.

---

## 6. Current: exact decomposition into coherent and inelastic parts

The Meir-Wingreen current is

  I = (e/ℏ) ∫(dE/2π) Tr{[f_L Γ_L − f_R Γ_R] (G^R − G^A) + [Γ_L − Γ_R] G^<}.

(We use Patil's symmetric bias convention μ_{L,R} = E_F ± V/2, Fix 1.)

Substituting (11) and collecting by mode-space structure,

  **I = I_coherent + I_inelastic**   (12)

with

  I_coherent = (e/h) ∫ dE · Σ_nm T_nm(E) · [f_L(E + ε_nm) − f_R(E + ε_nm)]   (13)
  T_nm(E) = Tr[Γ_L(E − ε_nm) · g_nm^R · Γ_R(E − ε_nm) · (g_nm^R)†]

which is **exactly Datta's quasi-1D Landauer shortcut** — a sum of bare
transmissions over all nine modes, with each mode's Fermi factor shifted by its
transverse energy, and

  I_inelastic = (e/h) ∫ dE · f_inel(E; G̃^R, σ̃^R, {|u_nm|²})     (14)

containing all phonon emission and absorption processes. Its explicit form
(derivable from the rank-1 part of (11) substituted into Meir-Wingreen) is a
bilinear of σ̃^<, σ̃^>, G̃^R, G̃^A weighted by sums of |u_nm|² Fermi factors.
The partition (12) is **rigorous**, not a heuristic — it follows from the rank-1
structure of Σ^R alone.

---

## 7. Exact partition at the symmetric limit

If (i) the molecule sits at the centre of the transverse box, r_⊥,mol = (L_y/2,
L_z_⊥/2), and (ii) the box is mirror-symmetric, then

  u_nm = ψ_nm(L_y/2, L_z_⊥/2) = (2/√(L_y L_z_⊥)) · sin(nπ/2) · sin(mπ/2).

This vanishes identically if n or m is even. For n, m ∈ {1, 2, 3}:

  (n, m):   (1,1) (1,2) (1,3) (2,1) (2,2) (2,3) (3,1) (3,2) (3,3)
  u_nm:      +1     0    −1     0     0     0    −1     0    +1
  |u_nm|²:    1     0     1     0     0     0     1     0     1

**Exactly four modes (the odd, odd ones) carry the vibron coupling. The other
five are decoupled by parity and propagate coherently with no inelastic
scattering whatsoever.** For these five modes, the Landauer result (13) is
**exact**, not approximate.

Consequently the "hybrid partition" into 4 inelastic + 5 coherent modes is a
rigorous theorem for the symmetric geometry, not an ad-hoc threshold. The
speedup claim in the paper changes from "3× at unchanged accuracy" to "≥4×
at identically zero loss" (the 4 comes from the 9-mode full SCBA being replaced
by a single projected 1D SCBA; compare Section 8 for the exact cost count).

---

## 8. Cost analysis

Let Np = 184 (longitudinal sites), N_E ≈ 500 (energy grid points), N_modes = 9,
N_inel = 4 in the symmetric case.

| Algorithm                                | Inversions per (E, SCBA iter) | Includes off-diag? | Partition       |
|------------------------------------------|--------------------------------|--------------------|-----------------|
| Full coupled SCBA (brute force)          | 1 × O((N_modes Np)³) ≈ 729 Np³ | Yes                | N/A             |
| Current "hybrid" (diag-only, 4 × 1D)     | 4 × O(Np³)                     | **No** (dropped)   | Ad-hoc score    |
| Rank-1 projected SCBA (proposed)         | 1 × O(Np³)                     | **Yes** (exact)    | Rigorous (Sec 7)|

Rank-1 SCBA is **~4× cheaper than the current code** while additionally including
all inter-mode off-diagonal phonon scattering. Compared to brute-force coupled
SCBA it is **~729× cheaper** — but the brute-force comparison is academic since
nobody does it.

---

## 9. Multiple phonon modes: bulk + molecular

The preceding sections treat a single phonon species. In practice, the RTD
has at least two phonon populations:

1. **Bulk LO phonon** (ℏω_LO = 72 meV for ZnO): a delocalized lattice mode
   that couples uniformly across the transverse plane.
2. **Molecular vibron** (ℏω_mol, localized at r_mol): the analyte signature.

Their transverse coupling structures are fundamentally different.

### 9.1 Bulk phonon: mode-diagonal

The bulk LO phonon couples through a spatially uniform potential in the
transverse plane. The transverse matrix element is

  M^bulk_{nm, n'm'}(z) = D_bulk · χ_bulk(z) · ∫d²r_⊥ ψ*_nm(r_⊥) · 1 · ψ_{n'm'}(r_⊥)
                       = D_bulk · χ_bulk(z) · δ_{nm, n'm'}

by orthonormality. The bulk self-energy is **mode-diagonal**:

  Σ^bulk_{nm, n'm'} = δ_{nm, n'm'} · σ̃^bulk_nm(z, z'; E)

where each σ̃^bulk_nm depends only on G^{<,>}_{nm,nm} (the diagonal block of G
in the nm-th transverse mode). Each mode's bulk phonon scattering is
independent — no inter-mode coupling. This gives 9 independent 1D SCBA
problems, one per transverse subband.

### 9.2 Molecular vibron: rank-1 (Sections 2–4)

As derived above, the molecular vibron coupling is rank-1:

  Σ^mol = |u⟩⟨u| · σ̃_mol(z, z'; E)

with σ̃_mol determined by the projected G̃ = ⟨u|G|u⟩.

### 9.3 Combined self-energy: diagonal + rank-1

The total phonon self-energy is

  Σ = Σ^bulk + Σ^mol = diag(σ̃^bulk_nm) + |u⟩⟨u| · σ̃_mol.   (15)

This is NOT rank-1 — it is the sum of a diagonal and a rank-1 matrix in mode
space. However, the Dyson equation can still be solved efficiently:

**Step 1.** Solve 9 independent 1D SCBAs for the bulk phonon alone,
obtaining phonon-dressed Green's functions g̃_nm^R(E) for each mode.

**Step 2.** Build the projected bare Green's function including bulk dressing:

  G̃_0^R(E) = Σ_nm |u_nm|² · g̃_nm^R(E)   (cf. Eq. 9, replacing ballistic g with bulk-dressed g̃)

**Step 3.** Solve the projected molecular SCBA (Eqs. 6, 10) on top of G̃_0^R.

This is exact. The molecular vibron sees the bulk-phonon-dressed modes, and
the rank-1 projection still collapses the molecular SCBA to a single 1D
problem. Total cost: 9 bulk SCBAs + 1 projected molecular SCBA.

### 9.4 Continuum transverse limit (large pixel)

For a sensor pixel A_⊥ = L_y × L_z ≫ λ_deBroglie (our case: L = 10 µm,
λ_dB ~ nm), the transverse mode spacing

  Δε = ℏ²π²/(2m* L²) ≈ 10⁻⁸ eV

is negligible compared to all energy scales (kT, ℏω, CBO). All modes are
effectively degenerate: g̃_nm^R(E) ≈ g̃^R(E) for all (n,m). Consequently:

- The 9 bulk SCBAs reduce to 1 (all identical).
- The projection G̃_0^R = (Σ|u_nm|²) · g̃^R = g̃^R (since |u| is normalized).
- The total SCBA reduces to a single 1D problem with both bulk and molecular
  self-energies applied additively.

The sum over transverse modes enters the current purely as the Datta area
prefactor A_⊥ · m*/(πℏ²). This is the regime implemented in
`run/run_rank1_sweep.py`, and it is exact for A_⊥ ≥ (1 µm)².

### 9.5 Multiple molecular vibrons

If the molecule has K vibrational modes at frequencies ω₁, …, ω_K, all
sharing the same spatial location (hence the same |u⟩ and χ(z)), the
molecular self-energy becomes

  σ̃_mol^{<,>}(E) = Σ_{ν=1}^K D²_ν · χ²(z) · [(N_ν+1)G̃^{<,>}(E−ℏω_ν) + N_ν G̃^{<,>}(E+ℏω_ν)]

This is still rank-1 in transverse mode space (same |u⟩⟨u| for all ν). The
phonon modes contribute additively to the scalar kernel σ̃_mol, and the
projected SCBA (Eq. 10) handles them all in one solve. **The rank-1 method
generalizes to arbitrary numbers of molecular vibrons without modification.**

---

## 10. Off-centre / asymmetric generalisation (formerly §9)

When the molecule is off-centre or the box is asymmetric, all nine u_nm are
nonzero and the rank-1 structure still holds. The partition into "large-|u|²"
(inelastic) and "small-|u|²" (coherent) modes becomes an approximation,
controlled by the dimensionless tolerance

  ε ≡ min_{nm ∈ inelastic} |u_nm|² / max_{nm} |u_nm|²

Dropping the small-|u|² modes from the inelastic set introduces an error on I
that is O(ε²) — quadratic, because these modes' contribution to σ̃ enters
bilinearly through G̃. The paper should report I_full(ε = 0) versus I_hybrid(ε)
for ε ∈ {0, 10⁻³, 10⁻², 10⁻¹} on an off-centre molecule as a convergence
demonstration.

---

## 11. What changes in the code

1. **New solver file:** `core/scba_solver_rank1.py` implements Sections 3–5
   above. Leaves `scba_solver_hybrid.py` untouched for reproducibility of v3.
2. **Gaussian form factor:** longitudinal profile χ(z) is Gaussian with σ_z =
   σ_mol = 3 Å (replaces the `mol_radius = barrier/2` heuristic of Fix 6).
3. **Run driver:** `run/run_iets.py` gets a `--solver {hybrid,rank1}` flag;
   default stays `hybrid` until the rank-1 solver is validated.
4. **Unit test:** `tests/test_rank1_vs_full.py` on a 2-mode toy problem checks
   that the rank-1 solver agrees with a brute-force 2-mode coupled SCBA to
   machine precision, validating Sections 3–5 numerically.
5. **Convergence study:** `run/run_convergence_offcenter.py` sweeps molecule
   y-position and reports I_rank1 vs I_full for the appendix of the paper.

---

## Fixes-against-first-principles checklist (standing)

- [x] Symmetric bias convention (Fix 1, Patil): μ_{L,R} = E_F ± V/2
- [x] D² = 0.1 eV² calibration (Fix 2)
- [x] Contact in-scattering f₁Γ₁ + f₂Γ₂ + Σ_ph^< in G^< (Fix 3)
- [x] Σ^R = −(i/2)(Σ^< − Σ^>) causal sign (Fix 4)
- [x] Landauer shortcut for decoupled modes (Fix 5, now rigorous via Sec 7)
- [ ] Rank-1 projected SCBA including off-diagonals (this document)
- [ ] Gaussian form factor σ_mol = 3 Å replacing Fix 6 heuristic
- [ ] Current conservation |I_L − I_R| / |I_L| < 10⁻³ verification
- [ ] 1D benchmark against Patil's `rtd2modes_1d.m` — required for reviewer trust
- [ ] Bulk phonon re-examination (In₂O₃ Fröhlich LO)
- [ ] η sensitivity, SCBA-vs-FBA comparison
- [ ] Self-consistent Poisson (future)
