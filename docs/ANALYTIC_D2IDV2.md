# Analytic d²I/dV² — SISPAD Tier-1 item α

**Status:** 2026-04-08, in progress.
**Problem file:** `run/run_iets.py:594-595` takes two successive `np.gradient`
over an 81-point V-grid — roughly 4×10⁴ noise amplification. Current IETS
plots are mostly noise; the signal at V = ℏω/e is buried.
**Goal:** replace numerical differentiation with a closed-form
d²I/dV² derived by differentiating the Meir–Wingreen current *symbolically*
w.r.t. V, evaluated at whatever sparse V-grid the user requests. This file
is the derivation; the implementation lives in `core/iets_analytic.py`.

## 1. Bias convention

Patil / SISPAD symmetric bias:

    μ_L(V) = E_F + V/2,   μ_R(V) = E_F − V/2.       (1)

Fermi functions f_α(E, V) = 1 / (1 + exp[(E − μ_α(V)) / kT]).

Chain-rule identities (used throughout):

    ∂f_L/∂V = −(1/2) ∂f_L/∂E,    ∂f_R/∂V = +(1/2) ∂f_R/∂E.      (2)
    ∂²f_L/∂V² = (1/4) ∂²f_L/∂E², ∂²f_R/∂V² = (1/4) ∂²f_R/∂E².    (3)

(All mixed partials commute; f_α depends on V only through μ_α(V) which is
linear, so higher V-derivatives just pick up (1/2)^n and the corresponding
E-derivative with an alternating sign per contact.)

Analytic E-derivatives of the Fermi function, with x ≡ (E − μ_α)/kT:

    f     = 1/(1 + e^x)
    f'    = −(1/kT)·f(1−f)                                      (4)
    f''   = +(1/(kT)²)·f(1−f)(1−2f)                             (5)

## 2. Meir–Wingreen current

For a two-contact device with retarded and advanced Green's functions G^R,
G^A and contact broadenings Γ_L, Γ_R, the current flowing into the left
contact is (Meir–Wingreen, 1992):

    I_L(V) = (e/h) ∫ dE Tr{ Γ_L · [ f_L · A − i G^< ] },        (6)

with A(E,V) = i(G^R − G^A) the spectral function and G^< the lesser
Green's function. The symmetric form used in SISPAD work is

    I(V) = (e/h) ∫ dE Tr{[f_L Γ_L − f_R Γ_R] A + i (Γ_L − Γ_R) G^<}.  (7)

We split this into a **coherent (Landauer)** and **inelastic (Keldysh)**
part — the decomposition is rigorous in the rank-1 projected SCBA (see
`docs/METHOD_DERIVATION.md` §6):

    I(V) = I_coh(V) + I_inel(V),                                (8)

    I_coh(V)  = (e/h) ∫ dE T(E) · [f_L − f_R],                  (9a)
    I_inel(V) = (e/h) ∫ dE Q_inel[{Σ^{<,>}_ph, G^{R,A,<}}; f_L, f_R](E). (9b)

    T(E) = Tr[Γ_L G₀^R Γ_R G₀^A]  is the elastic transmission function,  (9c)

where G₀ refers to the purely elastic (ballistic or potential-only) Green's
function — Σ_ph is peeled into the inelastic bookkeeper Q_inel. Under the
first Born approximation (FBA) Σ_ph is evaluated once at V = 0, so its
V-dependence enters only through (f_L, f_R) in the Keldysh integrands of
Q_inel. Under self-consistent Born (SCBA), Σ_ph is itself a functional of
G^<(V), and a full chain rule is needed — see §5.

## 3. Analytic d²I/dV² — coherent path

Differentiating (9a) twice w.r.t. V and using (3):

    d²I_coh/dV² = (e/h) ∫ dE T(E) · ∂²(f_L − f_R)/∂V²
                = (e/h) ∫ dE T(E) · K_therm(E, V),             (10)

with the **thermal d² kernel**

    K_therm(E, V) ≡ (1/4) [∂²f_L/∂E² − ∂²f_R/∂E²]
                  = (1/(4(kT)²)) [f_L(1−f_L)(1−2f_L)
                                  − f_R(1−f_R)(1−2f_R)].        (11)

This kernel is *exact*, antisymmetric about E = E_F (for V > 0), and strongly
peaked near E = μ_L and E = μ_R with opposite signs. It vanishes identically
at V = 0 and at T → 0 except at the Fermi step.

**Numerical noise comparison.** For 81 V-points and η = 20 meV, the
double `np.gradient` path has relative noise ≈ (δV)⁻² · ε_machine ≈
10⁻³/(5e−3)² ≈ 4×10⁴, which swamps the signal. Replacing (10) removes the
δV⁻² factor entirely — the only numerical error is the E-quadrature on
T(E), which is well-resolved at 500 E-points.

## 4. IETS from the coherent path alone — sanity check

The coherent kernel in (11) is *not* what produces IETS peaks at
V = ℏω/e. IETS peaks live in I_inel. In a single-barrier resonant diode,
I_coh dominates the current but d²I_coh/dV² only has features at the
transmission resonances in T(E) — the main NDR peak and its shoulder. A
sanity check of the implementation is therefore:

    d²I_coh/dV²  has peaks at V = 2·(E_res − E_F)  (factor 2 from symmetric bias),

and the IETS-specific peaks at V = ℏω/e appear only once I_inel(V) is added.

This is useful for validating (10) before (12) is in place: on any molecule
with D₀ = 0 (baseline run), the analytic d²I/dV² must agree with the
numerical d²I/dV² (up to the numerical noise floor) because only I_coh
contributes.

## 5. Analytic d²I/dV² — inelastic path *(follow-up)*

At FBA, Σ_ph is evaluated at V = 0, so it's an explicit constant. The
inelastic current is a bilinear in (G^<, G^>, Σ^{<,>}_ph) times Γ_L − Γ_R
(Eq. 7), and G^{<,>} depend on V only through (f_L, f_R) via

    G^< = G^R [f_L Γ_L + f_R Γ_R + Σ^<_ph] G^A,
    G^> = G^R [(1−f_L) Γ_L + (1−f_R) Γ_R + Σ^>_ph] G^A.          (12)

The V-derivatives therefore pick up (f_L, f_R) factors through the chain
rule — no divergences, no numerical differentiation. The resulting analytic
expression is a sum of integrals of the form

    ∫ dE g_a(E) · f_L^{(n_L)}(E, V) · f_R^{(n_R)}(E, V),        (13)

with n_L + n_R ≤ 2 and g_a(E) a V-independent kernel built from (Γ_L, Γ_R,
G₀^R, G₀^A, Σ^<_ph, Σ^>_ph). Each of these kernels is computable once per
E-grid after the rank-1 projected solve converges.

**Status / blocker.** Implementing (13) requires the Keldysh rank-1
wrapper — the layer that builds σ̃^{<,>} from G̃^{<,>} on the 1D
longitudinal grid and folds in the contact in-scattering terms. That
wrapper is next on the v4 critical path. Until it lands, this file
implements (10) and (11) only; the inelastic analytic path is laid out in
pseudocode at the end of this document as a design artifact and a
TODO-magnet for the next session.

Under SCBA (as opposed to FBA), Σ_ph is a functional of G^<(V), so the
chain rule picks up a term ∂Σ_ph/∂V that must be computed implicitly via
the self-consistent fixed point (∂G/∂V = ∂G/∂G · ∂Σ/∂V · ∂G/∂V + ∂G/∂V|_f ).
The SISPAD primary-stack run uses `--max-iter 10 --scba-mix 0.4` for Mol_A
only (Tier-1 β); for FBA the inelastic chain rule truncates at (13).

## 6. Implementation plan (core/iets_analytic.py)

Minimum surface:

    def fermi(E, mu, kT):                          # f(E, μ)
    def fermi_deriv_E(E, mu, kT):                  # ∂f/∂E
    def fermi_second_deriv_E(E, mu, kT):           # ∂²f/∂E²
    def thermal_d2_kernel(E, V, kT, E_F=0.0):      # K_therm from (11)
    def analytic_d2idv2_coherent(T_of_E, E_grid,   # (e/h) ∫ K·T dE from (10)
                                 V_grid, kT, E_F=0.0):

All scalar-in-float, broadcast in V for vectorised evaluation. Validated in
`tests/test_analytic_d2idv2.py` against:

1. `f` / `f'` / `f''` analytic identities against `scipy.misc.derivative`-style
    finite differences at the 1e-8 level.
2. The coherent d²I/dV² against a dense (5000-point) finite-difference
    d²/dV² on a toy T(E) = Lorentzian — agreement to 1e-6.
3. The V = 0 limit — K_therm(E, 0) ≡ 0 for all E.
4. Charge conservation: ∫ dE K_therm(E, V) = 0 (kernel has zero mass).

The inelastic analytic path lives in §5 as pseudocode and a TODO.

## 7. Pseudocode for the inelastic path *(blocked on Keldysh rank-1 wrapper)*

    def analytic_d2idv2_inelastic(rank1_state, V_grid, kT, E_F):
        # rank1_state = output of rank1_scba_loop: G̃^R, Σ̃^R, G̃^<, Σ̃^<,
        # plus the per-mode Γ_{L,R} projected with |u_nm|² weights.

        for V in V_grid:
            mu_L = E_F + V/2
            mu_R = E_F - V/2

            for E in E_grid:
                # Evaluate (13) term by term:
                #   ∂²I/∂V² |_V = Σ_{n_L + n_R ≤ 2} g_{n_L,n_R}(E) ·
                #                  f_L^{(n_L)}(E,V) · f_R^{(n_R)}(E,V)
                # with g_{n_L,n_R}(E) precomputed from
                #   Tr[Γ_L G̃^R (Γ_L - Γ_R) G̃^A · σ̃^{<,>}_ph · G̃^A]
                # etc. See Ryndyk Ch. 5 for the combinatorics.
                ...

        return d2I_inel_of_V

This is deliberately a sketch. Filling it in needs:
  - The rank-1 Keldysh wrapper exposing σ̃^{<,>} and the per-mode |u_nm|²
    Γ projections.
  - A symbolic differentiation pass (or hand derivation) of (13) to
    enumerate the (n_L, n_R) terms and the corresponding g kernels.
  - Validation against the numerical FBA d²I/dV² on Baseline and Mol_A at
    V-grid ≥ 500 pts (once numerical noise is tolerable).

Deferred to the same sprint that lands `core/scba_rank1_keldysh.py`.

## References

- Meir & Wingreen, PRL 68, 2512 (1992) — Eq. (6).
- Ryndyk, *Theory of Quantum Transport at Nanoscale*, Springer 2016, Ch. 5 —
  FBA/SCBA and IETS kernels.
- Patil et al., `papers/rtd2modes_1d.m` — 1D symmetric-bias reference.
- `docs/METHOD_DERIVATION.md` — rank-1 projected SCBA and the rigorous
  coherent/inelastic decomposition (Eq. 13, 14).
