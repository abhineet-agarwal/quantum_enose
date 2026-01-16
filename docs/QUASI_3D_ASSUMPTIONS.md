# Quasi-3D Mode Summation: Assumptions and Approximations

## Overview

This document lists all assumptions, approximations, and idealizations in the quasi-3D transverse mode summation implementation for NEGF-based RTD transport.

Understanding these assumptions is critical for:
- Interpreting simulation results correctly
- Identifying when results may be inaccurate
- Planning future improvements

**Last Updated**: 2026-01-13
**Implementation**: Phase 2 (Ballistic 3D with transverse modes)

---

## 1. Inter-Mode Coupling

### Assumption: No Inter-Mode Scattering
**Statement**: Phonon scattering only occurs within each transverse mode (n,m). No coupling between different modes (n,m) ↔ (n',m').

**Mathematical Expression**:
```
Σ_S[nm, n'm'] = 0  for (n,m) ≠ (n',m')
Σ_S is diagonal in mode space
```

**Physical Justification**:
- Bulk acoustic/optical phonons approximately conserve transverse momentum k_⊥
- For small-angle scattering, Δk_⊥ << π/Ly
- Elastic scattering dominates in quantum wells

**When Valid**:
- ✓ Bulk phonon scattering
- ✓ Devices with transverse dimensions >> phonon wavelength
- ⚠️ May break for molecular vibrations if molecule size ~ Ly

**Impact**:
- Modes evolve independently in SCBA iteration
- Computational cost: O(Nm) instead of O(Nm²)

**Code Location**: `core/scba_solver.py:474-493`

**Improvement Path**: Implement off-diagonal coupling matrix D[nm, n'm']

---

## 2. Contact Properties

### Assumption: Mode-Independent Contact Self-Energies
**Statement**: All transverse modes couple identically to the left and right contacts.

**Mathematical Expression**:
```
Σ₁(E) same for all modes
Σ₂(E) same for all modes
Γ₁_nm = Γ₁ for all (n,m)
```

**Physical Justification**:
- Valid if contact regions are much larger than transverse dimensions
- Uniform doping and composition in contact regions
- No transverse electric field variation

**When Valid**:
- ✓ Uniform planar contacts with area >> Ly × Lz
- ✓ Ohmic contacts with uniform doping
- ✗ Structured contacts (e.g., nanowire arrays)
- ✗ Non-uniform contact barriers

**Impact**:
- All modes have same injection/collection efficiency
- Simplifies current calculation significantly

**Code Location**:
- `core/green_functions.py:349` - Sigma1, Sigma2 passed as 2D
- `core/scba_solver.py:407-411` - Gamma arrays are 3D not 4D

**Improvement Path**: Implement mode-dependent contact self-energies

---

## 3. Hamiltonian Structure

### Assumption: Separable Longitudinal-Transverse Hamiltonian
**Statement**: The total Hamiltonian cleanly separates into longitudinal and transverse parts.

**Mathematical Expression**:
```
H_total = H_long(x) ⊗ I_trans + I_long ⊗ ε_nm

where:
  H_long(x) = kinetic + potential along x
  ε_nm = transverse mode energy
```

**Physical Justification**:
- Device has uniform cross-section (no x-dependence in Ly, Lz)
- Potential V(x) is uniform in transverse plane
- No transverse electric fields

**When Valid**:
- ✓ Rectangular quantum wire/waveguide
- ✓ Planar RTD with uniform layers
- ✗ Tapered structures
- ✗ Devices with transverse potential variation

**Impact**:
- Each mode sees same longitudinal device, shifted in energy
- Enables independent mode calculations

**Code Location**: `core/green_functions.py:395`

**Limitation**: Cannot model transverse band bending or confinement variation

---

## 4. Phonon Coupling

### Assumption A: Bulk Phonons - Mode-Independent Coupling
**Statement**: Bulk phonon coupling strength D is same for all transverse modes.

**Mathematical Expression**:
```
D_bulk[nm] = D_bulk  (scalar, not matrix)
```

**Physical Justification**:
- Bulk phonons are material property
- Deformation potential coupling nearly uniform
- Long-wavelength phonons don't resolve transverse structure

**When Valid**:
- ✓ LO/LA phonon scattering
- ✓ Wavelength λ_phonon >> Ly
- ⚠️ Approximate for confined acoustic phonons

**Impact**: Minor - bulk scattering affects all modes equally

**Code Location**: `run/run_single_molecule.py:100-109`

---

### Assumption B: Molecular Vibrations - Mode-Independent Coupling ⚠️ CRITICAL ISSUE

**Current Implementation**: Molecular coupling D_molecule is scalar → all modes couple equally

**What Should Happen**:
```
D_nm(molecule) = D₀ × |ψ_nm(y_mol, z_mol)|²
```

**Physical Reality**:
- Molecule at position (y_mol, z_mol) couples to mode (n,m) with strength ∝ |ψ_nm|²
- Modes with nodes at molecule location → NO coupling
- Modes with maxima at molecule location → STRONG coupling

**Example** (Molecule at center Ly/2, Lz/2):
```
Mode (1,1): ψ₁₁(Ly/2, Lz/2) = max  → D₁₁ = D₀ (strong)
Mode (1,2): ψ₁₂(Ly/2, Lz/2) = 0    → D₁₂ = 0 (NO coupling)
Mode (2,1): ψ₂₁(Ly/2, Lz/2) = 0    → D₂₁ = 0 (NO coupling)
Mode (2,2): ψ₂₂(Ly/2, Lz/2) = max  → D₂₂ = D₀ (strong)
```

**Impact**:
- **CRITICAL for molecular sensing**
- Without this, all modes scatter equally → wrong physics
- This is the "spatial locality" principle from explain.md

**Status**: ⚠️ **NEEDS TO BE FIXED**

**Code Location**: `run/run_single_molecule.py:118-133`

**Fix Priority**: **IMMEDIATE** - This is essential for correct molecular sensing

---

## 5. Boundary Conditions

### Assumption: Hard-Wall Transverse Confinement
**Statement**: Transverse modes use particle-in-a-box wavefunctions with hard walls.

**Mathematical Expression**:
```
ψ_nm(y,z) = (2/√(Ly·Lz)) sin(nπy/Ly) sin(mπz/Lz)
ψ_nm(0,z) = ψ_nm(Ly,z) = 0  (hard walls)
ψ_nm(y,0) = ψ_nm(y,Lz) = 0
```

**Physical Justification**:
- Infinite potential barriers at device edges
- Abrupt sidewall definition
- No leakage into surrounding material

**When Valid**:
- ✓ Etched semiconductor structures
- ✓ Strong dielectric/air interfaces
- ⚠️ Approximate for graded edges
- ✗ Smooth confinement potentials

**Impact**:
- Mode energies: ε_nm = (ℏ²π²/2m*)[(n²/Ly²) + (m²/Lz²)]
- Wavefunction nodes and antinodes at specific locations

**Code Location**: `core/hamiltonian.py:268-291`

**Alternative**: Harmonic oscillator (parabolic confinement) or numerical solution

---

## 6. Effective Mass

### Assumption: Single Effective Mass for All Modes
**Statement**: All transverse modes and all energies use same m*.

**Mathematical Expression**:
```
m*(E, n, m) = m* = constant
```

**Physical Justification**:
- Parabolic band approximation near zone center
- Energy range small compared to band gap
- Transverse confinement doesn't affect band structure

**When Valid**:
- ✓ Small energy range (< 0.5 eV from band edge)
- ✓ Direct-gap semiconductors (GaAs, In₂O₃)
- ⚠️ Approximate for wide energy ranges
- ✗ Non-parabolic bands, high energies

**Impact**:
- Mode energies might be slightly wrong for high-energy modes
- Minor effect on transmission

**Code Location**: `core/transverse_modes.py:73-75`

**Improvement Path**: Energy-dependent m*(E) from k·p theory

---

## 7. Interface Effects

### Assumption: No Mode Conversion at Interfaces
**Statement**: Electron maintains same transverse mode (n,m) when crossing interfaces/barriers.

**Mathematical Expression**:
```
Transmission: T_nm(E) independent
No T_nm→n'm' conversion
```

**Physical Justification**:
- Smooth interfaces preserve transverse momentum
- Barrier width >> Ly (no transverse scattering)
- No interface roughness

**When Valid**:
- ✓ Atomically smooth interfaces
- ✓ Gradual composition changes
- ⚠️ Approximate for real interfaces
- ✗ Rough interfaces, defects

**Impact**:
- Overestimates transmission for higher modes
- Each mode evolves independently

**Code Location**: Implicit in mode-independent Green's function calculation

**Reality**: Some mode mixing occurs at real interfaces

---

## 8. Molecular Properties

### Assumption A: Point Molecule (No Spatial Extent)
**Statement**: Molecule treated as point perturbation at position (x_mol, y_mol, z_mol).

**Physical Reality**:
- Molecules have finite size (benzene ring ~ 0.28 nm diameter)
- Vibrational modes distributed over molecular volume
- Coupling depends on overlap integral

**When Valid**:
- ✓ Molecule size << transverse dimensions (Ly ~ 1 µm >> 0.3 nm)
- ✓ Approximate for small molecules
- ⚠️ Less accurate for large molecules

**Impact**:
- Simplified spatial coupling
- Cannot model molecular shape effects

**Code Location**: `core/self_energy.py:get_molecule_location()`

---

### Assumption B: No Molecular Orientation
**Statement**: Molecule treated as isotropic perturbation.

**Physical Reality**:
- Benzene ring lies in plane → anisotropic
- C-H stretches have direction
- Different modes couple differently to device axes

**When Valid**:
- ⚠️ Rough approximation
- Better for spherical molecules

**Impact**:
- Cannot distinguish orientation effects
- Averages over orientations

**Code Location**: Molecule position is 1D (x only)

**Improvement Path**: Include 3D molecular structure and orientation

---

### Assumption C: Single Molecule
**Statement**: Only one molecule present in device.

**Physical Reality**:
- Real sensors may have multiple molecules
- Surface coverage effects
- Molecule-molecule interactions

**When Valid**:
- ✓ Low concentration / single-molecule limit
- ✗ High coverage

**Impact**: Cannot model concentration effects

---

## 9. Temperature Effects

### Assumption: Uniform Temperature
**Statement**: Entire device at single temperature T.

**Mathematical Expression**:
```
f_FD(E, μ, T) uses single T
n_BE(ℏω, T) uses single T
```

**Physical Justification**:
- Good thermal contact
- Small device (fast heat dissipation)
- Low current (negligible Joule heating)

**When Valid**:
- ✓ Equilibrium or near-equilibrium
- ✓ Low bias voltages
- ⚠️ Approximate at high bias

**Impact**: Cannot model thermal gradients

**Code Location**: `core/scba_solver.py:386`

---

## 10. Computational Approximations

### Assumption A: All-Inelastic Mode Treatment
**Statement**: All transverse modes receive full SCBA treatment (inscattering + outscattering self-energies).

**Alternative (Hybrid)**: Select top M modes for SCBA, rest coherent

**Physical Reality**:
- Many modes contribute negligibly:
  - High-energy modes: f_FD(ε_nm) ≈ 0 (thermal suppression)
  - Modes with nodes at molecule: |ψ_nm|² ≈ 0 (spatial suppression)
- Only ~10-20% of modes contribute significantly

**When Valid**:
- ✓ Always accurate (no approximation)
- ✗ Computationally expensive

**Impact**:
- ~9x slower than hybrid approach (for 9 modes)
- explain.md reports 7-8x speedup with hybrid, 90-95% accuracy

**Code Location**: `core/scba_solver.py:474-493`

**Status**: Accurate but slow - candidate for optimization

---

### Assumption B: Energy Grid Discretization
**Statement**: Energy integral approximated by finite grid.

**Mathematical Expression**:
```
∫ dE f(E) ≈ Σᵢ f(Eᵢ) ΔE
```

**When Valid**:
- ✓ Sufficient grid points (NE > 100)
- ✓ dE smaller than feature width
- ⚠️ May miss narrow resonances

**Impact**:
- Numerical integration error
- Trade-off: accuracy vs speed

**Code Location**: `run/run_single_molecule.py:241-248`

**Typical**: NE = 200, dE ~ 10 meV

---

## 11. Self-Consistency

### Assumption: SCBA Loop Convergence
**Statement**: SCBA iteration converges to self-consistent solution.

**Reality**:
- May not converge for strong coupling (D > 0.2 eV)
- Mixing parameter α may need tuning
- Sometimes oscillates instead of converging

**When Valid**:
- ✓ Weak-to-moderate coupling (D < 0.15 eV)
- ✓ Proper mixing (α ~ 0.3-0.5)
- ⚠️ May fail for strong coupling

**Impact**:
- Non-converged results are approximate
- May need more iterations or different mixing

**Code Location**: `core/scba_solver.py:495-514`

**Monitor**: Convergence flag and residual

---

## Summary Table: Assumption Validity

| Assumption | Validity | Impact | Fix Priority |
|------------|----------|--------|--------------|
| **No inter-mode scattering** | ✓ Good | Low | Low |
| **Mode-independent contacts** | ✓ OK | Low | Low |
| **Separable Hamiltonian** | ✓ Good | Low | Low |
| **Bulk phonon mode-independent** | ✓ Good | Low | Low |
| **⚠️ Molecular coupling mode-independent** | ✗ **WRONG** | **CRITICAL** | **IMMEDIATE** |
| **Hard-wall boundaries** | ✓ OK | Low-Medium | Medium |
| **Single effective mass** | ✓ OK | Low | Low |
| **No mode conversion** | ⚠️ Approximate | Medium | Medium |
| **Point molecule** | ✓ OK | Low | Future |
| **No molecular orientation** | ⚠️ Simplification | Medium | Future |
| **Uniform temperature** | ✓ Good | Low | Low |
| **All-inelastic treatment** | ✓ Accurate | None (slow) | Medium (speed) |
| **Energy grid discretization** | ✓ OK | Low | Low |

---

## Critical Issue: Mode-Dependent Molecular Coupling

### Current Implementation ⚠️
```python
# All modes get same coupling
phonon_modes.append({
    'coupling': D_molecule,  # Scalar - WRONG!
    'is_local': True
})
```

### Required Fix ✅
```python
# Mode-dependent coupling based on wavefunction overlap
D_nm = np.zeros(Nm)
for im in range(Nm):
    n, m = transverse_modes.modes[im]
    psi = transverse_modes.get_wavefunction(n, m, y_mol, z_mol)
    D_nm[im] = D₀ × psi**2

phonon_modes.append({
    'coupling': D_nm,  # Array[Nm] - CORRECT!
    'is_local': True
})
```

**This is the spatial selectivity that makes quasi-3D meaningful!**

---

## Conclusion

Most assumptions are reasonable approximations for:
- Uniform rectangular devices
- Weak-to-moderate coupling
- Small molecules at low concentration
- Smooth interfaces

**The one critical exception** is mode-dependent molecular coupling, which must be fixed for correct molecular sensing physics.

**Future improvements** (in order of priority):
1. ✅ **IMMEDIATE**: Mode-dependent molecular coupling (D_nm array)
2. Hybrid mode selection (7-8x speedup)
3. Mode mixing at interfaces (more accurate transmission)
4. Molecular orientation (anisotropic coupling)
5. Energy-dependent effective mass (better high-energy physics)

---

**References**:
- explain.md: Project requirements and mode selection strategy
- QUASI_3D_IMPLEMENTATION_COMPLETE.md: Implementation details
- Datta (2000): NEGF formalism with transverse modes

*Document Version: 1.0*
*Date: 2026-01-13*

