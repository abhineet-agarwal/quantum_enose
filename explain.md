# Research Progress Summary: Quantum Biomimetic Electronic Nose

## Project Overview
PhD research at IIT Bombay focused on developing a quantum biomimetic electronic nose based on resonant tunneling diodes (RTDs) and inelastic electron tunneling spectroscopy (IETS) for room-temperature odorant detection.

## Core Technical Achievement

### Validated 3D Hybrid NEGF-SCBA Simulation Framework
Successfully implemented and validated a sophisticated simulation framework combining:
- **3D transverse mode treatment**: Hard-wall quantization in y-z plane
- **1D longitudinal discretization**: Finite-difference method on discrete lattice (spacing a = 0.3 nm)
- **Hybrid coherent/inelastic approach**: Energy-based mode selection criterion
- **Self-consistent Born approximation (SCBA)**: Local phonon scattering with spatial projection operators

**Key Performance Metrics:**
- Peak current density: 5.83×10¹⁰ A/m² at 0.684V
- SCBA convergence: 100% success rate, typically 2-6 iterations
- Computational speedup: 7-8x with hybrid approach (vs full inelastic treatment)
- Accuracy maintained: 90-95% relative to full calculation

### Physical Validation
Framework produces proper RTD characteristics:
- Resonant tunneling peak in forward bias
- Negative differential resistance region
- Satellite peaks from phonon-assisted tunneling
- Proper contact self-energies (surface Green's functions)

## Material Selection

### Comprehensive Analytical Baseline Study
Evaluated 8 RTD material stacks using transfer matrix methods with established formulas (Tsu-Esaki, Schulman, Coon-Liu):

**Primary Candidate: In₂O₃/HfO₂**
- Quality factor Q = 87.9
- Barrier height ΔEc = 0.70 eV
- Well-validated from Intel/Purdue collaboration papers
- BEOL-compatible processing (<400°C)

**Alternative candidates analyzed:**
- IGZO/Al₂O₃, ZnO/Al₂O₃, SnO₂/Al₂O₃
- Trade-offs between Q-factor, barrier height, and BEOL compatibility

## Current Implementation Status

### Working Simulation Components

1. **Hamiltonian Construction**
   - Variable-mass tight-binding with position-dependent hopping
   - Transverse mode energies: εₙₘ = (ℏ²π²/2m*)[(n²/Wᵧ²) + (m²/Wᵧ²)]
   - Proper treatment of band-edge offsets at interfaces

2. **Contact Modeling**
   - Analytic surface self-energies for semi-infinite 1D chains
   - Σₗ/ᵣ(1,1) = -t exp(±ika) for open boundary conditions
   - Energy-dependent broadening: Γ = i(Σ - Σ†)

3. **Coherent Transport** 
   - Caroli formula: T(E) = Tr{ΓₗGΓᵣG†}
   - Proper integration with Fermi-Dirac distributions
   - Validated against ballistic conductor quantization

4. **Inelastic SCBA Loop**
   ```
   For each energy E:
     1. Calculate G(E) with current Σⁱⁿ/ᵒᵘᵗ
     2. Compute A(E), n(E), p(E) matrices
     3. Apply energy shifts for phonon coupling
     4. Project to molecular sites: P·Σ·P
     5. Update Σⁱⁿ/ᵒᵘᵗ with mixing parameter
     6. Check convergence: ||Σₙₑᵥ - Σₒₗ𝒹||₁
   ```

5. **Mode Selection Strategy**
   - Rank by thermal-weighted spatial overlap: |ψₙₘ(ymol, zmol)|² × fFD(εₙₘ)
   - Select top M modes for full inelastic treatment
   - Coherent sum over remaining modes
   - Exploits spatial locality: modes with nodes at molecule → negligible coupling

### Code Structure (Python/NumPy)
- **negf_hybrid_selective_inelastic.py**: Core physics library (~500 lines)
  - Transverse mode generators
  - Hamiltonian builders with local perturbations
  - Contact self-energies
  - SCBA solver with local projection
  - Hybrid driver with mode selection

- **run_simulation.py**: Batch runner (~150 lines)
  - Problem definition (device stack, molecule placement)
  - Bias sweep with self-consistent iteration
  - Post-processing: dI/dV, d²I/dV² calculation
  - CSV export and plotting

**Key Implementation Details:**
- Molecule placement: Central transverse position (Wᵧ/2, Wᵧ/2)
- Longitudinal localization: Central sites with neighbor_radius parameter
- Coupling decay: 1/(r+1) for neighbors
- Phonon parameters: Debye strength D, mode energies ℏω

## Next Phase: Molecular Vibrational Coupling

### Implementation Strategy (Modular Approach)

**Batch 1: Phonon Mode Extension**
- Modify SCBA to handle multiple phonon branches (bulk + molecular)
- Implement proper energy indexing for multiple ℏω values
- Test with 2-3 modes: 1 bulk + 1-2 molecular

**Batch 2: Spatial Localization**
- Develop projection operators for molecular sites
- Implement neighbor radius with decay function
- Validate localization vs broadening trade-off

**Batch 3: Molecular Database**
- Parse Pandey 2021 reference for all 20 odorants
- Extract vibrational frequencies (IR/Raman data)
- Organize by perceptual class:
  - Aromatic (benzaldehyde, styrene, etc.)
  - Roasted Coffee (pyridine derivatives)
  - Moth-ball (naphthalene, camphor)
  - Fruity (esters, ethyl butyrate)
  - Musk (macrocyclic compounds)
  - Garlicky (sulfur-containing)

**Batch 4: Selective Detection**
- Systematic sweep: vary molecular ℏω from 0.05-0.25 eV
- Compute IETS spectra: d²I/dV²(V)
- Demonstrate peak position ↔ vibrational energy mapping
- Validate selectivity across perceptual classes

### Critical Physics Considerations

1. **Mode Energy Calculation**
   - Must use proper quantum confinement scaling
   - Incorrect formula can cause order-of-magnitude errors
   - Cross-validate with analytical expressions

2. **Thermal Population**
   - Essential for mode ranking: nBE(ℏω) = 1/(exp(ℏω/kBT) - 1)
   - High-energy modes contribute negligibly despite strong coupling
   - Validates hybrid approach at T = 300K

3. **SCBA Convergence**
   - Mixing parameter α = 0.5 typically optimal
   - May need adjustment for strong coupling (D > 0.2 eV)
   - Monitor iteration count and residual norm

4. **Material Parameters**
   - Bulk phonon: ℏω ≈ 35 meV (LO-phonon in oxides)
   - Molecular vibrations: ℏω = 50-250 meV (400-2000 cm⁻¹)
   - Coupling strength: D ~ 0.05-0.15 eV (fitting parameter)

## Key Insights from Literature

### Patil et al. 2018 (Nature Scientific Reports)
- Demonstrated RTD-based IETS for explosive detection
- Key result: Quantum confinement in well provides energy filtering at 300K
- Avoided thermal broadening problem of traditional MIM IETS (cryogenic requirement)
- Important finding: Only emitter barrier effective for sensing (not well or collector)
  → Suggests asymmetric design with wider emitter barrier

### Pandey et al. 2021
- Established vibration-based odor classification for 20 molecules
- Organized by perceptual classes (not just chemical structure)
- Provides vibrational database for benchmark testing

### Datta Tutorial (Superlattices & Microstructures 2000)
- Clear exposition of NEGF for device modeling
- Self-energy concept for open boundaries
- Self-consistent Schrödinger-Poisson framework
- Emphasis on transport + electrostatics interplay

## Validation Methodology

### Current Level (Complete)
✅ Analytical baselines for coherent transport
✅ Proper I-V with resonant peak and NDR
✅ SCBA convergence for bulk phonons
✅ Mode selection energy-based ranking
✅ Spatial localization of phonon self-energies

### Next Level (In Progress)
🔄 Molecular vibration coupling
🔄 IETS peak resolution (d²I/dV²)
🔄 Selectivity demonstration across odorants
🔄 Multi-mode interference effects

### Future Validation
⏳ Comparison with experimental RTD-IETS data (if available)
⏳ Benchmark against NEMO1D/NEMO3D results
⏳ Sensitivity analysis: minimum detectable concentration
⏳ Temperature dependence: verify 300K operation

## Technical Challenges Addressed

1. **Computational Efficiency**
   - Problem: Full 3D NEGF-SCBA is prohibitively expensive
   - Solution: Hybrid mode selection based on physics (spatial locality)
   - Result: 7-8x speedup with <10% accuracy loss

2. **SCBA Stability**
   - Problem: Self-consistency loops can diverge for strong coupling
   - Solution: Proper mixing (α = 0.5), initialization from coherent result
   - Result: 100% convergence rate in test cases

3. **Mode Energy Scaling**
   - Problem: Incorrect quantization formulas give wrong transverse energies
   - Solution: Careful derivation from effective mass equation
   - Result: Validated against analytical expressions

4. **Material Parameter Uncertainty**
   - Problem: Limited data for oxide semiconductor interfaces
   - Solution: Leverage Intel/Purdue collaboration papers, IEEE VLSI symposium
   - Result: Identified In₂O₃/HfO₂ as well-characterized candidate

## Collaboration & Guidance

**Prof. Swaroop Ganguly (Advisor)**
- Critical physics insights on mode selection criteria
- Guidance on transport theory validation
- Emphasis on self-consistent field methods
- Connection to experimental fabrication constraints

**Literature Resources**
- Supriyo Datta's transport theory tutorials (transverse mode summation)
- NEMO documentation (numerical methods)
- Material parameter databases (band offsets, BEOL compatibility)

## Documentation & Reproducibility

- **Code**: Python, ~40 lines per example (tutorial style)
- **Platform**: Laptop-compatible (100×100 Hamiltonian matrices)
- **Visualization**: Matplotlib for I-V, IETS, band diagrams
- **Data Management**: CSV export for post-processing

## Upcoming Milestones

### Short-term (Current Phase)
1. Complete molecular vibrational coupling implementation
2. Generate IETS spectra for test molecules
3. Demonstrate peak position vs ℏω linearity

### Medium-term (Next 3-6 Months)
1. Full 20-odorant database integration
2. Multi-class discrimination validation
3. Optimize device geometry (barrier thickness, well width)

### Long-term (Paper Preparation)
1. Comprehensive performance metrics (sensitivity, selectivity, speed)
2. Comparison with existing e-nose technologies
3. Roadmap for experimental realization
4. SISPAD 2025 conference submission

## Key Formulas Reference

**Transverse Quantization:**
```
ε_nm = (ℏ²π²/2m*) × [(n²/W_y²) + (m²/W_z²)]
```

**Contact Self-Energy:**
```
Σ_L(1,1) = -t exp(ik_1 a)  where  E = E_c + U_1 + 2t(1 - cos(k_1 a))
```

**SCBA Inscattering:**
```
Σ_in(E) = Σ_η D × [(n_BE + 1)·G^n(E - ℏω) + n_BE·G^n(E + ℏω)]
```
With spatial projection: `Σ_in → P · Σ_in · P`

**Current from NEGF:**
```
I = (q²/h) × ∫dE × Tr[Γ_R · G^n - Σ_R^in · G^p]
```

**Quality Factor:**
```
Q = ΔE/Γ  where  ΔE = resonance energy spacing, Γ = level broadening
```

---

## Summary Statement for LLM Context

This research has successfully developed and validated a 3D hybrid NEGF-SCBA framework for RTD-based quantum sensing, achieving 7-8x computational speedup through physics-informed mode selection while maintaining 90-95% accuracy. The framework properly reproduces RTD I-V characteristics with resonant peaks and NDR, demonstrates 100% SCBA convergence for bulk phonon scattering, and implements spatial localization of molecular vibrations through projection operators. Material analysis identified In₂O₃/HfO₂ as the optimal BEOL-compatible stack (Q=87.9, ΔEc=0.70eV).

**Current status**: Core simulation validated; now extending to molecular vibrational coupling for odorant discrimination.

**Next steps**: Implement multi-mode phonon handling, build molecular database from Pandey 2021, demonstrate IETS peak selectivity across 20 odorants organized by perceptual classes.

**Technical foundation**: Modular Python code (~500 lines), laptop-executable, self-consistent Schrödinger-Poisson framework, proper open boundary conditions via analytic self-energies, validated against analytical baselines and literature results.