# Quantum Electronic Nose: Technical Documentation

## Comprehensive Guide for Paper Writing

This document provides a complete technical reference for the quantum electronic nose (e-nose) simulation framework. It covers the physics, mathematical formalism, implementation details, configurable parameters, and research capabilities.

---

## Table of Contents

1. [Overview](#1-overview)
2. [Physical System](#2-physical-system)
3. [Theoretical Framework](#3-theoretical-framework)
4. [Mathematical Formulation](#4-mathematical-formulation)
5. [Numerical Implementation](#5-numerical-implementation)
6. [Configurable Parameters](#6-configurable-parameters)
7. [Material Systems](#7-material-systems)
8. [Molecular Database](#8-molecular-database)
9. [Output and Analysis](#9-output-and-analysis)
10. [Research Capabilities](#10-research-capabilities)
11. [Physical Constants](#11-physical-constants)
12. [Code Architecture](#12-code-architecture)

---

## 1. Overview

### 1.1 Concept

The quantum e-nose is a gas sensing device based on **Inelastic Electron Tunneling Spectroscopy (IETS)** through a **Resonant Tunneling Diode (RTD)**. When an odorant molecule adsorbs onto the device surface, its vibrational modes couple to tunneling electrons, creating characteristic signatures in the second derivative of the current-voltage curve (d²I/dV²). These signatures serve as molecular "fingerprints" enabling selective gas identification.

### 1.2 Key Physics

- **Resonant tunneling**: Electrons tunnel through quantum-confined states in a double-barrier heterostructure
- **Inelastic scattering**: Electrons exchange energy with molecular vibrations (phonons)
- **IETS signature**: Vibrational modes appear as peaks in d²I/dV² at bias V ≈ ℏω/e
- **Spatial selectivity**: Molecule position affects coupling strength to tunneling electrons
- **Quasi-3D transport**: Multiple transverse conduction channels contribute to current

### 1.3 Simulation Approach

The framework uses the **Non-Equilibrium Green's Function (NEGF)** formalism with **Self-Consistent Born Approximation (SCBA)** to capture:
- Coherent quantum tunneling through barriers
- Inelastic electron-phonon scattering
- Non-equilibrium carrier distributions
- Mode-resolved transport in quasi-3D

---

## 2. Physical System

### 2.1 Device Structure

The RTD consists of a layered heterostructure:

```
┌─────────────────────────────────────────────────────────┐
│                    COLLECTOR CONTACT                     │
│                  (n-doped semiconductor)                 │
├─────────────────────────────────────────────────────────┤
│                  COLLECTOR BARRIER                       │
│                    (wide-gap insulator)                  │
├─────────────────────────────────────────────────────────┤
│                    QUANTUM WELL                          │
│                  (narrow-gap semiconductor)              │
├─────────────────────────────────────────────────────────┤
│                   EMITTER BARRIER                        │←── Molecule adsorbs here
│                    (wide-gap insulator)                  │
├─────────────────────────────────────────────────────────┤
│                    EMITTER CONTACT                       │
│                  (n-doped semiconductor)                 │
└─────────────────────────────────────────────────────────┘
```

### 2.2 Material Systems

**III-V Semiconductors (Legacy/Reference):**
- GaAs/AlAs: Well-characterized RTDs
- InGaAs/InAlAs: Lower effective mass, higher mobility

**Oxide Semiconductors (BEOL-Compatible):**
- In₂O₃/Al₂O₃: Primary candidate, ALD process
- IGZO/Al₂O₃: Amorphous, excellent uniformity
- ZnO/Al₂O₃: Mature process technology
- SnO₂/Al₂O₃: Backup candidate

### 2.3 Quantum Well States

The ground state energy in a quantum well of width L:

$$E_1 = \frac{\pi^2 \hbar^2}{2m^* L^2}$$

This determines the resonance bias voltage:

$$V_{\text{res}} \approx 2E_1/e$$

### 2.4 Molecule Position

The odorant molecule can adsorb at different locations:
1. **Emitter barrier**: Highest sensitivity (tunneling injection region)
2. **Quantum well**: Moderate sensitivity (resonant state region)
3. **Collector barrier**: Lower sensitivity (tunneling exit region)

---

## 3. Theoretical Framework

### 3.1 NEGF Formalism

The Non-Equilibrium Green's Function approach provides a rigorous framework for quantum transport including dissipation.

**Key quantities:**

| Symbol | Name | Physical Meaning |
|--------|------|------------------|
| G^R | Retarded Green's function | Causal response function |
| G^A | Advanced Green's function | (G^R)† |
| G^< | Lesser Green's function | Electron correlation |
| G^> | Greater Green's function | Hole correlation |
| Σ | Self-energy | Interaction effects |
| A | Spectral function | Density of states |
| Γ | Broadening function | Level width |

### 3.2 Self-Consistent Born Approximation

SCBA treats electron-phonon interactions to lowest order in perturbation theory while maintaining self-consistency:

1. **Initialize**: Σ_S = 0
2. **Compute Green's function**: G^R = [(E+iη)I - H - Σ_1 - Σ_2 - Σ_S]^(-1)
3. **Compute correlations**: n = G^< , p = G^>
4. **Update phonon self-energy**: Σ_S^{in/out} from n, p
5. **Mix and iterate** until convergence

### 3.3 Quasi-3D Transport

For devices with transverse confinement, electrons occupy discrete transverse modes:

$$E_{\text{total}} = E_{\text{longitudinal}} + \varepsilon_{nm}$$

where the transverse mode energy is:

$$\varepsilon_{nm} = \frac{\hbar^2\pi^2}{2m^*}\left[\frac{n^2}{L_y^2} + \frac{m^2}{L_z^2}\right]$$

Each mode contributes independently to total current:

$$I = \sum_{n,m} I_{nm}$$

---

## 4. Mathematical Formulation

### 4.1 Hamiltonian

**Tight-binding representation** on a 1D grid with spacing a:

$$H = \sum_i \epsilon_i |i\rangle\langle i| - \sum_i t_i (|i\rangle\langle i+1| + |i+1\rangle\langle i|)$$

**On-site energy:**
$$\epsilon_i = E_c(x_i) + U_{\text{ext}}(x_i) + t_{i-1} + t_i$$

**Hopping parameter** (position-dependent effective mass):
$$t_i = \frac{\hbar^2}{2m^*_{\text{avg}} a^2}$$

where $m^*_{\text{avg}} = (m^*_i + m^*_{i+1})/2$

### 4.2 Contact Self-Energy

For semi-infinite 1D leads with uniform hopping t and on-site energy U:

**Propagating regime** (|E - U| < 2t):
$$\Sigma = t \cdot e^{-ika}$$
where $\cos(ka) = (E - U)/(2t)$

**Evanescent regime** (|E - U| > 2t):
$$\Sigma = t \cdot e^{-\kappa a}$$
where $\cosh(\kappa a) = |E - U|/(2t)$

**Broadening function:**
$$\Gamma = i(\Sigma - \Sigma^\dagger)$$

### 4.3 Green's Function

**Retarded Green's function:**
$$G^R(E) = [(E + i\eta)I - H - \Sigma_1 - \Sigma_2 - \Sigma_S]^{-1}$$

**Spectral function:**
$$A(E) = i[G^R - G^A] = i[G - G^\dagger]$$

**Contact-resolved spectral functions:**
$$A_1 = G\Gamma_1 G^\dagger, \quad A_2 = G\Gamma_2 G^\dagger$$

### 4.4 Correlation Functions

**Electron correlation (occupied states):**
$$G^< = G^R \Sigma^{in} G^A$$

**Hole correlation (empty states):**
$$G^> = G^R \Sigma^{out} G^A$$

**Relation to spectral function:**
$$A = G^> - G^< = G^> + G^<$$ (with proper normalization)

**Total inscattering:**
$$\Sigma^{in} = f_1 \Gamma_1 + f_2 \Gamma_2 + \Sigma_S^{in}$$

where $f_{1,2}$ are Fermi-Dirac distributions at contacts 1 and 2.

### 4.5 Phonon Self-Energy

For electron-phonon coupling strength D and phonon energy ℏω:

**Inscattering (electrons gaining states):**
$$\Sigma_S^{in}(E) = D^2 [(n_B + 1) \cdot n(E + \hbar\omega) + n_B \cdot n(E - \hbar\omega)]$$

**Outscattering (electrons losing states):**
$$\Sigma_S^{out}(E) = D^2 [n_B \cdot p(E + \hbar\omega) + (n_B + 1) \cdot p(E - \hbar\omega)]$$

**Bose-Einstein distribution:**
$$n_B = \frac{1}{e^{\hbar\omega/k_BT} - 1}$$

**Local phonon modes** use a projection operator P to restrict coupling to specific sites:
$$\Sigma_S^{\text{local}} = P \cdot \Sigma_S^{\text{bulk}} \cdot P$$

### 4.6 Mode-Dependent Coupling

For quasi-3D transport, molecular coupling depends on transverse wavefunction overlap:

$$D_{nm} = D_{\text{base}} \times \frac{|\psi_{nm}(y_{\text{mol}}, z_{\text{mol}})|^2}{\max_{nm}|\psi_{nm}|^2}$$

where the transverse wavefunction is:

$$\psi_{nm}(y,z) = \frac{2}{\sqrt{L_y L_z}} \sin\left(\frac{n\pi y}{L_y}\right) \sin\left(\frac{m\pi z}{L_z}\right)$$

### 4.7 Current Calculation

**Landauer-Büttiker formula with NEGF:**
$$I = \frac{2e^2}{h} \int dE \, T(E) [f_1(E) - f_2(E)]$$

**Transmission function:**
$$T(E) = \text{Tr}[\Gamma_1 G \Gamma_2 G^\dagger]$$

**Mode-resolved current (quasi-3D):**
$$I = \frac{2e^2}{h} \sum_{nm} \int dE \, T_{nm}(E) [f_1(E) - f_2(E)]$$

### 4.8 IETS Signal

**Differential conductance:**
$$G(V) = \frac{dI}{dV}$$

**IETS spectrum:**
$$\frac{d^2I}{dV^2}$$

Peaks in d²I/dV² appear at bias voltages corresponding to vibrational energies:
$$eV = \hbar\omega_{\text{vib}}$$

---

## 5. Numerical Implementation

### 5.1 Discretization

**Grid spacing:** a = 0.1 - 0.3 nm (typically)

**Number of points:** $N_p = L_{\text{total}}/a$ (typically 50-100)

**Energy grid:** Uniform spacing, typically 30-200 points

**Bias sweep:** Uniform spacing, typically 5-50 points

### 5.2 SCBA Iteration

```
Algorithm: Self-Consistent Born Approximation
───────────────────────────────────────────────
Input: H, Σ₁(E), Σ₂(E), phonon modes, μ₁, μ₂, T
Output: G, Σ_S^{in}, Σ_S^{out}, I-V curve

1. Initialize Σ_S = 0
2. For iteration = 1 to max_iter:
   a. For each energy E:
      - Compute G(E) = [(E+iη)I - H - Σ₁ - Σ₂ - Σ_S]⁻¹
      - Compute A(E) = i[G - G†]
      - Compute n(E) = G · Σ^{in} · G†
      - Compute p(E) = A - n
   b. For each phonon mode:
      - Compute Σ_S^{in,out} from n, p
   c. Mix: Σ_S ← α·Σ_S^{new} + (1-α)·Σ_S^{old}
   d. Check convergence: |ΔΣ_S| < tolerance
3. Compute current from converged G
```

### 5.3 Convergence Parameters

| Parameter | Symbol | Typical Value | Effect |
|-----------|--------|---------------|--------|
| Mixing factor | α | 0.3 | Stability vs. speed |
| Tolerance | tol | 10⁻² - 10⁻⁵ | Accuracy |
| Max iterations | max_iter | 30-100 | Computation time |
| Broadening | η | 10⁻⁴ eV | Numerical stability |

### 5.4 Hybrid Mode Selection

For speedup in quasi-3D calculations:
- **Inelastic modes**: Top N modes by importance metric (full SCBA)
- **Coherent modes**: Remaining modes (ballistic transport)

**Importance metric:**
$$\text{importance}_{nm} = \alpha \cdot w_{\text{coupling}} + (1-\alpha) \cdot f_{\text{FD}}(\varepsilon_{nm})$$

where:
- $w_{\text{coupling}} = |\psi_{nm}(y_{\text{mol}}, z_{\text{mol}})|^2 / \max|\psi|^2$
- $f_{\text{FD}} = 1/(1 + e^{\varepsilon_{nm}/k_BT})$
- α = 0.3 (typical)

**Speedup:** ~(Total modes)/(Inelastic modes)

---

## 6. Configurable Parameters

### 6.1 Device Parameters

| Parameter | Variable | Default | Units | Description |
|-----------|----------|---------|-------|-------------|
| Device name | `device_name` | "GaAs_AlAs_symmetric" | - | Selects from device library |
| Grid spacing | `grid_spacing` | 0.12×10⁻⁹ | m | Spatial discretization |

### 6.2 Molecule Parameters

| Parameter | Variable | Default | Units | Description |
|-----------|----------|---------|-------|-------------|
| Molecule name | `molecule_name` | "Benzene" | - | Selects from molecular database |
| Coupling scale | `molecular_coupling_scale` | 1.0 | - | Multiplier for all molecular couplings |

### 6.3 Phonon Parameters

| Parameter | Variable | Default | Units | Description |
|-----------|----------|---------|-------|-------------|
| Bulk phonon energy | `bulk_phonon_energy` | 0.036 | eV | Material LO phonon |
| Bulk coupling | `bulk_phonon_coupling` | 0.010 | eV | Bulk electron-phonon coupling |

### 6.4 Energy Grid

| Parameter | Variable | Default | Units | Description |
|-----------|----------|---------|-------|-------------|
| E minimum | `E_min` | -0.3 | eV | Lower bound of energy grid |
| E maximum | `E_max` | 1.5 | eV | Upper bound of energy grid |
| E points | `E_points` | 200 | - | Number of energy points |

### 6.5 Bias Sweep

| Parameter | Variable | Default | Units | Description |
|-----------|----------|---------|-------|-------------|
| V minimum | `V_min` | 0.0 | V | Starting bias |
| V maximum | `V_max` | 0.5 | V | Ending bias |
| V points | `V_points` | 26 | - | Number of bias points |

### 6.6 Temperature

| Parameter | Variable | Default | Units | Description |
|-----------|----------|---------|-------|-------------|
| Temperature | `temperature` | 300 | K | Operating temperature |

### 6.7 SCBA Parameters

| Parameter | Variable | Default | Units | Description |
|-----------|----------|---------|-------|-------------|
| Max iterations | `scba_max_iter` | 50 | - | Maximum SCBA iterations |
| Tolerance | `scba_tolerance` | 10⁻⁴ | - | Convergence criterion |
| Mixing | `scba_mixing` | 0.3 | - | Linear mixing parameter |

### 6.8 Quasi-3D Parameters

| Parameter | Variable | Default | Units | Description |
|-----------|----------|---------|-------|-------------|
| Enable multimode | `use_multimode` | False | - | Enable quasi-3D transport |
| Transverse width y | `Ly` | 1.0×10⁻⁶ | m | Device width in y direction |
| Transverse width z | `Lz` | 1.0×10⁻⁶ | m | Device width in z direction |
| Max mode index y | `n_max` | 3 | - | Highest mode in y |
| Max mode index z | `m_max` | 3 | - | Highest mode in z |
| Transverse mass | `m_trans` | None (device m*) | kg | Effective mass for transverse motion |

### 6.9 Hybrid Mode Selection

| Parameter | Variable | Default | Units | Description |
|-----------|----------|---------|-------|-------------|
| Enable hybrid | `use_hybrid` | False | - | Enable hybrid mode selection |
| Inelastic modes | `n_inelastic` | 4 | - | Number of modes with full SCBA |

---

## 7. Material Systems

### 7.1 III-V Semiconductors

| Material | m*/m₀ | Ec (eV) | ε_r | ℏω_LO (meV) |
|----------|-------|---------|-----|-------------|
| GaAs | 0.067 | 0.0 | 12.9 | 36 |
| AlAs | 0.15 | 0.57 | 10.1 | 50 |
| InGaAs | 0.041 | -0.15 | 13.9 | 30 |
| InAlAs | 0.075 | 0.52 | 12.2 | 43 |

### 7.2 Oxide Semiconductors (BEOL-Compatible)

| Material | m*/m₀ | Ec (eV) | ε_r | ℏω_LO (meV) | Process T (°C) |
|----------|-------|---------|-----|-------------|----------------|
| In₂O₃ | 0.30 | 0.0 | 9.0 | 70 | 350 |
| IGZO | 0.34 | 0.0 | 10.0 | 60 | 300 |
| ZnO | 0.28 | 0.0 | 8.5 | 72 | 350 |
| SnO₂ | 0.25 | 0.0 | 9.0 | 75 | 300 |

### 7.3 Barrier Materials

| Material | m*/m₀ | Ec (eV) | ε_r | ℏω_LO (meV) |
|----------|-------|---------|-----|-------------|
| Al₂O₃ | 0.45 | 2.8 | 9.0 | 100 |
| HfO₂ | 0.50 | 2.5 | 25.0 | 95 |
| SiO₂ | 0.50 | 3.2 | 3.9 | 120 |

### 7.4 Available Device Configurations

| Device Name | Well | Barrier | Notes |
|-------------|------|---------|-------|
| GaAs_AlAs_symmetric | GaAs 2nm | AlAs 1.5nm | Validated baseline |
| GaAs_AlAs_asymmetric | GaAs 2nm | AlAs 3.0/1.0nm | Optimized sensing |
| GaAs_AlAs_thin | GaAs 2.5nm | AlAs 1.0nm | High current |
| GaAs_AlAs_wide_well | GaAs 4.0nm | AlAs 1.5nm | Low resonance |
| InGaAs_InAlAs_symmetric | InGaAs 2nm | InAlAs 1.5nm | Low m* system |
| In2O3_Al2O3_symmetric | In₂O₃ 2.5nm | Al₂O₃ 2.0nm | BEOL primary |
| IGZO_Al2O3_symmetric | IGZO 2.5nm | Al₂O₃ 2.0nm | Amorphous |
| ZnO_Al2O3_symmetric | ZnO 2.5nm | Al₂O₃ 2.0nm | Mature process |

---

## 8. Molecular Database

### 8.1 Perceptual Classes

Based on Pandey et al. (Scientific Reports 2021):

| Class | Molecules | Characteristic |
|-------|-----------|----------------|
| Aromatic | Benzene, Anthracene, Thiophene | Strong/weak benzene-like |
| Roasted Coffee | Furan, Furan_Methanethiol | Coffee aroma |
| Moth-ball | Naphthalene, Tetralin, Fluorene | PAH-like |
| Fruity | Heptanal, N-Amyl_Butyrate, Gamma_Octalactone | Citrus/pear/coconut |
| Musk | Civetone, Moxalone, Galaxolide, Helvetolide | Animal/synthetic musk |
| Garlicky | Benzyl_Mercaptan, Allyl_Thiol, Dimethyl_Sulfide, Diallyl_Disulfide, Allicin | Sulfurous/garlic |

### 8.2 Benzene Example

| Mode | Wavenumber (cm⁻¹) | Energy (meV) | Coupling (meV) |
|------|-------------------|--------------|----------------|
| 1 | 399 | 49.5 | 5.0 |
| 2 | 637 | 79.0 | 8.0 |
| 3 | 1084 | 134.4 | 10.0 |
| 4 | 1485 | 184.1 | 8.0 |
| 5 | 3189 | 395.4 | 5.0 |

### 8.3 Conversion Factors

- **Wavenumber to energy:** E (meV) = ν (cm⁻¹) × 0.123984
- **Energy to bias:** V ≈ E/e (for IETS peak position)

---

## 9. Output and Analysis

### 9.1 Primary Outputs

| Output | Variable | Units | Description |
|--------|----------|-------|-------------|
| Voltage array | `V_array` | V | Bias sweep points |
| Current array | `I_array` | A | Current at each bias |
| Conductance | `dIdV` | S | First derivative |
| IETS spectrum | `d2IdV2` | S/V | Second derivative |
| Green's function | `G` or `G_nm` | - | Complex, energy-resolved |
| Transmission | `T` | - | Energy-resolved |

### 9.2 IETS Analysis

**Peak extraction:**
- Identify local maxima in d²I/dV²
- Extract peak positions (voltage → energy)
- Extract peak intensities

**Fingerprint comparison:**
- Euclidean distance
- Correlation coefficient
- Cosine similarity

**Classification:**
- Nearest-neighbor based on fingerprint distance
- Per-class accuracy metrics

### 9.3 Selectivity Metrics

**Selectivity score:** 1 - (mean off-diagonal similarity)

Higher selectivity = better discrimination between different molecules

---

## 10. Research Capabilities

### 10.1 What Can Be Studied

**Device Physics:**
- Resonant tunneling characteristics (peak position, width, amplitude)
- NDR (Negative Differential Resistance) regime
- Temperature dependence of transport
- Transverse mode contributions

**Molecular Sensing:**
- IETS signatures for different molecules
- Peak positions vs. vibrational energies
- Molecular fingerprint discrimination
- Classification accuracy by perceptual class

**Design Optimization:**
- Barrier thickness (tunneling rate)
- Well width (resonance energy)
- Material selection (effective mass, barrier height)
- Device asymmetry (position sensitivity)

**Spatial Selectivity:**
- Molecule position effects
- Coupling strength vs. location
- Multi-position sensing arrays

**Quasi-3D Effects:**
- Mode-resolved current contributions
- Transverse confinement effects
- Mode-dependent molecular coupling

### 10.2 Parameter Studies

**Easily varied:**
- Temperature (0-500 K)
- Bias range (0-several V)
- Molecule type (20+ in database)
- Device geometry (layer thicknesses)
- Grid resolution (accuracy vs. speed)

**Requires code modification:**
- New material properties
- New molecules (add to database)
- Custom device structures

### 10.3 Typical Studies

1. **I-V characterization**: Sweep bias, plot I-V and find resonance
2. **IETS measurement**: High-resolution bias sweep, extract d²I/dV²
3. **Temperature study**: Multiple temperatures, observe thermal broadening
4. **Molecule comparison**: Same device, different molecules
5. **Device comparison**: Different devices, same molecule
6. **Position study**: Vary molecule location, compare sensitivity

---

## 11. Physical Constants

| Constant | Symbol | Value | Units |
|----------|--------|-------|-------|
| Electron charge | e | 1.602176634×10⁻¹⁹ | C |
| Reduced Planck | ℏ | 1.054571817×10⁻³⁴ | J·s |
| Planck constant | h | 6.62607015×10⁻³⁴ | J·s |
| Electron mass | m₀ | 9.10938356×10⁻³¹ | kg |
| Boltzmann | k_B | 1.380649×10⁻²³ | J/K |
| Boltzmann (eV) | k_B | 8.617333262145×10⁻⁵ | eV/K |
| Permittivity | ε₀ | 8.854187817×10⁻¹² | F/m |

**Thermal energy at 300 K:** k_B T = 25.85 meV

**Conductance quantum:** G₀ = 2e²/h = 7.748×10⁻⁵ S

---

## 12. Code Architecture

### 12.1 Directory Structure

```
quantum_enose/
├── config/
│   ├── device_library.py      # Device definitions and material properties
│   └── molecular_database.py  # Molecular vibrational data (Pandey 2021)
├── core/
│   ├── hamiltonian.py         # Discretization, Hamiltonian construction
│   ├── green_functions.py     # NEGF solver, spectral functions
│   ├── self_energy.py         # Contact and phonon self-energies
│   ├── scba_solver.py         # SCBA iteration, current calculation
│   ├── transverse_modes.py    # Quasi-3D mode infrastructure
│   ├── hybrid_modes.py        # Hybrid mode selection
│   └── scba_solver_hybrid.py  # Hybrid SCBA implementation
├── run/
│   ├── run_single_molecule.py # Single molecule simulation pipeline
│   └── run_batch.py           # Batch screening of multiple molecules
├── analysis/
│   └── iets_analysis.py       # Fingerprint extraction, classification
└── [simulation scripts]       # Various run scripts
```

### 12.2 Key Functions

**Device Setup:**
- `discretize_device(device_config, grid_spacing)` → grid dict
- `build_hamiltonian(grid, U_external)` → H, t arrays

**Transport:**
- `retarded_greens_function(E, H, Σ₁, Σ₂, Σ_S, η)` → G matrix
- `spectral_function(G)` → A matrix
- `broadening_function(Σ)` → Γ matrix

**SCBA:**
- `scba_iteration(...)` → result dict (1D)
- `scba_iteration_multimode(...)` → result dict (quasi-3D)
- `scba_iteration_hybrid(...)` → result dict (hybrid quasi-3D)

**Current:**
- `compute_current(result, Γ₁, Γ₂, E, μ₁, μ₂, T)` → I, I(E)
- `compute_current_multimode(result, E, μ₁, μ₂, T)` → I, I(E), I_{nm}(E)

**IETS:**
- `compute_iets(V, I)` → dI/dV, d²I/dV²

### 12.3 Typical Workflow

```python
# 1. Load configuration
config = SimulationConfig()
config.device_name = "GaAs_AlAs_symmetric"
config.molecule_name = "Benzene"

# 2. Setup device
device = get_device(config.device_name)
grid = discretize_device(device, config.grid_spacing)
H, t = build_hamiltonian(grid)

# 3. Setup molecule
molecule = get_molecule(config.molecule_name)
phonon_modes = build_phonon_modes(molecule, config, grid, device)

# 4. Define contact self-energies
def Sigma1(E): return contact_self_energy_matrix(E, grid, t, 'left')
def Sigma2(E): return contact_self_energy_matrix(E, grid, t, 'right')

# 5. Bias sweep
for V in V_array:
    mu1, mu2 = +V/2, -V/2
    result = scba_iteration(E_array, H, Sigma1, Sigma2, mu1, mu2,
                           temperature, phonon_modes, grid)
    I, _ = compute_current(result, Gamma1, Gamma2, E_array, mu1, mu2, T)

# 6. Compute IETS
dIdV, d2IdV2 = compute_iets(V_array, I_array)
```

---

## References

1. **NEGF Formalism:** Datta, S. "Electronic Transport in Mesoscopic Systems" (Cambridge, 1995)
2. **SCBA:** Lake et al., "Single and multiband modeling of quantum electron transport through layered semiconductor devices," J. Appl. Phys. 81, 7845 (1997)
3. **RTD Gas Sensing:** Patil et al., "Ultra-sensitive room-temperature gas sensing using resonant tunneling diodes," Sensors and Actuators B: Chemical (2018)
4. **Molecular Database:** Pandey et al., "Vibration-based biomimetic odor classification," Scientific Reports 11, 11119 (2021)
5. **Oxide RTDs:** Intel/Berkeley VLSI 2021, "BEOL-compatible oxide RTDs"

---

## Appendix: Quick Reference

### Energy Conversions
- 1 meV = 11.6 K (thermal)
- 1 meV = 8.07 cm⁻¹ (wavenumber)
- 1 eV = 1.602×10⁻¹⁹ J

### Typical Values
- GaAs m* = 0.067 m₀
- Room temperature k_BT = 26 meV
- AlAs barrier height = 0.57 eV
- GaAs LO phonon = 36 meV

### Computational Cost (approximate)
- 1D SCBA: ~1 second per bias point
- Quasi-3D (9 modes): ~30 seconds per bias point
- Hybrid (4 inelastic): ~15 seconds per bias point

---

*Document generated for quantum e-nose paper writing. Last updated: January 2026*
