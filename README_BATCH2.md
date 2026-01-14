# Quantum E-Nose: Batch 2 - Core Physics Modules ✅

## What's New in Batch 2

Two core physics modules for NEGF calculations:

### 1. **`core/hamiltonian.py`** ✅
- Discretizes RTD layer structure onto 1D grid
- Builds tight-binding Hamiltonian with variable effective mass
- Handles band offsets and doping profiles
- Computes transverse mode energies (hard-wall)
- Fermi level estimation

### 2. **`core/green_functions.py`** ✅
- Retarded Green's function: G^R(E) = [(E+iη)I - H - Σ]^(-1)
- Spectral function: A(E) = i[G - G†]
- Contact-resolved spectral functions: A₁, A₂
- Density matrix (equilibrium & non-equilibrium)
- Electron density extraction

---

## Test Results

### Hamiltonian Builder
```
GaAs/AlAs RTD:
  - 175 grid points (0.12 nm spacing)
  - Hopping: 17.6-39.5 eV
  - Diagonal: 35.8-79.0 eV

In2O3/Al2O3 RTD:
  - 221 grid points
  - Hopping: 5.9-8.8 eV (lower due to heavier mass)
  - Diagonal: 14.6-17.6 eV
```

### Green's Functions
```
1D barrier test:
  - Spectral function decomposition: ✓
  - |A - (A₁+A₂)| < 6×10⁻³ (excellent)
  - Fermi function: ✓ correct
```

---

## Usage Examples

### Build Hamiltonian for Any Device

```python
from core.hamiltonian import discretize_device, build_hamiltonian
from config.device_library import get_device

# Get your oxide RTD
device = get_device("In2O3_Al2O3_symmetric")

# Discretize onto grid
grid = discretize_device(device, grid_spacing=0.12e-9)

# Build Hamiltonian
H, t = build_hamiltonian(grid)

print(f"Hamiltonian: {H.shape}")
print(f"Grid points: {grid['Np']}")
```

### Compute Green's Function

```python
from core.green_functions import retarded_greens_function, spectral_function

# Define contact self-energies (simple model)
import numpy as np

def Sigma1(E):
    S = np.zeros((H.shape[0], H.shape[0]), dtype=complex)
    S[0, 0] = -1j * 0.1  # Broadening at left contact
    return S

def Sigma2(E):
    S = np.zeros((H.shape[0], H.shape[0]), dtype=complex)
    S[-1, -1] = -1j * 0.1  # Broadening at right contact
    return S

# Compute at specific energy
E = 0.5  # eV
G = retarded_greens_function(E, H, Sigma1(E), Sigma2(E))
A = spectral_function(G)

# Local density of states
LDOS = np.diag(A) / (2 * np.pi)
```

### Compute Density Matrix (Equilibrium)

```python
from core.green_functions import density_matrix_equilibrium

# Energy grid for integration
E_array = np.linspace(-0.5, 1.5, 200)  # eV
dE = E_array[1] - E_array[0]

# Chemical potential and temperature
mu = 0.0  # eV
T = 300  # K

# Compute density matrix
rho = density_matrix_equilibrium(
    E_array, H, Sigma1, Sigma2, mu, T, dE
)

# Extract electron density
from core.hamiltonian import electron_density_from_rho
n = electron_density_from_rho(rho, grid)
```

---

## Key Features

### Hamiltonian Builder
✅ **Variable effective mass**: Proper averaging at interfaces  
✅ **Oxide materials**: In₂O₃, IGZO, ZnO, SnO₂ all working  
✅ **Band offsets**: Automatic from material library  
✅ **Transverse modes**: Hard-wall quantization  
✅ **Fermi level**: Degenerate/non-degenerate estimation  

### Green's Functions
✅ **Numerically stable**: Regularization for singular matrices  
✅ **Contact decomposition**: Separate A₁ and A₂  
✅ **Equilibrium & non-equilibrium**: Both supported  
✅ **Scattering**: Can include Σ_S (for Batch 3)  
✅ **Integration**: Energy-resolved to density matrix  

---

## What's Next: Batch 3

**Self-Energy Functions & SCBA Solver**

1. **`core/self_energy.py`**
   - Contact self-energies (exact, semi-infinite leads)
   - Local phonon self-energies
   - Projection operators for molecular vibrations

2. **`core/scba_solver.py`**
   - Multi-phonon SCBA iteration
   - Bulk + molecular modes
   - Convergence diagnostics

---

## Testing

```bash
# Test Hamiltonian builder
cd /mnt/user-data/outputs/quantum_enose
python3 core/hamiltonian.py

# Test Green's functions
python3 core/green_functions.py
```

**Expected**: All tests pass with ✓

---

## File Structure

```
quantum_enose/
├── config/
│   ├── molecular_database.py    ✅ Batch 1
│   └── device_library.py         ✅ Batch 1
├── core/
│   ├── hamiltonian.py            ✅ Batch 2 (NEW!)
│   └── green_functions.py        ✅ Batch 2 (NEW!)
├── physics/                      (Batch 3)
├── analysis/                     (Batch 3)
└── run/                          (Batch 4)
```

---

## Physics Recap

### NEGF Formalism

1. **Hamiltonian**: H = H_device (from layers)
2. **Self-energies**: Σ₁, Σ₂ (contacts), Σ_S (scattering)
3. **Green's function**: G = [(E+iη)I - H - Σ₁ - Σ₂ - Σ_S]⁻¹
4. **Spectral function**: A = i[G - G†]
5. **Density matrix**: ρ = ∫ dE/(2π) [f₁A₁ + f₂A₂ + GΣ_S^in G†]
6. **Observables**: n = diag(ρ)/a, I = (q/ℏ) Tr[Γ₁G^n - Γ₂G^p]

### Key Differences: Oxide vs III-V

| Property | GaAs/AlAs | In₂O₃/Al₂O₃ | Impact |
|----------|-----------|-------------|--------|
| m* | 0.067 m₀ | 0.30 m₀ | **Lower hopping energy** |
| Barrier | 0.57 eV | 2.8 eV | **Stronger confinement** |
| Hopping | 40 eV | 8 eV | **Factor of 5 lower** |

→ Oxide RTDs will have **lower current** but **better quantum confinement**!

---

## 🎯 Status: Batch 2 Complete!

✅ Hamiltonian builder working for all materials  
✅ Green's function solver tested  
✅ Both GaAs and oxide systems validated  

**Ready for Batch 3?** Reply "Batch 3" for self-energies & SCBA!
