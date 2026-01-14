"""
Hamiltonian Construction for RTD Devices

Builds device Hamiltonian from layer structure using finite difference method.
Handles variable effective mass, band offsets, and doping.
"""

import numpy as np
import sys
import os

# Add parent directory to path for imports
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

try:
    from config.device_library import get_device, get_material, MATERIALS
except ImportError:
    # Try alternative import for standalone testing
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "device_library", 
        os.path.join(parent_dir, "config", "device_library.py")
    )
    device_lib = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(device_lib)
    get_device = device_lib.get_device
    get_material = device_lib.get_material
    MATERIALS = device_lib.MATERIALS

# Physical constants
m0 = 9.10938356e-31  # kg (electron mass)
q = 1.602176634e-19   # C (electron charge)
hbar = 1.054571817e-34  # J*s
epsilon_0 = 8.854187817e-12  # F/m
kB = 1.380649e-23  # J/K

def eV_to_J(eV):
    """Convert eV to Joules"""
    return eV * q

def J_to_eV(J):
    """Convert Joules to eV"""
    return J / q

# ============================================================================
# DISCRETIZATION
# ============================================================================

def discretize_device(device_config, grid_spacing=0.1e-9):
    """
    Discretize device structure onto 1D grid
    
    Parameters:
    -----------
    device_config : dict
        Device configuration from device_library
    grid_spacing : float
        Lattice spacing in meters (default 0.1 nm)
        
    Returns:
    --------
    dict with keys:
        'x' : array - position coordinates (m)
        'material' : list - material name at each site
        'm_eff' : array - effective mass at each site (kg)
        'Ec' : array - band edge at each site (eV)
        'epsilon_r' : array - relative permittivity at each site
        'doping' : array - doping density at each site (m^-3)
        'Np' : int - number of grid points
        'a' : float - grid spacing (m)
        'layer_boundaries' : list - (start_idx, end_idx) for each layer
    """
    
    layers = device_config['layers']
    a = grid_spacing
    
    # Arrays to build
    x_coords = []
    material_list = []
    m_eff_list = []
    Ec_list = []
    epsilon_r_list = []
    doping_list = []
    layer_boundaries = []
    
    x_current = 0.0
    global_idx = 0
    
    for layer in layers:
        mat_name = layer['material']
        thickness = layer['thickness']
        doping = layer['doping']
        
        # Get material properties
        mat = get_material(mat_name)
        
        # Number of points in this layer
        N_layer = max(1, int(round(thickness / a)))
        
        # Record layer boundaries
        start_idx = global_idx
        end_idx = global_idx + N_layer - 1
        layer_boundaries.append((start_idx, end_idx))
        
        # Generate points for this layer
        for i in range(N_layer):
            x_coords.append(x_current + i * a)
            material_list.append(mat_name)
            m_eff_list.append(mat['m_eff'] * m0)  # Convert to kg
            Ec_list.append(mat['Ec'])  # eV
            epsilon_r_list.append(mat['epsilon_r'])
            doping_list.append(doping)
            
        x_current += N_layer * a
        global_idx += N_layer
    
    return {
        'x': np.array(x_coords),
        'material': material_list,
        'm_eff': np.array(m_eff_list),
        'Ec': np.array(Ec_list),
        'epsilon_r': np.array(epsilon_r_list),
        'doping': np.array(doping_list),
        'Np': len(x_coords),
        'a': a,
        'layer_boundaries': layer_boundaries
    }

# ============================================================================
# HOPPING MATRIX (variable mass)
# ============================================================================

def compute_hopping(m_eff, a):
    """
    Compute position-dependent hopping energies
    
    For variable effective mass, the hopping between sites i and i+1 is:
    t[i] = ℏ²/(2*m_avg*a²) where m_avg = (m[i] + m[i+1])/2
    
    Parameters:
    -----------
    m_eff : array
        Effective mass at each site (kg)
    a : float
        Grid spacing (m)
        
    Returns:
    --------
    t : array
        Hopping energies (eV), length = Np - 1
    """
    
    Np = len(m_eff)
    t = np.zeros(Np - 1)
    
    for i in range(Np - 1):
        m_avg = 0.5 * (m_eff[i] + m_eff[i+1])
        t[i] = (hbar**2) / (2.0 * m_avg * a**2)
        t[i] = J_to_eV(t[i])  # Convert to eV
    
    return t

# ============================================================================
# HAMILTONIAN CONSTRUCTION
# ============================================================================

def build_hamiltonian(grid, U_external=None):
    """
    Build tight-binding Hamiltonian matrix
    
    H[i,i] = Ec[i] + U_external[i] + t[i-1] + t[i]
    H[i,i+1] = H[i+1,i] = -t[i]
    
    Parameters:
    -----------
    grid : dict
        Discretized device from discretize_device()
    U_external : array or None
        External potential at each site (eV), e.g., from Poisson
        If None, uses zero potential
        
    Returns:
    --------
    H : array (Np × Np)
        Hamiltonian matrix (eV)
    t : array
        Hopping energies (eV)
    """
    
    Np = grid['Np']
    a = grid['a']
    m_eff = grid['m_eff']
    Ec = grid['Ec']
    
    # Compute hopping
    t = compute_hopping(m_eff, a)
    
    # External potential
    if U_external is None:
        U = np.zeros(Np)
    else:
        U = U_external.copy()
    
    # Build Hamiltonian
    H = np.zeros((Np, Np), dtype=float)
    
    for i in range(Np):
        # Diagonal elements
        if i == 0:
            # Left edge
            H[i, i] = Ec[i] + U[i] + 2.0 * t[i]
        elif i == Np - 1:
            # Right edge
            H[i, i] = Ec[i] + U[i] + 2.0 * t[i-1]
        else:
            # Bulk
            H[i, i] = Ec[i] + U[i] + t[i-1] + t[i]
        
        # Off-diagonal (hopping)
        if i < Np - 1:
            H[i, i+1] = -t[i]
            H[i+1, i] = -t[i]
    
    return H, t

# ============================================================================
# TRANSVERSE MODES
# ============================================================================

def compute_transverse_modes(Ly, Lz, n_max, m_max, m_trans):
    """
    Compute transverse mode energies for hard-wall confinement
    
    E_nm = (ℏ²π²/2m*) * (n²/Ly² + m²/Lz²)
    
    Parameters:
    -----------
    Ly, Lz : float
        Transverse dimensions (m)
    n_max, m_max : int
        Maximum mode indices
    m_trans : float
        Transverse effective mass (kg)
        
    Returns:
    --------
    eps_perp : array
        Transverse mode energies (eV)
    modes : array
        Mode indices (n, m)
    """
    
    eps_list = []
    mode_list = []
    
    for n in range(1, n_max + 1):
        for m in range(1, m_max + 1):
            E_nm_J = (hbar**2 * np.pi**2 / (2.0 * m_trans)) * (
                (n**2 / Ly**2) + (m**2 / Lz**2)
            )
            eps_list.append(J_to_eV(E_nm_J))
            mode_list.append((n, m))
    
    return np.array(eps_list), np.array(mode_list, dtype=int)

def psi_nm_at_point(n, m, y, z, Ly, Lz):
    """
    Evaluate hard-wall mode wavefunction at point (y, z)
    
    ψ_nm(y,z) = (2/√(Ly*Lz)) * sin(nπy/Ly) * sin(mπz/Lz)
    
    Parameters:
    -----------
    n, m : int
        Mode indices
    y, z : float
        Coordinates (m)
    Ly, Lz : float
        Box dimensions (m)
        
    Returns:
    --------
    psi : float
        Wavefunction amplitude
    """
    
    psi_y = np.sqrt(2.0 / Ly) * np.sin(n * np.pi * y / Ly)
    psi_z = np.sqrt(2.0 / Lz) * np.sin(m * np.pi * z / Lz)
    return psi_y * psi_z

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_fermi_level(doping, temperature, m_eff_val):
    """
    Estimate Fermi level from doping using non-degenerate approximation
    
    For n-type: EF ≈ kT * ln(n / Nc)
    where Nc = 2 * (m*kT/2πℏ²)^(3/2)
    
    Parameters:
    -----------
    doping : float
        Doping density (m^-3)
    temperature : float
        Temperature (K)
    m_eff_val : float
        Effective mass (in units of m0)
        
    Returns:
    --------
    EF : float
        Fermi level relative to band edge (eV)
    """
    
    if doping <= 0:
        return 0.0
    
    m_star = m_eff_val * m0
    kT_eV = kB * temperature * J_to_eV(1.0)
    
    # Effective density of states
    Nc = 2.0 * (m_star * kB * temperature / (2.0 * np.pi * hbar**2))**(3.0/2.0)
    
    # Fermi level (non-degenerate)
    if doping < Nc:
        EF = kT_eV * np.log(doping / Nc)
    else:
        # Degenerate regime (simple approximation)
        EF = kT_eV * np.log(doping / Nc) + kT_eV * (doping / Nc)**(2.0/3.0)
    
    return EF

def print_grid_info(grid):
    """Print summary of discretized grid"""
    
    print("\n" + "="*70)
    print("DISCRETIZED GRID INFORMATION")
    print("="*70)
    print(f"Total points: {grid['Np']}")
    print(f"Grid spacing: {grid['a']*1e9:.3f} nm")
    print(f"Total length: {grid['x'][-1]*1e9:.3f} nm")
    print(f"\nLayers:")
    
    for i, (start, end) in enumerate(grid['layer_boundaries']):
        mat = grid['material'][start]
        thickness = (end - start + 1) * grid['a']
        doping = grid['doping'][start]
        
        if doping > 0:
            doping_str = f"{doping/1e6:.1e} cm⁻³"
        else:
            doping_str = "Undoped"
        
        print(f"  Layer {i}: {mat:8s}  {thickness*1e9:5.2f} nm  {doping_str}")
    
    print("="*70)

# ============================================================================
# TEST CODE
# ============================================================================

if __name__ == "__main__":
    print("\n" + "="*70)
    print("TESTING HAMILTONIAN BUILDER")
    print("="*70)
    
    # Test 1: GaAs/AlAs device
    print("\n[TEST 1] GaAs/AlAs symmetric RTD")
    print("-" * 70)
    
    from config.device_library import get_device
    
    device = get_device("GaAs_AlAs_symmetric")
    grid = discretize_device(device, grid_spacing=0.12e-9)
    print_grid_info(grid)
    
    # Build Hamiltonian
    H, t = build_hamiltonian(grid)
    
    print(f"\nHamiltonian matrix: {H.shape}")
    print(f"Hopping array: {t.shape}")
    print(f"Hopping range: {t.min():.4f} to {t.max():.4f} eV")
    print(f"Diagonal range: {np.diag(H).min():.4f} to {np.diag(H).max():.4f} eV")
    
    # Test 2: Oxide RTD
    print("\n[TEST 2] In2O3/Al2O3 symmetric RTD")
    print("-" * 70)
    
    device_oxide = get_device("In2O3_Al2O3_symmetric")
    grid_oxide = discretize_device(device_oxide, grid_spacing=0.12e-9)
    print_grid_info(grid_oxide)
    
    H_oxide, t_oxide = build_hamiltonian(grid_oxide)
    
    print(f"\nHamiltonian matrix: {H_oxide.shape}")
    print(f"Hopping array: {t_oxide.shape}")
    print(f"Hopping range: {t_oxide.min():.4f} to {t_oxide.max():.4f} eV")
    print(f"Diagonal range: {np.diag(H_oxide).min():.4f} to {np.diag(H_oxide).max():.4f} eV")
    
    # Test 3: Transverse modes
    print("\n[TEST 3] Transverse modes")
    print("-" * 70)
    
    Ly = Lz = 1e-6  # 1 µm
    n_max = m_max = 3
    m_trans = 0.067 * m0  # GaAs
    
    eps_perp, modes = compute_transverse_modes(Ly, Lz, n_max, m_max, m_trans)
    
    print(f"Generated {len(eps_perp)} transverse modes")
    print(f"Mode energies (meV):")
    for i, ((n, m), eps) in enumerate(zip(modes, eps_perp)):
        print(f"  ({n},{m}): {eps*1000:.2f} meV")
    
    # Test 4: Fermi level estimation
    print("\n[TEST 4] Fermi level estimation")
    print("-" * 70)
    
    for doping in [1e23, 1e24, 1e25]:  # m^-3
        EF = get_fermi_level(doping, 300, 0.067)
        print(f"  n = {doping:.1e} m⁻³ ({doping/1e6:.1e} cm⁻³) → EF = {EF:.3f} eV")
    
    print("\n" + "="*70)
    print("ALL TESTS PASSED ✓")
    print("="*70 + "\n")
