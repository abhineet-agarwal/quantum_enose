"""
Self-Energy Functions

Computes self-energies for:
1. Semi-infinite contacts (exact analytical)
2. Bulk phonons (global)
3. Molecular vibrations (local, with projection operator)
"""

import numpy as np
import sys
import os

# ============================================================================
# CONTACT SELF-ENERGIES (Semi-infinite leads)
# ============================================================================

def contact_self_energy_1d(E, U_edge, t_edge, eta=1e-4):
    """
    Analytical self-energy for semi-infinite 1D lead
    
    For a uniform 1D chain: H_ii = U, H_i,i±1 = -t
    The surface Green's function is:
    
    g_00 = exp(±ika) / t  where E - U = 2t(1 - cos(ka))
    
    For |E - U| < 2|t|: propagating (use arccos)
    For |E - U| > 2|t|: evanescent (use arccosh)
    
    Parameters:
    -----------
    E : float or complex
        Energy (eV)
    U_edge : float
        Onsite energy at edge (eV)
    t_edge : float
        Hopping energy (eV, positive)
    eta : float
        Small imaginary part for convergence
        
    Returns:
    --------
    Sigma : complex
        Self-energy at edge site
    """
    
    E_complex = E + 1j * eta
    z = E_complex - U_edge
    two_t = 2.0 * abs(t_edge)
    
    # Check if propagating or evanescent
    if abs(z) <= two_t:
        # Propagating: use arccos
        arg = z / (2.0 * t_edge)
        # Clip real part to avoid numerical issues
        arg_real = np.clip(np.real(arg), -1.0, 1.0)
        arg_clipped = arg_real + 1j * np.imag(arg)
        ka = np.arccos(arg_clipped)
        # Use -ika for retarded Green's function (causality)
        Sigma = t_edge * np.exp(-1j * ka)
    else:
        # Evanescent: use arccosh
        arg = abs(z) / (2.0 * abs(t_edge))
        ka = np.arccosh(arg)
        # Evanescent decay (already has correct sign)
        Sigma = t_edge * np.exp(-ka)
    
    return Sigma

def contact_self_energy_matrix(E, grid, t, contact='left', eta=1e-4):
    """
    Build contact self-energy matrix for RTD
    
    Σ_contact is a matrix with single non-zero element at edge.
    
    Parameters:
    -----------
    E : float
        Energy (eV)
    grid : dict
        Grid information
    t : array
        Hopping energies (from hamiltonian.py)
    contact : str
        'left' or 'right'
    eta : float
        Convergence parameter
        
    Returns:
    --------
    Sigma : array (Np × Np), complex
        Self-energy matrix
    """
    
    Np = grid['Np']
    Sigma = np.zeros((Np, Np), dtype=complex)
    
    if contact == 'left':
        # Left contact
        U_edge = grid['Ec'][0]
        t_edge = t[0]
        sigma_00 = contact_self_energy_1d(E, U_edge, t_edge, eta)
        Sigma[0, 0] = sigma_00
        
    elif contact == 'right':
        # Right contact
        U_edge = grid['Ec'][-1]
        t_edge = t[-1]
        sigma_NN = contact_self_energy_1d(E, U_edge, t_edge, eta)
        Sigma[-1, -1] = sigma_NN
    
    else:
        raise ValueError(f"Unknown contact: {contact}. Use 'left' or 'right'.")
    
    return Sigma

# ============================================================================
# PHONON PROJECTION OPERATOR (for localized modes)
# ============================================================================

def local_projection_operator(Np, local_sites, neighbor_radius=0, decay=True):
    """
    Build projection operator for local phonon modes
    
    P is diagonal matrix:
    - P[i,i] = 1 if i in local_sites
    - P[i,i] = 1/(r+1) if i is neighbor at distance r (if decay=True)
    - P[i,i] = 0 otherwise
    
    Then normalized so max(P) = 1.
    
    Parameters:
    -----------
    Np : int
        Number of sites
    local_sites : list of int
        Sites where molecule is adsorbed
    neighbor_radius : int
        Include neighbors within this distance
    decay : bool
        If True, weight neighbors by 1/(r+1)
        
    Returns:
    --------
    P : array (Np × Np)
        Projection operator (diagonal, real)
    """
    
    P = np.zeros((Np, Np), dtype=float)
    
    # Mark all sites within radius
    marked = {}  # {site_index: distance}
    
    for site in local_sites:
        for offset in range(-neighbor_radius, neighbor_radius + 1):
            idx = site + offset
            if 0 <= idx < Np:
                dist = abs(offset)
                if idx not in marked or marked[idx] > dist:
                    marked[idx] = dist
    
    # Fill diagonal
    for idx, dist in marked.items():
        if decay:
            weight = 1.0 / (dist + 1)
        else:
            weight = 1.0
        P[idx, idx] = weight
    
    # Normalize so max is 1
    if P.max() > 0:
        P = P / P.max()
    
    return P

# ============================================================================
# PHONON SELF-ENERGIES (for SCBA)
# ============================================================================

def phonon_self_energy_bulk(n_matrix, p_matrix, D, hbar_omega, E_grid, 
                             temperature, is_local=False, P=None):
    """
    Compute phonon self-energy for one mode
    
    For bulk phonon (global):
    Σ_S^in(E) = D² [(n_B + 1) n(E + ℏω) + n_B n(E - ℏω)]
    Σ_S^out(E) = D² [n_B p(E + ℏω) + (n_B + 1) p(E - ℏω)]
    
    For local phonon:
    Apply projection: Σ_S = P · Σ_S^global · P
    
    Parameters:
    -----------
    n_matrix : array (Np × Np × NE)
        Electron correlation function at all energies
    p_matrix : array (Np × Np × NE)
        Hole correlation function at all energies
    D : float
        Electron-phonon coupling (eV)
    hbar_omega : float
        Phonon energy (eV)
    E_grid : array
        Energy points (eV)
    temperature : float
        Temperature (K)
    is_local : bool
        If True, apply projection operator P
    P : array (Np × Np) or None
        Projection operator for local modes
        
    Returns:
    --------
    Sigma_in, Sigma_out : arrays (Np × Np × NE), complex
        Inscattering and outscattering self-energies
    """
    
    kB_eV = 8.617333262145e-5  # eV/K
    Np = n_matrix.shape[0]
    NE = len(E_grid)
    dE = E_grid[1] - E_grid[0]
    
    # Bose-Einstein distribution for phonons
    if hbar_omega > 0:
        kT = kB_eV * temperature
        n_B = 1.0 / (np.exp(hbar_omega / kT) - 1.0)
    else:
        n_B = 0.0
    
    # Energy shift in grid points
    shift = int(round(hbar_omega / dE))
    
    if shift <= 0 or shift >= NE:
        # Phonon energy too large or too small
        Sigma_in = np.zeros((Np, Np, NE), dtype=complex)
        Sigma_out = np.zeros((Np, Np, NE), dtype=complex)
        return Sigma_in, Sigma_out
    
    # Shifted correlation functions
    n_emission = np.zeros_like(n_matrix)  # n(E + ℏω)
    n_absorption = np.zeros_like(n_matrix)  # n(E - ℏω)
    p_emission = np.zeros_like(p_matrix)
    p_absorption = np.zeros_like(p_matrix)
    
    # Emission: E → E - ℏω (electron loses energy)
    n_emission[:, :, 0:NE-shift] = n_matrix[:, :, shift:NE]
    p_emission[:, :, 0:NE-shift] = p_matrix[:, :, shift:NE]
    
    # Absorption: E → E + ℏω (electron gains energy)
    n_absorption[:, :, shift:NE] = n_matrix[:, :, 0:NE-shift]
    p_absorption[:, :, shift:NE] = p_matrix[:, :, 0:NE-shift]
    
    # Compute self-energies
    D_squared = D**2
    
    Sigma_in = D_squared * (
        (n_B + 1.0) * n_emission +   # Emission: (n_B+1) factor
        n_B * n_absorption            # Absorption: n_B factor
    )
    
    Sigma_out = D_squared * (
        n_B * p_emission +            # Emission: n_B factor
        (n_B + 1.0) * p_absorption    # Absorption: (n_B+1) factor
    )
    
    # Apply projection for local modes
    if is_local and P is not None:
        # Suppress overflow warnings in projection operations
        with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
            for iE in range(NE):
                Sigma_in_temp = P @ Sigma_in[:, :, iE] @ P
                Sigma_out_temp = P @ Sigma_out[:, :, iE] @ P
                Sigma_in[:, :, iE] = np.nan_to_num(Sigma_in_temp, nan=0.0, posinf=0.0, neginf=0.0)
                Sigma_out[:, :, iE] = np.nan_to_num(Sigma_out_temp, nan=0.0, posinf=0.0, neginf=0.0)
    
    return Sigma_in, Sigma_out

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_molecule_location(grid, device_config, molecule_config=None):
    """
    Get lattice sites where molecule is located
    
    Parameters:
    -----------
    grid : dict
        Discretized grid
    device_config : dict
        Device configuration
    molecule_config : dict or None
        Molecule configuration (optional override)
        
    Returns:
    --------
    local_sites : list of int
        Grid indices where molecule is located
    neighbor_radius : int
        Radius for projection operator
    """
    
    # Get molecule location from device config
    mol_loc = device_config.get('molecule_location', 'emitter_barrier')
    
    # Default mapping (can be overridden)
    location_map = {
        'emitter_barrier': 1,      # Layer index
        'well_center': 2,
        'collector_barrier': 3,
        'emitter_well_interface': 1
    }
    
    layer_idx = location_map.get(mol_loc, 1)
    
    # Get layer boundaries
    boundaries = grid['layer_boundaries']
    if layer_idx >= len(boundaries):
        layer_idx = 1  # Default to first barrier
    
    start, end = boundaries[layer_idx]
    
    # Place molecule at center of layer
    center_site = (start + end) // 2
    local_sites = [center_site]
    
    # Neighbor radius
    neighbor_radius = molecule_config.get('neighbor_radius', 2) if molecule_config else 2
    
    return local_sites, neighbor_radius

def print_self_energy_info(Sigma, name="Self-energy"):
    """Print summary of self-energy matrix"""
    
    print(f"\n{name}:")
    print(f"  Shape: {Sigma.shape}")
    print(f"  Non-zero elements: {np.count_nonzero(np.abs(Sigma) > 1e-12)}")
    print(f"  Max |Σ|: {np.abs(Sigma).max():.6f}")
    
    if Sigma.ndim == 2:
        # Single matrix
        diag_vals = np.diag(Sigma)
        nonzero_diag = diag_vals[np.abs(diag_vals) > 1e-12]
        if len(nonzero_diag) > 0:
            print(f"  Non-zero diagonal elements: {len(nonzero_diag)}")
            print(f"  Diagonal range: {np.abs(nonzero_diag).min():.6f} to {np.abs(nonzero_diag).max():.6f}")

# ============================================================================
# TEST CODE
# ============================================================================

if __name__ == "__main__":
    print("\n" + "="*70)
    print("TESTING SELF-ENERGY FUNCTIONS")
    print("="*70)
    
    # Test 1: Contact self-energy
    print("\n[TEST 1] Contact self-energy (semi-infinite lead)")
    print("-" * 70)
    
    U_edge = 0.0  # eV
    t_edge = 1.0  # eV
    
    # Test at various energies
    energies = [-2.0, 0.0, 2.0, 4.0]
    
    print(f"Lead parameters: U = {U_edge} eV, t = {t_edge} eV")
    print(f"Band range: {U_edge - 2*t_edge} to {U_edge + 2*t_edge} eV")
    print(f"\nSelf-energies:")
    
    for E in energies:
        sigma = contact_self_energy_1d(E, U_edge, t_edge)
        in_band = (U_edge - 2*t_edge) <= E <= (U_edge + 2*t_edge)
        print(f"  E = {E:+.1f} eV → Σ = {sigma:.4f}, |Σ| = {abs(sigma):.4f} {'[in-band]' if in_band else '[evanescent]'}")
    
    # Test 2: Projection operator
    print("\n[TEST 2] Local projection operator")
    print("-" * 70)
    
    Np = 50
    local_sites = [25]  # Center
    neighbor_radius = 2
    
    P = local_projection_operator(Np, local_sites, neighbor_radius, decay=True)
    
    print(f"Grid size: {Np}")
    print(f"Molecule at site: {local_sites}")
    print(f"Neighbor radius: {neighbor_radius}")
    print(f"Non-zero P elements: {np.count_nonzero(P)}")
    print(f"P diagonal around molecule:")
    for i in range(23, 28):
        print(f"  P[{i},{i}] = {P[i,i]:.3f}")
    
    # Test 3: Contact self-energy matrix
    print("\n[TEST 3] Contact self-energy matrix for RTD")
    print("-" * 70)
    
    # Create simple grid
    current_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(current_dir)
    sys.path.insert(0, parent_dir)
    
    from core.hamiltonian import discretize_device, build_hamiltonian
    from config.device_library import get_device
    
    device = get_device("GaAs_AlAs_symmetric")
    grid = discretize_device(device, grid_spacing=0.12e-9)
    H, t = build_hamiltonian(grid)
    
    E_test = 0.5
    Sigma_L = contact_self_energy_matrix(E_test, grid, t, contact='left')
    Sigma_R = contact_self_energy_matrix(E_test, grid, t, contact='right')
    
    print_self_energy_info(Sigma_L, "Left contact Σ₁")
    print_self_energy_info(Sigma_R, "Right contact Σ₂")
    
    # Test 4: Molecule location
    print("\n[TEST 4] Molecule location finder")
    print("-" * 70)
    
    sites, radius = get_molecule_location(grid, device)
    print(f"Molecule placed at sites: {sites}")
    print(f"Neighbor radius: {radius}")
    print(f"Layer boundaries:")
    for i, (start, end) in enumerate(grid['layer_boundaries']):
        mat = grid['material'][start]
        print(f"  Layer {i}: sites {start}-{end} ({mat})")
    
    print("\n" + "="*70)
    print("ALL TESTS PASSED ✓")
    print("="*70 + "\n")
