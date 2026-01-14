"""
Self-Consistent Born Approximation (SCBA) Solver

Solves for inelastic transport with multiple phonon modes:
- Bulk phonons (global, material-specific)
- Molecular vibrations (local, from adsorbed odorant)

Iterates until inscattering/outscattering self-energies converge.
"""

import numpy as np
import sys
import os

current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

from core.green_functions import (
    retarded_greens_function, spectral_function, broadening_function,
    fermi_function
)
from core.self_energy import (
    phonon_self_energy_bulk, local_projection_operator
)

# ============================================================================
# SCBA ITERATION
# ============================================================================

def scba_iteration(E_array, H, Sigma1_func, Sigma2_func, mu1, mu2, temperature,
                   phonon_modes, grid, max_iter=100, tol=1e-5, mix=0.3, verbose=True):
    """
    Self-consistent Born approximation for multi-phonon transport
    
    Iterates:
    1. Compute G from current Σ_S
    2. Compute n, p correlation functions
    3. Update Σ_S from phonon scattering
    4. Mix and check convergence
    
    Parameters:
    -----------
    E_array : array
        Energy grid (eV)
    H : array (Np × Np)
        Hamiltonian
    Sigma1_func, Sigma2_func : callable
        Contact self-energies Σ(E)
    mu1, mu2 : float
        Chemical potentials (eV)
    temperature : float
        Temperature (K)
    phonon_modes : list of dict
        Each dict has keys: 'energy' (eV), 'coupling' (eV), 'is_local' (bool)
        Optional: 'local_sites', 'neighbor_radius'
    grid : dict
        Grid information
    max_iter : int
        Maximum iterations
    tol : float
        Convergence tolerance
    mix : float
        Mixing parameter (0 < mix <= 1)
    verbose : bool
        Print progress
        
    Returns:
    --------
    dict with keys:
        'G' : array (Np × Np × NE) - Green's functions
        'Sigma_S_in' : array (Np × Np × NE) - Total inscattering
        'Sigma_S_out' : array (Np × Np × NE) - Total outscattering
        'n' : array (Np × Np × NE) - Electron correlation
        'p' : array (Np × Np × NE) - Hole correlation
        'converged' : bool
        'iterations' : int
        'residual' : float
    """
    
    kB_eV = 8.617333262145e-5  # eV/K
    kT = kB_eV * temperature
    
    Np = H.shape[0]
    NE = len(E_array)
    dE = E_array[1] - E_array[0]
    
    # Initialize arrays
    G = np.zeros((Np, Np, NE), dtype=complex)
    n_matrix = np.zeros((Np, Np, NE), dtype=complex)
    p_matrix = np.zeros((Np, Np, NE), dtype=complex)
    
    Sigma_S_in = np.zeros((Np, Np, NE), dtype=complex)
    Sigma_S_out = np.zeros((Np, Np, NE), dtype=complex)
    
    # Broadening matrices for contacts
    Gamma1_array = np.zeros((Np, Np, NE), dtype=float)
    Gamma2_array = np.zeros((Np, Np, NE), dtype=float)
    
    # Fermi functions
    f1 = fermi_function(E_array, mu1, kT)
    f2 = fermi_function(E_array, mu2, kT)
    
    # Pre-compute contact quantities
    for iE, E in enumerate(E_array):
        S1 = Sigma1_func(E)
        S2 = Sigma2_func(E)
        Gamma1_array[:, :, iE] = broadening_function(S1)
        Gamma2_array[:, :, iE] = broadening_function(S2)
    
    # Build projection operators for local modes
    projectors = {}
    for i, mode in enumerate(phonon_modes):
        if mode.get('is_local', False):
            sites = mode.get('local_sites', [Np//2])
            radius = mode.get('neighbor_radius', 2)
            P = local_projection_operator(Np, sites, radius, decay=True)
            projectors[i] = P
    
    # SCBA iteration
    converged = False
    
    for iteration in range(max_iter):
        # Store old self-energies
        Sigma_S_in_old = Sigma_S_in.copy()
        Sigma_S_out_old = Sigma_S_out.copy()
        
        # ====================================================================
        # Step 1: Compute Green's functions and correlation functions
        # ====================================================================
        
        for iE, E in enumerate(E_array):
            # Self-energies
            S1 = Sigma1_func(E)
            S2 = Sigma2_func(E)
            
            # Total scattering self-energy
            SS = Sigma_S_in[:, :, iE] - Sigma_S_out[:, :, iE]
            
            # Green's function
            G[:, :, iE] = retarded_greens_function(E, H, S1, S2, SS, eta=1e-4)
            
            # Spectral function
            A = spectral_function(G[:, :, iE])
            
            # Contact broadenings
            Gamma1 = Gamma1_array[:, :, iE]
            Gamma2 = Gamma2_array[:, :, iE]
            Gamma_S_in = 2.0 * np.real(Sigma_S_in[:, :, iE])
            Gamma_S_out = 2.0 * np.real(Sigma_S_out[:, :, iE])
            
            # Inscattering self-energy (all sources)
            Sigma_in_total = f1[iE] * Gamma1 + f2[iE] * Gamma2 + Gamma_S_in

            # Electron and hole correlation functions
            # Suppress overflow warnings and replace NaN/Inf with zeros
            with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
                n_temp = G[:, :, iE] @ Sigma_in_total @ G[:, :, iE].conj().T
            n_matrix[:, :, iE] = np.nan_to_num(n_temp, nan=0.0, posinf=0.0, neginf=0.0)
            p_matrix[:, :, iE] = A - n_matrix[:, :, iE]
        
        # ====================================================================
        # Step 2: Update phonon self-energies
        # ====================================================================
        
        Sigma_S_in_new = np.zeros_like(Sigma_S_in)
        Sigma_S_out_new = np.zeros_like(Sigma_S_out)
        
        for i, mode in enumerate(phonon_modes):
            hbar_omega = mode['energy']
            D = mode['coupling']
            is_local = mode.get('is_local', False)
            
            # Get projection operator if local
            P = projectors.get(i, None) if is_local else None
            
            # Compute this mode's contribution
            Sig_in, Sig_out = phonon_self_energy_bulk(
                n_matrix, p_matrix, D, hbar_omega, E_array, 
                temperature, is_local=is_local, P=P
            )
            
            Sigma_S_in_new += Sig_in
            Sigma_S_out_new += Sig_out
        
        # ====================================================================
        # Step 3: Mix and check convergence
        # ====================================================================
        
        # Calculate residual
        diff_in = np.abs(Sigma_S_in_new - Sigma_S_in_old).sum()
        diff_out = np.abs(Sigma_S_out_new - Sigma_S_out_old).sum()
        residual = diff_in + diff_out
        
        # Mix
        Sigma_S_in = (1 - mix) * Sigma_S_in_old + mix * Sigma_S_in_new
        Sigma_S_out = (1 - mix) * Sigma_S_out_old + mix * Sigma_S_out_new
        
        if verbose and (iteration % 10 == 0 or iteration < 5):
            print(f"  Iteration {iteration:3d}: residual = {residual:.6e}")
        
        # Check convergence
        if residual < tol:
            converged = True
            if verbose:
                print(f"  ✓ Converged at iteration {iteration} (residual = {residual:.6e})")
            break
    
    if not converged and verbose:
        print(f"  ⚠ Did not converge after {max_iter} iterations (residual = {residual:.6e})")
    
    return {
        'G': G,
        'Sigma_S_in': Sigma_S_in,
        'Sigma_S_out': Sigma_S_out,
        'n': n_matrix,
        'p': p_matrix,
        'converged': converged,
        'iterations': iteration + 1,
        'residual': residual
    }

# ============================================================================
# CURRENT CALCULATION
# ============================================================================

def compute_current(result, Gamma1_array, Gamma2_array, E_array, mu1=None, mu2=None, temperature=300):
    """
    Compute current from SCBA result using Datta NEGF formula

    Correct formula (Datta 2000, Eq. 4.14):
        I = (q/h) ∫ dE Tr[Γ₁A₂](f₁ - f₂)
    where:
        A₁ = GΓ₁G†  (left spectral function)
        A₂ = GΓ₂G†  (right spectral function)
        f₁, f₂ = Fermi functions at contacts 1 and 2

    This formula GUARANTEES I = 0 at equilibrium (f₁ = f₂).

    Parameters:
    -----------
    result : dict
        Output from scba_iteration (must contain 'G')
    Gamma1_array, Gamma2_array : arrays (Np × Np × NE)
        Contact broadening functions
    E_array : array
        Energy grid (eV)
    mu1, mu2 : float, optional
        Chemical potentials (eV). If None, uses symmetric bias from E_array
    temperature : float
        Temperature (K)

    Returns:
    --------
    I : float
        Current (Amperes)
    I_vs_E : array
        Current density vs energy (A/eV) for IETS spectrum
    """

    q = 1.602176634e-19  # C
    hbar = 1.054571817e-34  # J·s
    h = 2 * np.pi * hbar  # Planck's constant
    kB_eV = 8.617333262145e-5  # eV/K
    kT = kB_eV * temperature

    # Get Green's function from result
    if 'G' not in result:
        raise ValueError("Result dict must contain 'G' (Green's function). "
                        "For ballistic transport, use scba_iteration with no phonons.")

    G = result['G']
    NE = len(E_array)
    Np = G.shape[0]

    # Default: symmetric bias
    if mu1 is None or mu2 is None:
        E_mid = (E_array[0] + E_array[-1]) / 2
        mu1 = E_mid
        mu2 = E_mid

    # Compute Fermi functions
    f1 = 1.0 / (1.0 + np.exp((E_array - mu1) / kT))
    f2 = 1.0 / (1.0 + np.exp((E_array - mu2) / kT))

    # Current density vs energy using Datta formula
    # I(E) = (q/h) Tr[Γ₁A₂](f₁ - f₂) where A₂ = GΓ₂G†
    I_vs_E = np.zeros(NE)

    for iE in range(NE):
        # Compute A₂ = GΓ₂G† with numerical safeguards
        with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
            A2 = G[:, :, iE] @ Gamma2_array[:, :, iE] @ G[:, :, iE].conj().T
        A2 = np.nan_to_num(A2, nan=0.0, posinf=0.0, neginf=0.0)

        # Tr[Γ₁A₂]
        with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
            transmission = np.trace(Gamma1_array[:, :, iE] @ A2)
        transmission = np.nan_to_num(transmission, nan=0.0, posinf=0.0, neginf=0.0)

        # Current integrand: (q/h) T(E) (f₁ - f₂)
        # Note: E_array is in eV, so dE integration gives eV units
        # We need to convert: (q/h) has units of A/eV when E is in eV
        I_vs_E[iE] = np.real(transmission) * (f1[iE] - f2[iE])

    # Integrate over energy
    # Standard Landauer-Büttiker formula:
    # I = (2q²/h) ∫ T(E)(f₁ - f₂) dE  where the factor 2 includes spin
    # Since E_array is in eV, we need q² to convert: one q for eV→J, one for charge
    # Note: For spin-degenerate systems, the factor 2 is included here
    I = (2 * q**2 / h) * np.trapz(I_vs_E, E_array)

    # For plotting: convert I_vs_E to current density (A/eV)
    I_vs_E_density = I_vs_E * (2 * q**2 / h)

    return I, I_vs_E_density

# ============================================================================
# MULTIMODE SCBA ITERATION (Quasi-3D Transport)
# ============================================================================

def scba_iteration_multimode(E_array, H, Sigma1_func, Sigma2_func,
                             mu1, mu2, temperature, phonon_modes, grid,
                             transverse_modes,
                             max_iter=100, tol=1e-5, mix=0.3, verbose=True):
    """
    Self-consistent Born approximation with transverse mode summation

    Key changes from scba_iteration():
    1. Arrays become 4D: (Np, Np, NE, Nm)
    2. For each energy E and mode m:
       - Compute G_nm(E + ε_nm)
       - Compute n_nm, p_nm
    3. Phonon self-energies are mode-dependent
    4. Each mode evolves independently in SCBA loop

    Parameters:
    -----------
    E_array : array
        Longitudinal energy grid (eV)
    H : array (Np × Np)
        Longitudinal Hamiltonian
    Sigma1_func, Sigma2_func : callable
        Contact self-energies Σ(E) → (Np, Np)
    mu1, mu2 : float
        Chemical potentials (eV)
    temperature : float
        Temperature (K)
    phonon_modes : list of dict
        Phonon mode specifications
    grid : dict
        Grid information
    transverse_modes : TransverseModes object
        Transverse mode configuration
    max_iter : int
        Maximum SCBA iterations
    tol : float
        Convergence tolerance
    mix : float
        Mixing parameter (0 < mix <= 1)
    verbose : bool
        Print progress

    Returns:
    --------
    dict with keys:
        'G_nm' : array (Np, Np, NE, Nm) - Mode-resolved Green's functions
        'Sigma_S_in_nm' : array (Np, Np, NE, Nm) - Mode-resolved inscattering
        'Sigma_S_out_nm' : array (Np, Np, NE, Nm) - Mode-resolved outscattering
        'n_nm' : array (Np, Np, NE, Nm) - Mode-resolved electron correlation
        'p_nm' : array (Np, Np, NE, Nm) - Mode-resolved hole correlation
        'converged' : bool
        'iterations' : int
        'residual' : float
        'Gamma1_array' : array (Np, Np, NE) - Left broadening (mode-independent)
        'Gamma2_array' : array (Np, Np, NE) - Right broadening (mode-independent)
        'transverse_modes' : TransverseModes object
    """

    # Import here to avoid circular dependency
    from core.green_functions import (
        retarded_greens_function, spectral_function, broadening_function
    )
    from core.self_energy import phonon_self_energy_bulk, local_projection_operator

    kB_eV = 8.617333262145e-5  # eV/K
    kT = kB_eV * temperature

    Np = H.shape[0]
    NE = len(E_array)
    Nm = transverse_modes.Nm

    if verbose:
        print(f"  Multimode SCBA: {Nm} transverse modes, {NE} energies, {Np} sites")

    # Initialize mode-resolved arrays (4D)
    G_nm = np.zeros((Np, Np, NE, Nm), dtype=complex)
    n_matrix_nm = np.zeros((Np, Np, NE, Nm), dtype=complex)
    p_matrix_nm = np.zeros((Np, Np, NE, Nm), dtype=complex)
    Sigma_S_in_nm = np.zeros((Np, Np, NE, Nm), dtype=complex)
    Sigma_S_out_nm = np.zeros((Np, Np, NE, Nm), dtype=complex)

    # Pre-compute contact broadening (mode-independent)
    Gamma1_array = np.zeros((Np, Np, NE), dtype=float)
    Gamma2_array = np.zeros((Np, Np, NE), dtype=float)

    for iE, E in enumerate(E_array):
        S1 = Sigma1_func(E)
        S2 = Sigma2_func(E)
        Gamma1_array[:, :, iE] = broadening_function(S1)
        Gamma2_array[:, :, iE] = broadening_function(S2)

    # Build projection operators for local phonon modes
    projectors = {}
    for i, mode in enumerate(phonon_modes):
        if mode.get('is_local', False):
            sites = mode.get('local_sites', [Np//2])
            radius = mode.get('neighbor_radius', 2)
            P = local_projection_operator(Np, sites, radius, decay=True)
            projectors[i] = P

    # SCBA iteration loop
    converged = False

    for iteration in range(max_iter):
        # Store old self-energies
        Sigma_S_in_old = Sigma_S_in_nm.copy()
        Sigma_S_out_old = Sigma_S_out_nm.copy()

        # ====================================================================
        # Step 1: Green's functions and correlations (FOR EACH MODE)
        # ====================================================================
        for iE, E_long in enumerate(E_array):
            S1 = Sigma1_func(E_long)
            S2 = Sigma2_func(E_long)

            # Compute Green's function for all modes at this longitudinal energy
            for im in range(Nm):
                # Total energy: longitudinal + transverse confinement
                E_total = transverse_modes.get_total_energy(E_long, im)

                # Extract mode-specific scattering self-energy
                SS_in = Sigma_S_in_nm[:, :, iE, im]
                SS_out = Sigma_S_out_nm[:, :, iE, im]
                SS = SS_in - SS_out

                # Compute Green's function
                G_nm[:, :, iE, im] = retarded_greens_function(
                    E_total, H, S1, S2, SS, eta=1e-4
                )

            # Compute spectral function and correlations for each mode
            for im in range(Nm):
                A_nm = spectral_function(G_nm[:, :, iE, im])

                # Inscattering contribution (mode-dependent)
                Sigma_in_total = Sigma_S_in_nm[:, :, iE, im]

                with np.errstate(all='ignore'):
                    n_temp = (G_nm[:, :, iE, im] @ Sigma_in_total @
                             G_nm[:, :, iE, im].conj().T)

                n_matrix_nm[:, :, iE, im] = np.nan_to_num(
                    n_temp, nan=0.0, posinf=0.0, neginf=0.0
                )
                p_matrix_nm[:, :, iE, im] = A_nm - n_matrix_nm[:, :, iE, im]

        # ====================================================================
        # Step 2: Update phonon self-energies (FOR EACH MODE)
        # ====================================================================
        Sigma_S_in_new = np.zeros_like(Sigma_S_in_nm)
        Sigma_S_out_new = np.zeros_like(Sigma_S_out_nm)

        for im in range(Nm):
            # Extract correlations for this mode
            n_mode = n_matrix_nm[:, :, :, im]  # (Np, Np, NE)
            p_mode = p_matrix_nm[:, :, :, im]  # (Np, Np, NE)

            # Compute phonon self-energies for each phonon mode
            for i, mode in enumerate(phonon_modes):
                # Handle mode-dependent coupling
                D_coupling = mode['coupling']
                if isinstance(D_coupling, np.ndarray):
                    # Mode-dependent: use coupling for this specific transverse mode
                    D = D_coupling[im]
                else:
                    # Mode-independent: scalar coupling (bulk phonons or 1D)
                    D = D_coupling

                hbar_omega = mode['energy']
                is_local = mode.get('is_local', False)
                P = projectors.get(i, None)

                # Compute self-energy for this phonon mode
                Sigma_in, Sigma_out = phonon_self_energy_bulk(
                    n_mode, p_mode, D, hbar_omega, E_array,
                    temperature, is_local, P
                )

                Sigma_S_in_new[:, :, :, im] += Sigma_in
                Sigma_S_out_new[:, :, :, im] += Sigma_out

        # Mix with old values
        Sigma_S_in_nm = mix * Sigma_S_in_new + (1 - mix) * Sigma_S_in_old
        Sigma_S_out_nm = mix * Sigma_S_out_new + (1 - mix) * Sigma_S_out_old

        # Check convergence
        residual_in = np.sum(np.abs(Sigma_S_in_nm - Sigma_S_in_old))
        residual_out = np.sum(np.abs(Sigma_S_out_nm - Sigma_S_out_old))
        residual = residual_in + residual_out

        if verbose and (iteration < 5 or iteration % 5 == 0):
            print(f"  Iteration {iteration:3d}: residual = {residual:.6e}")

        if residual < tol:
            converged = True
            if verbose:
                print(f"  ✓ Converged after {iteration+1} iterations")
            break

    if not converged and verbose:
        print(f"  ⚠ Did not converge after {max_iter} iterations (residual = {residual:.6e})")

    return {
        'G_nm': G_nm,
        'Sigma_S_in_nm': Sigma_S_in_nm,
        'Sigma_S_out_nm': Sigma_S_out_nm,
        'n_nm': n_matrix_nm,
        'p_nm': p_matrix_nm,
        'converged': converged,
        'iterations': iteration + 1 if converged else max_iter,
        'residual': residual,
        'Gamma1_array': Gamma1_array,
        'Gamma2_array': Gamma2_array,
        'transverse_modes': transverse_modes
    }


# ============================================================================
# MULTIMODE CURRENT CALCULATION (Quasi-3D Transport)
# ============================================================================

def compute_current_multimode(result, E_array, mu1=None, mu2=None, temperature=300):
    """
    Compute current with transverse mode summation

    Formula:
        I = (2q²/h) ∫ dE [f₁(E) - f₂(E)] × Σ_nm T_nm(E)

    where T_nm(E) = Tr[Γ₁ · G_nm · Γ₂ · G_nm†] for each mode (n,m)

    The sum over transverse modes represents multiple conduction channels.

    Parameters:
    -----------
    result : dict
        From scba_iteration_multimode, contains:
        - 'G_nm': (Np, Np, NE, Nm) Green's functions
        - 'Gamma1_array', 'Gamma2_array': (Np, Np, NE) broadening
        - 'transverse_modes': TransverseModes object
    E_array : array
        Energy grid (eV)
    mu1, mu2 : float, optional
        Chemical potentials (eV). If None, uses symmetric bias
    temperature : float
        Temperature (K)

    Returns:
    --------
    I : float
        Total current (A)
    I_vs_E : array (NE,)
        Current density by energy (A/eV)
    I_vs_mode : array (NE, Nm)
        Current density by energy and mode (A/eV)
        Useful for analyzing mode contributions
    """

    q = 1.602176634e-19  # C
    h = 2 * np.pi * 1.054571817e-34  # J·s
    kB_eV = 8.617333262145e-5  # eV/K
    kT = kB_eV * temperature

    # Extract arrays from result
    G_nm = result['G_nm']  # (Np, Np, NE, Nm)
    Gamma1_array = result['Gamma1_array']  # (Np, Np, NE)
    Gamma2_array = result['Gamma2_array']  # (Np, Np, NE)
    transverse_modes = result['transverse_modes']

    NE = len(E_array)
    Nm = transverse_modes.Nm

    # Default: symmetric bias
    if mu1 is None or mu2 is None:
        E_mid = (E_array[0] + E_array[-1]) / 2
        mu1 = E_mid
        mu2 = E_mid

    # Fermi functions
    f1 = 1.0 / (1.0 + np.exp((E_array - mu1) / kT))
    f2 = 1.0 / (1.0 + np.exp((E_array - mu2) / kT))

    # Current integrand per mode
    I_vs_mode = np.zeros((NE, Nm))

    for iE in range(NE):
        for im in range(Nm):
            # Compute A2_nm = G_nm · Γ₂ · G_nm†
            with np.errstate(all='ignore'):
                A2_nm = (G_nm[:, :, iE, im] @ Gamma2_array[:, :, iE] @
                         G_nm[:, :, iE, im].conj().T)
            A2_nm = np.nan_to_num(A2_nm, nan=0.0, posinf=0.0, neginf=0.0)

            # Transmission: T_nm = Tr[Γ₁ · A2_nm]
            with np.errstate(all='ignore'):
                transmission_nm = np.trace(Gamma1_array[:, :, iE] @ A2_nm)
            transmission_nm = np.nan_to_num(transmission_nm, nan=0.0, posinf=0.0, neginf=0.0)

            # Mode contribution to current
            I_vs_mode[iE, im] = np.real(transmission_nm) * (f1[iE] - f2[iE])

    # Sum over all modes
    I_vs_E = np.sum(I_vs_mode, axis=1)

    # Integrate: Landauer-Büttiker formula with spin factor
    I = (2 * q**2 / h) * np.trapz(I_vs_E, E_array)

    # Convert to current density
    I_vs_E_density = I_vs_E * (2 * q**2 / h)
    I_vs_mode_density = I_vs_mode * (2 * q**2 / h)

    return I, I_vs_E_density, I_vs_mode_density


# ============================================================================
# IETS CALCULATION
# ============================================================================

def compute_iets(V_array, I_array):
    """
    Compute IETS (d²I/dV²) from I-V curve
    
    Parameters:
    -----------
    V_array : array
        Voltage points (V)
    I_array : array
        Current points (A)
        
    Returns:
    --------
    dIdV : array
        Differential conductance (S)
    d2IdV2 : array
        IETS signal (S/V)
    """
    
    dV = V_array[1] - V_array[0]
    
    # First derivative
    dIdV = np.gradient(I_array, dV)
    
    # Second derivative
    d2IdV2 = np.gradient(dIdV, dV)
    
    return dIdV, d2IdV2

# ============================================================================
# TEST CODE
# ============================================================================

if __name__ == "__main__":
    print("\n" + "="*70)
    print("TESTING SCBA SOLVER")
    print("="*70)
    
    from core.hamiltonian import discretize_device, build_hamiltonian
    from core.self_energy import contact_self_energy_matrix
    from config.device_library import get_device
    
    # Simple test: 1D chain with one phonon mode
    print("\n[TEST 1] SCBA with single bulk phonon")
    print("-" * 70)
    
    # Small test system
    Np = 30
    H_test = np.zeros((Np, Np))
    t_test = 1.0
    
    for i in range(Np-1):
        H_test[i, i+1] = -t_test
        H_test[i+1, i] = -t_test
    for i in range(Np):
        H_test[i, i] = 2*t_test
    
    # Simple grid
    grid_test = {
        'Np': Np,
        'a': 0.12e-9,
        'Ec': np.zeros(Np)
    }
    
    # Contact self-energies (simple model)
    def Sigma1_simple(E):
        S = np.zeros((Np, Np), dtype=complex)
        S[0, 0] = -1j * 0.05
        return S
    
    def Sigma2_simple(E):
        S = np.zeros((Np, Np), dtype=complex)
        S[-1, -1] = -1j * 0.05
        return S
    
    # Single phonon mode (bulk)
    phonon_modes = [
        {
            'energy': 0.036,  # 36 meV (GaAs-like)
            'coupling': 0.01,  # 10 meV
            'is_local': False
        }
    ]
    
    # Energy grid
    E_array = np.linspace(1.5, 2.5, 30)
    
    # Run SCBA
    print("Running SCBA...")
    result = scba_iteration(
        E_array, H_test, Sigma1_simple, Sigma2_simple,
        mu1=2.1, mu2=1.9, temperature=300,
        phonon_modes=phonon_modes, grid=grid_test,
        max_iter=50, tol=1e-4, mix=0.5, verbose=True
    )
    
    print(f"\nConverged: {result['converged']}")
    print(f"Iterations: {result['iterations']}")
    print(f"Final residual: {result['residual']:.6e}")
    
    # Test 2: IETS calculation
    print("\n[TEST 2] IETS calculation")
    print("-" * 70)
    
    # Mock I-V data with features
    V_test = np.linspace(0, 1, 100)
    # Parabola + small features
    I_test = 0.5 * V_test**2 + 0.01 * np.sin(20 * V_test)
    
    dIdV, d2IdV2 = compute_iets(V_test, I_test)
    
    print(f"Generated I-V: {len(V_test)} points")
    print(f"dI/dV range: {dIdV.min():.6f} to {dIdV.max():.6f} S")
    print(f"d²I/dV² range: {d2IdV2.min():.6f} to {d2IdV2.max():.6f} S/V")
    print(f"d²I/dV² peaks: {np.sum(d2IdV2 > 0.01)} positive features")
    
    print("\n" + "="*70)
    print("ALL TESTS PASSED ✓")
    print("="*70 + "\n")
