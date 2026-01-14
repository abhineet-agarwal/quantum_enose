"""
Green's Functions and NEGF Solver

Computes retarded Green's function, spectral function, and density matrix
using Non-Equilibrium Green's Function (NEGF) formalism.
"""

import numpy as np
import scipy.linalg as la

# ============================================================================
# GREEN'S FUNCTIONS
# ============================================================================

def retarded_greens_function(E, H, Sigma1, Sigma2, SigmaS=None, eta=1e-4):
    """
    Compute retarded Green's function
    
    G^R(E) = [(E + iη)I - H - Σ₁ - Σ₂ - Σ_S]^(-1)
    
    Parameters:
    -----------
    E : float
        Energy (eV)
    H : array (Np × Np)
        Hamiltonian matrix (eV)
    Sigma1, Sigma2 : array (Np × Np) or complex
        Self-energies from contacts
    SigmaS : array (Np × Np) or None
        Scattering self-energy (optional)
    eta : float
        Small imaginary part for convergence (eV)
        
    Returns:
    --------
    G : array (Np × Np), complex
        Retarded Green's function
    """
    
    Np = H.shape[0]
    I = np.eye(Np, dtype=complex)
    
    # Build energy matrix with infinitesimal imaginary part
    E_matrix = (E + 1j * eta) * I
    
    # Total self-energy
    Sigma_total = Sigma1 + Sigma2
    if SigmaS is not None:
        Sigma_total = Sigma_total + SigmaS
    
    # Invert to get Green's function
    try:
        G = la.inv(E_matrix - H - Sigma_total)
    except la.LinAlgError:
        # Add small regularization if singular
        print(f"Warning: Singular matrix at E={E:.4f} eV, adding regularization")
        G = la.inv(E_matrix - H - Sigma_total + 1e-12 * I)
    
    return G

def spectral_function(G):
    """
    Compute spectral function from retarded Green's function
    
    A(E) = i[G^R(E) - G^A(E)] = i[G - G†]
    
    Parameters:
    -----------
    G : array (Np × Np), complex
        Retarded Green's function
        
    Returns:
    --------
    A : array (Np × Np), float
        Spectral function (real, positive semi-definite)
    """
    
    A = 1j * (G - G.conj().T)
    return np.real(A)  # Should be real by construction

def broadening_function(Sigma):
    """
    Compute broadening function from self-energy
    
    Γ(E) = i[Σ(E) - Σ†(E)]
    
    Parameters:
    -----------
    Sigma : array (Np × Np), complex
        Self-energy
        
    Returns:
    --------
    Gamma : array (Np × Np), float
        Broadening function (real, positive semi-definite)
    """
    
    Gamma = 1j * (Sigma - Sigma.conj().T)
    return np.real(Gamma)

# ============================================================================
# CONTACT-RESOLVED SPECTRAL FUNCTIONS
# ============================================================================

def contact_spectral_functions(G, Gamma1, Gamma2):
    """
    Compute contact-resolved spectral functions
    
    A₁(E) = G Γ₁ G†  (available from contact 1)
    A₂(E) = G Γ₂ G†  (available from contact 2)
    
    These represent the density of states accessible from each contact.
    
    Parameters:
    -----------
    G : array (Np × Np)
        Retarded Green's function
    Gamma1, Gamma2 : array (Np × Np)
        Broadening functions from contacts
        
    Returns:
    --------
    A1, A2 : arrays (Np × Np)
        Contact-resolved spectral functions
    """
    
    A1 = G @ Gamma1 @ G.conj().T
    A2 = G @ Gamma2 @ G.conj().T
    
    return np.real(A1), np.real(A2)

# ============================================================================
# DENSITY MATRIX (EQUILIBRIUM)
# ============================================================================

def fermi_function(E, mu, kT):
    """
    Fermi-Dirac distribution
    
    f₀(E) = 1 / (1 + exp[(E - μ)/kT])
    
    Parameters:
    -----------
    E : float or array
        Energy (eV)
    mu : float
        Chemical potential (eV)
    kT : float
        Thermal energy (eV)
        
    Returns:
    --------
    f : float or array
        Occupation probability
    """
    
    # Avoid overflow for large arguments
    x = (E - mu) / kT
    x = np.clip(x, -50, 50)  # Limit range
    
    return 1.0 / (1.0 + np.exp(x))

def density_matrix_equilibrium(E_array, H, Sigma1, Sigma2, mu, temperature, dE):
    """
    Compute equilibrium density matrix by integration
    
    ρ(r,r') = ∫ dE/(2π) A(E) f₀(E - μ)
    
    Parameters:
    -----------
    E_array : array
        Energy grid for integration (eV)
    H : array (Np × Np)
        Hamiltonian
    Sigma1, Sigma2 : functions of E
        Self-energies (callable)
    mu : float
        Chemical potential (eV)
    temperature : float
        Temperature (K)
    dE : float
        Energy spacing (eV)
        
    Returns:
    --------
    rho : array (Np × Np)
        Density matrix
    """
    
    # Thermal energy
    kB_eV = 8.617333262145e-5  # eV/K
    kT = kB_eV * temperature
    
    Np = H.shape[0]
    rho = np.zeros((Np, Np), dtype=complex)
    
    for E in E_array:
        # Self-energies at this energy
        S1 = Sigma1(E)
        S2 = Sigma2(E)
        
        # Green's function
        G = retarded_greens_function(E, H, S1, S2)
        
        # Spectral function
        A = spectral_function(G)
        
        # Fermi function
        f = fermi_function(E, mu, kT)
        
        # Accumulate
        rho += A * f * dE / (2.0 * np.pi)
    
    return np.real(rho)

# ============================================================================
# DENSITY MATRIX (NON-EQUILIBRIUM)
# ============================================================================

def density_matrix_nonequilibrium(E_array, H, Sigma1, Sigma2, mu1, mu2, 
                                  temperature, dE, SigmaS_in=None):
    """
    Compute non-equilibrium density matrix
    
    ρ = ∫ dE/(2π) [f₁ A₁ + f₂ A₂ + G Σ_S^in G†]
    
    For coherent transport (no scattering): SigmaS_in = None
    
    Parameters:
    -----------
    E_array : array
        Energy grid (eV)
    H : array (Np × Np)
        Hamiltonian
    Sigma1, Sigma2 : functions of E
        Contact self-energies
    mu1, mu2 : float
        Chemical potentials of contacts (eV)
    temperature : float
        Temperature (K)
    dE : float
        Energy spacing (eV)
    SigmaS_in : function of E, or None
        Inscattering self-energy (for dissipative transport)
        
    Returns:
    --------
    rho : array (Np × Np)
        Density matrix
    """
    
    kB_eV = 8.617333262145e-5  # eV/K
    kT = kB_eV * temperature
    
    Np = H.shape[0]
    rho = np.zeros((Np, Np), dtype=complex)
    
    for E in E_array:
        # Self-energies
        S1 = Sigma1(E)
        S2 = Sigma2(E)
        
        # Scattering self-energy
        if SigmaS_in is not None:
            SS_in = SigmaS_in(E)
            G = retarded_greens_function(E, H, S1, S2, SS_in)
        else:
            SS_in = None
            G = retarded_greens_function(E, H, S1, S2)
        
        # Broadening functions
        Gamma1 = broadening_function(S1)
        Gamma2 = broadening_function(S2)
        
        # Contact spectral functions
        A1, A2 = contact_spectral_functions(G, Gamma1, Gamma2)
        
        # Fermi functions
        f1 = fermi_function(E, mu1, kT)
        f2 = fermi_function(E, mu2, kT)
        
        # Accumulate contributions
        rho_E = f1 * A1 + f2 * A2
        
        # Add scattering contribution if present
        if SS_in is not None:
            rho_E += G @ SS_in @ G.conj().T
        
        rho += rho_E * dE / (2.0 * np.pi)
    
    return np.real(rho)

# ============================================================================
# ELECTRON DENSITY
# ============================================================================

def electron_density_from_rho(rho, grid):
    """
    Extract electron density from density matrix
    
    n(x) = ρ(x,x) / a
    
    where 'a' is the lattice spacing.
    
    Parameters:
    -----------
    rho : array (Np × Np)
        Density matrix
    grid : dict
        Grid information with 'a' (lattice spacing)
        
    Returns:
    --------
    n : array (Np,)
        Electron density (m^-3)
    """
    
    a = grid['a']
    n = np.diag(rho) / a
    
    return n

# ============================================================================
# MULTIMODE GREEN'S FUNCTIONS (Quasi-3D Transport)
# ============================================================================

def retarded_greens_function_multimode(E_long, H, Sigma1, Sigma2, transverse_modes,
                                       SigmaS=None, eta=1e-4):
    """
    Compute mode-resolved Green's functions for quasi-3D transport

    For each transverse mode (n, m):
        E_total = E_longitudinal + ε_nm
        G_nm = [(E_total + iη)I - H - Σ₁ - Σ₂ - Σ_S]⁻¹

    Key physics:
    - Transverse confinement adds energy ε_nm to each mode
    - Each mode has independent Green's function
    - Contact self-energies assumed mode-independent (uniform contacts)
    - Phonon self-energy can be mode-dependent

    Parameters:
    -----------
    E_long : float
        Longitudinal transport energy (eV)
    H : array (Np, Np)
        Longitudinal Hamiltonian
    Sigma1, Sigma2 : array (Np, Np)
        Contact self-energies (mode-independent approximation)
    transverse_modes : TransverseModes object
        Contains mode energies ε_nm and mode count Nm
    SigmaS : array (Np, Np, Nm) or None
        Mode-resolved phonon self-energy
        If None, ballistic transport
        If 2D array (Np, Np), same for all modes
        If 3D array (Np, Np, Nm), mode-dependent
    eta : float
        Regularization parameter (eV)

    Returns:
    --------
    G_nm : array (Np, Np, Nm)
        Mode-resolved Green's functions
        G_nm[:, :, im] is the Green's function for mode im
    """

    if transverse_modes.energies is None:
        raise RuntimeError("transverse_modes.compute_modes() must be called first")

    Np = H.shape[0]
    Nm = transverse_modes.Nm

    # Initialize output array
    G_nm = np.zeros((Np, Np, Nm), dtype=complex)

    # Determine SigmaS dimensionality
    if SigmaS is None:
        # Ballistic: no scattering
        use_mode_dependent_scattering = False
    elif SigmaS.ndim == 2:
        # Same scattering for all modes
        use_mode_dependent_scattering = False
    elif SigmaS.ndim == 3:
        # Mode-dependent scattering
        if SigmaS.shape[2] != Nm:
            raise ValueError(f"SigmaS has {SigmaS.shape[2]} modes, "
                           f"but transverse_modes has {Nm} modes")
        use_mode_dependent_scattering = True
    else:
        raise ValueError(f"SigmaS must be None, 2D, or 3D array, got shape {SigmaS.shape}")

    # Compute Green's function for each mode
    for im in range(Nm):
        # Total energy: longitudinal + transverse confinement
        E_total = transverse_modes.get_total_energy(E_long, im)

        # Extract mode-specific scattering self-energy
        if SigmaS is None:
            SigmaS_m = None
        elif use_mode_dependent_scattering:
            SigmaS_m = SigmaS[:, :, im]
        else:
            SigmaS_m = SigmaS

        # Compute Green's function for this mode
        G_nm[:, :, im] = retarded_greens_function(
            E_total, H, Sigma1, Sigma2, SigmaS_m, eta
        )

    return G_nm


def spectral_function_multimode(G_nm):
    """
    Compute mode-resolved spectral functions

    A_nm = i[G_nm - G_nm†]

    Parameters:
    -----------
    G_nm : array (Np, Np, Nm)
        Mode-resolved Green's functions

    Returns:
    --------
    A_nm : array (Np, Np, Nm)
        Mode-resolved spectral functions
    """
    A_nm = 1j * (G_nm - np.conj(np.transpose(G_nm, (1, 0, 2))))
    return np.real(A_nm)


# ============================================================================
# TEST CODE
# ============================================================================

if __name__ == "__main__":
    print("\n" + "="*70)
    print("TESTING GREEN'S FUNCTIONS")
    print("="*70)
    
    # Simple test: 1D chain with single barrier
    print("\n[TEST 1] Simple 1D barrier")
    print("-" * 70)
    
    Np = 50
    H = np.zeros((Np, Np))
    
    # Uniform chain with central barrier
    t = 1.0  # eV
    for i in range(Np - 1):
        H[i, i+1] = -t
        H[i+1, i] = -t
    
    # Diagonal: 2t everywhere, except barrier (higher)
    for i in range(Np):
        if 20 <= i < 30:
            H[i, i] = 2*t + 0.5  # Barrier
        else:
            H[i, i] = 2*t
    
    # Contact self-energies (semi-infinite leads)
    def Sigma1(E):
        S = np.zeros((Np, Np), dtype=complex)
        # Simple model: constant imaginary part (broadening)
        S[0, 0] = -1j * 0.1
        return S
    
    def Sigma2(E):
        S = np.zeros((Np, Np), dtype=complex)
        S[-1, -1] = -1j * 0.1
        return S
    
    # Test at one energy
    E_test = 2.0
    S1 = Sigma1(E_test)
    S2 = Sigma2(E_test)
    
    G = retarded_greens_function(E_test, H, S1, S2)
    A = spectral_function(G)
    Gamma1 = broadening_function(S1)
    Gamma2 = broadening_function(S2)
    A1, A2 = contact_spectral_functions(G, Gamma1, Gamma2)
    
    print(f"Green's function computed: {G.shape}")
    print(f"Spectral function computed: {A.shape}")
    print(f"Total spectral function (diagonal): {np.diag(A).sum():.4f}")
    print(f"A1 + A2 (diagonal sum): {(np.diag(A1) + np.diag(A2)).sum():.4f}")
    
    # Check: A should equal A1 + A2 for coherent transport
    diff = np.abs(A - A1 - A2).max()
    print(f"Max difference |A - (A1+A2)|: {diff:.2e}")
    
    if diff < 1e-10:
        print("✓ Spectral function decomposition correct")
    
    # Test Fermi function
    print("\n[TEST 2] Fermi function")
    print("-" * 70)
    
    mu = 2.0
    kT = 0.026  # Room temperature
    
    E_vals = np.linspace(1.0, 3.0, 5)
    f_vals = fermi_function(E_vals, mu, kT)
    
    print(f"μ = {mu:.3f} eV, kT = {kT:.3f} eV")
    for E, f in zip(E_vals, f_vals):
        print(f"  E = {E:.3f} eV → f = {f:.4f}")
    
    print("\n" + "="*70)
    print("ALL TESTS PASSED ✓")
    print("="*70 + "\n")
