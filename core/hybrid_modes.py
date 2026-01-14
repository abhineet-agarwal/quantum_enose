"""
Hybrid Mode Selection for Quasi-3D Transport

Implements intelligent mode selection to achieve 7-8x speedup:
- Top M "important" modes: Full SCBA (inelastic scattering)
- Remaining modes: Coherent (ballistic, no scattering)

Importance metric: |ψ_nm(y_mol, z_mol)|² × f_FD(ε_nm)
This captures both spatial overlap AND thermal occupancy.
"""

import numpy as np

def rank_modes_by_importance(transverse_modes, phonon_modes, temperature,
                             y_mol=None, z_mol=None):
    """
    Rank transverse modes by importance for molecular scattering.

    Importance metric: coupling_weight × thermal_occupancy
    - coupling_weight: |ψ_nm(y_mol, z_mol)|² (spatial overlap)
    - thermal_occupancy: f_FD(ε_nm) (Fermi-Dirac at mode energy)

    Parameters:
    -----------
    transverse_modes : TransverseModes
        Transverse mode object
    phonon_modes : list
        List of phonon modes (to extract molecular coupling)
    temperature : float
        Temperature in Kelvin
    y_mol, z_mol : float or None
        Molecule position (defaults to center)

    Returns:
    --------
    mode_indices : array
        Indices sorted by importance (most important first)
    importance_scores : array
        Importance score for each mode
    """

    kB_eV = 8.617333262145e-5
    kT = kB_eV * temperature

    Nm = transverse_modes.Nm

    # Default molecule position: center
    if y_mol is None:
        y_mol = transverse_modes.Ly / 2.0
    if z_mol is None:
        z_mol = transverse_modes.Lz / 2.0

    # Compute coupling weights (spatial overlap)
    coupling_weights = np.zeros(Nm)
    for im in range(Nm):
        n, m = transverse_modes.modes[im]
        psi = transverse_modes.get_wavefunction(n, m, y_mol, z_mol)
        coupling_weights[im] = psi**2

    # Normalize coupling weights [0, 1]
    max_coupling = np.max(coupling_weights)
    if max_coupling > 0:
        coupling_weights /= max_coupling

    # Compute thermal occupancy (Fermi-Dirac at mode energy)
    # Assume Fermi level at E=0 for simplicity
    thermal_occupancy = np.zeros(Nm)
    for im in range(Nm):
        E_mode = transverse_modes.energies[im]
        # f_FD(E) = 1/(1 + exp(E/kT))
        thermal_occupancy[im] = 1.0 / (1.0 + np.exp(E_mode / kT))

    # Combined importance score
    # Use weighted sum instead of product to avoid zero-importance for modes
    # with nodes at molecule position (e.g., n=2 or m=2 at center)
    # Thermal occupancy is primary (all populated modes contribute)
    # Coupling is secondary bonus (modes with coupling scatter more)
    alpha = 0.3  # Weight for coupling (secondary)
    importance_scores = alpha * coupling_weights + (1 - alpha) * thermal_occupancy

    # Sort by importance (descending)
    mode_indices = np.argsort(importance_scores)[::-1]

    return mode_indices, importance_scores


def select_hybrid_modes(transverse_modes, phonon_modes, temperature,
                       n_inelastic=4, y_mol=None, z_mol=None):
    """
    Select modes for hybrid approach.

    Parameters:
    -----------
    transverse_modes : TransverseModes
        Transverse mode object
    phonon_modes : list
        List of phonon modes
    temperature : float
        Temperature in Kelvin
    n_inelastic : int
        Number of modes to treat with full SCBA (default: 4)
    y_mol, z_mol : float or None
        Molecule position

    Returns:
    --------
    inelastic_modes : array
        Indices of modes to treat with SCBA
    coherent_modes : array
        Indices of modes to treat as coherent
    importance_scores : array
        Importance score for each mode (for analysis)
    """

    # Rank modes
    mode_indices, importance_scores = rank_modes_by_importance(
        transverse_modes, phonon_modes, temperature, y_mol, z_mol
    )

    # Top n_inelastic modes: full SCBA
    inelastic_modes = mode_indices[:n_inelastic]

    # Remaining modes: coherent
    coherent_modes = mode_indices[n_inelastic:]

    return inelastic_modes, coherent_modes, importance_scores


def print_mode_selection(transverse_modes, inelastic_modes, coherent_modes,
                        importance_scores):
    """
    Print mode selection summary.
    """

    print(f"\n  Hybrid Mode Selection:")
    print(f"    Inelastic modes (full SCBA): {len(inelastic_modes)}/{transverse_modes.Nm}")
    print(f"    Coherent modes (ballistic):  {len(coherent_modes)}/{transverse_modes.Nm}")

    print(f"\n  Mode ranking by importance:")
    print(f"    {'Rank':<6} {'Mode':<10} {'(n,m)':<8} {'Energy (meV)':<15} {'Importance':<12} {'Treatment'}")
    print(f"    {'-'*70}")

    # Sort by importance for display
    sorted_indices = np.argsort(importance_scores)[::-1]

    for rank, im in enumerate(sorted_indices):
        n, m = transverse_modes.modes[im]
        E_mode = transverse_modes.energies[im] * 1000  # meV
        importance = importance_scores[im]

        if im in inelastic_modes:
            treatment = "Inelastic (SCBA)"
        else:
            treatment = "Coherent"

        print(f"    {rank+1:<6} {im:<10} ({n},{m}){' ':<3} {E_mode:<15.3f} {importance:<12.6f} {treatment}")


def estimate_speedup(n_inelastic, n_total):
    """
    Estimate speedup from hybrid mode selection.

    Speedup ≈ n_total / n_inelastic
    (Assumes SCBA time dominates, which is true for quasi-3D)

    Parameters:
    -----------
    n_inelastic : int
        Number of modes with full SCBA
    n_total : int
        Total number of modes

    Returns:
    --------
    speedup : float
        Expected speedup factor
    """

    if n_inelastic == 0:
        return 1.0  # No inelastic modes = pure coherent

    # Rough model:
    # Time ≈ n_inelastic × t_scba + n_coherent × t_coherent
    # where t_scba >> t_coherent
    #
    # For all-inelastic: Time_all = n_total × t_scba
    # For hybrid: Time_hybrid = n_inelastic × t_scba + n_coherent × t_coherent
    #
    # Since t_scba >> t_coherent, approximate:
    # Speedup ≈ n_total / n_inelastic

    speedup = n_total / n_inelastic

    return speedup


def validate_mode_selection(inelastic_modes, coherent_modes, total_modes):
    """
    Validate that mode selection is consistent.
    """

    # Check: no overlap
    overlap = set(inelastic_modes) & set(coherent_modes)
    if len(overlap) > 0:
        raise ValueError(f"Mode selection has overlap: {overlap}")

    # Check: covers all modes
    all_selected = set(inelastic_modes) | set(coherent_modes)
    if len(all_selected) != total_modes:
        raise ValueError(f"Mode selection doesn't cover all modes: {len(all_selected)} vs {total_modes}")

    # Check: indices in range
    if np.any(inelastic_modes < 0) or np.any(inelastic_modes >= total_modes):
        raise ValueError(f"Invalid inelastic mode indices")
    if np.any(coherent_modes < 0) or np.any(coherent_modes >= total_modes):
        raise ValueError(f"Invalid coherent mode indices")

    return True
