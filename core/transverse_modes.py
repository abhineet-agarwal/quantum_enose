"""
Transverse Mode Infrastructure for Quasi-3D Transport

Manages transverse mode energies and wavefunctions for quantum waveguide transport.
Implements hard-wall confinement in the transverse (y-z) plane.

Key Physics:
- Mode energies: ε_nm = (ℏ²π²/2m*)[(n²/Ly²) + (m²/Lz²)]
- Wavefunctions: ψ_nm(y,z) = (2/√(Ly*Lz)) sin(nπy/Ly) sin(mπz/Lz)
- Total energy: E_total = E_longitudinal + ε_nm
"""

import numpy as np
import sys
import os

# Import from hamiltonian module
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

from core.hamiltonian import compute_transverse_modes, psi_nm_at_point

# ============================================================================
# TRANSVERSE MODES CLASS
# ============================================================================

class TransverseModes:
    """
    Manages transverse mode energies and wavefunctions for quasi-3D transport

    Hard-wall confinement in y-z plane:
        ε_nm = (ℏ²π²/2m*) × [(n²/Ly²) + (m²/Lz²)]

    Where n, m = 1, 2, 3, ... are mode indices.

    Attributes:
    -----------
    Ly, Lz : float
        Transverse dimensions (m)
    n_max, m_max : int
        Maximum mode indices (generates n_max × m_max modes)
    m_trans : float
        Transverse effective mass (kg)
    energies : array
        Mode energies ε_nm in eV, shape (Nm,)
    modes : array
        Mode index pairs (n, m), shape (Nm, 2)
    Nm : int
        Total number of modes = n_max × m_max
    """

    def __init__(self, Ly, Lz, n_max, m_max, m_trans):
        """
        Initialize transverse mode configuration

        Parameters:
        -----------
        Ly : float
            Transverse width in y direction (m)
        Lz : float
            Transverse width in z direction (m)
        n_max : int
            Maximum mode index in y direction
        m_max : int
            Maximum mode index in z direction
        m_trans : float
            Transverse effective mass (kg)
        """
        self.Ly = Ly
        self.Lz = Lz
        self.n_max = n_max
        self.m_max = m_max
        self.m_trans = m_trans

        # Will be computed by compute_modes()
        self.energies = None  # ε_nm array (eV)
        self.modes = None     # (n, m) pairs
        self.Nm = None        # Total number of modes

        # Fermi function cache for mode ranking (future use)
        self._fermi_weights = None

    def compute_modes(self):
        """
        Compute all transverse mode energies and indices

        Uses hard-wall quantization formula:
            ε_nm = (ℏ²π²/2m*) × [(n²/Ly²) + (m²/Lz²)]

        Modes are ordered by energy (lowest first).
        """
        # Use existing function from hamiltonian.py
        self.energies, self.modes = compute_transverse_modes(
            self.Ly, self.Lz,
            self.n_max, self.m_max,
            self.m_trans
        )

        # Sort by energy (should already be sorted, but ensure)
        sort_idx = np.argsort(self.energies)
        self.energies = self.energies[sort_idx]
        self.modes = self.modes[sort_idx]

        self.Nm = len(self.energies)

    def get_wavefunction(self, n, m, y, z):
        """
        Evaluate transverse mode wavefunction at point (y, z)

        ψ_nm(y,z) = (2/√(Ly*Lz)) × sin(nπy/Ly) × sin(mπz/Lz)

        Parameters:
        -----------
        n, m : int
            Mode indices
        y, z : float
            Position in transverse plane (m)

        Returns:
        --------
        psi : float
            Wavefunction amplitude (m^-1)
        """
        return psi_nm_at_point(n, m, y, z, self.Ly, self.Lz)

    def get_wavefunction_by_index(self, mode_idx, y, z):
        """
        Evaluate wavefunction using linear mode index

        Parameters:
        -----------
        mode_idx : int
            Linear mode index (0 to Nm-1)
        y, z : float
            Position in transverse plane (m)

        Returns:
        --------
        psi : float
            Wavefunction amplitude (m^-1)
        """
        n, m = self.modes[mode_idx]
        return self.get_wavefunction(n, m, y, z)

    def get_total_energy(self, E_long, mode_idx):
        """
        Combine longitudinal and transverse energies

        E_total = E_longitudinal + ε_nm

        Parameters:
        -----------
        E_long : float
            Longitudinal transport energy (eV)
        mode_idx : int
            Mode index (0 to Nm-1)

        Returns:
        --------
        E_total : float
            Total energy including transverse confinement (eV)
        """
        if self.energies is None:
            raise RuntimeError("Must call compute_modes() first")

        return E_long + self.energies[mode_idx]

    def rank_modes(self, ymol, zmol, temperature):
        """
        Rank modes by thermal-weighted spatial overlap

        Metric: |ψ_nm(ymol, zmol)|² × f_FD(ε_nm, T)

        This is used for hybrid mode selection (coherent vs inelastic).
        Modes with high overlap and thermal population contribute most.

        Parameters:
        -----------
        ymol, zmol : float
            Molecule position in transverse plane (m)
        temperature : float
            Temperature (K)

        Returns:
        --------
        ranking : array (Nm,)
            Mode indices sorted by importance (highest first)
        weights : array (Nm,)
            Corresponding weight values
        """
        if self.energies is None:
            raise RuntimeError("Must call compute_modes() first")

        kB_eV = 8.617333262145e-5  # eV/K
        kT = kB_eV * temperature

        weights = np.zeros(self.Nm)

        for i in range(self.Nm):
            n, m = self.modes[i]

            # Spatial overlap: |ψ_nm(ymol, zmol)|²
            psi = self.get_wavefunction(n, m, ymol, zmol)
            spatial_weight = psi**2

            # Thermal population: Fermi-Dirac at mode energy
            # Assuming modes are above Fermi level, use f_FD(ε_nm)
            fermi_weight = 1.0 / (1.0 + np.exp(self.energies[i] / kT))

            weights[i] = spatial_weight * fermi_weight

        # Sort by weight (descending)
        ranking = np.argsort(weights)[::-1]

        return ranking, weights[ranking]

    def get_mode_info(self, mode_idx):
        """
        Get human-readable information for a mode

        Parameters:
        -----------
        mode_idx : int
            Mode index (0 to Nm-1)

        Returns:
        --------
        info : dict
            Contains 'n', 'm', 'energy_meV', 'energy_eV'
        """
        if self.energies is None:
            raise RuntimeError("Must call compute_modes() first")

        n, m = self.modes[mode_idx]
        return {
            'index': mode_idx,
            'n': int(n),
            'm': int(m),
            'energy_eV': float(self.energies[mode_idx]),
            'energy_meV': float(self.energies[mode_idx] * 1000)
        }

    def print_summary(self):
        """Print summary of transverse mode configuration"""
        if self.energies is None:
            print("TransverseModes (not computed yet)")
            print(f"  Configuration: n_max={self.n_max}, m_max={self.m_max}")
            print(f"  Dimensions: Ly={self.Ly*1e6:.2f} µm, Lz={self.Lz*1e6:.2f} µm")
            return

        print("TransverseModes Summary")
        print("=" * 60)
        print(f"Transverse dimensions: Ly = {self.Ly*1e6:.2f} µm, Lz = {self.Lz*1e6:.2f} µm")
        print(f"Effective mass: m* = {self.m_trans/9.10938356e-31:.3f} m₀")
        print(f"Mode grid: n_max = {self.n_max}, m_max = {self.m_max}")
        print(f"Total modes: {self.Nm}")
        print()
        print("Lowest 5 modes:")
        print("  Mode    (n,m)   Energy (meV)")
        print("-" * 35)
        for i in range(min(5, self.Nm)):
            n, m = self.modes[i]
            print(f"  {i:3d}     ({n},{m})     {self.energies[i]*1000:6.2f}")

        if self.Nm > 5:
            print("  ...")
            i = self.Nm - 1
            n, m = self.modes[i]
            print(f"  {i:3d}     ({n},{m})     {self.energies[i]*1000:6.2f}")
        print("=" * 60)


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def create_default_modes(device, n_max=3, m_max=3):
    """
    Create TransverseModes with default parameters from device

    Assumes square cross-section of 1 µm × 1 µm.

    Parameters:
    -----------
    device : dict
        Device specification with 'm_eff' (relative to m₀)
    n_max, m_max : int
        Maximum mode indices

    Returns:
    --------
    transverse_modes : TransverseModes
        Configured and computed mode object
    """
    m0 = 9.10938356e-31  # kg
    m_trans = device['m_eff'] * m0

    # Default: 1 µm × 1 µm cross-section
    Ly = 1.0e-6  # m
    Lz = 1.0e-6  # m

    modes = TransverseModes(Ly, Lz, n_max, m_max, m_trans)
    modes.compute_modes()

    return modes


# ============================================================================
# TESTING CODE
# ============================================================================

if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("TEST: TransverseModes Class")
    print("=" * 70 + "\n")

    # Test 1: GaAs modes
    print("Test 1: GaAs transverse modes (m* = 0.067 m₀)")
    print("-" * 70)
    m0 = 9.10938356e-31  # kg
    m_GaAs = 0.067 * m0

    modes = TransverseModes(Ly=1e-6, Lz=1e-6, n_max=3, m_max=3, m_trans=m_GaAs)
    modes.compute_modes()
    modes.print_summary()

    # Test 2: Wavefunction evaluation
    print("\nTest 2: Wavefunction at center (Ly/2, Lz/2)")
    print("-" * 70)
    y_center = modes.Ly / 2
    z_center = modes.Lz / 2

    for i in range(min(3, modes.Nm)):
        n, m = modes.modes[i]
        psi = modes.get_wavefunction(n, m, y_center, z_center)
        print(f"  Mode {i} (n={n}, m={m}): ψ = {psi:.3e} m⁻¹")

    # Test 3: Total energy calculation
    print("\nTest 3: Total energy (E_long = 0.5 eV)")
    print("-" * 70)
    E_long = 0.5  # eV

    for i in range(min(3, modes.Nm)):
        E_total = modes.get_total_energy(E_long, i)
        print(f"  Mode {i}: E_total = {E_total:.4f} eV (ε_nm = {modes.energies[i]*1000:.2f} meV)")

    # Test 4: Mode ranking
    print("\nTest 4: Mode ranking (molecule at center, T=300K)")
    print("-" * 70)
    ranking, weights = modes.rank_modes(y_center, z_center, 300)

    print("  Rank  Mode  (n,m)   Weight")
    for rank in range(min(5, modes.Nm)):
        idx = ranking[rank]
        n, m = modes.modes[idx]
        print(f"  {rank+1:4d}  {idx:4d}  ({n},{m})   {weights[rank]:.3e}")

    print("\n" + "=" * 70)
    print("✓ All tests completed successfully!")
    print("=" * 70 + "\n")
