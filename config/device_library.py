"""
Device Library: Pre-defined RTD structures

Provides commonly used RTD configurations for quantum e-nose simulations.
Users can select devices by name or define custom structures.

Band Offset References:
- IGZO/HfO2: Hays et al., Appl. Phys. Rev. 4, 021301 (2017)
- ZnO/Al2O3, ZnO/HfO2: PMC5438334, XPS measurements
- GaAs/AlAs: Standard 60:40 rule, validated
- Oxide systems: Extrapolated from XPS and IPE measurements
"""

import numpy as np

# ============================================================================
# CONDUCTION BAND OFFSETS (Literature Values)
# ============================================================================
# Key: (well_material, barrier_material) -> CBO in eV
# CBO = Ec_barrier - Ec_well (positive means barrier is higher)

BAND_OFFSETS = {
    # III-V semiconductors (well-established)
    ("GaAs", "AlAs"): 0.57,       # 60:40 rule, validated
    ("InGaAs", "InAlAs"): 0.67,   # Higher offset for InGaAs wells

    # Oxide semiconductors with Al2O3 barrier
    # Al2O3 bandgap ~7.0 eV, electron affinity ~1.0 eV
    ("In2O3", "Al2O3"): 2.8,      # In2O3 Eg~3.6 eV, χ~4.3 eV
    ("IGZO", "Al2O3"): 2.9,       # IGZO Eg~3.2 eV, slightly lower χ
    ("ZnO", "Al2O3"): 3.0,        # ZnO Eg~3.3 eV, χ~4.5 eV (XPS measured)
    ("SnO2", "Al2O3"): 2.9,       # SnO2 Eg~3.6 eV, χ~4.5 eV

    # Oxide semiconductors with HfO2 barrier
    # HfO2 bandgap ~5.6 eV, electron affinity ~2.5 eV
    ("In2O3", "HfO2"): 2.0,       # Lower offset than Al2O3
    ("IGZO", "HfO2"): 2.3,        # Hays et al. APR 2017: 2.26-2.39 eV
    ("ZnO", "HfO2"): 2.2,         # XPS measured (PMC5438334)
    ("SnO2", "HfO2"): 2.2,        # Estimated from similar oxides

    # Oxide semiconductors with SiO2 barrier
    # SiO2 bandgap ~9.0 eV, electron affinity ~0.9 eV
    ("In2O3", "SiO2"): 3.4,       # Largest offset
    ("IGZO", "SiO2"): 3.5,        # From Hays et al. review
    ("ZnO", "SiO2"): 3.5,         # Estimated
    ("SnO2", "SiO2"): 3.5,        # Estimated

    # ==========================================================================
    # LOW-BARRIER OXIDE COMBINATIONS (IDEAL FOR IETS - similar to Patil 0.6 eV)
    # ==========================================================================
    # References:
    # - ACS Appl. Electron. Mater. 2020: κ-Ga2O3/In2O3 CBO = 0.45 eV
    # - In2O3/Ga2O3 CBO = 0.07 eV (type-I alignment)
    # - Scientific Reports 2017: ZnO/MgZnO tunable CBO

    # In2O3/Ga2O3 system (PRIMARY CANDIDATE FOR IETS)
    ("In2O3", "Ga2O3"): 0.07,           # β-Ga2O3 barrier, very low CBO
    ("In2O3", "Ga2O3_kappa"): 0.45,     # κ-Ga2O3 barrier, ideal for IETS!

    # ZnO/MgZnO system (tunable barrier)
    ("ZnO", "MgZnO"): 0.47,             # For Mg0.3Zn0.7O (CBO = 1.57x)

    # Cross-combinations
    ("IGZO", "Ga2O3"): 0.10,            # Estimated from similar systems
    ("ZnO", "Ga2O3"): 0.61,             # From ZnO/α-Ga2O3 XPS measurements
}

def get_band_offset(well_material, barrier_material):
    """
    Get conduction band offset for a specific well/barrier pair.

    Parameters:
    -----------
    well_material : str
        Name of the well (channel) material
    barrier_material : str
        Name of the barrier material

    Returns:
    --------
    float : Conduction band offset in eV (Ec_barrier - Ec_well)
    """
    key = (well_material, barrier_material)
    if key in BAND_OFFSETS:
        return BAND_OFFSETS[key]

    # Fallback: use the Ec values from MATERIALS
    well = MATERIALS.get(well_material, {})
    barrier = MATERIALS.get(barrier_material, {})

    if 'Ec' in well and 'Ec' in barrier:
        return barrier['Ec'] - well['Ec']

    raise ValueError(f"No band offset defined for {well_material}/{barrier_material}")

# ============================================================================
# MATERIAL PROPERTIES
# ============================================================================

MATERIALS = {
    # ========== III-V SEMICONDUCTORS (legacy reference) ==========
    "GaAs": {
        "m_eff": 0.067,  # Effective mass (in units of m0)
        "Ec": 0.0,       # Conduction band edge (eV) - reference
        "Eg": 1.42,      # Band gap (eV)
        "chi": 4.07,     # Electron affinity (eV)
        "epsilon_r": 12.9,  # Relative permittivity
        "phonon_energy_meV": 36.0,  # LO phonon energy
        "description": "Gallium Arsenide"
    },
    "AlAs": {
        "m_eff": 0.15,
        "Ec": 0.57,      # CBO relative to GaAs (60:40 rule)
        "Eg": 2.16,      # Band gap (eV) - indirect
        "chi": 3.50,     # Electron affinity (eV)
        "epsilon_r": 10.1,
        "phonon_energy_meV": 50.0,
        "description": "Aluminum Arsenide (barrier)"
    },
    "InGaAs": {
        "m_eff": 0.041,  # In0.53Ga0.47As lattice-matched to InP
        "Ec": -0.15,     # Below GaAs CB
        "Eg": 0.74,      # Band gap (eV)
        "chi": 4.50,     # Electron affinity (eV)
        "epsilon_r": 13.9,
        "phonon_energy_meV": 30.0,
        "description": "Indium Gallium Arsenide"
    },
    "InAlAs": {
        "m_eff": 0.075,
        "Ec": 0.52,      # Relative to InGaAs: 0.67 eV offset
        "Eg": 1.46,      # Band gap (eV)
        "chi": 4.05,     # Electron affinity (eV)
        "epsilon_r": 12.2,
        "phonon_energy_meV": 43.0,
        "description": "Indium Aluminum Arsenide (barrier)"
    },

    # ========== OXIDE SEMICONDUCTORS (BEOL-compatible, <400°C) ==========
    "In2O3": {
        "m_eff": 0.30,   # Effective mass (Presley et al. JAP 2004)
        "Ec": 0.0,       # Reference level for oxide system
        "Eg": 3.6,       # Direct band gap (eV)
        "chi": 4.3,      # Electron affinity (eV)
        "epsilon_r": 9.0,
        "phonon_energy_meV": 70.0,  # Optical phonon (Walsh et al. PRB 2009)
        "description": "Indium Oxide (BEOL-compatible)",
        "process_temp": 350,
        "notes": "ALD/sputtering <400°C, excellent BEOL compatibility"
    },
    "SnO2": {
        "m_eff": 0.25,   # Effective mass (Godinho et al. JPC 2009)
        "Ec": 0.0,       # Reference (adjusted for each heterostructure)
        "Eg": 3.6,       # Direct band gap (eV)
        "chi": 4.5,      # Electron affinity (eV)
        "epsilon_r": 9.0,
        "phonon_energy_meV": 75.0,
        "description": "Tin Oxide (BEOL-compatible)",
        "process_temp": 300,
        "notes": "Sputtering/ALD <400°C, good stability"
    },
    "IGZO": {
        "m_eff": 0.34,   # InGaZnO (Nomura et al. Nature 2004)
        "Ec": 0.0,       # Reference for IGZO system
        "Eg": 3.2,       # Band gap (eV) - amorphous
        "chi": 4.2,      # Electron affinity (eV)
        "epsilon_r": 10.0,
        "phonon_energy_meV": 60.0,
        "description": "Indium Gallium Zinc Oxide (amorphous)",
        "process_temp": 300,
        "notes": "Sputtering <350°C, amorphous, excellent uniformity"
    },
    "ZnO": {
        "m_eff": 0.28,   # Effective mass (Look et al. SSC 1998)
        "Ec": 0.0,       # Reference for ZnO system
        "Eg": 3.3,       # Direct band gap (eV)
        "chi": 4.5,      # Electron affinity (eV)
        "epsilon_r": 8.5,
        "phonon_energy_meV": 72.0,  # LO phonon
        "description": "Zinc Oxide (BEOL-compatible)",
        "process_temp": 350,
        "notes": "ALD/sputtering <400°C, mature process"
    },

    # ========== OXIDE BARRIERS (BEOL-compatible) ==========
    # NOTE: Ec values are reference only. Use get_band_offset() for
    # accurate well/barrier-specific offsets from literature.
    "Al2O3": {
        "m_eff": 0.45,   # High effective mass (insulator)
        "Ec": 2.9,       # Average CBO vs oxide wells (~2.8-3.0 eV)
        "Eg": 7.0,       # Band gap (eV) - wide gap insulator
        "chi": 1.0,      # Electron affinity (eV)
        "epsilon_r": 9.0,
        "phonon_energy_meV": 100.0,
        "description": "Aluminum Oxide (barrier, BEOL-compatible)",
        "process_temp": 300,
        "notes": "ALD <350°C, excellent barrier quality, low leakage"
    },
    "HfO2": {
        "m_eff": 0.50,
        "Ec": 2.2,       # Average CBO vs oxide wells (~2.0-2.3 eV)
        "Eg": 5.6,       # Band gap (eV) - high-k dielectric
        "chi": 2.5,      # Electron affinity (eV)
        "epsilon_r": 25.0,
        "phonon_energy_meV": 95.0,
        "description": "Hafnium Oxide (high-k barrier, BEOL-compatible)",
        "process_temp": 350,
        "notes": "ALD <400°C, higher permittivity than Al2O3"
    },
    "SiO2": {
        "m_eff": 0.50,
        "Ec": 3.5,       # Large CBO vs oxide wells (~3.4-3.5 eV)
        "Eg": 9.0,       # Band gap (eV) - very wide gap
        "chi": 0.9,      # Electron affinity (eV)
        "epsilon_r": 3.9,
        "phonon_energy_meV": 120.0,
        "description": "Silicon Dioxide (barrier, BEOL-compatible)",
        "process_temp": 300,
        "notes": "PECVD <400°C, standard CMOS dielectric"
    },

    # ========== LOW-BARRIER OXIDE SEMICONDUCTORS (for IETS) ==========
    "Ga2O3": {
        "m_eff": 0.28,   # Effective mass (β-Ga2O3, Peelaers & Van de Walle PRB 2015)
        "Ec": 0.07,      # CBO relative to In2O3 (from ACS Appl. Electron. Mater. 2020)
        "Eg": 4.9,       # Band gap (eV) - ultra-wide bandgap
        "chi": 4.0,      # Electron affinity (eV)
        "epsilon_r": 10.0,
        "phonon_energy_meV": 44.0,  # Dominant LO phonon
        "description": "Gallium Oxide (low-barrier for In2O3, BEOL-compatible)",
        "process_temp": 300,
        "notes": "ALD 100-400°C, sputtering 300°C, CBO~0.07-0.45 eV vs In2O3"
    },
    "Ga2O3_kappa": {
        "m_eff": 0.28,   # κ-phase Ga2O3
        "Ec": 0.45,      # Higher CBO for κ-phase/In2O3 interface
        "Eg": 4.9,       # Band gap (eV)
        "chi": 3.6,      # Slightly different electron affinity
        "epsilon_r": 10.0,
        "phonon_energy_meV": 44.0,
        "description": "κ-phase Gallium Oxide (0.45 eV barrier vs In2O3)",
        "process_temp": 350,
        "notes": "ALD/sputtering <400°C, CBO=0.45 eV ideal for IETS"
    },
    "MgZnO": {
        "m_eff": 0.30,   # Mg0.3Zn0.7O effective mass
        "Ec": 0.47,      # CBO = 1.57x for x=0.3 Mg content
        "Eg": 4.0,       # Band gap increases with Mg content
        "chi": 4.2,      # Electron affinity
        "epsilon_r": 8.5,
        "phonon_energy_meV": 70.0,
        "description": "MgZnO alloy barrier (tunable CBO)",
        "process_temp": 350,
        "notes": "ALD/sputtering <350°C, CBO=1.57x where x=Mg fraction"
    }
}

# Physical constants
m0 = 9.10938356e-31  # kg (electron mass)
q = 1.602176634e-19   # C (electron charge)
hbar = 1.054571817e-34  # J*s
epsilon_0 = 8.854187817e-12  # F/m

# ============================================================================
# DEVICE LIBRARY
# ============================================================================

DEVICES = {
    # ========== BASELINE DEVICE (from validation) ==========
    "GaAs_AlAs_symmetric": {
        "description": "Symmetric GaAs/AlAs RTD (validated baseline)",
        "layers": [
            {"material": "GaAs", "thickness": 8e-9, "doping": 1e24},  # Emitter
            {"material": "AlAs", "thickness": 1.5e-9, "doping": 0},   # Barrier 1
            {"material": "GaAs", "thickness": 2e-9, "doping": 0},     # Well
            {"material": "AlAs", "thickness": 1.5e-9, "doping": 0},   # Barrier 2
            {"material": "GaAs", "thickness": 8e-9, "doping": 1e24}   # Collector
        ],
        "transverse_size": (1e-6, 1e-6),  # (Ly, Lz) in meters
        "molecule_location": "emitter_barrier",  # Where to place molecule
        "notes": "Validated device from Phase 1, doping=10^18 cm^-3"
    },
    
    # ========== ASYMMETRIC DEVICE (optimized for sensing) ==========
    "GaAs_AlAs_asymmetric": {
        "description": "Asymmetric RTD with wide emitter barrier (optimized)",
        "layers": [
            {"material": "GaAs", "thickness": 8e-9, "doping": 1e24},
            {"material": "AlAs", "thickness": 3.0e-9, "doping": 0},   # WIDE emitter barrier
            {"material": "GaAs", "thickness": 2e-9, "doping": 0},
            {"material": "AlAs", "thickness": 1.0e-9, "doping": 0},   # THIN collector barrier
            {"material": "GaAs", "thickness": 8e-9, "doping": 1e24}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
        "notes": "Wide emitter barrier for larger sensing area (from Patil 2018)"
    },
    
    # ========== THIN BARRIER DEVICE (high current) ==========
    "GaAs_AlAs_thin": {
        "description": "Thin barrier RTD for high current operation",
        "layers": [
            {"material": "GaAs", "thickness": 10e-9, "doping": 5e23},
            {"material": "AlAs", "thickness": 1.0e-9, "doping": 0},
            {"material": "GaAs", "thickness": 2.5e-9, "doping": 0},
            {"material": "AlAs", "thickness": 1.0e-9, "doping": 0},
            {"material": "GaAs", "thickness": 10e-9, "doping": 5e23}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
        "notes": "Thinner barriers, higher current, doping=5×10^17 cm^-3"
    },
    
    # ========== WIDE WELL DEVICE (for lower energy resonance) ==========
    "GaAs_AlAs_wide_well": {
        "description": "Wide well RTD for lower resonance energy",
        "layers": [
            {"material": "GaAs", "thickness": 8e-9, "doping": 1e24},
            {"material": "AlAs", "thickness": 1.5e-9, "doping": 0},
            {"material": "GaAs", "thickness": 4.0e-9, "doping": 0},  # WIDE well
            {"material": "AlAs", "thickness": 1.5e-9, "doping": 0},
            {"material": "GaAs", "thickness": 8e-9, "doping": 1e24}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
        "notes": "Wider well lowers resonance energy"
    },
    
    # ========== INGAAS/INALAS SYSTEM (alternative material) ==========
    "InGaAs_InAlAs_symmetric": {
        "description": "InGaAs/InAlAs RTD (alternative material system)",
        "layers": [
            {"material": "InGaAs", "thickness": 8e-9, "doping": 1e24},
            {"material": "InAlAs", "thickness": 1.5e-9, "doping": 0},
            {"material": "InGaAs", "thickness": 2e-9, "doping": 0},
            {"material": "InAlAs", "thickness": 1.5e-9, "doping": 0},
            {"material": "InGaAs", "thickness": 8e-9, "doping": 1e24}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
        "notes": "Lower effective mass → higher tunneling probability"
    },
    
    # ========== LARGE AREA DEVICE (for practical sensor) ==========
    "GaAs_AlAs_large": {
        "description": "Large area RTD for practical sensing application",
        "layers": [
            {"material": "GaAs", "thickness": 8e-9, "doping": 1e24},
            {"material": "AlAs", "thickness": 2.0e-9, "doping": 0},
            {"material": "GaAs", "thickness": 2e-9, "doping": 0},
            {"material": "AlAs", "thickness": 2.0e-9, "doping": 0},
            {"material": "GaAs", "thickness": 8e-9, "doping": 1e24}
        ],
        "transverse_size": (5e-6, 5e-6),  # 5 µm × 5 µm
        "molecule_location": "emitter_barrier",
        "notes": "Larger area for more modes, realistic sensor size"
    },
    
    # ========== OXIDE RTD DEVICES (BEOL-compatible, <400°C) ==========
    
    "In2O3_Al2O3_symmetric": {
        "description": "In2O3/Al2O3 RTD (PRIMARY CANDIDATE - best BEOL evidence)",
        "layers": [
            {"material": "In2O3", "thickness": 10e-9, "doping": 1e25},  # 10^19 cm^-3
            {"material": "Al2O3", "thickness": 2.0e-9, "doping": 0},    # Barrier 1
            {"material": "In2O3", "thickness": 2.5e-9, "doping": 0},    # Well
            {"material": "Al2O3", "thickness": 2.0e-9, "doping": 0},    # Barrier 2
            {"material": "In2O3", "thickness": 10e-9, "doping": 1e25}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
        "notes": "ALD <350°C, Intel collaboration (VLSI 2021), excellent stability"
    },
    
    "In2O3_Al2O3_asymmetric": {
        "description": "In2O3/Al2O3 asymmetric RTD (optimized for sensing)",
        "layers": [
            {"material": "In2O3", "thickness": 10e-9, "doping": 1e25},
            {"material": "Al2O3", "thickness": 3.0e-9, "doping": 0},    # WIDE emitter
            {"material": "In2O3", "thickness": 2.5e-9, "doping": 0},
            {"material": "Al2O3", "thickness": 1.5e-9, "doping": 0},    # THIN collector
            {"material": "In2O3", "thickness": 10e-9, "doping": 1e25}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
        "notes": "Wide emitter barrier for large sensing area"
    },
    
    "IGZO_Al2O3_symmetric": {
        "description": "IGZO/Al2O3 RTD (amorphous, excellent uniformity)",
        "layers": [
            {"material": "IGZO", "thickness": 10e-9, "doping": 8e24},   # 8×10^18 cm^-3
            {"material": "Al2O3", "thickness": 2.0e-9, "doping": 0},
            {"material": "IGZO", "thickness": 2.5e-9, "doping": 0},
            {"material": "Al2O3", "thickness": 2.0e-9, "doping": 0},
            {"material": "IGZO", "thickness": 10e-9, "doping": 8e24}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
        "notes": "Sputtering <300°C, amorphous = excellent uniformity"
    },
    
    "ZnO_Al2O3_symmetric": {
        "description": "ZnO/Al2O3 RTD (mature ALD process)",
        "layers": [
            {"material": "ZnO", "thickness": 10e-9, "doping": 1e25},
            {"material": "Al2O3", "thickness": 2.0e-9, "doping": 0},
            {"material": "ZnO", "thickness": 2.5e-9, "doping": 0},
            {"material": "Al2O3", "thickness": 2.0e-9, "doping": 0},
            {"material": "ZnO", "thickness": 10e-9, "doping": 1e25}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
        "notes": "ALD <350°C, mature process, good electron mobility"
    },
    
    "SnO2_Al2O3_symmetric": {
        "description": "SnO2/Al2O3 RTD (backup candidate)",
        "layers": [
            {"material": "SnO2", "thickness": 10e-9, "doping": 1e25},
            {"material": "Al2O3", "thickness": 2.0e-9, "doping": 0},
            {"material": "SnO2", "thickness": 2.5e-9, "doping": 0},
            {"material": "Al2O3", "thickness": 2.0e-9, "doping": 0},
            {"material": "SnO2", "thickness": 10e-9, "doping": 1e25}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
        "notes": "Sputtering/ALD <300°C, good stability, lower mobility"
    },
    
    "In2O3_HfO2_symmetric": {
        "description": "In2O3/HfO2 RTD (high-k barrier alternative)",
        "layers": [
            {"material": "In2O3", "thickness": 10e-9, "doping": 1e25},
            {"material": "HfO2", "thickness": 2.0e-9, "doping": 0},
            {"material": "In2O3", "thickness": 2.5e-9, "doping": 0},
            {"material": "HfO2", "thickness": 2.0e-9, "doping": 0},
            {"material": "In2O3", "thickness": 10e-9, "doping": 1e25}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
        "notes": "High-k dielectric, may reduce leakage current"
    },

    # ==========================================================================
    # LOW-BARRIER RTD DEVICES FOR IETS (CBO ~ 0.5 eV like Patil et al.)
    # ==========================================================================
    # These devices use low-barrier oxide combinations suitable for room-temperature
    # IETS molecular detection, following Patil et al. 2018 (Sci. Rep.)

    "In2O3_Ga2O3_symmetric": {
        "description": "In2O3/Ga2O3 RTD (LOW BARRIER - ideal for IETS!)",
        "layers": [
            {"material": "In2O3", "thickness": 10e-9, "doping": 1e25},
            {"material": "Ga2O3", "thickness": 2.0e-9, "doping": 0},    # CBO = 0.07 eV
            {"material": "In2O3", "thickness": 3.0e-9, "doping": 0},    # 3nm well
            {"material": "Ga2O3", "thickness": 2.0e-9, "doping": 0},
            {"material": "In2O3", "thickness": 10e-9, "doping": 1e25}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
        "notes": "CBO=0.07eV, ALD <400°C, PRIMARY IETS CANDIDATE"
    },

    "In2O3_Ga2O3_1nm": {
        "description": "In2O3/Ga2O3 RTD with 1nm barriers",
        "layers": [
            {"material": "In2O3", "thickness": 10e-9, "doping": 1e25},
            {"material": "Ga2O3", "thickness": 1.0e-9, "doping": 0},    # Ultra-thin
            {"material": "In2O3", "thickness": 3.0e-9, "doping": 0},
            {"material": "Ga2O3", "thickness": 1.0e-9, "doping": 0},
            {"material": "In2O3", "thickness": 10e-9, "doping": 1e25}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
        "notes": "Ultra-thin barriers for maximum transmission"
    },

    "In2O3_Ga2O3_kappa_symmetric": {
        "description": "In2O3/κ-Ga2O3 RTD (CBO=0.45eV - matches Patil!)",
        "layers": [
            {"material": "In2O3", "thickness": 10e-9, "doping": 1e25},
            {"material": "Ga2O3_kappa", "thickness": 2.0e-9, "doping": 0},  # CBO = 0.45 eV
            {"material": "In2O3", "thickness": 3.0e-9, "doping": 0},
            {"material": "Ga2O3_kappa", "thickness": 2.0e-9, "doping": 0},
            {"material": "In2O3", "thickness": 10e-9, "doping": 1e25}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
        "notes": "CBO=0.45eV. SISPAD *comparison* stack (demoted 2026-04-07, see docs/STACK_DECISION.md); primary is ZnO_MgZnO_asymmetric. κ-phase metastability + single-source CBO + no demonstrated RTD are the reviewer flags."
    },

    "In2O3_Ga2O3_kappa_1nm": {
        "description": "In2O3/κ-Ga2O3 RTD with 1nm barriers",
        "layers": [
            {"material": "In2O3", "thickness": 10e-9, "doping": 1e25},
            {"material": "Ga2O3_kappa", "thickness": 1.0e-9, "doping": 0},
            {"material": "In2O3", "thickness": 3.0e-9, "doping": 0},
            {"material": "Ga2O3_kappa", "thickness": 1.0e-9, "doping": 0},
            {"material": "In2O3", "thickness": 10e-9, "doping": 1e25}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
        "notes": "CBO=0.45eV, ultra-thin barriers"
    },

    "In2O3_Ga2O3_kappa_asymmetric": {
        "description": "In2O3/κ-Ga2O3 asymmetric RTD (wide emitter for sensing)",
        "layers": [
            {"material": "In2O3", "thickness": 10e-9, "doping": 1e25},
            {"material": "Ga2O3_kappa", "thickness": 3.0e-9, "doping": 0},  # Wide emitter
            {"material": "In2O3", "thickness": 3.0e-9, "doping": 0},
            {"material": "Ga2O3_kappa", "thickness": 1.5e-9, "doping": 0},  # Thin collector
            {"material": "In2O3", "thickness": 10e-9, "doping": 1e25}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
        "notes": "Asymmetric Patil geometry. SISPAD *comparison* stack — the in-flight v3 cloud run uses this at 1 µm. Primary is ZnO_MgZnO_asymmetric (see docs/STACK_DECISION.md)."
    },

    "ZnO_MgZnO_symmetric": {
        "description": "ZnO/Mg0.3Zn0.7O RTD (tunable barrier, CBO~0.47eV) — SISPAD primary (symmetric variant)",
        "layers": [
            {"material": "ZnO", "thickness": 10e-9, "doping": 1e25},
            {"material": "MgZnO", "thickness": 2.0e-9, "doping": 0},    # CBO = 0.47 eV
            {"material": "ZnO", "thickness": 3.0e-9, "doping": 0},
            {"material": "MgZnO", "thickness": 2.0e-9, "doping": 0},
            {"material": "ZnO", "thickness": 10e-9, "doping": 1e25}
        ],
        "transverse_size": (10e-6, 10e-6),  # 10 µm × 10 µm = 10^-7 cm^2 sensor-pixel area (Datta A prefactor; see docs/STACK_DECISION.md)
        "molecule_location": "emitter_barrier",
        "notes": "SISPAD PRIMARY (sym variant). CBO=0.47eV for x_Mg=0.3 (multi-source, tunable as 1.57*x_Mg); ALD-matured process; demonstrated NDR in ZnO/Mg0.33Zn0.67O RTD by Tampo et al., IEEE NANO 2019. m*=0.28 m0 (ZnO well), LO phonon ~72 meV. Swapped in for In2O3/κ-Ga2O3 2026-04-07."
    },

    "ZnO_MgZnO_asymmetric": {
        "description": "ZnO/Mg0.3Zn0.7O RTD, Patil-style asymmetric (3 nm emitter / 1.5 nm collector) — SISPAD primary",
        "layers": [
            {"material": "ZnO", "thickness": 10e-9, "doping": 1e25},
            {"material": "MgZnO", "thickness": 3.0e-9, "doping": 0},    # Wide emitter barrier (sensing area)
            {"material": "ZnO", "thickness": 3.0e-9, "doping": 0},      # Quantum well
            {"material": "MgZnO", "thickness": 1.5e-9, "doping": 0},    # Thin collector (good extraction)
            {"material": "ZnO", "thickness": 10e-9, "doping": 1e25}
        ],
        "transverse_size": (10e-6, 10e-6),  # 10 µm × 10 µm = 10^-7 cm^2 sensor-pixel area (Datta A prefactor; see docs/STACK_DECISION.md)
        "molecule_location": "emitter_barrier",
        "notes": "SISPAD PRIMARY (asym variant). Patil's 3/1.5 nm asymmetry mirrored onto ZnO/Mg0.3Zn0.7O: wide emitter barrier maximises molecule-resonance overlap, thin collector keeps transmission. Physics numbers match In2O3/κ-Ga2O3 closely (CBO 0.47 vs 0.45 eV, m* 0.28 vs 0.30, ℏω 72 vs 70 meV) so results transfer. Gains over the κ-Ga2O3 stack: stable wurtzite phase (vs metastable κ), decades-old ALD maturity, tunable CBO via x_Mg, and an *experimentally demonstrated* NDR device. Added 2026-04-07, paired with launch_cloud_v4.sh and --solver rank1."
    },

    "ZnO_MgZnO_patil_like": {
        "description": "ZnO/Mg0.3Zn0.7O RTD, Patil-like geometry (2nm/5nm/2nm symmetric)",
        "layers": [
            {"material": "ZnO", "thickness": 10e-9, "doping": 1e25},
            {"material": "MgZnO", "thickness": 2.0e-9, "doping": 0},    # CBO = 0.47 eV
            {"material": "ZnO", "thickness": 5.0e-9, "doping": 0},      # 5nm well (Patil used 5nm)
            {"material": "MgZnO", "thickness": 2.0e-9, "doping": 0},
            {"material": "ZnO", "thickness": 10e-9, "doping": 1e25}
        ],
        "transverse_size": (10e-6, 10e-6),
        "molecule_location": "emitter_barrier",
        "notes": "Patil-like thicknesses (2nm/5nm/2nm symmetric) on ZnO/MgZnO stack. Patil's GaAs/AlGaAs had CBO=0.23eV, m*=0.067, V_res~0.336V. This stack has CBO=0.47eV, m*=0.28 — heavier mass and higher barrier shift V_res. The wider 5nm well lowers E1 and should bring NDR into the 0-500mV range. Added 2026-04-09."
    },

    "ZnO_Ga2O3_symmetric": {
        "description": "ZnO/Ga2O3 RTD (CBO~0.61eV - closest to Patil)",
        "layers": [
            {"material": "ZnO", "thickness": 10e-9, "doping": 1e25},
            {"material": "Ga2O3", "thickness": 2.0e-9, "doping": 0},    # CBO = 0.61 eV
            {"material": "ZnO", "thickness": 3.0e-9, "doping": 0},
            {"material": "Ga2O3", "thickness": 2.0e-9, "doping": 0},
            {"material": "ZnO", "thickness": 10e-9, "doping": 1e25}
        ],
        "transverse_size": (1e-6, 1e-6),
        "molecule_location": "emitter_barrier",
        "notes": "CBO=0.61eV very close to Patil 0.6eV!"
    }
}

# ============================================================================
# MOLECULE LOCATION DEFINITIONS
# ============================================================================

MOLECULE_LOCATIONS = {
    "emitter_barrier": {
        "description": "On emitter barrier (most sensitive - from Patil 2018)",
        "layer_index": 1,  # Index in layers list (0-indexed)
        "position_fraction": 0.5,  # Middle of barrier
        "neighbor_radius": 2  # Lattice sites
    },
    "well_center": {
        "description": "In quantum well center",
        "layer_index": 2,
        "position_fraction": 0.5,
        "neighbor_radius": 1
    },
    "collector_barrier": {
        "description": "On collector barrier (less sensitive)",
        "layer_index": 3,
        "position_fraction": 0.5,
        "neighbor_radius": 2
    },
    "emitter_well_interface": {
        "description": "At emitter barrier / well interface",
        "layer_index": 1,
        "position_fraction": 1.0,  # End of barrier
        "neighbor_radius": 2
    }
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_device(name):
    """Get device configuration by name"""
    if name not in DEVICES:
        raise ValueError(f"Device '{name}' not found.\nAvailable: {list(DEVICES.keys())}")
    return DEVICES[name]

def get_material(name):
    """Get material properties by name"""
    if name not in MATERIALS:
        raise ValueError(f"Material '{name}' not found.\nAvailable: {list(MATERIALS.keys())}")
    return MATERIALS[name]

def list_all_devices():
    """List all available devices"""
    return list(DEVICES.keys())

def list_all_materials():
    """List all available materials"""
    return list(MATERIALS.keys())

def print_device_info(name):
    """Print detailed device information"""
    device = get_device(name)
    
    print(f"\n{'='*70}")
    print(f"DEVICE: {name}")
    print(f"{'='*70}")
    print(f"Description: {device['description']}")
    print(f"Transverse size: {device['transverse_size'][0]*1e6:.2f} µm × {device['transverse_size'][1]*1e6:.2f} µm")
    print(f"Molecule location: {device['molecule_location']}")
    print(f"\nLayer structure:")
    print(f"  {'Layer':<5} {'Material':<10} {'Thickness (nm)':<15} {'Doping (cm⁻³)':<15}")
    print(f"  {'-'*50}")
    
    total_thickness = 0
    for i, layer in enumerate(device['layers']):
        thick_nm = layer['thickness'] * 1e9
        doping_cm3 = layer['doping'] / 1e6  # Convert m^-3 to cm^-3
        total_thickness += thick_nm
        
        if doping_cm3 > 0:
            doping_str = f"{doping_cm3:.1e}"
        else:
            doping_str = "Undoped"
        
        mat = layer['material']
        mat_props = get_material(mat)
        print(f"  {i:<5} {mat:<10} {thick_nm:<15.2f} {doping_str:<15}")
    
    print(f"  {'-'*50}")
    print(f"  Total thickness: {total_thickness:.2f} nm")
    
    if 'notes' in device:
        print(f"\nNotes: {device['notes']}")
    print(f"{'='*70}\n")

def print_material_info(name):
    """Print material properties"""
    mat = get_material(name)
    
    print(f"\n{'='*70}")
    print(f"MATERIAL: {name}")
    print(f"{'='*70}")
    print(f"Description: {mat['description']}")
    print(f"Effective mass: {mat['m_eff']:.3f} m₀")
    print(f"Conduction band edge: {mat['Ec']:.3f} eV")
    print(f"Relative permittivity: {mat['epsilon_r']:.1f}")
    print(f"Phonon energy: {mat['phonon_energy_meV']:.1f} meV")
    print(f"{'='*70}\n")

def print_library_summary():
    """Print summary of device library"""
    print("="*70)
    print("DEVICE LIBRARY SUMMARY")
    print("="*70)
    print(f"Available devices: {len(DEVICES)}")
    print(f"Available materials: {len(MATERIALS)}")
    
    print("\nDevices:")
    for name, device in DEVICES.items():
        print(f"  {name:<30} - {device['description']}")
    
    print("\nMaterials:")
    for name, mat in MATERIALS.items():
        print(f"  {name:<15} - {mat['description']}")
    print("="*70)

def create_custom_device(layers, transverse_size=(1e-6, 1e-6), 
                        molecule_location="emitter_barrier", 
                        description="Custom device"):
    """
    Create a custom device configuration
    
    Parameters:
    -----------
    layers : list of dict
        Each dict has keys: 'material', 'thickness', 'doping'
    transverse_size : tuple
        (Ly, Lz) in meters
    molecule_location : str
        One of the predefined locations
    description : str
        Device description
        
    Returns:
    --------
    dict : Device configuration
    """
    # Validate layers
    for layer in layers:
        if 'material' not in layer or 'thickness' not in layer or 'doping' not in layer:
            raise ValueError("Each layer must have 'material', 'thickness', and 'doping'")
        if layer['material'] not in MATERIALS:
            raise ValueError(f"Unknown material: {layer['material']}")
    
    # Create device dict
    device = {
        "description": description,
        "layers": layers,
        "transverse_size": transverse_size,
        "molecule_location": molecule_location,
        "notes": "Custom user-defined device"
    }
    
    return device

# ============================================================================
# TEST CODE
# ============================================================================

if __name__ == "__main__":
    # Print library summary
    print_library_summary()
    
    # Print device details
    print_device_info("GaAs_AlAs_symmetric")
    print_device_info("GaAs_AlAs_asymmetric")
    
    # Print material info
    print_material_info("GaAs")
    print_material_info("AlAs")
    
    # Test custom device creation
    print("\n" + "="*70)
    print("EXAMPLE: Creating custom device")
    print("="*70)
    
    custom = create_custom_device(
        layers=[
            {"material": "GaAs", "thickness": 10e-9, "doping": 1e24},
            {"material": "AlAs", "thickness": 2.5e-9, "doping": 0},
            {"material": "GaAs", "thickness": 3e-9, "doping": 0},
            {"material": "AlAs", "thickness": 2.5e-9, "doping": 0},
            {"material": "GaAs", "thickness": 10e-9, "doping": 1e24}
        ],
        description="My custom RTD"
    )
    
    print(f"Custom device created: {custom['description']}")
    print(f"Total layers: {len(custom['layers'])}")
    print(f"Transverse size: {custom['transverse_size']}")
    print("="*70)
