"""
Molecular Database from Pandey et al. Scientific Reports 2021
"Vibration-based biomimetic odor classification"

Contains vibrational data for 20 odorant molecules across 6 perceptual classes.
"""

import numpy as np

def cm1_to_meV(cm1):
    """Convert wavenumber (cm^-1) to energy (meV)"""
    return cm1 * 0.123984

# ============================================================================
# MOLECULAR DATABASE
# ============================================================================

MOLECULES = {
    # ========== BASELINE (NO MOLECULE) ==========
    "Baseline": {
        "class": "None",
        "subclass": None,
        "modes_meV": [],
        "modes_cm1": [],
        "coupling_meV": [],
        "formula": "None",
        "description": "No odorant molecule (baseline simulation)"
    },
    
    # ========== AROMATIC CLASS ==========
    "Benzene": {
        "class": "Aromatic",
        "subclass": "Strong",
        "modes_cm1": [399, 637, 1084, 1485, 3189],
        "modes_meV": [49.5, 79.0, 134.4, 184.1, 395.4],
        "coupling_meV": [5.0, 8.0, 10.0, 8.0, 5.0],
        "formula": "C6H6",
        "description": "Strong aromatic odor"
    },
    
    "Anthracene": {
        "class": "Aromatic",
        "subclass": "Weak",
        "modes_cm1": [615, 826, 1147, 1446, 3206],
        "modes_meV": [76.3, 102.4, 142.2, 179.3, 397.5],
        "coupling_meV": [5.0, 8.0, 10.0, 8.0, 5.0],
        "formula": "C14H10",
        "description": "Weak aromatic odor"
    },
    
    "Thiophene": {
        "class": "Aromatic",
        "subclass": "Weak",
        "modes_cm1": [538, 758, 1043, 1424, 3258],
        "modes_meV": [66.7, 94.0, 129.3, 176.6, 403.9],
        "coupling_meV": [5.0, 8.0, 10.0, 8.0, 5.0],
        "formula": "C4H4S",
        "description": "Weak aromatic odor with sulfur"
    },
    
    # ========== ROASTED COFFEE CLASS ==========
    "Furan": {
        "class": "Roasted Coffee",
        "subclass": None,
        "modes_cm1": [633, 816, 1011, 1498, 3311],
        "modes_meV": [78.5, 101.2, 125.4, 185.7, 410.6],
        "coupling_meV": [5.0, 8.0, 10.0, 8.0, 5.0],
        "formula": "C4H4O",
        "description": "Roasted coffee aroma"
    },
    
    "Furan_Methanethiol": {
        "class": "Roasted Coffee",
        "subclass": None,
        "modes_cm1": [157, 684, 1025, 1526, 2650, 3251],
        "modes_meV": [19.5, 84.8, 127.1, 189.2, 328.6, 403.1],
        "coupling_meV": [3.0, 5.0, 8.0, 10.0, 7.0, 5.0],
        "formula": "C5H6OS",
        "description": "Roasted coffee with meaty note"
    },
    
    # ========== MOTH-BALL CLASS ==========
    "Naphthalene": {
        "class": "Moth-ball",
        "subclass": None,
        "modes_cm1": [151, 441, 842, 1150, 1515, 3213],
        "modes_meV": [18.7, 54.7, 104.4, 142.6, 187.8, 398.4],
        "coupling_meV": [3.0, 5.0, 8.0, 10.0, 8.0, 5.0],
        "formula": "C10H8",
        "description": "Classic moth-ball odor"
    },
    
    "Tetralin": {
        "class": "Moth-ball",
        "subclass": None,
        "modes_cm1": [118, 468, 827, 1147, 1471, 3080],
        "modes_meV": [14.6, 58.0, 102.5, 142.2, 182.4, 381.9],
        "coupling_meV": [3.0, 5.0, 8.0, 10.0, 8.0, 5.0],
        "formula": "C10H12",
        "description": "Moth-ball with petroleum note"
    },
    
    "Fluorene": {
        "class": "Moth-ball",
        "subclass": None,
        "modes_cm1": [137, 463, 871, 1126, 1456, 3193],
        "modes_meV": [17.0, 57.4, 108.0, 139.6, 180.5, 395.9],
        "coupling_meV": [3.0, 5.0, 8.0, 10.0, 8.0, 5.0],
        "formula": "C13H10",
        "description": "Moth-ball odor"
    },
    
    # ========== FRUITY CLASS ==========
    "Heptanal": {
        "class": "Fruity",
        "subclass": None,
        "modes_cm1": [79, 883, 1287, 1817, 3046],
        "modes_meV": [9.8, 109.5, 159.6, 225.3, 377.7],
        "coupling_meV": [3.0, 8.0, 10.0, 8.0, 5.0],
        "formula": "C7H14O",
        "description": "Citrus fruity odor"
    },
    
    "N-Amyl_Butyrate": {
        "class": "Fruity",
        "subclass": None,
        "modes_cm1": [118, 846, 1157, 1386, 1780, 3066],
        "modes_meV": [14.6, 104.9, 143.5, 171.9, 220.7, 380.2],
        "coupling_meV": [3.0, 5.0, 8.0, 10.0, 8.0, 5.0],
        "formula": "C9H18O2",
        "description": "Pear/apple fruity odor"
    },
    
    "Gamma_Octalactone": {
        "class": "Fruity",
        "subclass": None,
        "modes_cm1": [112, 820, 1263, 1816, 3066],
        "modes_meV": [13.9, 101.7, 156.6, 225.2, 380.2],
        "coupling_meV": [3.0, 8.0, 10.0, 8.0, 5.0],
        "formula": "C8H14O2",
        "description": "Coconut fruity odor"
    },
    
    # ========== MUSK CLASS ==========
    "Civetone": {
        "class": "Musk",
        "subclass": None,
        "modes_cm1": [240, 1044, 1377, 3009],
        "modes_meV": [29.8, 129.5, 170.7, 373.1],
        "coupling_meV": [5.0, 8.0, 10.0, 5.0],
        "formula": "C17H30O",
        "description": "Animal musk odor"
    },
    
    "Moxalone": {
        "class": "Musk",
        "subclass": None,
        "modes_cm1": [292, 993, 1370, 3029],
        "modes_meV": [36.2, 123.1, 169.9, 375.6],
        "coupling_meV": [5.0, 8.0, 10.0, 5.0],
        "formula": "C12H16",
        "description": "Synthetic musk"
    },
    
    "Galaxolide": {
        "class": "Musk",
        "subclass": None,
        "modes_cm1": [318, 989, 1373, 3035],
        "modes_meV": [39.4, 122.6, 170.2, 376.3],
        "coupling_meV": [5.0, 8.0, 10.0, 5.0],
        "formula": "C18H26O",
        "description": "Polycyclic musk"
    },
    
    "Helvetolide": {
        "class": "Musk",
        "subclass": None,
        "modes_cm1": [280, 999, 1359, 3029],
        "modes_meV": [34.7, 123.9, 168.5, 375.6],
        "coupling_meV": [5.0, 8.0, 10.0, 5.0],
        "formula": "C17H26O",
        "description": "Macrocyclic musk"
    },
    
    # ========== GARLICKY CLASS ==========
    "Benzyl_Mercaptan": {
        "class": "Garlicky",
        "subclass": "Artificial",
        "modes_cm1": [260, 881, 1091, 1454, 2585, 3193],
        "modes_meV": [32.2, 109.2, 135.3, 180.3, 320.5, 395.9],
        "coupling_meV": [5.0, 8.0, 10.0, 8.0, 7.0, 5.0],
        "formula": "C7H8S",
        "description": "Artificial sulfurous-garlic"
    },
    
    "Allyl_Thiol": {
        "class": "Garlicky",
        "subclass": "Artificial",
        "modes_cm1": [272, 927, 1313, 1696, 2602, 3086],
        "modes_meV": [33.7, 114.9, 162.8, 210.3, 322.6, 382.7],
        "coupling_meV": [5.0, 8.0, 10.0, 8.0, 7.0, 5.0],
        "formula": "C3H6S",
        "description": "Artificial garlic odor"
    },
    
    "Dimethyl_Sulfide": {
        "class": "Garlicky",
        "subclass": "Artificial",
        "modes_cm1": [226, 885, 1371, 3051],
        "modes_meV": [28.0, 109.7, 170.0, 378.3],
        "coupling_meV": [5.0, 8.0, 10.0, 5.0],
        "formula": "C2H6S",
        "description": "Sulfurous odor"
    },
    
    "Diallyl_Disulfide": {
        "class": "Garlicky",
        "subclass": "Natural",
        "modes_cm1": [110, 439, 922, 1254, 1724, 3070],
        "modes_meV": [13.6, 54.4, 114.3, 155.5, 213.8, 380.7],
        "coupling_meV": [3.0, 5.0, 8.0, 10.0, 8.0, 5.0],
        "formula": "C6H10S2",
        "description": "Natural garlic (from Allium)"
    },
    
    "Allicin": {
        "class": "Garlicky",
        "subclass": "Natural",
        "modes_cm1": [120, 382, 958, 1313, 1692, 3109],
        "modes_meV": [14.9, 47.4, 118.8, 162.8, 209.8, 385.5],
        "coupling_meV": [3.0, 5.0, 8.0, 10.0, 8.0, 5.0],
        "formula": "C6H10OS2",
        "description": "Natural garlic (main component)"
    }
}

# ============================================================================
# PERCEPTUAL CLASS ORGANIZATION
# ============================================================================

PERCEPTUAL_CLASSES = {
    "Aromatic": ["Benzene", "Anthracene", "Thiophene"],
    "Roasted Coffee": ["Furan", "Furan_Methanethiol"],
    "Moth-ball": ["Naphthalene", "Tetralin", "Fluorene"],
    "Fruity": ["Heptanal", "N-Amyl_Butyrate", "Gamma_Octalactone"],
    "Musk": ["Civetone", "Moxalone", "Galaxolide", "Helvetolide"],
    "Garlicky": ["Benzyl_Mercaptan", "Allyl_Thiol", "Dimethyl_Sulfide", 
                 "Diallyl_Disulfide", "Allicin"]
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_molecule(name):
    """Get molecule data by name"""
    if name not in MOLECULES:
        raise ValueError(f"Molecule '{name}' not found in database.\nAvailable: {list(MOLECULES.keys())}")
    return MOLECULES[name]

def get_class_molecules(perceptual_class):
    """Get all molecules in a perceptual class"""
    if perceptual_class not in PERCEPTUAL_CLASSES:
        raise ValueError(f"Class '{perceptual_class}' not found.\nAvailable: {list(PERCEPTUAL_CLASSES.keys())}")
    return PERCEPTUAL_CLASSES[perceptual_class]

def list_all_molecules():
    """List all available molecules"""
    return list(MOLECULES.keys())

def list_all_classes():
    """List all perceptual classes"""
    return list(PERCEPTUAL_CLASSES.keys())

def print_database_summary():
    """Print summary of molecular database"""
    print("="*70)
    print("MOLECULAR DATABASE SUMMARY (Pandey et al. 2021)")
    print("="*70)
    print(f"Total molecules: {len(MOLECULES)}")
    print(f"Perceptual classes: {len(PERCEPTUAL_CLASSES)}")
    print("\nMolecules per class:")
    for pclass, mols in PERCEPTUAL_CLASSES.items():
        print(f"  {pclass:20s}: {len(mols)} molecules")
    
    print("\nVibrational energy ranges:")
    all_modes = []
    for mol_data in MOLECULES.values():
        all_modes.extend(mol_data["modes_meV"])
    
    if all_modes:
        print(f"  Min: {min(all_modes):.1f} meV")
        print(f"  Max: {max(all_modes):.1f} meV")
        print(f"  Mean: {np.mean(all_modes):.1f} meV")
    print("="*70)

def print_molecule_info(name):
    """Print detailed info for a specific molecule"""
    mol = get_molecule(name)
    print(f"\n{'='*70}")
    print(f"MOLECULE: {name}")
    print(f"{'='*70}")
    print(f"Formula: {mol['formula']}")
    print(f"Class: {mol['class']}", end="")
    if mol['subclass']:
        print(f" ({mol['subclass']})")
    else:
        print()
    print(f"Description: {mol['description']}")
    
    if mol['modes_meV']:
        print(f"\nVibrational modes:")
        print(f"  {'Mode':<6} {'cm⁻¹':<8} {'meV':<8} {'Coupling (meV)':<15}")
        print(f"  {'-'*40}")
        for i, (cm1, mev, coup) in enumerate(zip(
            mol['modes_cm1'], mol['modes_meV'], mol['coupling_meV']
        ), 1):
            print(f"  {i:<6} {cm1:<8.0f} {mev:<8.1f} {coup:<15.1f}")
    else:
        print("\n(No vibrational modes - baseline)")
    print(f"{'='*70}\n")

# ============================================================================
# TEST CODE
# ============================================================================

if __name__ == "__main__":
    # Print database summary
    print_database_summary()
    
    # Example: Get Benzene data
    print("\n" + "="*70)
    print("EXAMPLE: Getting Benzene data")
    print("="*70)
    benzene = get_molecule("Benzene")
    print(f"Modes (meV): {benzene['modes_meV']}")
    print(f"Coupling (meV): {benzene['coupling_meV']}")
    
    # Print detailed info
    print_molecule_info("Benzene")
    print_molecule_info("Furan")
    print_molecule_info("Baseline")
    
    # List all aromatics
    print("\n" + "="*70)
    print("All Aromatic molecules:")
    aromatics = get_class_molecules("Aromatic")
    for mol in aromatics:
        print(f"  - {mol}")
    print("="*70)
