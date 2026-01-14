"""
Test configuration modules thoroughly
"""
import sys
sys.path.append('.')

from config.molecular_database import get_molecule, list_all_molecules, get_class_molecules
from config.device_library import get_device, list_all_devices, get_material

print("=" * 70)
print("TESTING CONFIGURATION MODULES")
print("=" * 70)

# Test 1: Molecular database
print("\n[TEST 1] Molecular Database")
print("-" * 70)

try:
    # Test individual molecule
    benzene = get_molecule("Benzene")
    print(f"✓ Loaded Benzene: {benzene['formula']}")
    print(f"  - Vibrational modes: {len(benzene['modes_meV'])}")
    print(f"  - Coupling strengths: {benzene['coupling_meV']}")

    # Test all molecules
    all_mols = list_all_molecules()
    print(f"✓ Total molecules: {len(all_mols)}")

    # Test by class
    aromatic = get_class_molecules("Aromatic")
    print(f"✓ Aromatic class: {len(aromatic)} molecules")
    for mol in aromatic:
        print(f"    - {mol}")

    # Test invalid molecule
    try:
        invalid = get_molecule("InvalidMolecule")
        print("✗ Should have raised error for invalid molecule")
    except ValueError as e:
        print(f"✓ Correctly caught invalid molecule: {e}")

except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback
    traceback.print_exc()

# Test 2: Device library
print("\n[TEST 2] Device Library")
print("-" * 70)

try:
    # Test individual device
    device = get_device("In2O3_Al2O3_symmetric")
    print(f"✓ Loaded device: In2O3_Al2O3_symmetric")
    print(f"  - Description: {device['description']}")
    print(f"  - Number of layers: {len(device['layers'])}")
    print(f"  - Transverse size: {device['transverse_size']}")

    # Test all devices
    all_devs = list_all_devices()
    print(f"✓ Total devices: {len(all_devs)}")

    # Test material properties
    in2o3 = get_material("In2O3")
    print(f"✓ In2O3 properties:")
    print(f"    - Effective mass: {in2o3['m_eff']} m0")
    print(f"    - Band offset: {in2o3['Ec']} eV")

    # Test invalid device
    try:
        invalid = get_device("InvalidDevice")
        print("✗ Should have raised error for invalid device")
    except ValueError as e:
        print(f"✓ Correctly caught invalid device: {e}")

except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback
    traceback.print_exc()

# Test 3: Data consistency
print("\n[TEST 3] Data Consistency Checks")
print("-" * 70)

try:
    all_molecules = list_all_molecules()

    # Check each molecule has required fields
    for mol_name in all_molecules:
        mol = get_molecule(mol_name)
        assert 'formula' in mol, f"{mol_name} missing formula"
        assert 'modes_meV' in mol, f"{mol_name} missing modes_meV"
        assert 'coupling_meV' in mol, f"{mol_name} missing coupling_meV"
        assert len(mol['modes_meV']) == len(mol['coupling_meV']), \
            f"{mol_name} mode/coupling mismatch"

    print(f"✓ All {len(all_molecules)} molecules have consistent data")

    # Check each device has required fields
    all_devices = list_all_devices()
    for dev_name in all_devices:
        dev = get_device(dev_name)
        assert 'description' in dev, f"{dev_name} missing description"
        assert 'layers' in dev, f"{dev_name} missing layers"
        assert 'transverse_size' in dev, f"{dev_name} missing transverse_size"

        # Check all layer materials exist
        for layer in dev['layers']:
            get_material(layer['material'])

    print(f"✓ All {len(all_devices)} devices have consistent data")

except Exception as e:
    print(f"✗ FAILED: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 70)
print("CONFIGURATION TESTS COMPLETE")
print("=" * 70)
