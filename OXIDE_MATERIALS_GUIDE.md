# Oxide RTD Material Systems - Quick Reference

## Your Target Systems (BEOL-compatible <400°C)

### 1. **In₂O₃/Al₂O₃** (PRIMARY CANDIDATE) ⭐
- **Why**: Best BEOL evidence (Intel/Purdue VLSI 2021)
- **Process**: ALD <350°C
- **Effective mass**: 0.30 m₀ (In₂O₃), 0.45 m₀ (Al₂O₃)
- **Band offset**: 2.8 eV
- **Doping**: 10¹⁹ cm⁻³
- **Device**: `In2O3_Al2O3_symmetric`, `In2O3_Al2O3_asymmetric`
- **Status**: Strong literature support, mature ALD process

### 2. **IGZO/Al₂O₃** (BACKUP)
- **Why**: Amorphous = excellent uniformity
- **Process**: Sputtering <300°C
- **Effective mass**: 0.34 m₀ (IGZO)
- **Band offset**: 2.8 eV (assuming similar to In₂O₃)
- **Doping**: 8×10¹⁸ cm⁻³
- **Device**: `IGZO_Al2O3_symmetric`
- **Status**: TFT technology mature, RTD application novel

### 3. **ZnO/Al₂O₃** (BACKUP)
- **Why**: Mature ALD process, good mobility
- **Process**: ALD <350°C
- **Effective mass**: 0.28 m₀ (ZnO)
- **Band offset**: 2.8 eV
- **Doping**: 10¹⁹ cm⁻³
- **Device**: `ZnO_Al2O3_symmetric`
- **Status**: Widely studied, very mature

### 4. **SnO₂/Al₂O₃** (BACKUP)
- **Why**: Low-temp sputtering, good stability
- **Process**: Sputtering/ALD <300°C
- **Effective mass**: 0.25 m₀ (SnO₂)
- **Band offset**: 2.8 eV
- **Doping**: 10¹⁹ cm⁻³
- **Device**: `SnO2_Al2O3_symmetric`
- **Status**: Lower mobility, but very stable

---

## Material Property Summary

| Material | m*/m₀ | Ec (eV) | εᵣ | ℏω (meV) | Process T |
|----------|-------|---------|-----|----------|-----------|
| **In₂O₃** | 0.30 | 0.0 | 9.0 | 70 | 350°C |
| **SnO₂**  | 0.25 | 0.0 | 9.0 | 75 | 300°C |
| **IGZO**  | 0.34 | 0.0 | 10.0 | 60 | 300°C |
| **ZnO**   | 0.28 | 0.0 | 8.5 | 72 | 350°C |
| **Al₂O₃** | 0.45 | 2.8 | 9.0 | 100 | 300°C |
| **HfO₂**  | 0.50 | 2.5 | 25.0 | 95 | 350°C |

---

## Usage in Code

```python
from config.device_library import get_device, print_device_info

# Get your primary device
device = get_device("In2O3_Al2O3_symmetric")

# Print details
print_device_info("In2O3_Al2O3_symmetric")

# Access properties
print(device['layers'])
print(device['transverse_size'])
print(device['molecule_location'])
```

---

## Key Advantages of Oxide RTDs

✅ **BEOL-compatible**: All processes <400°C  
✅ **CMOS integration**: Can be fabricated on top of Si logic  
✅ **No III-V**: Simpler, cheaper, more manufacturable  
✅ **Mature processes**: ALD/sputtering well-established  
✅ **Scalable**: Standard semiconductor manufacturing  

---

## References for Material Parameters

1. **In₂O₃**: Presley et al. JAP 2004 (m*), Walsh et al. PRB 2009 (phonon)
2. **SnO₂**: Godinho et al. JPC 2009
3. **IGZO**: Nomura et al. Nature 2004
4. **ZnO**: Look et al. SSC 1998
5. **Al₂O₃**: Robertson & Wallace, MSE R 2015 (band offsets)

---

## Simulation Priority

1. **First**: Validate with GaAs/AlAs (known system)
2. **Second**: Run In₂O₃/Al₂O₃ (primary candidate)
3. **Third**: Compare all 4 oxide systems
4. **Fourth**: Optimize barrier thickness/well width for each

---

## Expected Differences from GaAs/AlAs

| Property | GaAs/AlAs | Oxide RTDs | Impact |
|----------|-----------|------------|--------|
| Effective mass | 0.067 | 0.25-0.34 | **Lower tunneling current** |
| Band offset | 0.57 eV | 2.8 eV | **Better confinement** |
| Process temp | 500-800°C | <400°C | **BEOL-compatible** |
| Phonon energy | 36 meV | 60-100 meV | **Different IETS spectrum** |
| Mobility | Very high | Moderate | **Lower peak current** |

**Net result**: Oxide RTDs will have:
- Lower overall current (but still measurable)
- Better quantum confinement (sharper peaks)
- Different background phonon spectrum
- MORE MANUFACTURABLE! 🎉
