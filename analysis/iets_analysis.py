"""
Analysis Tools

Analyze IETS spectra for selectivity, clustering, and classification.
"""

import numpy as np
import sys
import os
import csv

# Setup path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

from config.molecular_database import get_molecule, PERCEPTUAL_CLASSES

# ============================================================================
# IETS FINGERPRINT ANALYSIS
# ============================================================================

def extract_iets_fingerprint(V_array, d2IdV2, n_features=10):
    """
    Extract IETS fingerprint (peak positions and intensities)
    
    Parameters:
    -----------
    V_array : array
        Voltage points
    d2IdV2 : array
        IETS spectrum
    n_features : int
        Number of features to extract
        
    Returns:
    --------
    fingerprint : dict
        Peak positions, intensities, and feature vector
    """
    
    # Normalize
    d2IdV2_norm = d2IdV2 / (np.abs(d2IdV2).max() + 1e-12)
    
    # Find peaks (simple: local maxima above threshold)
    threshold = 0.1
    peaks = []
    
    for i in range(1, len(d2IdV2_norm)-1):
        if d2IdV2_norm[i] > threshold:
            if d2IdV2_norm[i] > d2IdV2_norm[i-1] and d2IdV2_norm[i] > d2IdV2_norm[i+1]:
                peaks.append({
                    'voltage': V_array[i],
                    'intensity': d2IdV2_norm[i],
                    'index': i
                })
    
    # Sort by intensity
    peaks = sorted(peaks, key=lambda p: p['intensity'], reverse=True)
    
    # Take top n_features
    peaks = peaks[:n_features]
    
    # Create feature vector (fixed size)
    feature_vector = np.zeros(n_features * 2)  # positions + intensities
    
    for i, peak in enumerate(peaks):
        feature_vector[i] = peak['voltage']
        feature_vector[n_features + i] = peak['intensity']
    
    return {
        'peaks': peaks,
        'n_peaks': len(peaks),
        'feature_vector': feature_vector
    }

# ============================================================================
# SIMILARITY METRICS
# ============================================================================

def compute_similarity(fp1, fp2, method='euclidean'):
    """
    Compute similarity between two IETS fingerprints
    
    Parameters:
    -----------
    fp1, fp2 : dict
        IETS fingerprints from extract_iets_fingerprint
    method : str
        'euclidean', 'correlation', 'cosine'
        
    Returns:
    --------
    similarity : float
        Similarity score (higher = more similar)
    distance : float
        Distance metric (lower = more similar)
    """
    
    v1 = fp1['feature_vector']
    v2 = fp2['feature_vector']
    
    if method == 'euclidean':
        distance = np.linalg.norm(v1 - v2)
        similarity = 1.0 / (1.0 + distance)
        
    elif method == 'correlation':
        if len(v1) > 0 and len(v2) > 0:
            correlation = np.corrcoef(v1, v2)[0, 1]
            similarity = (correlation + 1.0) / 2.0  # Scale to [0, 1]
            distance = 1.0 - correlation
        else:
            similarity = 0.0
            distance = 2.0
            
    elif method == 'cosine':
        norm1 = np.linalg.norm(v1)
        norm2 = np.linalg.norm(v2)
        if norm1 > 0 and norm2 > 0:
            similarity = np.dot(v1, v2) / (norm1 * norm2)
            distance = 1.0 - similarity
        else:
            similarity = 0.0
            distance = 1.0
    
    else:
        raise ValueError(f"Unknown method: {method}")
    
    return similarity, distance

# ============================================================================
# SELECTIVITY MATRIX
# ============================================================================

def build_selectivity_matrix(batch_results, method='euclidean'):
    """
    Build selectivity matrix showing pairwise similarities
    
    Parameters:
    -----------
    batch_results : dict
        Results from batch screening
    method : str
        Similarity method
        
    Returns:
    --------
    matrix : dict
        Selectivity matrix and metadata
    """
    
    # Extract successful simulations
    sims = [s for s in batch_results['simulations'] if 'error' not in s]
    
    if len(sims) == 0:
        print("No successful simulations to analyze!")
        return None
    
    n_molecules = len(sims)
    molecule_names = [s['molecule'] for s in sims]
    
    # Extract fingerprints
    fingerprints = []
    for sim in sims:
        V = sim['results']['V_array']
        d2IdV2 = sim['results']['d2IdV2']
        fp = extract_iets_fingerprint(V, d2IdV2)
        fingerprints.append(fp)
    
    # Compute pairwise similarities
    similarity_matrix = np.zeros((n_molecules, n_molecules))
    distance_matrix = np.zeros((n_molecules, n_molecules))
    
    for i in range(n_molecules):
        for j in range(n_molecules):
            sim, dist = compute_similarity(
                fingerprints[i], fingerprints[j], method
            )
            similarity_matrix[i, j] = sim
            distance_matrix[i, j] = dist
    
    # Compute selectivity score (off-diagonal mean)
    mask = ~np.eye(n_molecules, dtype=bool)
    mean_similarity = similarity_matrix[mask].mean()
    selectivity = 1.0 - mean_similarity
    
    return {
        'molecule_names': molecule_names,
        'fingerprints': fingerprints,
        'similarity_matrix': similarity_matrix,
        'distance_matrix': distance_matrix,
        'selectivity': selectivity,
        'method': method
    }

# ============================================================================
# CLASSIFICATION ANALYSIS
# ============================================================================

def analyze_classification(selectivity_matrix):
    """
    Analyze classification accuracy by perceptual class
    
    Parameters:
    -----------
    selectivity_matrix : dict
        From build_selectivity_matrix
        
    Returns:
    --------
    classification_report : dict
        Accuracy metrics
    """
    
    molecules = selectivity_matrix['molecule_names']
    sim_matrix = selectivity_matrix['similarity_matrix']
    
    # Get true classes
    true_classes = []
    for mol_name in molecules:
        mol = get_molecule(mol_name)
        true_classes.append(mol['perceptual_class'])
    
    # Predict class for each molecule (nearest neighbor excluding self)
    predictions = []
    
    for i, mol_name in enumerate(molecules):
        # Get similarities to all others
        similarities = sim_matrix[i, :].copy()
        similarities[i] = -1  # Exclude self
        
        # Find most similar
        most_similar_idx = similarities.argmax()
        predicted_class = true_classes[most_similar_idx]
        predictions.append(predicted_class)
    
    # Compute accuracy
    correct = sum(1 for t, p in zip(true_classes, predictions) if t == p)
    accuracy = correct / len(molecules)
    
    # Per-class accuracy
    class_accuracy = {}
    for class_name in set(true_classes):
        class_indices = [i for i, c in enumerate(true_classes) if c == class_name]
        class_correct = sum(1 for i in class_indices if predictions[i] == class_name)
        class_accuracy[class_name] = class_correct / len(class_indices) if class_indices else 0.0
    
    return {
        'overall_accuracy': accuracy,
        'class_accuracy': class_accuracy,
        'true_classes': true_classes,
        'predictions': predictions,
        'correct': correct,
        'total': len(molecules)
    }

# ============================================================================
# SAVE ANALYSIS
# ============================================================================

def save_selectivity_matrix(selectivity_matrix, output_dir, filename="selectivity_matrix.csv"):
    """Save selectivity matrix to CSV"""
    
    filepath = os.path.join(output_dir, filename)
    
    molecules = selectivity_matrix['molecule_names']
    sim_matrix = selectivity_matrix['similarity_matrix']
    
    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Header
        writer.writerow([''] + molecules)
        
        # Data
        for i, mol in enumerate(molecules):
            row = [mol] + [f'{sim_matrix[i,j]:.4f}' for j in range(len(molecules))]
            writer.writerow(row)
        
        # Summary
        writer.writerow([])
        writer.writerow(['Overall selectivity:', f"{selectivity_matrix['selectivity']:.4f}"])
    
    print(f"Selectivity matrix saved to: {filepath}")
    return filepath

def save_classification_report(classification_report, selectivity_matrix, 
                               output_dir, filename="classification_report.csv"):
    """Save classification report to CSV"""
    
    filepath = os.path.join(output_dir, filename)
    
    with open(filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Overall metrics
        writer.writerow(['Overall Accuracy:', f"{classification_report['overall_accuracy']:.4f}"])
        writer.writerow(['Correct:', classification_report['correct']])
        writer.writerow(['Total:', classification_report['total']])
        writer.writerow([])
        
        # Per-class accuracy
        writer.writerow(['Per-Class Accuracy'])
        writer.writerow(['Class', 'Accuracy'])
        for class_name, acc in classification_report['class_accuracy'].items():
            writer.writerow([class_name, f'{acc:.4f}'])
        writer.writerow([])
        
        # Predictions
        writer.writerow(['Molecule', 'True Class', 'Predicted Class', 'Correct'])
        for mol, true_c, pred_c in zip(
            selectivity_matrix['molecule_names'],
            classification_report['true_classes'],
            classification_report['predictions']
        ):
            correct = '✓' if true_c == pred_c else '✗'
            writer.writerow([mol, true_c, pred_c, correct])
    
    print(f"Classification report saved to: {filepath}")
    return filepath

# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    print("\n" + "="*70)
    print("IETS ANALYSIS TOOLS - DEMO")
    print("="*70)
    
    # Demo: Synthetic data
    print("\n[DEMO] Creating synthetic IETS spectra...")
    
    # Create 3 synthetic molecules with different peak patterns
    V_test = np.linspace(0, 0.5, 100)
    
    # Molecule A: peaks at 0.1, 0.3 V
    d2IdV2_A = np.exp(-((V_test - 0.1)/0.02)**2) + 0.5*np.exp(-((V_test - 0.3)/0.02)**2)
    
    # Molecule B: peaks at 0.15, 0.35 V (similar to A)
    d2IdV2_B = np.exp(-((V_test - 0.15)/0.02)**2) + 0.5*np.exp(-((V_test - 0.35)/0.02)**2)
    
    # Molecule C: peaks at 0.2, 0.4 V (different)
    d2IdV2_C = np.exp(-((V_test - 0.2)/0.02)**2) + 0.8*np.exp(-((V_test - 0.4)/0.02)**2)
    
    # Extract fingerprints
    print("\n[DEMO] Extracting fingerprints...")
    fp_A = extract_iets_fingerprint(V_test, d2IdV2_A)
    fp_B = extract_iets_fingerprint(V_test, d2IdV2_B)
    fp_C = extract_iets_fingerprint(V_test, d2IdV2_C)
    
    print(f"  Molecule A: {fp_A['n_peaks']} peaks")
    print(f"  Molecule B: {fp_B['n_peaks']} peaks")
    print(f"  Molecule C: {fp_C['n_peaks']} peaks")
    
    # Compute similarities
    print("\n[DEMO] Computing similarities...")
    sim_AB, dist_AB = compute_similarity(fp_A, fp_B, method='euclidean')
    sim_AC, dist_AC = compute_similarity(fp_A, fp_C, method='euclidean')
    sim_BC, dist_BC = compute_similarity(fp_B, fp_C, method='euclidean')
    
    print(f"  A-B similarity: {sim_AB:.3f} (distance: {dist_AB:.3f})")
    print(f"  A-C similarity: {sim_AC:.3f} (distance: {dist_AC:.3f})")
    print(f"  B-C similarity: {sim_BC:.3f} (distance: {dist_BC:.3f})")
    print(f"  → A and B are more similar (as expected)")
    
    print("\n" + "="*70)
    print("DEMO COMPLETE ✓")
    print("="*70 + "\n")
