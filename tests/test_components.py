#!/usr/bin/env python
"""Check PCA component properties."""

import pickle
import numpy as np
from pathlib import Path

decomp_file = Path('output/pca_components/decomposition_2017-02-01_GR_F367_2017021_components.pkl')

with open(decomp_file, 'rb') as f:
    decomp = pickle.load(f)

if isinstance(decomp, dict):
    components = decomp.get('components')
    variance = decomp.get('explained_variance_ratio')
    print("Decomposition format: Dictionary")
else:
    components = decomp.components_
    variance = decomp.explained_variance_ratio_
    print("Decomposition format: sklearn PCA object")

print(f"\nComponent shape: {components.shape}")
print(f"Number of components: {len(components)}")
print(f"Spectral length: {components.shape[1]}")

print("\nComponent norms (L2):")
for i, comp in enumerate(components):
    norm = np.linalg.norm(comp)
    print(f"  Component {i}: norm={norm:.4f}, min={comp.min():.4f}, max={comp.max():.4f}, mean={comp.mean():.6f}")

print("\nExplained variance ratio:")
for i, var in enumerate(variance):
    print(f"  Component {i}: {var*100:.2f}%")

print(f"\nTotal variance: {sum(variance)*100:.2f}%")

# Check orthogonality
print("\nOrthogonality check (dot products between different components):")
for i in range(len(components)):
    for j in range(i+1, min(i+2, len(components))):  # Just check adjacent
        dot = np.dot(components[i], components[j])
        print(f"  C{i} · C{j} = {dot:.6f}")
