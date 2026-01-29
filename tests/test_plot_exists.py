#!/usr/bin/env python
"""Analyze one of the generated plots to see if correction is visible."""

import matplotlib.image as mpimg
import numpy as np
from pathlib import Path

plot_file = Path('/home/diskB/sofiaobsdata/m51/m51_central_fits/pca_plots/pca_correction_UNKNOWN_17562_LFAH_PX00_S.png')

if plot_file.exists():
    img = mpimg.imread(plot_file)
    print(f"Plot loaded: {plot_file.name}")
    print(f"Image shape: {img.shape}")
    print(f"Image dtype: {img.dtype}")
    print(f"Pixel value range: min={img.min():.4f}, max={img.max():.4f}, mean={img.mean():.4f}")
    print(f"\nPlot file size: {plot_file.stat().st_size / 1024:.1f} KB")
    print("✓ Plot exists and has content")
else:
    print(f"Plot not found: {plot_file}")
