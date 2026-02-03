#!/bin/bash
# Activate the cygrid environment (with cygrid for gridding, no OpenCV)
# Usage: source activate_cygrid.sh
source ~/.venv/cygrid/bin/activate
echo "✓ Activated ~/.venv/cygrid (Python 3.12, with cygrid)"
echo "  Available: cygrid for gridding (e.g., via map_integrated)"
echo "  ⚠ pca_correct NOT available in this environment (requires OpenCV)"
