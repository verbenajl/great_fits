#!/usr/bin/env python3
"""
Detailed inspection of spectral axis parameters from table columns.
"""

from pathlib import Path
from astropy.io import fits
import numpy as np
import sys

# Path to the input FITS file
fits_file = Path("/home/diskB/sofiaobsdata/m51/m51_central_region/kalibrated_all/m51_central_tile_singlehdu.fits")

if not fits_file.exists():
    print(f"Error: FITS file not found: {fits_file}")
    sys.exit(1)

# Open the FITS file
hdul = fits.open(fits_file)
table = hdul[1].data
header = hdul[1].header

print("=" * 80)
print("SPECTRAL AXIS PARAMETERS FROM TABLE COLUMNS")
print("=" * 80)

# Key spectral columns
print("\n1. RESTFREQ (Rest Frequency)")
if 'RESTFREQ' in table.names:
    restfreq = table['RESTFREQ']
    print(f"   Shape: {restfreq.shape}")
    print(f"   Data type: {restfreq.dtype}")
    print(f"   Unique values: {len(np.unique(restfreq))}")
    print(f"   Min: {np.min(restfreq):.10e} Hz")
    print(f"   Max: {np.max(restfreq):.10e} Hz")
    print(f"   Mean: {np.mean(restfreq):.10e} Hz")
    # Print first few values
    print(f"   First 5 values: {restfreq[:5]}")

print("\n2. LOFREQ (Line-of-sight / Local Oscillator Frequency)")
if 'LOFREQ' in table.names:
    lofreq = table['LOFREQ']
    print(f"   Shape: {lofreq.shape}")
    print(f"   Data type: {lofreq.dtype}")
    print(f"   Unique values: {len(np.unique(lofreq))}")
    print(f"   Min: {np.min(lofreq):.10e} Hz")
    print(f"   Max: {np.max(lofreq):.10e} Hz")
    print(f"   Mean: {np.mean(lofreq):.10e} Hz")
    print(f"   First 5 values: {lofreq[:5]}")

print("\n3. DELTAV (Delta Velocity - velocity channel spacing)")
if 'DELTAV' in table.names:
    deltav = table['DELTAV']
    print(f"   Shape: {deltav.shape}")
    print(f"   Data type: {deltav.dtype}")
    print(f"   Unique values: {len(np.unique(deltav))}")
    print(f"   Min: {np.min(deltav):.6f} m/s")
    print(f"   Max: {np.max(deltav):.6f} m/s")
    print(f"   Mean: {np.mean(deltav):.6f} m/s")
    print(f"   First 5 values: {deltav[:5]}")

print("\n4. FOFFSET (Frequency Offset)")
if 'FOFFSET' in table.names:
    foffset = table['FOFFSET']
    print(f"   Shape: {foffset.shape}")
    print(f"   Data type: {foffset.dtype}")
    print(f"   Unique values: {len(np.unique(foffset))}")
    print(f"   Min: {np.min(foffset):.10e} Hz")
    print(f"   Max: {np.max(foffset):.10e} Hz")
    print(f"   Mean: {np.mean(foffset):.10e} Hz")
    print(f"   First 5 values: {foffset[:5]}")

print("\n5. VELDEF (Velocity Definition)")
if 'VELDEF' in table.names:
    veldef = table['VELDEF']
    print(f"   Shape: {veldef.shape}")
    print(f"   Data type: {veldef.dtype}")
    print(f"   Unique values: {np.unique(veldef)}")

print("\n6. VELOCITY (Velocity)")
if 'VELOCITY' in table.names:
    velocity = table['VELOCITY']
    print(f"   Shape: {velocity.shape}")
    print(f"   Data type: {velocity.dtype}")
    print(f"   Unique values: {len(np.unique(velocity))}")
    print(f"   Min: {np.min(velocity):.6f} m/s")
    print(f"   Max: {np.max(velocity):.6f} m/s")
    print(f"   Mean: {np.mean(velocity):.6f} m/s")
    print(f"   First 5 values: {velocity[:5]}")

# Check SPECTRUM shape and relate to spectral channels
print("\n" + "=" * 80)
print("SPECTRUM ARRAY ANALYSIS")
print("=" * 80)
spectrum = table['SPECTRUM']
nobs, nchans = spectrum.shape
print(f"\nSPECTRUM shape: {nobs} observations × {nchans} channels")
print(f"This means we have {nchans} spectral channels")

# Check if DELTAV * nchans makes sense
if 'DELTAV' in table.names and 'RESTFREQ' in table.names:
    deltav_sample = table['DELTAV'][0]
    restfreq_sample = table['RESTFREQ'][0]
    print(f"\nFrom first observation:")
    print(f"  DELTAV: {deltav_sample:.6f} m/s (velocity per channel)")
    print(f"  RESTFREQ: {restfreq_sample:.10e} Hz")
    print(f"  Total velocity span: {deltav_sample * nchans:.6f} m/s = {deltav_sample * nchans / 1000:.2f} km/s")
    
    # Speed of light
    c = 299792458.0  # m/s
    print(f"\n  This corresponds to a frequency span of: {(deltav_sample * nchans / c) * restfreq_sample:.10e} Hz")

# Check LOFREQ vs RESTFREQ
if 'LOFREQ' in table.names and 'RESTFREQ' in table.names:
    lofreq = table['LOFREQ']
    restfreq = table['RESTFREQ']
    foffset = table['FOFFSET'] if 'FOFFSET' in table.names else None
    
    print("\n" + "=" * 80)
    print("FREQUENCY RELATIONSHIP")
    print("=" * 80)
    print(f"\nRESTFREQ - LOFREQ (rest - LO):")
    diff = restfreq - lofreq
    print(f"  Min: {np.min(diff):.10e} Hz")
    print(f"  Max: {np.max(diff):.10e} Hz")
    print(f"  Mean: {np.mean(diff):.10e} Hz")
    
    if foffset is not None:
        print(f"\nFOFFSET vs (RESTFREQ - LOFREQ):")
        print(f"  Are they equal? {np.allclose(foffset, diff)}")

# Check header keywords related to spectral axis
print("\n" + "=" * 80)
print("BINARY TABLE HEADER SPECTRAL KEYWORDS")
print("=" * 80)
spectral_keys = ['DELTAV', 'FOFFSET', 'RESTFREQ', 'LOFREQ', 'VELDEF', 'IMAGFREQ']
for key in spectral_keys:
    if key in header:
        comment = header.comments[key] if key in header.comments else ''
        print(f"{key:15s} = {str(header[key]):40s} / {comment}")

hdul.close()
print("\nDone!")
