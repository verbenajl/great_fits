#!/usr/bin/env python3
"""
Investigate the reference channel for velocity axis.
Look for CRPIX-like keywords or other clues about which channel has the reference velocity.
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
print("INVESTIGATING REFERENCE CHANNEL FOR VELOCITY AXIS")
print("=" * 80)

# Extract spectral parameters
spectrum = table['SPECTRUM']
restfreq = table['RESTFREQ'][0]
lofreq = table['LOFREQ'][0]
deltav = table['DELTAV'][0]
velocity = table['VELOCITY'][0]
nchans = spectrum.shape[1]

print(f"\nBasic Parameters:")
print(f"  Number of channels: {nchans}")
print(f"  RESTFREQ: {restfreq:.10e} Hz")
print(f"  LOFREQ (first obs): {lofreq:.10e} Hz")
print(f"  DELTAV: {deltav:.2f} m/s ({deltav/1000:.4f} km/s)")
print(f"  VELOCITY: {velocity:.2f} m/s ({velocity/1000:.2f} km/s)")

# Look for all header keywords that might relate to reference pixel
print("\n" + "=" * 80)
print("HEADER KEYWORDS POTENTIALLY RELATED TO REFERENCE CHANNEL")
print("=" * 80)

reference_keywords = [
    # Standard WCS keywords
    'CRPIX1', 'CRVAL1', 'CDELT1',  # Usually for spectral axis
    'CRPIX2', 'CRVAL2', 'CDELT2',  # Usually for spatial
    'CRPIX3', 'CRVAL3', 'CDELT3',
    # Frequency related
    'CRFREQ', 'CRFREQ1',
    # Alternative names
    'FREQ-REF', 'FREQREF', 'FREQCEN', 'FREQ',
    # Channel/pixel related
    'CHAN', 'CHANNEL', 'PIXREF', 'REFCHAN',
    # Velocity reference
    'VREF', 'VELOREF', 'VELO', 'VELO-REF',
    # Frequency offset
    'FOFFSET', 'LOFREQ', 'IFFREQ',
]

print("\nFrom PRIMARY HDU header:")
if len(hdul) > 0 and hdul[0].header:
    primary_header = hdul[0].header
    for key in reference_keywords:
        if key in primary_header:
            val = primary_header[key]
            comment = primary_header.comments[key] if key in primary_header.comments else ''
            print(f"  {key:15s} = {str(val):30s} / {comment}")

print("\nFrom BINARY TABLE HDU header:")
for key in reference_keywords:
    if key in header:
        val = header[key]
        comment = header.comments[key] if key in header.comments else ''
        print(f"  {key:15s} = {str(val):30s} / {comment}")

# Check all header keywords (comprehensive)
print("\n" + "=" * 80)
print("ALL HEADER KEYWORDS IN BINARY TABLE")
print("=" * 80)
print("\nAll header cards in BINARY TABLE:")
for i, card in enumerate(header):
    comment = header.comments[card] if card in header.comments else ''
    val = str(header[card])
    # Only show if it's short enough and meaningful
    if len(val) < 50 and card not in ['COMMENT', 'HISTORY']:
        print(f"  {card:12s} = {val:30s} / {comment}")

# Look for LOFREQ and IMAGFREQ relationship
print("\n" + "=" * 80)
print("FREQUENCY ANALYSIS")
print("=" * 80)

print(f"\nRESTFREQ: {restfreq:.10e} Hz")
print(f"LOFREQ:   {lofreq:.10e} Hz")
print(f"Difference (REST - LO): {restfreq - lofreq:.10e} Hz")

if 'IMAGFREQ' in table.names:
    imagfreq = table['IMAGFREQ'][0]
    print(f"IMAGFREQ: {imagfreq:.10e} Hz")

# Analyze frequency offset
print(f"\nFOFFSET: {table['FOFFSET'][0]:.10e} Hz")

# Try to understand the spectral setup
print("\n" + "=" * 80)
print("SPECTRAL SETUP ANALYSIS")
print("=" * 80)

# Speed of light
c = 299792458.0  # m/s

# Frequency range from RESTFREQ and LOFREQ
freq_diff = restfreq - lofreq
print(f"\nFrequency offset (RESTFREQ - LOFREQ): {freq_diff:.10e} Hz")

# Convert to velocity
v_offset = (freq_diff / restfreq) * c
print(f"Velocity offset: {v_offset:.2f} m/s ({v_offset/1000:.2f} km/s)")

# Check against VELOCITY
print(f"\nVELOCITY column value: {velocity:.2f} m/s ({velocity/1000:.2f} km/s)")
print(f"Calculated velocity offset: {v_offset:.2f} m/s ({v_offset/1000:.2f} km/s)")
print(f"Difference: {abs(velocity - v_offset):.2f} m/s")

# Try to find reference channel
print("\n" + "=" * 80)
print("INFERRING REFERENCE CHANNEL")
print("=" * 80)

# The VELOCITY might be at the center of the band, or at a specific channel
# Let's check: if VELOCITY is the velocity of a specific channel
# and we know DELTAV, can we infer which channel?

# Hypothesis: VELOCITY might be at the center (channel = nchans // 2)
center_channel = nchans // 2
print(f"\nHypothesis 1: VELOCITY at channel center ({center_channel}):")
print(f"  If v(center) = {velocity:.2f} m/s")
print(f"  Then v(0) = {velocity - center_channel * deltav:.2f} m/s")
print(f"              = {(velocity - center_channel * deltav)/1000:.2f} km/s")

# Hypothesis: VELOCITY might be at the reference frequency channel
# If LOFREQ is the local oscillator and we're observing at RESTFREQ
# The frequency offset tells us which channel we're at
print(f"\nHypothesis 2: VELOCITY at frequency-matched channel:")
print(f"  Frequency offset corresponds to velocity offset")
print(f"  v_offset = {v_offset:.2f} m/s = {v_offset/1000:.2f} km/s")
print(f"  If VELOCITY = reference at that frequency...")
print(f"  Channel with v_offset: {v_offset / deltav:.1f} channels from reference")

# Hypothesis: VELOCITY could be at channel 0
print(f"\nHypothesis 3: VELOCITY at channel 0:")
print(f"  Then velocity span would be {velocity:.2f} to {velocity + nchans * deltav:.2f} m/s")
print(f"                             = {velocity/1000:.2f} to {(velocity + nchans * deltav)/1000:.2f} km/s")

# Hypothesis: VELOCITY could be at first valid channel (accounting for USB/LSB)
print(f"\nHypothesis 4: VELOCITY might relate to observing mode (USB/LSB):")
print(f"  Check LOFREQ variation across observations...")

# Look at LOFREQ distribution
lofreqs = table['LOFREQ']
print(f"\n  LOFREQ range:")
print(f"    Min: {np.min(lofreqs):.10e} Hz")
print(f"    Max: {np.max(lofreqs):.10e} Hz")
print(f"    Unique values: {len(np.unique(lofreqs))}")

# Check if there's a relationship between LOFREQ and channel coverage
# The frequency range covered by DELTAV * NCHANS
freq_span = deltav * nchans / c * restfreq
print(f"\n  Frequency span (from DELTAV × NCHANS):")
print(f"    {freq_span:.10e} Hz = {freq_span / 1e9:.3f} GHz")

# Try to find channel center based on RESTFREQ as center
print(f"\nHypothesis 5: RESTFREQ at CHANNEL CENTER")
center_ch = (nchans - 1) / 2.0
print(f"  Channel center: {center_ch:.1f}")
print(f"  If RESTFREQ is observed at center frequency:")
print(f"    Then VELOCITY should map to channel {center_ch:.1f}")
print(f"    Reference pixel (CRPIX): {center_ch + 1:.1f} (FITS convention)")
v_at_center = velocity
v_at_ch0 = v_at_center - center_ch * deltav
print(f"    v(center={center_ch:.0f}): {v_at_center:.2f} m/s = {v_at_center/1000:.2f} km/s")
print(f"    v(ch=0): {v_at_ch0:.2f} m/s = {v_at_ch0/1000:.2f} km/s")

hdul.close()
