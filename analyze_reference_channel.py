#!/usr/bin/env python3
"""
Infer the reference channel from frequency information.

The relationship between frequency and velocity is:
    v = c * (f_rest - f_obs) / f_rest   [optical definition]
or
    v = c * (f_obs - f_rest) / f_rest   [radio definition]

Where:
    c = speed of light
    f_rest = rest frequency
    f_obs = observed frequency
    v = radial velocity

We can use LOFREQ + frequency offset to get the observed frequency at each channel,
then calculate the velocity.
"""

from pathlib import Path
from astropy.io import fits
import numpy as np

fits_file = Path("/home/diskB/sofiaobsdata/m51/m51_central_region/kalibrated_all/m51_central_tile_singlehdu.fits")

if not fits_file.exists():
    print(f"Error: FITS file not found: {fits_file}")
    exit(1)

c = 299792458.0  # Speed of light in m/s

with fits.open(fits_file) as hdul:
    table = hdul[1].data
    header = hdul[1].header
    
    print("=" * 80)
    print("SPECTRAL AXIS ANALYSIS: FREQUENCY TO VELOCITY CONVERSION")
    print("=" * 80)
    
    # Get first observation
    restfreq = table['RESTFREQ'][0]
    lofreq = table['LOFREQ'][0]
    deltav = table['DELTAV'][0]
    velocity = table['VELOCITY'][0]
    foffset = table['FOFFSET'][0]
    spectrum = table['SPECTRUM'][0]
    
    nchans = len(spectrum)
    
    print(f"\nFirst observation parameters:")
    print(f"  RESTFREQ: {restfreq:.10e} Hz ({restfreq/1e12:.6f} THz)")
    print(f"  LOFREQ:   {lofreq:.10e} Hz")
    print(f"  FOFFSET:  {foffset:.10e} Hz")
    print(f"  DELTAV:   {deltav:.2f} m/s")
    print(f"  VELOCITY: {velocity:.2f} m/s ({velocity/1000:.2f} km/s)")
    print(f"  SPECTRUM shape: {spectrum.shape}")
    print(f"  Number of channels: {nchans}")
    
    # Calculate observed frequency (LO + offset)
    # In a typical receiver setup, the LO frequency is the local oscillator
    # The observed frequency at channel i would be related to the channel offset
    
    print(f"\n" + "=" * 80)
    print("FREQUENCY-VELOCITY RELATIONSHIP")
    print("=" * 80)
    
    # The VELOCITY column likely indicates the velocity at channel 0
    # But let's verify by checking if DELTAV maps correctly
    
    print(f"\nAssuming VELOCITY is the reference velocity at some reference channel:")
    print(f"  Reference velocity: {velocity:.2f} m/s = {velocity/1000:.2f} km/s")
    print(f"  Velocity step: {deltav:.2f} m/s = {deltav/1000:.4f} km/s per channel")
    
    # Calculate velocities for all channels
    # If VELOCITY is at channel i_ref, then:
    # v(i) = VELOCITY + (i - i_ref) * DELTAV
    
    # Let's assume VELOCITY is at channel 0 (standard assumption)
    print(f"\nAssuming VELOCITY applies to channel 0:")
    velocities_ch0 = velocity + np.arange(nchans) * deltav
    print(f"  Channel 0: {velocities_ch0[0]:.2f} m/s = {velocities_ch0[0]/1000:.2f} km/s")
    print(f"  Channel {nchans//2}: {velocities_ch0[nchans//2]:.2f} m/s = {velocities_ch0[nchans//2]/1000:.2f} km/s")
    print(f"  Channel {nchans-1}: {velocities_ch0[nchans-1]:.2f} m/s = {velocities_ch0[nchans-1]/1000:.2f} km/s")
    
    # But the user said 470 km/s is not at channel 0
    # Let's see if there's a pattern with frequency
    
    print(f"\n" + "=" * 80)
    print("ALTERNATIVE: VELOCITY AT MIDDLE CHANNEL?")
    print("=" * 80)
    
    # Try assuming VELOCITY is at the middle channel
    mid_chan = nchans // 2
    velocities_mid = velocity + (np.arange(nchans) - mid_chan) * deltav
    print(f"\nAssuming VELOCITY applies to channel {mid_chan} (middle):")
    print(f"  Channel 0: {velocities_mid[0]:.2f} m/s = {velocities_mid[0]/1000:.2f} km/s")
    print(f"  Channel {mid_chan}: {velocities_mid[mid_chan]:.2f} m/s = {velocities_mid[mid_chan]/1000:.2f} km/s")
    print(f"  Channel {nchans-1}: {velocities_mid[nchans-1]:.2f} m/s = {velocities_mid[nchans-1]/1000:.2f} km/s")
    
    # Try assuming VELOCITY is at channel nchans//2 or somewhere else
    print(f"\n" + "=" * 80)
    print("TESTING DIFFERENT REFERENCE CHANNELS")
    print("=" * 80)
    
    # Test a range of possible reference channels
    test_channels = [0, nchans//4, nchans//2, 3*nchans//4, nchans-1]
    
    for ref_ch in test_channels:
        v_at_0 = velocity - ref_ch * deltav
        v_at_last = velocity + (nchans - 1 - ref_ch) * deltav
        print(f"\nIf VELOCITY={velocity:.0f} m/s is at channel {ref_ch}:")
        print(f"  Channel 0: {v_at_0:.2f} m/s = {v_at_0/1000:.2f} km/s")
        print(f"  Channel {ref_ch}: {velocity:.2f} m/s = {velocity/1000:.2f} km/s")
        print(f"  Channel {nchans-1}: {v_at_last:.2f} m/s = {v_at_last/1000:.2f} km/s")
    
    print(f"\n" + "=" * 80)
    print("CHECKING: WHICH CHANNEL HAS VELOCITY=470 km/s?")
    print("=" * 80)
    
    # The user observed that 470 km/s must be at some channel, not channel 0
    # Let's find which channel would give 470 km/s if it's the velocity range
    target_v = 470000.0  # 470 km/s in m/s
    
    # If v(i) = velocity + (i - ref_ch) * deltav = target_v
    # Then: ref_ch = (velocity - target_v) / deltav + i
    
    # Assuming velocity is the ref velocity and we want to find channel where v = 470 km/s
    for test_ref in [0, nchans//2, nchans//4]:
        for test_i in range(0, nchans, 100):
            v_at_i = velocity + (test_i - test_ref) * deltav
            if abs(v_at_i - target_v) < 1000:  # Within 1 km/s
                print(f"Channel {test_i} ≈ 470 km/s if reference is channel {test_ref}")
