#!/usr/bin/env python3
"""
Visualize the velocity axis structure with CRPIX1 reference pixel.
This demonstrates how velocity maps to channel index.
"""

import numpy as np
import matplotlib.pyplot as plt

# M51 Spectral Parameters
CRVAL3 = 470000  # m/s (reference velocity)
CRPIX3 = 506     # FITS 1-indexed (≈506)
CDELT3 = 500     # m/s per channel
NAXIS3 = 1264    # total channels

# Convert CRPIX3 to 0-indexed
ref_channel_0indexed = CRPIX3 - 1  # ≈505

# Calculate velocity for each channel
channels = np.arange(NAXIS3)
velocities = CRVAL3 + (channels - ref_channel_0indexed) * CDELT3

# Convert to km/s for plotting
velocities_kms = velocities / 1000

# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

# Plot 1: Full velocity axis
ax1.plot(channels, velocities_kms, 'b-', linewidth=2, label='Velocity axis')
ax1.axvline(ref_channel_0indexed, color='r', linestyle='--', linewidth=2, label=f'Reference channel {ref_channel_0indexed:.0f}')
ax1.axhline(CRVAL3/1000, color='r', linestyle='--', linewidth=1, alpha=0.5)
ax1.scatter([ref_channel_0indexed], [CRVAL3/1000], color='r', s=100, zorder=5, label=f'Reference point (v={CRVAL3/1000:.0f} km/s)')

# Mark key channels
key_channels = [0, 500, 505, 1000, 1263]
for ch in key_channels:
    if ch < NAXIS3:
        v = CRVAL3 + (ch - ref_channel_0indexed) * CDELT3
        ax1.scatter([ch], [v/1000], color='green', s=50, zorder=4)
        ax1.annotate(f'Ch {ch}\n{v/1000:.0f} km/s', xy=(ch, v/1000), 
                    xytext=(10, 10), textcoords='offset points', fontsize=8,
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.5))

ax1.set_xlabel('Channel Index (0-based)', fontsize=12)
ax1.set_ylabel('Velocity (km/s)', fontsize=12)
ax1.set_title('M51 Spectral Axis: Velocity vs Channel Index\n(CRPIX1=506, CDELT3=500 m/s)', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=10, loc='upper left')
ax1.set_xlim(-50, 1313)

# Plot 2: Zoomed in around reference channel
zoom_range = 50
zoom_start = max(0, ref_channel_0indexed - zoom_range)
zoom_end = min(NAXIS3, ref_channel_0indexed + zoom_range)

zoom_channels = channels[zoom_start:zoom_end+1]
zoom_velocities = velocities_kms[zoom_start:zoom_end+1]

ax2.plot(zoom_channels, zoom_velocities, 'b-', linewidth=2, marker='o', markersize=4)
ax2.axvline(ref_channel_0indexed, color='r', linestyle='--', linewidth=2, label=f'Reference channel {ref_channel_0indexed:.0f}')
ax2.axhline(CRVAL3/1000, color='r', linestyle='--', linewidth=1, alpha=0.5)
ax2.scatter([ref_channel_0indexed], [CRVAL3/1000], color='r', s=200, marker='*', zorder=5, 
           label=f'CRVAL3 = {CRVAL3/1000:.0f} km/s')

# Annotate every 5 channels in zoom view
for ch in zoom_channels[::5]:
    v = CRVAL3 + (ch - ref_channel_0indexed) * CDELT3
    ax2.annotate(f'{v/1000:.1f}', xy=(ch, v/1000), 
                xytext=(0, 5), textcoords='offset points', fontsize=8, ha='center')

ax2.set_xlabel('Channel Index (0-based)', fontsize=12)
ax2.set_ylabel('Velocity (km/s)', fontsize=12)
ax2.set_title(f'Zoomed: Reference Region (Channels {zoom_start}-{zoom_end})', fontsize=12, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.legend(fontsize=10)

plt.tight_layout()
plt.savefig('/home/verbena/software/oi_zeigt/velocity_axis_visualization.png', dpi=150, bbox_inches='tight')
print("✓ Saved visualization to velocity_axis_visualization.png")

# Print detailed channel information
print("\n" + "="*70)
print("VELOCITY AXIS STRUCTURE - M51 DATASET")
print("="*70)
print(f"\nWCS Parameters:")
print(f"  CRVAL3 (reference velocity):     {CRVAL3:>10} m/s = {CRVAL3/1000:>7.1f} km/s")
print(f"  CRPIX3 (reference pixel):        {CRPIX3:>10.0f} (FITS 1-indexed)")
print(f"  CRPIX3 (0-indexed):              {ref_channel_0indexed:>10.1f}")
print(f"  CDELT3 (velocity step):          {CDELT3:>10} m/s = {CDELT3/1000:>7.3f} km/s/ch")
print(f"  NAXIS3 (total channels):         {NAXIS3:>10}")

print(f"\nKey Channel Velocities:")
print(f"  {'Channel':<10} {'Velocity (m/s)':<20} {'Velocity (km/s)':<20}")
print(f"  {'-'*10} {'-'*20} {'-'*20}")

for ch in [0, 100, 200, 300, 400, 500, 505, 506, 600, 700, 800, 900, 1000, 1100, 1200, 1263]:
    v = CRVAL3 + (ch - ref_channel_0indexed) * CDELT3
    marker = " ← REFERENCE" if ch == int(ref_channel_0indexed) else ""
    print(f"  {ch:<10} {v:<20.0f} {v/1000:<20.2f}{marker}")

print(f"\nVelocity Range:")
v_first = CRVAL3 + (0 - ref_channel_0indexed) * CDELT3
v_last = CRVAL3 + (NAXIS3-1 - ref_channel_0indexed) * CDELT3
print(f"  First channel (0):       {v_first/1000:>7.2f} km/s")
print(f"  Last channel ({NAXIS3-1}):      {v_last/1000:>7.2f} km/s")
print(f"  Total span:              {(v_last - v_first)/1000:>7.2f} km/s")

print("\n" + "="*70)
