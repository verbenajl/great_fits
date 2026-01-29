from astropy.io import fits

hdul = fits.open('/home/diskB/data_soft/tests/SKYCHOPDIFF/cycle6.fits')
data = hdul[1].data

# Access the spectrum you saw in position "Row 42":
spectrum_index = 14595
my_spectrum = data[spectrum_index]['SPECTRUM']
my_object = data[spectrum_index]['OBJECT']

# Look at the values directly
print(my_spectrum)
print(f"Shape: {my_spectrum.shape}")
print(f"Min: {my_spectrum.min()}, Max: {my_spectrum.max()}")

# Find blank values
import numpy as np
blank_mask = my_spectrum == 0.0  # or whatever blank value you identify
print(f"Blank channels at indices: {np.where(blank_mask)[0]}")
print(f"Blank values: {my_spectrum[blank_mask]}")
