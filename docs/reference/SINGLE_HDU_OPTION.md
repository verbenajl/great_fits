# Single HDU Option - Implementation Summary

## Overview
Added a `--single-hdu` flag to the `combine_fits` command that allows you to merge all spectral data from multiple FITS files into a single binary table HDU, instead of keeping them as separate extension HDUs.

## Usage

### Default behavior (Multiple HDUs):
```bash
combine_fits --input list.txt --output combined.fits
```
- Keeps all extension HDUs from all files as separate HDUs
- Result: Primary HDU + multiple extension HDUs

### Single HDU mode:
```bash
combine_fits --input list.txt --output combined.fits --single-hdu
```
- Merges all spectral data into a single binary table
- Result: Primary HDU + one large SPECTRA binary table with all combined data

## Key Differences

| Aspect | Multiple HDUs | Single HDU |
|--------|---------------|-----------|
| File structure | Primary HDU + N extension HDUs | Primary HDU + 1 SPECTRA HDU |
| Data organization | Organized by source file | All data merged together |
| File size | Slightly larger (headers for each HDU) | Slightly smaller |
| Use case | Keep original structure | Flatten structure for processing |

## Python API

You can also use the functions directly:

```python
from oi_zeigt.basic_io import combine_fits_from_list, combine_fits_files

# Multiple HDUs (default)
combined = combine_fits_from_list("list.txt", "combined.fits")

# Single HDU
combined = combine_fits_from_list("list.txt", "combined_single.fits", single_hdu=True)

# Or with the lower-level function
files = ["/path/to/file1.fits", "/path/to/file2.fits"]
combined = combine_fits_files(files, single_hdu=True)
combined.writeto("combined_single.fits", overwrite=True)
```

## Implementation Details

The `single_hdu` parameter was added to:
1. **`combine_fits_files()`** in `basic_io.py` - Core combining logic
2. **`combine_fits_from_list()`** in `basic_io.py` - File list reader
3. **`combine_fits()` CLI command** in `cli.py` - Command-line interface

When `single_hdu=True`:
- Reads all FITS files
- Extracts all extension HDUs (binary tables)
- Converts them to Astropy Table objects
- Concatenates all tables together
- Creates a single binary table HDU named "SPECTRA"
- Combines with the primary HDU

## Example Workflow

```bash
# Create list of files
echo "/path/to/file1.fits" > my_files.txt
echo "/path/to/file2.fits" >> my_files.txt
echo "/path/to/file3.fits" >> my_files.txt

# Combine with single HDU
combine_fits --input my_files.txt --output combined_single.fits --single-hdu

# Check result
print_oifits_info --fits combined_single.fits
```

The output will show only 2 HDUs:
- HDU 0: Primary HDU
- HDU 1: SPECTRA (binary table with all combined data)
