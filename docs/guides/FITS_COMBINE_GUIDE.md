# FITS File Combination Guide

## Overview

The `combine_fits` command allows you to combine multiple FITS files into a single FITS file. This is useful when you have data split across multiple files that need to be merged together.

## Usage

```bash
combine_fits --input list.txt --output combined.fits
```

### Parameters

- `--input`: (Required) Path to a text file containing a list of FITS file paths, one per line
- `--output`: (Required) Path where the combined FITS file will be written

## Input File Format

The input file is a simple text file where each line contains the full path to a FITS file:

```
# lines starting with # are comments (optional)
/path/to/file1.fits
/path/to/file2.fits
/path/to/file3.fits

# You can have blank lines
/path/to/file4.fits
```

### Features of the input file:
- One FITS file path per line
- Empty lines are automatically ignored
- Lines starting with `#` are treated as comments and ignored
- Absolute paths are recommended for clarity

## Example Usage

### Example 1: Simple combination

```bash
combine_fits --input my_files.txt --output combined.fits
```

### Example 2: With absolute paths

Create `my_files.txt`:
```
/home/diskB/data_soft/tests/SKYCHOPDIFF/cycle6_part1.fits
/home/diskB/data_soft/tests/SKYCHOPDIFF/cycle6_part2.fits
/home/diskB/data_soft/tests/SKYCHOPDIFF/cycle6_part3.fits
```

Then run:
```bash
combine_fits --input my_files.txt --output /home/diskB/data_soft/tests/SKYCHOPDIFF/combined.fits
```

## How It Works

1. The command reads the input list file
2. It opens each FITS file in sequence
3. It combines all HDUs (Header/Data Units) into a single HDUList:
   - The primary HDU from the first file is preserved
   - All extension HDUs from all files are appended in order
4. The combined HDUList is written to the output file

## Python API

You can also use the combining functionality directly in Python:

```python
from oi_zeigt.basic_io import combine_fits_from_list

# Combine FITS files from a list file
combined_hdul = combine_fits_from_list("list.txt", "combined.fits")
combined_hdul.info()
```

Or if you prefer to work with a list directly:

```python
from oi_zeigt.basic_io import combine_fits_files

files = [
    "/path/to/file1.fits",
    "/path/to/file2.fits",
    "/path/to/file3.fits"
]

combined_hdul = combine_fits_files(files)
combined_hdul.writeto("combined.fits", overwrite=True)
```

## Error Handling

The command will report errors if:
- The input list file doesn't exist
- Any FITS file in the list doesn't exist
- The input list is empty (no valid file paths)
- There are permission issues writing the output file

## Sample Input File

See `sample_fits_list.txt` for a template that you can copy and modify.
