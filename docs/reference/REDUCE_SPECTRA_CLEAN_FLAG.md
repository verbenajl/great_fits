# reduce_spectra --clean Flag Implementation

## Overview

Added a `--clean` flag to the `reduce_spectra` command that automatically reads the clean FITS file from the configuration file (`[output].clean_fits`) instead of the default input file.

## Usage

### Default behavior (reads from `[input].fits_file`):
```bash
reduce_spectra --config config.toml --baseline
```

### With `--clean` flag (reads from `[output].clean_fits`):
```bash
reduce_spectra --config config.toml --clean --baseline
```

### Override with explicit file (reads specified file):
```bash
reduce_spectra --fits my_clean_file.fits --baseline
```

## Priority Order

The command respects the following priority for input file selection:

1. **`--fits` flag** (highest priority) - If specified, uses this file
2. **`--clean` flag** - If specified, uses `[output].clean_fits` from config
3. **Default** (lowest priority) - Uses `[input].fits_file` from config

## Combined Usage Example

A typical workflow using the `--clean` flag:

```bash
# Step 1: Filter the data (creates clean_fits and rejected_fits)
filter_fits --config config.toml

# Step 2: Apply reduction to the clean data
reduce_spectra --config config.toml --clean --baseline --unblank

# Step 3: Average the reduced spectra
average --config config.toml --reduced
```

## Error Handling

If you use the `--clean` flag but the `[output].clean_fits` path is not defined in your config, you'll get a clear error message:

```
Error: --clean flag specified but [output].clean_fits not defined in config
```

## Configuration Requirement

Your `config.toml` must have the `clean_fits` path defined in the `[output]` section:

```toml
[output]
clean_fits = "/path/to/clean_data.fits"
rejected_fits = "/path/to/rejected_data.fits"
reduced_fits = "/path/to/reduced_data.fits"
averaged_fits = "/path/to/averaged_data.fits"
```

## Examples

### Example 1: Reduce only the clean M51 data
```bash
reduce_spectra --config config.toml --clean --baseline --baseline-order 2
```

### Example 2: Reduce clean data with multiple methods
```bash
reduce_spectra --config config.toml --clean --unblank --baseline --smooth
```

### Example 3: Reduce and specify output file
```bash
reduce_spectra --config config.toml --clean --baseline --output my_reduced.fits
```

## Benefits

✅ **Convenience**: No need to remember or specify the clean FITS file path
✅ **Pipeline integration**: Makes it easy to create automated reduction pipelines
✅ **Configuration-driven**: Respects your config.toml settings
✅ **Flexible**: Still allows manual file specification with `--fits` if needed
