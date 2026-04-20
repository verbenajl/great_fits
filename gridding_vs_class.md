# Gridding comparison: oi_zeigt vs GILDAS-CLASS

## Background

The collapsed moment-0 maps produced by `oi_zeigt collapse_cube` look noisier than
equivalent maps produced by GILDAS-CLASS `xy_map`.  This document records the
root-cause analysis so we can revisit and align the two approaches.

Source investigated: `~/gildas-src-apr25a/packages/class/lib/map/xymap.f90`
and `convolve.f90`.

---

## Key parameter comparison

| Parameter         | CLASS default                          | oi_zeigt default          |
|-------------------|----------------------------------------|---------------------------|
| Pixel size        | `reso / 2  ≈  beam / 2`  (Nyquist)    | `beam / 3`                |
| Kernel FWHM       | `beam / 3`                             | `beam`  (3× wider)        |
| Kernel support    | `3 × FWHM  =  beam`                    | `3σ  ≈  1.27 × beam`      |
| Kernel function   | `exp(-(r/σ)²)` (Gaussian via `ctype_exponential`, parm(3)=2) | `exp(-(r/σ)²)` (Gaussian via cygrid) |
| LUT               | Tabulated at 1/100th pixel for speed   | HEALPix-based (cygrid)    |

For a 14.1″ beam these become:

| Parameter      | CLASS     | oi_zeigt  |
|----------------|-----------|-----------|
| Pixel size     | 7.1″      | 4.7″      |
| Kernel FWHM    | 4.7″      | 14.1″     |
| Kernel σ       | 2.8″      | 6.0″      |
| Support radius | 14.1″     | 18.0″     |

---

## Where this is set in the CLASS source

**`xymap.f90`, lines 940 / 953** — pixel size default:
```fortran
map%cell(1) = -map%reso(1)/2.0   ! Nyquist sampling: pixel = resolution/2
```

**`xymap.f90`, lines 973–986** — kernel FWHM default:
```fortran
! Default: beam/3.0
fwhm(1) = map%beam/3.0
map%reso(1) = sqrt(map%beam**2+fwhm(1)**2)   ! effective reso ≈ 1.054*beam
```

**`xymap.f90`, lines 1002–1009** — kernel parameters passed to convolution:
```fortran
map%conv%x%parm(1) = map%support(1)/abs(map%cell(1))          ! support in pixels
map%conv%x%parm(2) = fwhm(1)/(2*sqrt(log(2.0)))/abs(map%cell(1))  ! sigma in pixels
map%conv%x%parm(3) = 2                                         ! Gaussian exponent
map%conv%x%ctype   = ctype_exponential  ! exp(-(r/sigma)^parm3)
```

**`convolve.f90`, lines 160–162** — the LUT lookup (1/100th pixel tabulation):
```fortran
ix   = nint(100.0*dx + conv%x%bias)
iy   = nint(100.0*dy + conv%y%bias)
resu = conv%x%buff(ix) * conv%y%buff(iy)
```
The LUT is a performance optimisation, not a conceptual difference.

---

## Why CLASS looks smoother

With pixel = beam/2, CLASS maps are **Nyquist-sampled** — no pixel is finer
than the resolution element warrants.  Our beam/3 pixels oversample the beam,
showing noise structure at scales smaller than the beam.

Additionally, CLASS's narrower kernel (FWHM = beam/3) places each spectrum
close to its true sky position with minimal spreading, while our wider kernel
(FWHM = beam) blurs neighbouring spectra together, potentially spreading
artefacts and residual baseline structure over a wider area.

In practice the two effects partially cancel, but the **pixel size is the
dominant visual difference**: more pixels means more noise instances visible
on screen.

---

## How to reproduce CLASS-like output today

Use the existing `--kernel-fwhm` and `--pixel-size-arcsec` CLI options to
match CLASS defaults:

```bash
oi-zeigt create_datacube \
    --fits reduced.fits \
    --beamsize 14.1 \
    --pixel-size-arcsec 7.05 \   # beam/2
    --kernel-fwhm 4.7 \          # beam/3
    --output cube_class_style.fits \
    --config config.toml
```

Then collapse normally:
```bash
oi-zeigt collapse_cube cube_class_style.fits \
    --snr-threshold 3 \
    -o moment0_class_style.fits
```

---

## Possible code changes to consider

1. **Change pixel default from beam/3 to beam/2.**
   Matches CLASS convention.  The cube will have ~55% fewer pixels per axis,
   which also reduces memory and I/O.  Set `pixel_size_arcsec` to `beamsize/2`
   in `get_gridding_params_from_config()`.

2. **Change kernel FWHM default from beam to beam/3.**
   Matches CLASS convention.  Produces sharper maps at native resolution.
   Update `create_datacube` so that when `kernel_fwhm_arcsec` is None the
   default passed to cygrid is `beamsize / 3` rather than `beamsize`.

3. **Add `--effective-beam` output flag.**
   When a non-zero kernel FWHM is used, write the correct effective resolution
   `sqrt(beam² + kernel_fwhm²)` into `BMAJ`/`BMIN` of the output FITS header.
   Currently `BMAJ = beamsize` regardless of kernel choice.

4. **Expose `--pixel-size-arcsec` as a first-class CLI option.**
   It is already wired internally but check it reaches `create_datacube_cmd`
   without ambiguity (the `--pixsize` alias exists but is not prominently
   documented).

---

## Notes

- The CLASS `LUT` (Look-Up Table) is only a speed trick: kernel values are
  pre-computed at 1/100th pixel intervals to avoid calling `exp()` in the
  inner loop.  cygrid uses a HEALPix-based approach for the same purpose.
  This is not the source of any visual difference.
- CLASS effective resolution after default gridding: `sqrt(beam² + (beam/3)²)
  ≈ 1.054 × beam` (a 5.4% broadening, noted in the source comments).
- The `--smooth` and `--snr-threshold` options in `collapse_cube` are useful
  post-gridding workarounds regardless of which kernel convention is used.
