# PCA Reduction Parameters

The pipeline splits into two steps — **decompose** (extract eigenspectra from the sky-reference data) and **correct** (fit and subtract components from science spectra). Parameters are grouped accordingly in the config file under `[common]`, `[decompose]`, and `[correct]`.

---

## Common parameters (affect both steps)

| Parameter | Default | Effect |
|---|---|---|
| `pca_source` | `SKYCHOPDIFF` | Source name of the spectra used to derive components. Use `SKYCHOPDIFF` (or `SKY-DIFF` for CLASS-pipeline data). Science spectra are always the correction target regardless of this setting. |
| `smoothing_kernel_size` | `3` | Boxcar kernel width (channels) applied to each component after decomposition. Smoothing suppresses high-frequency noise in the components. Set `false` to disable. Larger → smoother components, less sensitivity to narrow spectral features. |
| `rolling_noise_window` | `11` | Window for a rolling RMS estimate of each component, used to compute the *global structure* diagnostic score displayed in the component plots. Has **no effect** on the correction itself. |
| `decomposition_type` | `PCA` | Algorithm: `PCA`, `SparsePCA`, or `ICA`. PCA is the only well-tested option. |
| `do_scale` | `false` | Scale all input spectra by dividing by their standard deviation before decomposition. Rarely needed; default is fine. |
| `z_score_cutoff` | `false` | Integer. Before fitting, blank any channel whose absolute z-score exceeds this value (spike rejection). Use `10` as a conservative starting point if spike contamination is visible. |
| `line_window` | `false` | List of velocity pairs (km/s) that explicitly mark the emission line, e.g. `[450, 550]`. These channels are excluded from the component-fitting step so the line is not partially subtracted. Use when automatic line detection (`line_kernel_size`) misses the line or is unreliable. |
| `pca_exclude_range` | `false` | Velocity range physically removed from the spectrum before the entire PCA step. Channels in this range are filled back in (unchanged) after correction. Use for a known contaminated band that must not influence the decomposition or the correction fitting. |
| `tag` | auto | String label written to a tag file after decompose; the correct step reads it back to find the right pickle. Set explicitly in Snakemake/batch runs to guarantee matching. |

---

## Decompose parameters

| Parameter | Default | Effect |
|---|---|---|
| `number_components` | `5` | Number of eigenspectra to extract. More components → more sky structure captured, but also increasing risk of absorbing faint astrophysical signal. Typical range: 3–8. The code derives `number_components - 1` PCA components plus one mean spectrum. |
| `noise_cutoff` | `false` | Integer. When non-sky spectra are mixed into the decomposition input, drop any spectrum whose emission window exceeds `noise_cutoff × rms` in more than `channel_line_width` channels. SKY-DIFF/SKYCHOPDIFF spectra are always kept. **No effect in the new fits-based pipeline** (see `pca_param_todo.txt`). |
| `scramble` | `false` | Validation tool: randomises input spectra so components should be pure noise. Currently dead code in both pipelines. |

---

## Correct parameters

These control whether each PCA component is subtracted from each science spectrum.

### Line detection

| Parameter | Default | Effect |
|---|---|---|
| `line_kernel_size` | `51` (legacy) / `61` (M51 config) | Size of the Gaussian blur kernel (channels) used by the waterfall-image line detector. Larger → detects broader/fainter emission; smaller → more aggressive masking. If automatic detection flags too much baseline as line, reduce this value. |

### Component-level filtering (global)

| Parameter | Default | Effect |
|---|---|---|
| `variance_cutoff` | `false` | Explained-variance-ratio threshold. Components whose `explained_variance_ratio_ ≤ variance_cutoff` are never applied to any spectrum. The ratio is between 0 and 1; components are sorted highest to lowest, so the first component always has the largest ratio. **Higher cutoff = stricter** — only the strongest components pass. **Lower cutoff = more permissive** — only the very weakest components are dropped. Typical values if used: `0.001`–`0.01`. Set `false` to use all components regardless of explained variance. |
| `global_noise_ratio_cutoff` | `15` (M51) | Integer/float. A **global per-component** gate based on the noise ratio distribution. See below for what noise ratio means. The code fits a KDE to the noise ratio values of each component across all spectra and finds the center of the first (lowest) peak. If `center ≥ global_noise_ratio_cutoff`, the component is **rejected for all spectra** — its sky signal is too weak relative to the noise across the board. **Higher value = more permissive** (weak components still pass). **Lower value = stricter** (only components with a strong signal-to-noise pass). Set `false` or `0` to disable. Good starting value: `15`–`20`. |

### Per-spectrum filtering

The noise ratio for a given spectrum and component is:

```
noise_ratio = std(smoothed_spectrum, baseline channels) / std(coeff × component)
```

The denominator is the amplitude of the scaled component in this spectrum (how large the sky pattern actually is here). The numerator is the noise level of the spectrum.

- **High noise ratio** → the component is *small* relative to the noise. Subtracting it would add noise rather than remove a real sky feature. **Reject.**
- **Low noise ratio** → the component is *large* relative to the noise — a coherent sky signal worth removing. **Apply.**

The rejection condition is `noise_ratio > cutoff_value`, so **higher cutoff = more permissive** (weaker components still get applied), **lower cutoff = stricter**.

| Parameter | Default | Effect |
|---|---|---|
| `noise_ratio_cutoff` | `false` | Per-spectrum threshold. Spectra where `noise_ratio > noise_ratio_cutoff` are skipped for that component. If set to `false`/`0` while `global_noise_ratio_cutoff` is active, the cutoff is derived adaptively per component as `KDE_center + 3 × KDE_sigma` — recommended. If set to a fixed number, that scalar is used for all components equally, which is less flexible but useful when the KDE fit is unreliable (few spectra, multi-modal distributions). |

### Experimental

| Parameter | Default | Effect |
|---|---|---|
| `cut_coefficients` | `false` | Float (antenna-temperature units). The coefficient is the dot product of the component with the spectrum — it measures how much of that component's *shape* is present in this particular spectrum. Rejection condition: `abs(coeff) < cut_coefficients`. **Higher value = stricter** — only subtract a component if it has a strong projection onto this spectrum. **Lower value = more permissive.** This catches the case where a component passes the noise-ratio test globally but barely projects onto one particular spectrum, meaning the sky pattern is simply absent there. Typical value if enabled: `0.1` K. Keep disabled (`false` or `0`) unless noise-ratio filtering alone leaves residuals in spectra where a component has near-zero amplitude. |

---

## Accept / reject decision flow

```
For each component j and each science spectrum i:

1. variance_cutoff check (global)
   → if explained_variance_ratio[j] <= variance_cutoff: REJECT j for all spectra

2. global_noise_ratio_cutoff check (global, KDE peak)
   → if KDE_center[j] >= global_noise_ratio_cutoff: REJECT j for all spectra

3. per-spectrum noise_ratio check
   → if noise_ratio[i,j] > cutoff_value: REJECT j for spectrum i
      (cutoff_value = noise_ratio_cutoff if set, else adaptive KDE center + 3σ)

4. cut_coefficients check
   → if |coeff[i,j]| < cut_coefficients: REJECT j for spectrum i

If none of the above reject the component → APPLY: subtract coeff[j] × component[j] from spectrum i
```

---

## Practical starting points

- **Too much subtracted / line partly removed**: increase `line_kernel_size`, add an explicit `line_window`, or raise `global_noise_ratio_cutoff`.
- **Too little subtracted / baseline stripes persist**: lower `global_noise_ratio_cutoff`, reduce `number_components` (if weak components are noise), or lower `cutoff`.
- **Noisy components**: increase `smoothing_kernel_size` (3 → 5 → 9).
- **Spike contamination in input**: enable `z_score_cutoff = 10`.
