# Understanding Explained Variance Ratio in PCA

## Quick Answer

**Explained Variance Ratio** tells you what fraction of the *total variation* in your data is captured by the first N components.

- **High % (e.g., 99%)** = Spectra are **very similar to each other** (low variation)
- **Low % (e.g., 61%)** = Spectra are **very different from each other** (high variation)

## Detailed Explanation

### What is "Variance" in This Context?

In PCA, variance measures **how much the spectra differ from each other** across the dataset.

- If all 2268 spectra in a mission were identical, the variance would be 0%
- If spectra were all very different from each other, the variance would be high

### What Does Explained Variance Ratio Mean?

When we say "the 5 components explain 99% of variance", we mean:

> **99% of the differences between your spectra can be described using just these 5 principal component directions**

The remaining 1% would require additional components (6, 7, 8, ..., etc.) to capture.

### Interpretation

#### High Explained Variance (90-99%)
```
Example: 2016-05-12_GR_F at 99.73%
```

**What this means:**
- The spectra are **very coherent** and **similar to each other**
- They're dominated by a single pattern (the continuum level)
- Small variations are captured well by a few components
- The data is "clean" and well-behaved

**In practical terms:**
- The reference spectra form a tight cluster in high-dimensional space
- Most spectral variation follows predictable patterns
- Good candidate for using as a calibration reference

#### Moderate Explained Variance (70-85%)
```
Example: 2016-05-18_GR_F at 81.19%
```

**What this means:**
- The spectra have **notable differences** from each other
- More diversity in the reference data
- Still dominated by main patterns, but with more secondary variation
- The data has some noise or instrumental variation

**In practical terms:**
- Reference spectra are more spread out in feature space
- You need more components to capture the full behavior
- Still usable but with somewhat more uncertainty

#### Low Explained Variance (50-70%)
```
Example: 2016-05-19_GR_F at 61.54%
```

**What this means:**
- The spectra are **quite diverse** and **different from each other**
- Significant unexplained variation in the data
- No single dominant pattern
- Higher noise or instrumental effects

**In practical terms:**
- Reference spectra are scattered widely in feature space
- Would need 10+ components to reach 99%
- Less reliable as a calibration reference
- May indicate instrument problems during observation

## Your Mission Data Summary

| Mission | Spectra | Total Variance | Interpretation |
|---------|---------|-----------------|-----------------|
| 2016-05-12_GR_F | 2,268 | **99.73%** | Very clean, coherent data |
| 2017-02-01_GR_F | 13,986 | **99.76%** | Very clean, coherent data |
| 2016-05-18_GR_F | 6,426 | **81.19%** | Moderate variation, usable |
| 2016-05-19_GR_F | 3,024 | **61.54%** | High variation, noisier |

### Component Breakdown

#### Clean Missions (99%+)
The variance is dominated by Component 1:
```
2016-05-12: C1=99.41% (almost everything is Component 1)
2017-02-01: C1=99.24% (almost everything is Component 1)
```
This is expected and good - Component 1 captures the sky continuum level.

#### Noisier Missions
The variance is more distributed:
```
2016-05-18: C1=74.85%, C2=4.54%, ... (more variation spread across components)
2016-05-19: C1=58.10%, C2=1.85%, ... (even more spread out)
```
Components 2+ capture secondary effects (noise, instrumental variation).

## What This Means for Your Analysis

### For Spectral Correction

**Clean missions (99%+):**
- Simple correction: just subtract Component 1 (the sky continuum)
- High confidence that the correction is accurate
- Reference spectra are very reliable

**Noisier missions (61-81%):**
- May need to subtract Components 1 AND 2
- More uncertainty in correction
- Reference spectra have more noise mixed in
- Consider masking or filtering noisy spectra first

### For Using Components

You have mission-specific PCA bases:
- Each mission gets its own components
- Accounts for mission-specific instrumental effects
- Better correction than using a global basis
- 2016-05-19 may need special handling (more noisy)

## Key Takeaway

**Explained Variance Ratio is NOT about data quality in an absolute sense** - it's about **data coherence**.

- 99% = "Your spectra are very similar" = coherent reference set
- 61% = "Your spectra are very different" = incoherent reference set

For your use case (sky correction), you want high explained variance because it means your reference spectra form a reliable basis for describing the sky pattern. The 99%+ missions are ideal; the 61% mission is more problematic.
