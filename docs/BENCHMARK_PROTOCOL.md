# Benchmark protocol

**Protocol version:** `PROTOCOL-v1.0.0`  
> This protocol is locked before benchmark execution.  
> Changes after execution require a new protocol version and must not overwrite prior results.

## Scope and estimands

This protocol measures synthetic numerical behavior only. The primary estimand is
the absolute error between a predeclared focus-metric extremum and the known
relative focus distance, per scene/seed/condition. The statistical unit is an
independently generated scene realization, not a reconstruction plane or pixel.

## Configuration and preprocessing

Use the checked-in YAML configuration as the authority. Generate deterministic
point, extended, and textured objects from the recorded seed. Normalize only by
the configuration's specified deterministic rule. Propagate with
`wave_hologram.angular_spectrum`; record with `record_hologram`; estimate the
carrier with `estimate_carrier_2d`; reconstruct with `reconstruct_field`.
Apply no hidden denoising, crop, threshold, or manual selection.

## Metrics and focus rule

Evaluate the metrics in [`METRICS.md`](METRICS.md). The primary focus metric is
`realness`, minimized over the supplied distance sweep; `ac_realness` is a
predeclared secondary metric. Other metrics must declare their polarity in the
result. Focus is the coordinate of the global min/max according to that polarity.
Record boundary status; a boundary extremum is a failure, not a successful focus.
RMSE/PSNR/SSIM are computed on amplitude against the known reference reconstruction
after a predeclared global scale fit; do not use them as a physical-truth proxy.

## Noise and random seeds

Add zero-mean Gaussian noise to the hologram with sigma defined relative to the
noise-free hologram's RMS. Negative intensity clipping is enabled explicitly in the
checked-in configurations to model a non-negative camera readout. Each condition uses seeds 1001 through 1010 at minimum; the
nominal config uses 10 seeds. Scene seed, noise seed, and parameter configuration
are recorded independently. Repeated metrics on one seed are paired only for
method comparisons sharing that exact scene and noise realization.

## Statistics and decision thresholds

For each condition report N, mean, median, standard deviation, 95% confidence
interval, minimum, maximum, and failure fraction. Use percentile bootstrap intervals
(10,000 resamples, recorded bootstrap seed) when assumptions are unclear. Report
effect size and interval for comparisons; use paired bootstrap differences only
for matched scenes. Correct multiple mode/condition tests with Holm for confirmatory
tests; use Benjamini-Hochberg for exploratory tests. A p-value alone is not a
success criterion.

The minimum success rule is an interior focus extremum within one sweep step of
truth in at least 9/10 seeds in every predeclared cell. This is an engineering gate, not evidence of
generalization. Missing seeds, missing metrics, invalid checksums, absent commit
identity, or protocol mismatch invalidate the result bundle.

## Output inventory

Write machine-readable outputs before narrative output:

1. `per_seed_results.csv` and `metrics.csv` (long form).
2. `summary.csv` with all descriptive statistics and confidence limits.
3. `manifest.json`, `environment.json`, `config.yaml`, `code_commit.txt`,
   `summary.json`, and `checksums.sha256`.
4. Code-generated figures: object/hologram; amplitude/phase reconstruction;
   focus curve; focus error by noise; PSNR/SSIM by condition; and failure gallery.
   The real-data slope figure is permitted only from the locked real-data protocol.

Every figure caption must identify SYNTHETIC or REPORTED data, condition, metric,
uncertainty, units, and limitation. Do not use an axis range that exaggerates a
difference.

## Claim gates

Only condition-bounded claims with N, uncertainty, protocol version, and failure
rate may pass. Absolute-focus or product claims are forbidden. A result is
publishable only when every required seed, metric, provenance field, and checksum
is present and validated. The full claim vocabulary is in [`CLAIMS.md`](CLAIMS.md).
