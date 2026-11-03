# Scientific protocol

**Status:** protocol specification; the new factorial benchmark is UNRUN.
**Scope:** SYNTHETIC optical holography. Real-data evidence is handled separately.

## 1. Research question

**Primary question:** Under controlled synthetic off-axis digital-holography
conditions, how accurately and robustly does the implemented focus metric recover
relative focus distance?

**Secondary questions:**
- How does reconstruction quality vary with sampling, carrier frequency, noise, and object complexity?
- How does focus-metric stability vary around the true focus?
- What perturbations cause focus estimation to fail?
- Does the limited real-data microscope slope check reproduce under the locked protocol?

## 2. Hypotheses

- H1: For noiseless synthetic holograms with adequate sampling, the focus metric
  reaches an extremum at or near the true relative focus distance.
- H2: Increasing additive Gaussian noise degrades focus-metric sharpness and
  increases absolute focus-estimation error.
- H3: The real-data microscope slope check reproduces the previously reported
  positive focus trend only under the locked preprocessing and reconstruction
  protocol. The historical fit is negative under the currently documented stage
  convention; this sign discrepancy is unresolved, so H3 is not presently
  testable and no slope sign may be changed to force agreement.

## 3. Experimental design

Independent factors: grid resolution (256, 512, 1024); carrier (low, nominal,
high); object class (point, simple extended, textured); additive Gaussian sigma
(0, low, medium, high); and a reconstruction-distance sweep centered on known
truth. The exact numerical levels, sampling geometry, object definitions, and sweep
grid are versioned in the YAML configuration, not inferred from labels.

Dependent outcomes: focus-curve extremum and polarity, signed and absolute focus
error, within-tolerance success and failure rate, focus-curve sharpness, RMSE, PSNR,
SSIM, normalized cross-correlation, runtime, and peak process memory. A non-random
deterministic condition is run once for exact regression; every stochastic condition
uses at least 10 independent seeds.

Controlled variables include wavelength, pixel pitch, propagation model, reference
amplitude/phase policy, sideband method/window, object normalization, reconstruction
implementation, and boundary handling. Each random generator is initialized from
the recorded seed. No learned component is currently in the pipeline; therefore
there is no train/validation/test split. If a learned component is added, split by
scene family and seed before fitting or tuning.

Exclude only invalid numerical inputs, non-finite intermediate output, or a
predeclared sideband-separation violation. Preserve excluded rows and reasons.
Stop a run only on completion, resource failure, or invalid configuration; report
partial output as PARTIAL/FAILED, never as complete. Use Holm correction for
predeclared confirmatory comparisons and Benjamini-Hochberg for exploratory
families. Success means the metric-selected extremum is interior and the
predeclared focus tolerance is met; failure includes boundary extrema, missing
metrics, and tolerance failures. Do not select thresholds after observing results.

## 4. Preregistration

See [`BENCHMARK_PROTOCOL.md`](BENCHMARK_PROTOCOL.md), protocol `PROTOCOL-v1.0.0`.
No experiment output has been generated from that protocol yet.
