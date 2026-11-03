# WaveHolo

Reproducible off-axis digital-holography research pipeline for synthetic hologram
reconstruction and relative focus-distance estimation.

**Evidence status:** SYNTHETIC pipeline; limited real-data slope check REPORTED under a locked protocol.

**What this is:** A research prototype for evaluating holographic reconstruction and
relative focus estimation under controlled conditions.

**What this is not:** An absolute-focus validation, production autofocus system,
holographic display, microscope product, medical device, or validated biological-imaging tool.

## What this repository does

WaveHolo generates synthetic objects and off-axis holograms, estimates the carrier,
reconstructs complex fields with angular-spectrum propagation, and compares
distance-dependent focus metrics. It also retains a narrow external-data analysis:
530 microscope holograms at six stage positions.

## What this repository does not do

It does not establish absolute focus accuracy, generalization to real holograms,
clinical or biological utility, production readiness, or superiority over the field.
The browser demonstration is synthetic and is not measurement evidence.

## Scientific question

Under controlled synthetic off-axis digital-holography conditions, how accurately
and robustly does the implemented focus metric recover relative focus distance?

## Method overview

The implementation forms `H = |U + R|^2`, estimates the off-axis carrier, isolates
the object sideband, and numerically back-propagates across a distance sweep.
`realness` and `ac_realness` are minimized; each alternative metric has its own
declared polarity. See [the scientific protocol](docs/SCIENTIFIC_PROTOCOL.md).

### Related work

Digital-holography focus criteria are an established research area. WaveHolo's
metrics have not been compared against the external methods cited in the historical
studies. Earlier criterion families include amplitude-integral, Fourier-weighted,
gradient, Tamura, and sparsity approaches; see the original references and their
limits in [RESULTS_v0.1.3.md](RESULTS_v0.1.3.md).

## Primary metrics

Focus error and trend agreement are primary. Reconstruction RMSE, PSNR, SSIM, and
normalized cross-correlation are descriptive image metrics; none proves physical
correctness. Definitions and units are in [METRICS.md](docs/METRICS.md).

## Reproducible quickstart

```bash
python -m pip install -e ".[test]"
python examples/quickstart.py
make benchmark CONFIG=benchmarks/configs/waveholo_smoke.yaml
```

The benchmark command writes machine-readable outputs to `benchmark-output/` by default
or to the path passed as `OUTPUT=...`.
The smoke configuration is a small engineering check, not a basis for a scientific
performance claim. Its claim gate is expected to fail (fewer than 10 seeds and
boundary extrema); exit code 2 reports that failure while preserving the outputs.

## Benchmark results summary

**UNRUN:** the preregistered factorial benchmark has not been run. No result is
estimated or projected. Run `make benchmark CONFIG=benchmarks/configs/waveholo_nominal.yaml`
to produce a new, commit-identified result bundle.

Previously reported observations remain historical evidence: a clean synthetic
demonstration reported a 2 µm focus error in one scene; a separate real-data analysis
reported a negative slope of approximately 96 focus-µm per stage-µm across six
positions. The latter is not absolute-focus validation and is only partially
reproducible; it is a slope check, not an absolute-focus test. See
[REAL_DATA_PROTOCOL.md](docs/REAL_DATA_PROTOCOL.md).

Historical synthetic studies also recorded sharp-object bias of −50 to +64 µm from
sideband cross-talk and a noise study with 100 *peak* photons per pixel; those
finite explorations are not the new factorial benchmark and are not generalized to
other scenes or optics.

## Limitations and failure modes

Metric extrema can be multiple, shallow, or scene-dependent. Sampling, sideband
overlap, carrier, noise, object structure, and reconstruction distance affect results.
Negative cases—including phase-object failures and sideband cross-talk—are retained.
See [LIMITATIONS.md](docs/LIMITATIONS.md).

## Real-data status

REPORTED, limited microscope slope/trend check only; 530 holograms at 6 of 100 stage
positions. Independent rerun is PARTIALLY REPRODUCIBLE because the full external
dataset and some historical environment/protocol metadata are not bundled. The
check does not test absolute focus.

## Installation

Python 3.9+ is retained for compatibility with the current published package.
Install from a clone with `python -m pip install -e .`; development dependencies
are available with `python -m pip install -e ".[test]"`. No large external dataset
is downloaded by installation.

## Testing

Run `make test` for the maintained suite. CI runs unit and integration checks and a
fast smoke benchmark; the full factorial experiment is not a per-commit job.

## Citation

Cite the version identified in [`CITATION.cff`](CITATION.cff). The existing Zenodo
archive DOI is 10.5281/zenodo.23653721.

## Licence and provenance

AGPL-3.0-only. External datasets and assets retain their own licences, described in
[`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md). SYNTHETIC results are project
generated; REPORTED evidence must identify its external data source and protocol.
