# REVIEW.md — wave-hologram v0.1.0 + sensor test

Review date: 2026-10-24. All numbers below were measured by running the
code in this folder, not assumed. Truth throughout: object 5.00 mm from
the sensor, focus expected at -5.00 mm back-propagation.

## What is solid

- **Focus exactness.** Clean pipeline: measured -5.00 mm, error 0.000 mm,
  imag-energy 0.0094 vs median 0.4764 (50.8x contrast), recovery corr 0.951.
- **Tests.** 11/11 pass (`python3 -m unittest discover -s tests`), including
  round-trip propagation (corr > 0.999), focus exactness, the in-line
  twin-ambiguity proof, noise hold/break, shot noise, and the edge-window
  regression test added during this review (see below).
- **Physics.** Angular-spectrum propagation with evanescent cutoff and a
  cached transfer function; off-axis sideband isolation (the
  Leith-Upatnieks answer to Gabor's in-line twin problem); a focus metric
  (imaginary-energy fraction) that keys on field realness at focus rather
  than on brightness, which is why it survives noise that kills the image.
- **Documented failures.** The two earlier wrong-plane measurements
  (in-line + intensity variance at 0.00 mm; flat metric on smooth objects)
  are written up in DESIGN.md with causes and fixes, not hidden.

## What is weak

- **Synthetic only.** Every hologram so far was computed, not recorded.
  No real sensor image has been run through this pipeline yet.
- **Sideband truncation cap.** The binary object has sharp edges whose
  spectrum spills outside the ±48-bin sideband window, capping recovery
  at corr 0.951. The cap is stated, not rounded away, but it is a real
  ceiling on this configuration.
- **One geometry tested.** Object distance 5 mm, wavelength 532 nm, pixel
  pitch 6.5 µm, 256×256 grid, one binary object. Other distances, phase
  objects, and thicker samples are untested.
- **Noise models are partial.** Gaussian read noise and Poisson shot
  noise only. Vibration, thermal drift, partial coherence, and a
  non-planar reference wavefront are not modelled.
- **Sensor-test cases are not paired.** Each case draws fresh random
  numbers from one advancing RNG stream, so small contrast differences
  between cases are partly draw variation, not a ranking. The robust
  claims are the ones that hold across all cases (exact focus, corr floor).

## Sensor test: the closest-to-physical test available here

`sensor_test.py` (seed 11) degrades the hologram the way hardware would
and re-measures focus. Full table:

| case | best focus | error | contrast | corr |
|---|---:|---:|---:|---:|
| clean | -5.00 mm | 0.00 | 50.8x | 0.951 |
| vignette (Gaussian falloff) | -5.00 mm | 0.00 | 52.9x | 0.914 |
| dust (25 opaque disks) | -5.00 mm | 0.00 | 25.7x | 0.942 |
| gain (2% fixed-pattern) | -5.00 mm | 0.00 | 47.9x | 0.950 |
| quantize (12-bit ADC) | -5.00 mm | 0.00 | 50.8x | 0.951 |
| saturate (top 1% clipped) | -5.00 mm | 0.00 | 27.8x | 0.939 |
| carrier+2 (2-bin reference error) | -5.00 mm | 0.00 | 1.2x | 0.944 |
| combined (all + 5% read noise) | -5.00 mm | 0.00 | 34.6x | 0.912 |

**Reading: no artifact moves the focus plane.** Every case, including
the full combined stack, finds -5.00 mm exactly. What artifacts cost is
image quality (corr 0.912–0.951 across the board) and, in one case,
metric trust:

- Vignetting, dust, saturation each roughly halve the metric contrast
  but leave focus and image nearly intact. Quantization at 12 bits is
  invisible (identical numbers to clean). Fixed-pattern gain is negligible.
- **Carrier miscalibration is the one artifact that threatens the
  measurement itself** — not by moving the answer, but by flattening the
  metric to 1.2x. A 2-bin reference-angle error puts a residual phase
  ramp on the reconstructed field, so the field is never real even at
  focus and the dip nearly vanishes — while the recovered intensity
  (corr 0.944) is untouched, because intensity ignores a pure phase
  ramp. Under additional noise, a 1.2x dip would be the first thing to
  lose. Practical lesson: calibrate the reference angle; the pipeline
  tolerates everything else tested.
- The combined stack — vignetting + dust + gain + quantization + read
  noise together — still lands exactly on -5.00 mm at corr 0.912. That is
  the closest-to-physical result this folder can honestly claim today.

**Bug found and fixed during this review.** The first run of
`sensor_test.py` crashed on the carrier+2 case before reaching the
combined case:

```
ValueError: could not broadcast input array from shape (97,0) into shape (97,97)
```

Cause: with car=82 the assumed sideband window (columns -2..94 of the
shifted spectrum) runs off the array edge, and NumPy's negative slice
start wrapped it into an empty slice. Fix in `reconstruct_field`: the
window is now clipped to the array with frequency alignment preserved
(missing bins stay zero instead of wrapping in from the opposite side
of the spectrum), and a regression test guards it. Standard cases
(car=80) are bit-identical to before the fix; the full suite passes
11/11 after it.

## What a true physical test would need

A real recorded hologram plus its metadata, all three documented:

1. **Pixel pitch** of the sensor (sets the frequency grid).
2. **Wavelength** of the source (sets the propagation kernel).
3. **Known object distance** — the ground truth the measured focus is
   checked against. Without it, a reconstruction is a picture, not a test.

Plus, honestly: a characterization of the reference beam (angle, and
whether it is planar), because the sensor test shows reference-angle
error is the one miscalibration that weakens the measurement itself.

## Closest real-data option (search only, not tested)

The JoVE open-access database for digital lensless holographic
microscopy (DLHM): 11,760 experimental lensless holograms of bio and
non-bio samples, multi-wavelength (654/510/405 nm) and single-wavelength
(633 nm He-Ne) sets, with acquisition parameters described in the paper:
https://visualize.jove.com/38437255

Flagged honestly: this was found by web search only. Nothing was
downloaded, and nothing from it has been run through this pipeline. Two
caveats before anyone treats it as a drop-in test: DLHM is lensless and
in-line, so this package's off-axis sideband pipeline would not apply
directly (the twin ambiguity returns); and each recording's ground-truth
object distance would need to be confirmed from the dataset's metadata
before any focus claim could be checked against it. It is the closest
real data found, not a validated test.
