# DESIGN.md — wave-hologram

> **v0.1.3 note:** the imaginary-energy metric and the 'comb' described below are superseded. imag_energy depends on the arbitrary reference phase, and the comb is its lambda/2 oscillation aliased by 1 um sampling. See RESULTS_v0.1.3.md (sections 2-3). Kept unchanged as the historical record.


## What this is

A self-contained digital holography pipeline: record a hologram, reconstruct
the object by numerical back-propagation, and estimate the focus distance
instead of assuming it. Built from a Perplexity thread's focus-sweep export
snippet, corrected where the snippet's setup was under-determined.

## Physics

Propagation is the angular spectrum method. A field U(x, y; 0) propagated by
distance z becomes, in Fourier space, F(U) multiplied by

    H(fx, fy; z) = exp(i k z sqrt(1 - (λ fx)^2 - (λ fy)^2)),  k = 2π/λ,

with evanescent components (the square root's argument ≤ 0) set to zero.
Parameters match the thread: λ = 532 nm, pixel pitch 6.5 µm, 256×256 grid,
object 5 mm from the sensor. The transfer function depends only on z, so it
is cached per distance rather than recomputed per frame.

## Recording: why off-axis

The hologram is H = |U + R|^2 with a tilted reference R = exp(i 2π u0 x),
u0 chosen so the object sideband sits 80 FFT bins from DC. Reconstruction
crops that sideband (±48 bins) and shifts it to baseband, recovering U up
to a constant. This is the Leith-Upatnieks construction, and it matters:
with an in-line hologram (R = 1) the recorded field is real, so the
reconstruction intensity at +d equals that at -d exactly (tested in
`test_inline_twin_ambiguity`). The twin image makes the focus distance
ambiguous in sign. Off-axis recording removes the twin; there is no
ambiguity left to argue about.

## Focus metric: why imaginary energy

The object is a real amplitude mask. At true focus the reconstructed field
is real, so mean(Im^2)/mean(|.|^2) bottoms out. This beats intensity
statistics for two measured reasons: intensity variance is dominated by
hologram fringes at d = 0 (the first version measured 0.00 mm), and for
smooth objects intensity variance barely changes with defocus. The
imag-energy metric is structural, not amplitude-based: uniform noise
raises the floor everywhere without moving the dip, which is why focus
survives to 2x noise while image quality dies at 0.4x.

The focus curve oscillates (phase wrapping gives local dips), but the
global minimum is unambiguous: 0.0094 against ~0.06 for the next dip.

## Noise models

- Gaussian: additive read noise, sigma as a fraction of the hologram mean,
  clipped at zero (a sensor cannot record negative intensity).
- Poisson: shot noise. The hologram is scaled to a peak photon count,
  sampled, and scaled back. Fewer photons = noisier.

Both models act on the recorded hologram, before sideband isolation, so
the sideband crop's bandpass effect on broadband noise is honestly
included in the measurement.

## Measured results (2026-10-24, this machine)

| condition | best focus | contrast | corr |
|---|---:|---:|---:|
| clean | -5.00 mm | 50.8x | 0.951 |
| gaussian 0.20 | -5.00 mm | 8.7x | 0.943 |
| gaussian 0.40 | -5.00 mm | 3.2x | 0.919 |
| gaussian 0.80 | -5.00 mm | 1.6x | 0.803 |
| gaussian 2.00 | -5.00 mm | 1.1x | 0.348 |
| gaussian 4.00 | wanders | 1.0x | 0.091 |
| poisson 100 photons | -5.00 mm | 7.2x | 0.939 |
| poisson 10 photons | -5.00 mm | 1.7x | 0.808 |
| poisson 3 photons | -5.00 mm | 1.2x | 0.576 |
| poisson 1 photon | wanders | 1.1x | 0.293 |

Focus location holds to 2x gaussian noise and to 3 peak photons; image
quality stays above 0.9 through 0.4x gaussian noise, falls to 0.80 at
0.8x, and collapses past 2x. Past either break the reported plane is
noise (different draws wander to different planes). Finding the plane
is not the same as seeing the object.

## What was fixed, not hidden

1. v1 (in-line + intensity variance) measured 0.00 mm. Cause: fringe
   variance at d = 0 plus the twin ambiguity. Fix: off-axis + imag-energy.
2. v1's focus curve was flat for smooth objects. Fix: same metric change;
   the metric now keys on field realness, which is sharp by construction.
3. Sideband truncation caps recovery at corr 0.951 on the binary object
   (its sharp edges carry spectrum outside the ±48-bin window). The cap
   is stated, not rounded away.

## Honest limits

- Synthetic, noise-free baseline; real sensors add dust, non-uniform
  reference, and vibration. No real sensor image has been tested yet.
- The MP4 stores rendered views, not the complex field; playback speed
  is not propagation time.
- This is OPTICAL holography (wavefronts), not the membrane portfolio's
  associative holographic memory (W = sum |p><p|), and not the LED fan
  display. Same word, three different things.
