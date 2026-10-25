# REAL_DATA.md — first real hologram through the pipeline (2026-10-25)

## What was obtained

### 1. Real off-axis hologram — OBTAINED AND RUN (this is the headline)

- **File:** `fly-20x.bmp` (2048×2048, uint8, range 0–250), copied here
  from the public pyDHM repository (MIT),
  https://github.com/catrujilla/pyDHM, `data/compensation samples/`.
- **What it is:** an EXPERIMENTAL off-axis telecentric DHM hologram of a
  fly head, supplied by the 3D Imaging & Display Lab / Opto-Digital
  group to validate pyDHM (see that repo's README acknowledgment).
- **Documented parameters** (pyDHM `examples_compensation.py`, FRS/ERS/CFS
  examples for this exact file): wavelength **0.532 µm**, pixel size
  **2.4 µm × 2.4 µm**.
- **Carrier, measured from the hologram itself:** spectral peak at
  (−264, +290) bins from DC (392.2 bins away) — a diagonal tilt, so the
  package's x-only `reconstruct_field` does NOT apply directly. The
  script below generalises crop-and-shift to 2-D.
- **Script:** `real_fly_focus.py` (separate script; package code untouched,
  suite still passes). One adaptation forced by scale: the package's
  per-distance transfer cache OOM-killed the first run at 2048²
  (exit 137, 81 cached 67 MB transfer functions), so the script uses an
  uncached propagator with the identical formula.

### 2. JoVE DLHM database — metadata obtained, IMAGES NOT ON OSF (blocked)

- JoVE page https://visualize.jove.com/38437255 gives only the abstract.
- The paper (Appl. Opt. 63(7) B49, 2024) points to OSF
  https://doi.org/10.17605/OSF.IO/SYHT9 ("HQ-DLHM-OD-Database").
- OSF API listing is definitive: the node contains exactly ONE file,
  `ODP-DLHM-Database.csv` (680,229 bytes, downloaded here), and zero
  child components. No hologram images are hosted there.
- The CSV itself is real metadata: 11,826 rows (paper claims 11,760 —
  33 extra rows in each of MW/SW; minor discrepancy, noted not hidden),
  columns `path;Wavelength;L_value;z_value`; wavelengths 0.405 / 0.51 /
  0.633 / 0.654 µm; z_value 0.32–3.95 (units unstated in the CSV);
  paths are relative (`./MW-Dataset/405/...`) with no base URL and no
  pixel pitch anywhere in the file. A web search for a mirror of the
  actual `./MW-Dataset/` images found none.
- Even if the images surfaced, DLHM is lensless IN-LINE: the off-axis
  pipeline could not use it directly (twin ambiguity returns).
- Verdict: route blocked at "metadata only, images not publicly
  hosted at the cited location". Best next step for this route: contact
  the authors (UNAL-Optodigital) for the image hosting location, or use
  the paper's tables (paywalled at Optica) for pixel pitch — neither
  was pursued here (no author contact without Jean's word; no paywall).

### 3. Other routes checked (search only)

- **pyDHM siblings** (same repo, same provenance): `holonotele_fly-resize.png`
  (non-telecentric, 0.633 µm, 7 µm), `hologram UofM.png`, `hologram EAFIT.jpg`,
  `RBC1.tif`, `star-target.png` — all show clean off-axis carriers in
  their spectra. Not run; one real hologram was enough for the findings
  below, and none of them documents a ground-truth focus distance either.
- **Zenodo 10.5281/zenodo.23050910** (PLOS ONE off-axis paper): NUMERICAL
  (simulated) holograms — not real recordings, not pursued.
- **BioDIHM** (Univ. St-Étienne, labh-curien.univ-st-etienne.fr/public-data/BioDIHM/):
  real experimental IN-LINE holograms (beads, E. coli) — wrong geometry
  for this pipeline, noted as a fallback for in-line work only.
- **Hugging Face inline-digital-holography-v2**: SYNTHETIC in-line data.
  Not real, not off-axis. Not pursued.

## Measured results — fly hologram, sweep −10…+10 mm, dx = 2.4 µm

| metric | extremum | reading |
|---|---|---|
| imag-energy (package metric) | min +7.25 mm, 0.4308 vs median 0.4990 (1.16×) | No signal. Synthetic clean case was 50.8×. A fly is not a real amplitude object, so the field is never real at focus. |
| complex-field gradient sharpness (package `field_sharpness`) | max −1.75 mm, 7.5009 vs median 7.5003 (0.01% variation) | NOT A FOCUS METRIC. See proof below. The "maximum" is noise. |
| amplitude variance (diagnostic) | MIN +1.00 mm coarse, **+1.10 mm refined** (6.20e4 vs 6.93e4 at z=0) | Candidate focus. Shallow (10% dip), no ground truth to check against. |
| amplitude at z = 0 | crisp fly head (hairs, eye facets) | The recording is an image-plane hologram, so focus near 0 is expected from the geometry — an expectation, NOT a documented truth value. |
| amplitude at ±10 mm | corr vs z=0 drops to 0.49–0.69; +10 mm image visibly blurred/ringed | Propagation is working; the metrics, not the physics, are what fail. |

**No measured-vs-truth error can be claimed for this hologram.**
pyDHM documents no focus distance for it. The honest statement is:
reconstruction at z = 0 yields a crisp fly head; amplitude variance
suggests the true plane is ≈ +1.1 mm from the hologram plane; the two
package metrics do not identify any plane.

## Two findings that change the package

### A. `field_sharpness` is mathematically invariant — do not use it

Mean |∇U|² = ∫|k|²|Û(k)|² dk (Parseval). Angular-spectrum propagation
multiplies Û(k) by a unit-magnitude phase, so |Û(k)| is unchanged and
the metric is EXACTLY conserved (the observed 0.01% wobble is
finite-difference roundoff, and the evanescent cutoff removes nothing
here since λ·f_max = 0.11 < 1). The docstring claim that gradient
energy "peaks at true focus" is false for the complex field under this
propagator. The package docstring should be corrected; this script
does not edit the package, it only reports the falsification.
(Amplitude gradient is a different quantity — and on this phase-heavy
object it is HIGHER out of focus, because defocus converts phase into
amplitude contrast.)

### B. Focus metrics invert for phase objects

For a phase object, amplitude variance is MINIMAL at focus (phase
contributes nothing to amplitude there; defocus creates amplitude
texture). The synthetic pipeline's metrics all assume a real amplitude
object. On real biological samples the sign flips: look for a variance
MINIMUM, not a maximum — and even then the dip here is only 10%, so a
single global metric is a weak reed for thick samples.

## Caveats on the numbers

- Measured z scales as dx². If pyDHM's 2.4 µm is not the effective
  hologram-plane pitch (their commented propagation example uses
  5.4 µm for reasons never explained in the repo), every z above
  rescales by (dx/2.4 µm)² — e.g. ×5.06 at 5.4 µm. The +1.1 mm figure
  is conditional on dx = 2.4 µm.
- Sideband window ±200 bins truncates the object spectrum (same cap
  family as the synthetic 0.951 correlation ceiling); a wider window
  risks DC/twin overlap (DC is 392 bins from the carrier peak).
- FFT periodic-boundary ringing is visible at the frame edges in the
  ±10 mm reconstructions; it inflates variance away from focus and
  slightly shallows the +1.1 mm dip. The dip survives it but is not
  clean.
- Pixel pitch for the DLHM CSV is absent, and z_value units are
  unstated — no DLHM number is usable without the paper's tables.

## Best next route

1. **A calibration target with a documented distance.** The strongest
   next test is an off-axis hologram of a USAF target (or similar)
   recorded at a stated object–sensor distance — e.g. a lab recording
   where the distance is part of the experiment, not recovered by
   focusing. pyDHM's `star-target.png` / USAF phase-shifting samples
   are candidates but still lack a stated distance; a new recording
   with a ruler in the loop beats another archive hunt.
2. **Fix first, then re-measure:** correct the `field_sharpness`
   docstring, and add an amplitude-variance-minimum metric for phase
   objects (with its shallow-dip warning) before any further real-data
   claim.
3. **DLHM images:** only via the authors. Do not present the CSV as
   "the database obtained" — it is the index, not the data.

## Files in this directory

- `fly-20x.bmp` — the real hologram (pyDHM, MIT)
- `real_fly_focus.py` — the run script (separate from the package)
- `fly_amp_z0.png` — reconstruction at z = 0 (crisp fly head)
- `fly_amp_best.png` — reconstruction at the amp-variance minimum (+1.10 mm)
- `fly_amp_plus10mm.png`, `fly_amp_minus10mm.png` — defocused controls
- `fly_focus_curve.png` — imag-energy + sharpness curves (flat = the finding)
- `fly_spectrum_log.png` — log spectrum (DC + twin sidebands visible)
- `ODP-DLHM-Database.csv` — DLHM metadata index (NOT the hologram images)
