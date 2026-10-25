# DISCOVERY.md — synthetic discovery sweep (2026-10-24)

> **v0.1.3 note:** the imaginary-energy metric and the 'comb' described below are superseded. imag_energy depends on the arbitrary reference phase, and the comb is its lambda/2 oscillation aliased by 1 um sampling. See RESULTS_v0.1.3.md (sections 2-3). Kept unchanged as the historical record.


Run: `python3 discovery_sweep.py` (seed 23). Every number below was
measured by running that script, not assumed. The headline findings
change how the rest of this package should be read — in particular,
finding 5 re-qualifies every "exact focus" claim made so far.

## Finding 1 — distance: imag-energy holds 1–50 mm (on-grid truths)

| z_obj | imag best | err | contrast | pvar best | err | pvar contrast | corr |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1 mm | -1.00 | 0.00 | 45.7x | -1.00 | 0.00 | 2.4x | 0.945 |
| 2 mm | -2.00 | 0.00 | 18.5x | -2.00 | 0.00 | 2.2x | 0.939 |
| 5 mm | -5.00 | 0.00 | 50.8x | -5.00 | 0.00 | 3.4x | 0.951 |
| 10 mm | -10.00 | 0.00 | 50.1x | -10.00 | 0.00 | 5.1x | 0.951 |
| 20 mm | -20.00 | 0.00 | 52.1x | -20.00 | 0.00 | 7.4x | 0.957 |
| 50 mm | -50.00 | 0.00 | 61.0x | -50.00 | 0.00 | 15.6x | 0.955 |

No breakdown anywhere in 1–50 mm: no aliasing cliff, no metric
flattening, corr 0.94–0.96 throughout. Caveat from finding 5: these
truths are all on the 0.1 mm grid, so this is "holds when the grid
samples the central tooth", not "holds for any distance".

## Finding 2 — phase objects: BOTH metrics fail (clean failures)

Truth -5.00 mm throughout. Phase corr is measured at the
phase-variance plane, so a wrong imag argmin cannot fake it.

| object | imag best | imag err | imag contrast | pvar best | pvar err | pvar contrast | phase corr |
|---|---:|---:|---:|---:|---:|---:|---:|
| amplitude (control) | -5.00 | 0.00 | 50.8x | -5.00 | 0.00 | 3.4x | 0.951 |
| phase pi/4 | -8.30 | 3.30 | 19.5x | -10.00 | 5.00 | 1.3x | 0.001 |
| phase pi/2 | -8.30 | 3.30 | 9.6x | -10.00 | 5.00 | 1.4x | 0.085 |
| phase 3pi/4 | -6.70 | 1.70 | 13.5x | -10.00 | 5.00 | 1.2x | 0.191 |
| phase disk pi/2, fill .48 | -9.90 | 4.90 | 1.1x | -10.00 | 5.00 | 1.0x | 0.055 |
| mixed amp+phase | -9.90 | 4.90 | 3.3x | -10.00 | 5.00 | 1.0x | 0.020 |

Why imag-energy fails: it does not just shallow — defocus
oscillation dips BEAT the true-focus dip (pi/2: a dip at -8.30 with
contrast 9.6x wins over focus). The disk case is the predicted
inversion/flattening: fill 0.48 x sin^2(pi/2) meets the ~0.5 defocused
background, contrast 1.1x, no measurement at all.

Why phase variance fails: the uniform background phasor dominates
the intensity-weighted mean at EVERY plane (a plane wave propagates
unchanged), so the metric is flat to within ±0.02 across the whole
sweep (pi/2: 0.042–0.068) and its argmin is the sweep edge — noise,
not focus. The reconstructed phase at true focus is actually good
(phase corr 0.965 at -5.00 for pi/2); the metric just cannot find
the plane.

Also falsified along the way: gradient sharpness. It is
propagation-invariant by Parseval (measured 0.009720 at 0, -2.5, -5,
-7.5 mm) — any quadratic functional of the field is invariant the
same way, so a focus metric must be nonlinear. `field_sharpness`
stays in the module only as the documented dead end.

Other candidates probed, none a universal metric (all clean, 0.1 mm
grid): phasor-gradient energy is flat for phase objects; gradient
kurtosis is ELEVATED at focus for all types (amplitude 34.1 vs
median 17.7, phase 33.8 vs 20.8, mixed 44.2 vs 25.9) but its argmax
is unreliable (phase: sensor-plane edge; mixed: -5.50); intensity
variance is flat for amplitude (argmax -7.00, wrong) and mixed, and
only a shallow near-miss for pure phase (argmin -4.90, 0.10 off).
Phase-object autofocus remains OPEN. Best leads: intensity-variance
minimum for pure phase objects (±0.1 mm), and gradient kurtosis as
a focus-elevation signal in need of a better estimator.

Note for the mixed row: intensity corr 0.821 was measured at the
WRONG imag plane (-9.90) — defocused intensity still correlates
with the object, so corr alone never validates a focus plane.

## Finding 3 — new noise physics: shift is benign, tilt is malignant

Fringe period 3.2 px. Vibration uses 64-frame exposure averages.
When contrast is near 1x the argmin is noise: landing on truth then
is luck, not a measurement.

| case | imag best | err | imag contrast | pvar best | corr |
|---|---:|---:|---:|---:|---:|
| vibration shift 0.1 px | -5.00 | 0.00 | 53.3x | -5.00 | 0.951 |
| vibration shift 0.25 px | -5.00 | 0.00 | 64.1x | -5.00 | 0.952 |
| vibration shift 0.5 px | -5.00 | 0.00 | 116.5x | -5.00 | 0.953 |
| vibration shift 1.0 px | -5.00 | 0.00 | 400.8x | -5.00 | 0.951 |
| vibration shift 2.0 px | -5.00 | 0.00 | 1377.8x | -5.00 | 0.938 |
| vibration tilt 0.5 bins | -5.00 | 0.00 | 9.6x | -5.00 | 0.481 |
| vibration tilt 1.0 bins | -4.50 | 0.50 | 1.1x | -4.90 | 0.577 |
| vibration tilt 2.0 bins | -5.50 | 0.50 | 2.1x | -4.80 | 0.673 |
| vibration tilt 4.0 bins | -4.00 | 1.00 | 2.2x | -4.70 | 0.722 |
| visibility gamma=0.8 | -5.00 | 0.00 | 43.9x | -5.00 | 0.949 |
| visibility gamma=0.5 | -5.00 | 0.00 | 27.7x | -5.00 | 0.943 |
| visibility gamma=0.3 | -5.00 | 0.00 | 13.7x | -4.90 | 0.924 |
| visibility gamma=0.1 | -5.00 | 0.00 | 2.7x | -3.20 | 0.741 |
| visibility gamma=0.05 | -5.00 | 0.00 | 1.4x | -2.80 | 0.520 |
| visibility gamma=0 | -6.20 | 1.20 | 1.1x | -1.30 | 0.336 |
| blur sigma=0.5 px | -5.00 | 0.00 | 24.5x | -5.00 | 0.952 |
| blur sigma=1.0 px | -5.00 | 0.00 | 2.9x | -5.00 | 0.916 |
| blur sigma=1.5 px | -5.00 | 0.00 | 1.1x | -9.10 | 0.510 |
| blur sigma=2.0 px | -5.00 | 0.00 | 1.0x | +0.00 | 0.343 |

Reading:

- Sub-pixel shift vibration is not just benign, it IMPROVES the
  metric: contrast grows to 1377.8x at 2 px. Mechanism, verified:
  exposure averaging low-passes the sideband, and the imag-energy
  at focus falls from 0.00937 (clean) to 0.00018 (shift-averaged) —
  the blur removes the truncation ringing that sets the clean
  floor — while corr stays 0.938 because this object is mostly
  low-frequency. (Shift implementation sanity-checked: a 5 px
  shift correlates 1.000 at the +5 offset.)
- Reference-tilt jitter is the malignant vibration: averaging over
  tilt puts a Gaussian envelope on the interference term (width
  ~1/sigma_car), so the sideband is position-dependent, the window
  truncates it, and contrast collapses to <=2.2x for sigma >= 1
  bin with focus errors 0.5–1.0 mm. Corr is non-monotonic across
  tilt levels (0.48–0.72) because it is measured at a wandering
  argmin on a flat metric — do not over-read the ordering. The
  robust claim is the contrast collapse + focus error, not the
  exact corr sequence.
- Fringe-visibility threshold, two independent routes agreeing:
  focus holds to gamma ~ 0.1 (contrast 2.7x, corr 0.741) and breaks
  below it; blur holds to sigma = 1.0 px (contrast 2.9x, corr
  0.916) and breaks at 1.5 px (1.1x). A fringe transfer of ~0.1–
  0.15 is the floor: gamma = 0.1 directly, and blur sigma = 1.0
  attenuates the 3.2 px fringe by exp(-2 pi^2 sigma^2 / p^2)
  = 0.145. Below that floor the metric is flat and any plane it
  reports — including truth — is luck.
- Phase variance degrades EARLIER than imag-energy under every new
  noise type (tilt 0.5 bins: pvar contrast already 1.2x; gamma
  0.1: pvar err 1.80 mm while imag is still exact). For amplitude
  objects imag-energy remains the better metric.

## Finding 4 — carrier auto-calibration: the review's threat, closed

`estimate_carrier` locates the sideband peak in the hologram
spectrum (smoothed, DC blanked, centroid-refined) instead of
assuming the carrier.

| case | est car | imag best | err | contrast | corr |
|---|---:|---:|---:|---:|---:|
| true car 70 | 70 | -5.00 | 0.00 | 47.8x | 0.953 |
| true car 80 | 80 | -5.00 | 0.00 | 50.8x | 0.951 |
| true car 82 | 82 | -5.00 | 0.00 | 44.4x | 0.952 |
| true car 90 | 90 | -5.00 | 0.00 | 35.4x | 0.946 |
| car 82, assumed 80 | 80 | -5.00 | 0.00 | 1.2x | 0.945 |
| car 82, estimated | 82 | -5.00 | 0.00 | 44.4x | 0.952 |
| car 82 + gauss 0.20, est | 82 | -5.00 | 0.00 | 8.8x | 0.946 |

The review's collapsed case (1.2x) is fully recovered to 44.4x,
and estimation survives 0.20 gaussian noise (82 recovered, 8.8x,
corr 0.946). Practical consequence: the pipeline should never
assume the carrier again — estimate it from the data first. Note
the assumed-80 row still lands on -5.00: the carrier error never
MOVED focus, it only flattened the metric, which is why it is
dangerous under additional noise rather than visibly wrong.

## Finding 5 — the comb: every "exact focus" claim is grid luck

The imag-energy curve near focus is NOT a smooth dip. At 0.001 mm
steps around a -5.03 mm truth (reconstructed field):

| d (mm) | imag | d (mm) | imag | d (mm) | imag |
|---:|---:|---:|---:|---:|---:|
| -5.040 | 0.9084 | -5.033 | 0.5842 | -5.026 | 0.0242 |
| -5.039 | 0.2542 | -5.032 | 0.9866 | -5.025 | 0.3562 |
| -5.038 | 0.0627 | -5.031 | 0.4731 | -5.024 | 0.9589 |
| -5.037 | 0.6944 | -5.030 | 0.0102 | -5.023 | 0.6978 |
| -5.036 | 0.9602 | -5.029 | 0.4693 | -5.022 | 0.0644 |
| -5.035 | 0.3598 | -5.028 | 0.9862 | -5.021 | 0.2509 |
| -5.034 | 0.0233 | -5.027 | 0.5879 | -5.020 | 0.9064 |

Teeth every 0.004 mm, peaks ~0.99 between them, tooth depth
decaying with distance from truth (0.0102 central, 0.0233 at
±0.004, 0.063 at ±0.008). Confirmed on the TRUE field with no
hologram (symmetric comb, exactly 0.0000 at truth) — this is
physics of the metric, not a reconstruction artifact. Mechanism:
the transfer-function phase advances ~1.88 cycles per micron at
low frequencies (2 pi delta / lambda), so phase alignment
recurs quasi-periodically in propagation distance — a
Talbot-like longitudinal revival. Only exact focus aligns ALL
frequencies, so only the central tooth is deepest.

Consequences, all measured:

- Coarse 0.1 mm sweep, off-grid truth -5.03: argmin -5.30,
  err 0.27 mm (lands on the broad oscillation dip, whose minimum
  sits 0.32 mm from truth at 0.0158 — only 1.6x above the central
  tooth's 0.0102).
- Parabolic interpolation on the coarse grid makes it WORSE:
  -5.3413, err 0.31 mm. It fits the broad dip's shoulder, not the
  tooth. On a fine grid that resolves the tooth it adds nothing
  (0.01 mm grid: grid -5.03 exact, parabolic -5.0300 exact).
  Verdict: step size sharpens this minimum, interpolation does not.
- Noise does NOT smooth the comb: at gauss 0.20, 0.001 mm steps
  still give argmin -5.030 exactly (min 0.0554, median 0.4760).
- Therefore every "error 0.000 mm" claim in README/DESIGN/REVIEW
  holds only because the truth sat on the sweep grid. A blind
  sweep's honest accuracy is ~0.3 mm (broad-dip offset) unless
  the scan resolves the central tooth to sub-micron steps.
  Recommended procedure: coarse sweep to find the basin, then a
  sub-micron-step scan around each candidate dip, and report the
  tooth — never the envelope.

## Finding 6 — two planes: the bright plane wins, the dim one vanishes

Ring+cross at 5 mm, dots at 7.5 mm, 0.05 mm steps:

| | imag | pvar | corr vs that plane |
|---|---:|---:|---:|
| at plane1 (-5.0) | 0.1967 | 0.1872 | 0.838 |
| at plane2 (-7.5) | 0.7182 | 0.2209 | 0.425 |
| global imag argmin | -5.00 (median 0.480) | | |
| global pvar argmin | -5.25 (median 0.237) | | |

No local dip at -7.5 in either metric (imag 0.718 >> half-median
0.240; no pvar local minima below half-median at all). The
out-of-focus bright plane's complex background swamps the dim
plane's dip. Worse, oscillation dips at -6.55 and -3.40 (0.236,
0.237) nearly rival the true dip (0.197) — metric trust degrades
in 3D samples even for the plane that IS found. A global focus
metric cannot section a thick sample; per-region or
per-frequency metrics would be needed, and that is future work,
not claimed here.

## What changed in the code

- `wave_hologram.py`: added `estimate_carrier`,
  `phase_variance`, `parabolic_min`; `field_sharpness` kept with
  a falsification docstring.
- `discovery_sweep.py`: all six sections, reproducible (seed 23).
- `tests/test_wave_hologram.py`: 11 new tests (carrier
  estimation incl. under noise, contrast restoration, phase
  metric failures, comb teeth, coarse-miss vs fine-hit,
  sharpness invariance, two-plane no-dip). Suite: 22/22 pass.

## Consolidation (same night): first real hologram

See `real_data/REAL_DATA.md`. A real off-axis DHM hologram of a fly
head (pyDHM, MIT, λ 532 nm, 2.4 µm pitch) reconstructs to a crisp fly
head at z = 0 — propagation works on real data — but no focus metric
finds a trustworthy plane: imag-energy gives 1.16x contrast (a fly is
not a real amplitude object), gradient sharpness is flat (already
falsified above), and amplitude variance gives a shallow 10% dip at
+1.10 mm, with no ground-truth distance to check against. Two
consequences folded back into the package: `field_sharpness`'s
docstring already carries the falsification, and `amplitude_variance`
was added as a phase-object diagnostic only (minimum at focus, shallow
dip, never use alone), with a test (suite now 23/23). The JoVE DLHM
images turned out not to be hosted on OSF at all — only the index CSV
is there — so that route is closed unless the authors release them.
Best next route: a calibration-target hologram with a stated distance,
or a new recording with the distance logged.

## v0.1.2: calibration-target result (2026-10-25)

Cuenat et al. tiny-DHM dataset (Zenodo 6337535): amplitude-variance
basin vs logged hexapod Z gives slope -93.7 um per um (M^2 scaling,
true magnification ~9.68, not nominal 10x), residual RMS 0.047 mm.
Full tables in `real_data/dhm_tiny/FINDINGS.md`; recording protocol in
`real_data/RECORDING_PROTOCOL.md`. The 2.2 GB of extracted holograms
under `real_data/dhm_tiny/extracted/` are NOT in the release zip —
fetch them from the Zenodo record.
