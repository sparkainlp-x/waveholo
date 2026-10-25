# wave-hologram v0.1.3: what changed and what we now know

v0.1.2 measured focus with a metric that turned out to depend on an
arbitrary constant, the phase of the reference beam. v0.1.3 replaces it
with a metric that ignores that constant, estimates the carrier to a
fraction of a bin, refines focus off the grid, and runs the whole
pipeline on **530 real microscope holograms** with logged stage
positions. All numbers below come from `studies.py`
(`results/studies_output.txt`) and `real_data/dhm_v013.py`
(`real_data/dhm_v013/`) on this machine. The test suite has 43 tests,
all passing.

## 1. Headline: it works on a real microscope

The Cuenat et al. dataset (Zenodo 6337535, CC-BY-4.0, arXiv:2203.07772)
contains holograms from a Lyncée Tec DHM with a 10× objective,
λ = 674.99 nm, 5.86 µm pixels, imaging a phase target at hexapod
stage positions logged to the nanometre. We downloaded the first 2.2 GB
of the 48.6 GB archive: 2,030 holograms covering 6 stage positions
(Z = −50.0, −49.0, −40.1, −39.1, −38.1, −37.1 µm). Then we ran 100
holograms per position (30 for the last, partial one) through the
pipeline. The objective's phase curvature was removed with the new
`remove_quadratic_phase`, and focus was found blind on a −20…+4 mm
sweep.

| metric | slope (µm focus / µm stage) | R² | class-mean residual (µm stage) | single-hologram SD (µm stage) |
|---|---:|---:|---:|---:|
| **realness** (curvature removed) | **−96.4 ± 1.3** | 0.909 | **0.19** | 1.53 |
| ac_realness (curvature removed) | −108.4 ± 1.7 | 0.886 | 0.26 | 1.83 |
| amplitude_variance (raw) | −85.2 ± 1.2 | 0.900 | 0.65 | 1.46 |
| tamura (raw) | fails, sticks to sweep edge | 0.06 | — | — |

What this means:

- **The slope matches the optics.** A 10× objective magnifies depth by
  M² ≈ 100. Realness gives 96.4 ± 1.3. This is a measured-vs-truth
  test on hardware data, not a simulation.
- **Realness is the most linear metric.** Averaged per stage position,
  it sits 0.19 µm from a straight line, which is 3.4× tighter than
  amplitude variance (0.65 µm). The stage steps here are about 1 µm.
- **Single holograms scatter about 1.5 µm, and most of that is
  probably real motion.** Realness and amplitude variance work on
  different fields (curvature-removed vs raw) and use different
  physics. Yet their per-hologram deviations correlate at 0.67, and
  neither correlates with carrier jitter (−0.08, −0.03). A shared
  wobble seen by two independent probes is more likely the sample or
  stage than the metric. This is an inference: the dataset has no
  independent per-shot position, so it cannot be proven here.
- **The carrier moves between shots.** At a fixed stage position the
  carrier jitters by ±0.19 bins (rows) and ±0.13 bins (columns).
  v0.1.2's integer-bin carrier could not see this. The v0.1.3
  estimator tracks it per hologram.
- **Without curvature removal, realness is blind.** On the raw field it
  sits at 0.43–0.45 on every plane, because no field carrying about
  78 rad of lens curvature is ever real. After removal it drops to
  0.11 at focus, with 1.4× contrast.

Figures: `real_data/dhm_v013/dhm_focus_vs_stage.png` (all 530
holograms, four metrics) and `dhm_mean_curves.png` (class-mean focus
curves). Raw numbers: `dhm_per_hologram.csv`, `dhm_summary.txt`,
`class_mean_fit.txt`.

## 2. The v0.1.2 "comb" was an aliasing artefact (S2)

v0.1.2 reported sub-micron focus "teeth" every ~4 µm. They are not
optics. `imag_energy` has a fixed reference phase, so it oscillates
with the global phase of the propagated field, which advances by 2π
every wavelength. Squared, that gives a period of λ/2 = 0.266 µm.

- On a 25 nm grid, the measured minima spacing is 0.275 µm.
- On v0.1.2's 1 µm grid, a 0.266 µm period aliases to
  1/|1/0.266 − 4| = **4.16 µm**. v0.1.2 reported "teeth every 4 µm".
- Realness over the same 20 µm varies by 0.45%, which is smooth.

Figure: `results/s2_no_comb.png`.

## 3. Why imag_energy had to go: the reference piston (S1)

Every real setup has an unknown constant phase between the object and
reference beams. v0.1.2 always simulated it as 0. Sweeping it:

- `imag_energy` at the true focus ranges from 0.009 to 0.991 (it
  follows sin² of the piston). The focus error reaches **0.79 mm**.
- `realness = |Σ U²| / Σ |U|²` is exactly invariant: multiplying U by
  e^{iφ} only rotates Σ U². Its value at truth stays within
  0.0088–0.0106, and its worst focus error is 0.040 mm.

With a new random piston in every trial (S6), imag_energy lands within
0.1 mm of the truth in only 0–2 of 12 trials at *any* noise level,
including the cleanest. Realness and ac_realness land 12 of 12 up to
0.35× gaussian noise and down to 100 photons per pixel.

| noise | realness within 0.1 mm | imag_energy within 0.1 mm |
|---|---:|---:|
| gaussian 0.1–0.35 | 12/12 | 2/12 |
| gaussian 0.5 | 11/12 | 2/12 |
| gaussian 0.8 | 9/12 | 1/12 |
| gaussian 1.0 | 6/12 | 1/12 |
| gaussian 2.0 | 2/12 | 1/12 |
| poisson 100 photons | 12/12 | 2/12 |
| poisson 30 photons | 9/12 | 1/12 |

v0.1.2 claimed imag_energy holds "to 2× noise". That was true only for
piston = 0, a value no real instrument guarantees. Figures:
`results/s1_piston.png`, `results/s6_noise_random_piston.png`.

## 4. Sub-bin carrier estimation (S3)

The carrier is now located with a two-pass local DTFT around the
sideband peak, with 1-D and 2-D versions.

- Max error is 0.030 bins on a 1-D fractional sweep. Across 20 random
  2-D carriers the median error is 0.0034 bins and the max is 0.0114.
- Rounding a fractional carrier to an integer bin costs imag_energy up
  to 0.79 mm. Realness is much less sensitive (0.043 mm with integer
  rounding, 0.019 mm with sub-bin).

## 5. The remaining bias is sideband cross-talk, not the metric (S4)

On the true sensor field, the focus error is exactly 0. On a hologram,
the twin image and the |U|² term leak broadband tails into the sideband
window.

- Sharp binary object: −50 to +64 µm, depending on carrier, window and
  piston.
- Band-limited object (spectrum radius < 24 bins): ≤ 1.6 µm at every
  window tested.

So the residual error is set by how much spectrum overlaps between
sidebands, which is a recording-geometry problem. The fix is a higher
carrier, a band-limited object or a phase-shifted twin removal, not a
better metric. Figure: `results/s4_bias_vs_window.png`.

## 6. Off-grid focus (S5)

`find_focus` runs a coarse sweep and then a golden-section refinement
on the continuous metric, using about 125 metric evaluations. For truths
placed off-grid (−3.300, −5.030, −6.137, −7.770, −9.100 mm, clean and
at 0.2× noise), the errors are 2–49 µm. imag_energy on the bare 0.1 mm
grid errs by 0.04–0.63 mm.

## 7. Phase and mixed objects (S7)

There are 12 synthetic objects in 7 families, at a piston of 0.7 rad.
The count of objects found within 0.1 mm:

| metric | found |
|---|---:|
| **ac_realness** (realness of the mean-subtracted field) | **9/12** |
| realness | 3/12 |
| tamura | 3/12 |
| amplitude_variance | 2/12 |
| phase_variance | 2/12 |
| imag_energy | 0/12 |

ac_realness solves binary phase objects (π/2, π) and single-material
amplitude+phase objects, which defeated every v0.1.2 metric. It still
fails on smooth phase objects and on two materials with separate
supports. There, amplitude variance does better (about 0.26 mm), but
nothing is reliable. No single metric wins on every object class.

## 8. Depth maps for multi-plane samples (S8)

The test sample has a ring and cross at 5.0 mm and three dots at
7.5 mm. One global focus reports −5.09 mm, a single number for two
planes. The new `local_focus_map` (per-tile realness focus) recovers
both: dot tiles read 7.5–7.7 mm and ring tiles 5.0–5.1 mm, with one
mixed edge tile at 6.6 mm. Figure: `results/s8_two_plane_depth_map.png`.

## 9. Autofocus is frequency estimation (S9)

Σ|U_z|² is constant with distance (Parseval), and Σ U_z² = Σ_k F(k)
F(−k) e^{2i z k_z(k)}. So the realness curve depends only on the
**pair spectrum** P(k) = F(k)F(−k). Paraxially, it is the magnitude of
a 1-D Fourier transform of P over |f|², at "frequency" λz.

- Numerically verified: the max difference from direct propagation over
  101 planes is 2.2 × 10⁻¹⁶.
- The pair-spectrum route is 2.8× faster (0.29 s vs 0.80 s).
- Negative result: weighting P(k) by |P|^γ, by |f|^α, or high-passing
  it did not improve noise robustness. Flat weighting was best.

## 10. Real fly hologram, no ground truth (S10)

pyDHM `fly-20x.bmp`, 532 nm, 2.4 µm pixels:

| metric | focus (mm) | contrast |
|---|---:|---:|
| ac_realness | +1.057 | 1.108 |
| amplitude_variance | +1.165 | 1.195 |
| tamura | +1.341 | 1.045 |
| realness | −2.245 | 1.013 (flat) |
| imag_energy | −7.873 | 2.648 |

Three independent metrics agree within 0.3 mm near +1.1 mm. imag_energy
has the highest contrast, but at −7.9 mm, which is consistent with the
piston artefact above. This is agreement between metrics, not accuracy:
the fly has no logged truth. Figures: `results/s10_fly_refocus.png`,
`results/s10_fly_metrics.png`.

## Honest limits

- **Only 6 of 100 stage positions** were used: 2.2 GB of the 48.6 GB
  archive at about 1 MB/s. They span Z = −50 to −37 µm. Linearity over
  the full −50…+42 µm range is untested.
- **The absolute focus is not tested, only the slope.** The intercept
  absorbs where the DHM image plane sits. The three metrics' slopes
  differ (85, 96, 108) by more than their standard errors. Removing the
  quadratic phase acts like a thin lens, so the curvature-removed
  metrics live in a rescaled coordinate. "Which slope is the true M²"
  needs the objective's calibrated magnification, which we don't have.
- **The DHM target is periodic (9 µm).** Talbot self-images could
  create replica focus planes. The sweep window was fixed (−20…+4 mm)
  and replicas were not investigated.
- **Class index ±1.** Directory class_0 holds files named class_1_…
  Directory indices were used. If the filenames are right, every Z
  shifts by about +1 µm. The slope is unaffected.
- Synthetic studies use one 256² grid, one wavelength and up to 12
  seeds per cell. Treat the fractions as indicative, not as precise
  rates.

## Reproduce

```
python3 -m unittest discover -s tests          # 43 tests
python3 studies.py                             # S1–S10, ~8 min, writes results/
python3 real_data/dhm_v013.py --data <dir with class_* folders> --per-class 100
```

The DHM data comes from
https://zenodo.org/api/records/6337535/files/dataset.tar.gz/content. A
byte-range download of the first ~2.2 GB, extracted, gives the classes
used here.
