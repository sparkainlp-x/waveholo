# Validation: what is checked, and what is not

This page lists what the test suite compares with an independent reference
(a closed-form result, an exact symmetry, a conservation law, or a number
quoted in the README), and what has no such check. "Agreement" below is the
worst-case deviation measured when this page was written (Python 3.13,
numpy 2.5.4); the tests use looser bounds. `tests/test_physics_validation.py`
was also run locally on Python 3.9 to 3.14 (numpy 2.0.2 to 2.5.4, whatever each
local environment had installed), and the full suite runs in CI on 3.9 to 3.14
with the environment pinned in `uv.lock` (see Reproducibility below).

Line coverage (`pyproject.toml`, `fail_under = 85`) says which lines run, not
whether results are right. The tests listed here are the ones that check values.

## Units and conventions

- All distances and the wavelength are in **metres**; the pixel pitch `dx` is a
  scalar or a `(dy, dx)` pair in metres. Defaults: `DX = 6.5e-6`, `WL = 532e-9`,
  `N = 256`.
- `angular_spectrum(field, z)` multiplies the spectrum by `exp(+i kz z)`, with
  `kz = (2 pi / wl) sqrt(1 - (wl fx)^2 - (wl fy)^2)`. A positive `z` moves a
  wave with positive `fx` toward +x. Evanescent components are set to zero.
- Spatial frequencies, carriers and window sizes are in **FFT bins** (cycles
  per field of view), counted from DC.
- The field is **periodic**: the FFT grid wraps around, and no padding is added.

## Validated against a reference

All in `tests/test_physics_validation.py`.

| Check | Reference | Agreement |
|---|---|---|
| On-grid plane wave | `exp(i kz z)` phase, exact | 1.9e-15 (max abs) |
| Uniform field | free-space phase `exp(i k z)` | 0 (max abs, 16 x 16 field) |
| Two steps vs one | `P(z2) P(z1) = P(z1 + z2)` | 1.0e-11 |
| Energy, all waves propagating | `sum abs(U)^2` unchanged | 4e-16 (relative) |
| Realness of a real object | even in `z` (exact symmetry) | 5.6e-17 |
| Gaussian beam, `z` = 0.5, 1, 2 `zR` | `w(z) = w0 sqrt(1 + (z/zR)^2)`, and axial peak `(w0/w)^2` | width 7e-7, peak 9e-7 (relative) |
| Talbot self-imaging (cosine grating, 16 px period) | image repeats at `zT = wl / (1 - sqrt(1 - (wl/p)^2))`, half period shifted at `zT/2` | 5e-15 and 6e-15 (intensity) |
| Talbot distance vs the paraxial `2 p^2 / wl` | 40.661 mm exact vs paraxial | differ by 6.5e-6 (relative) |
| Circular aperture, axial intensity | Fresnel `4 sin^2(pi N_F / 2)`, `N_F = a^2 / (wl z)`; `a` from the pixel count (19.91 px) | within 0.022 absolute at `N_F` = 0.5 to 3 (values of 0 to 4); 2.3e-5 at `N_F` = 2 |
| Lateral shift of a tilted packet | `z tan(theta)`, `sin(theta) = wl fx` | 2e-4 px at most, including a shift of 177 px that wraps around the grid |
| Wrapped packet, back-propagated | recovers the input (unitary on the grid) | 4.5e-16 |
| Hologram of one tone | Fourier coefficients: DC `1 + a^2`, two fringe terms of size `a`, nothing else | below 1e-15 |
| Sideband demodulation of a tone | returns `U exp(-i piston)` | 7e-15 |
| Fringe frequency above Nyquist | records exactly the same hologram as its alias | 0 |
| Carrier at or above Nyquist, or at 0 | `ValueError` | exact |
| README quickstart | focus -4.998 mm, carrier within 0.03 bins, correlation 0.956 | matches to the printed digits |
| Repeatability in one process | same arrays, bit for bit | exact |

The Fresnel and Gaussian rows are paraxial results compared with a non-paraxial
propagator. The Fresnel comparison depends on the pixelated aperture radius, so
its tolerance is much looser than the other rows.

Older tests check, among other things: the round trip, `imag_energy` oscillating
with period wl/2 in `z`, realness invariance to the reference piston, the
pair-spectrum identity against direct propagation, sub-bin carrier estimation,
and seeded noise cases (`tests/test_wave_hologram.py`, `tests/test_v015_coverage.py`).

## Reproduced in this audit

- The recorded S6 noise result for `realness` (`results/studies_output.txt`):
  Gaussian 0.35 and Poisson 100 photons, 12 seeded trials each. Reproduced:
  median error 0.015 mm and 0.014 mm, both 12/12 within 0.1 mm.
- The game table in the README matches the data stored in `docs/game/index.html`,
  and `scripts/make_game_frames.py --check` (run in CI) confirms the page equals
  a regeneration from fixed seeds.

## Input validation (`tests/test_input_validation.py`)

Public functions now check their numeric arguments and raise `ValueError`
(bad value, shape or range) or `TypeError` (not a real number, or a complex
hologram) with a message that names the argument. Previously the same
inputs gave a `ZeroDivisionError`, a bare numpy message, an `OverflowError`, or
a silent NaN-filled or wrong result.

| Argument | Rule |
|---|---|
| `dx` (all propagation paths, `cropped_pitch`, `record_hologram`) | positive finite scalar, or a `(dy, dx)` pair of them |
| `wl` | positive finite scalar |
| `z`, entries of `sweep` | finite real scalars; `sweep` is 1-D, at least 2 entries for `focus_sweep`, at least 1 for `find_focus` and `local_focus_map` |
| fields (`field`, `field0`, `sensor_field`) | 2-D, non-empty, finite |
| holograms | 2-D, non-empty, finite, real (a complex array with zero imaginary part is accepted) |
| `car` | finite real, `0 < car < n/2` (unchanged range) |
| `carrier_y` | finite real, `abs(carrier_y) < ny/2` |
| `ref_amp` | positive finite; `ref_phase` finite real |
| `carrier`, `hint` | pair of finite reals |
| `win` | finite scalar or `(wy, wx)` pair (existing range checks kept) |
| `dc_radius` | finite, non-negative; `hint_radius` positive |
| `tol`, `max_iter` | positive; `max_iter` an integer; `golden_min` also needs `a <= b` |
| `local_focus_map`: `tile`, `min_energy_frac` | integer in `1..min(ny, nx)`; value in `[0, 1]` |
| `remove_quadratic_phase`: `iters`, field size | positive integer; at least 3 x 3 |
| `run`: `z_obj`, `sweep_min`, `sweep_max`, `frames`, `fps` | finite; `z_obj`, `fps` positive; `frames` an integer (existing `frames >= 2` kept); nothing is created on disk before these pass |
| `make_object(n)` | integer; `n < 122` raises (the cross does not fit); `122 <= n < 160` returns the clipped object and emits a `UserWarning` |

Unchanged on purpose: the metric functions (`realness`, `tamura`, and the rest)
still return a number for a NaN or empty array, and `realness` returns 0.5 for a
field with zero or non-finite energy; the fields they receive from the
propagators are already checked. `find_focus` accepts a one-point sweep.
`estimate_carrier_2d(hint_radius=...)` is ignored when `hint` is not given.

Valid input behaves as before: the 126 earlier tests, the README quickstart
(-4.998 mm), the game data check (`scripts/make_game_frames.py --check`) and the
golden tests below all pass unchanged.

## Sideband separation rule (`tests/test_sideband_separation.py`)

For an object field with spectrum confined to `|kx|, |ky| <= B` bins, an x-carrier
of `car` bins and a window of half-width `B`:

- the zero-order terms occupy `|k| <= 2B` (`|U|^2` is the autocorrelation of the
  spectrum of `U`);
- the wanted sideband occupies `car - B .. car + B`;
- the two are disjoint when `car >= 3B + 1` (at `car = 3B` they share one bin);
- the sideband stays below Nyquist when `car <= n/2 - B - 1`.

Measured with exactly band-limited complex fields, fixed seeds, `n = 128`,
integer carriers, relative L2 error of the reconstructed field against
`U * exp(-i * piston)`:

| B | car = 3B - 1 | car = 3B | car = 3B + 1 | car = n/2 - B - 1 | car = n/2 - B |
|---|---|---|---|---|---|
| 4 | 3.2e-2 | 2.4e-2 | 3.8e-15 | 4.6e-15 | 3.3e-1 |
| 6 | 2.1e-2 | 1.6e-2 | 4.8e-15 | 4.9e-15 | 3.3e-1 |
| 10 | 1.3e-2 | 6.7e-3 | 5.5e-15 | 5.5e-15 | 2.2e-1 |

For B = 8 (seed 2026): 6.1 at `car` = 4 or 8 (the twin image overlaps), 8.2e-2 at
12, 6.1e-2 at 16, 3.7e-2 at 20, 1.7e-2 at 23, 7.0e-3 at 24, at most 1.9e-14 from 25
(= 3B + 1) through 55, and 1.9e-1 at 56 (= n/2 - B). The error decreases with
`car` over 12, 16, 20, 23, 24, 25 for each of seeds 0 to 9 (tested for 0 to 4).
The window matters as much as the object: at `car = 25`, `win` = 7 gives
4.6e-1 and `win` = 9 gives 7.8e-3.

Focus estimate, for real positive band-limited amplitude objects (B = 8,
propagated 5 mm, `find_focus` on the default sweep, `win` = 8), four seeds: the
error is below 1e-12 um at `car` = 25, 40 and 55, and between 1.4 and 187 um at
`car` = 12, 16, 20 and 24 (not monotonic in `car`; the leak is a random
interference term).

What this does not show: it is exact only for exactly band-limited objects and
integer carriers. Real objects are not band-limited (sharp edges), so the rule
gives a design guideline, not a pass/fail limit. The default demo geometry
(`CAR_DEFAULT` 80, `WIN_DEFAULT` 48, `N` 256) cannot satisfy it for an object that
fills the window (`145 <= car <= 79` has no solution); the demo object is not
band-limited to 48 bins, and its measured residual bias is in the README and
`RESULTS_v0.1.3.md`. Noise, fractional carriers and 2-D carriers are not
covered by this rule's tests. Fractional and 2-D carriers are covered by the
estimator tests in `tests/test_wave_hologram.py`.

## Reproducibility

### Pinned environments

- **Matrix jobs (Python 3.9 to 3.14):** `uv.lock` pins every package, including
  the optional extras, for each Python version (for example numpy 2.0.2 on 3.9,
  2.2.6 on 3.10, 2.4.6 on 3.11, 2.5.4 on 3.12 to 3.14). CI runs
  `uv sync --locked`, which fails if `uv.lock` does not match `pyproject.toml`.
  A version bump in `pyproject.toml` therefore needs `uv lock` in the same change.
- **Pyodide/WASM job (container, Python 3.10):** `requirements-ci.txt`, generated
  from `requirements-ci.in` with
  `uv pip compile requirements-ci.in --python-version 3.10 --python-platform x86_64-manylinux_2_28 -o requirements-ci.txt`.
  It pins `playwright==1.44.0` to match the browsers in the image
  `mcr.microsoft.com/playwright/python:v1.44.0-jammy`.
  Previously the job installed the newest Playwright (1.64.1 in the last
  run on `main`); the image's browsers did not match it, and the browser test was
  skipped with "Playwright browser binary unavailable", so the job passed
  without running the test.
- Not pinned: the build backend (`setuptools>=77`, fetched in pip's build
  isolation), the runner image (`ubuntu-latest`), and uv itself.
- Regenerate after changing dependencies: `uv lock`, then the `uv pip compile`
  command above.

### What is claimed about bit-identical results

Bit-identical results are claimed only for repeated runs in one process, and
within one pinned environment on one machine. Measured on one x86-64 machine
across Python 3.9 to 3.14 and numpy 2.0.2 to 2.5.4: the raw bytes of the sensor
field and of the hologram in the README quickstart were identical in all seven
environments tested (SHA-256 prefixes 4875f93b7c67 and f4cb780a12c6). The
realness curve differed by at most 1.1e-16 between numpy versions, and `value`
and `contrast` from `find_focus` differ in the last digit between numpy 2.0.2 and
2.5.4. No claim is made for other CPUs, operating systems, or other BLAS/FFT
builds, which this audit did not test.

### Golden tests (`tests/test_reproducibility.py`)

- Scalars (focus, carrier, five curve points) are compared with a tolerance
  (1e-9 for values, 1e-12 m for the focus distance), not exactly, because of the
  last-digit differences above.
- Arrays are compared by SHA-256 of the intensity rounded to 1e-6 and stored as
  little-endian int64. The hashes were identical in the six Python versions tested (3.9 to 3.14) and in a
  seventh environment (3.13.5, numpy 2.5.4).
  The nearest value to a rounding boundary is 3.4e-6 quanta away for the clean
  hologram and 1.3e-5 for the noisy one (seed 1000), against a floating-point
  spread of about 1e-10 quanta. If a platform changes a hash, the tolerance tests
  show whether the physics changed.
- Regenerate after an intentional change to the pipeline (snippet below):

```python
import hashlib, numpy as np, wave_hologram as W
obj = W.make_object(); sensor = W.angular_spectrum(obj.astype(complex), 5e-3)
holo = W.record_hologram(sensor, car=80.37, ref_phase=1.2)
print(hashlib.sha256(np.round(holo * 1e6).astype('<i8').tobytes()).hexdigest())
```

## Not validated

- **Real microscope data.** The result "530 holograms, 6 stage positions, slope
  -96.4 +/- 1.3 um/um" comes from `real_data/dhm_v013.py` on the Zenodo 6337535
  dataset (a 48.6 GB archive; the script uses the first ~2 GB). It is not in
  the repository, not run in CI and not re-run here. It is a slope check; the
  absolute focus distance on real data is not validated.
- **Bias from sideband cross-talk.** Quoted from the recorded studies
  (synthetic objects only); not re-run here and not asserted by a test beyond the
  documented 30 um tolerance in the existing tests.
- **Physics beyond scalar, paraxial-quality diffraction.** There is no vector or
  polarisation model, no sensor model (fill factor, pixel response, saturation),
  and no comparison with a measured diffraction pattern.
- **Boundary handling.** The propagator is periodic, so a distance with
  `z tan(theta) > N dx / 2` wraps. For content up to the Nyquist frequency this
  is `z > N dx^2 / wl` (20.3 mm for the default grid). Padding is the caller's
  job; there is no test that a padded workflow gives the unwrapped result.
- **Off-axis separation, beyond the tested case.** The rule above is tested
  with exactly band-limited objects, integer x-carriers and no noise. It is not
  tested with noise, fractional or 2-D carriers, or objects that are not
  band-limited (sharp edges), so it is a design guideline for real data.
- **Metrics on invalid arrays.** `realness`, `tamura` and the other metric
  functions still return a number for NaN or empty arrays (`realness` returns
  0.5 for zero or non-finite energy). The propagators that feed them reject such
  input.
- **Phase objects and two-material scenes.** The README lists these as unsolved;
  the tests record the limit rather than a fix.
- **Browser sandbox.** `docs/sandbox/` shows an illustrative curve; the browser
  test checks that it runs, not that it measures anything. It is skipped in the
  unit-test jobs (no browser binary). The container job runs it with
  Playwright 1.44.0, pinned to match the image; before that pin the job skipped
  it.
- **Bit-identical results on other CPUs or operating systems.** Not tested; see
  Reproducibility. Not pinned: the build backend, the runner image and uv itself.
