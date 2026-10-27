# Validation follow-up: inputs, sideband separation, reproducibility

This page adds to `docs/validation.md` (pull request #8). It is a separate file
so that it does not conflict with that pull request or with the README and
CHANGELOG edits in #4 and #5. After #8 is merged, the three sections below can be
folded into `docs/validation.md` and the matching lines removed from its
"Not validated" list (input checks on `dx` and `wl`, off-axis separation,
dependencies).

Numbers were measured with Python 3.13 and numpy 2.5.4 unless stated, from the
tests named in each section.

## 1. Input validation (`tests/test_input_validation.py`)

Public functions now check their numeric arguments and raise `ValueError`
(bad value, shape or range) or `TypeError` (not a real number, or a complex
hologram) with a message that names the argument. Before this change the same
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

## 2. Sideband separation rule (`tests/test_sideband_separation.py`)

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

## 3. Reproducibility

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
  Before this change the job installed the newest Playwright (1.64.1 in the last
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
- Regenerate after an intentional change to the pipeline:

```python
import hashlib, numpy as np, wave_hologram as W
obj = W.make_object(); sensor = W.angular_spectrum(obj.astype(complex), 5e-3)
holo = W.record_hologram(sensor, car=80.37, ref_phase=1.2)
print(hashlib.sha256(np.round(holo * 1e6).astype('<i8').tobytes()).hexdigest())
```

## Still not validated

- Real microscope data, absolute focus on real data, and the cross-talk bias
  figures (see `docs/validation.md`).
- Vector and sensor effects; padding workflows; the metrics on NaN input
  (returned as numbers, see section 1).
- The separation rule with noise, with fractional or 2-D carriers, and for
  objects that are not exactly band-limited.
- Bit-identical results on other CPUs or operating systems.
- Whether the Pyodide browser test passes with Playwright 1.44.0 in the container
  is reported in the pull request (CI), not assumed here.
