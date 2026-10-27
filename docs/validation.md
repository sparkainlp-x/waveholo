# Validation: what is checked, and what is not

This page lists what the test suite compares with an independent reference
(a closed-form result, an exact symmetry, a conservation law, or a number
quoted in the README), and what has no such check. "Agreement" below is the
worst-case deviation measured when this page was written (Python 3.13,
numpy 2.5.4); the tests use looser bounds. `tests/test_physics_validation.py`
was also run locally on Python 3.9 to 3.14 (numpy 2.0.2 to 2.5.4, whatever each
local environment had installed), and the full suite runs in CI on 3.9 to 3.14.

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
- **Off-axis separation.** Only the Nyquist alias of the fringe frequency is
  tested. Whether the sideband separates from the zero-order term for a given
  object bandwidth and carrier is not tested as a rule.
- **Input checks on `dx` and `wl`.** A negative pixel pitch or wavelength is
  accepted without an error; zero raises `ZeroDivisionError`. No test pins this
  behaviour.
- **Phase objects and two-material scenes.** The README lists these as unsolved;
  the tests record the limit rather than a fix.
- **Browser sandbox.** `docs/sandbox/` shows an illustrative curve; the browser
  test checks that it runs, not that it measures anything. It is skipped when
  the browser binary is missing from the unit-test job.
- **Dependencies.** `pyproject.toml` gives lower bounds only (`numpy>=1.24.0`)
  and there is no lock file, so CI installs the newest releases each run.
  Floating-point results are not guaranteed to be bit-identical across numpy
  versions or platforms; the tests use tolerances.
