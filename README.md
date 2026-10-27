# WaveHolo (`wave-hologram`)

[![CI](https://github.com/sparkainlp-x/waveholo/actions/workflows/ci.yml/badge.svg)](https://github.com/sparkainlp-x/waveholo/actions/workflows/ci.yml)

Real digital holography, end to end: record an off-axis hologram,
reconstruct the object by numerical back-propagation, and **estimate**
the focus distance from the complex wavefield instead of assuming it.
The estimate carries a documented error (see "Honest limits"). As of
v0.1.3 it has been **checked** against logged stage positions on 530
real microscope holograms taken at 6 of 100 stage positions. That is a
slope check, not an absolute-focus test.

The name is WaveHolo; the Python module, package and Zenodo record are
called `wave-hologram` / `wave_hologram`.

The released v0.1.5 (Zenodo DOI 10.5281/zenodo.23653721) added a vectorized
pair-spectrum solver for realness sweeps, a browser-only Pyodide/NetCDF
research sandbox at [`docs/sandbox/`](docs/sandbox/), and a real-data benchmark
harness at [`real_data/live_dhm_benchmark.py`](real_data/live_dhm_benchmark.py).
Changes after the v0.1.5 tag are listed under "Unreleased" in the
[CHANGELOG](CHANGELOG.md).

## Install

`wave_hologram.py` is a single flat module. **No release has been published
on PyPI**: do not expect `pip install wave-hologram` to give you this project.
(`.github/workflows/publish.yml` exists but has never been run for a
release.) Use it from a clone:

```
git clone https://github.com/sparkainlp-x/waveholo.git
cd waveholo
python3 -m venv .venv && . .venv/bin/activate
pip install -e .                 # numpy + matplotlib
pip install -e ".[dhm]"          # adds Pillow, for real_data/dhm_v013.py
pip install -e ".[notebook]"     # adds jupyter, to run examples/waveholo_demo.ipynb
```

The tests pass on Python 3.9 to 3.14 (checked locally and in CI for 3.9 to
3.14). Plain `python3 wave_hologram.py` also works without installing.

## Quickstart

```python
import numpy as np
import wave_hologram as W

obj = W.make_object()                                    # synthetic ring, cross, dots
sensor = W.angular_spectrum(obj.astype(complex), 5.0e-3) # 5 mm to the sensor
holo = W.record_hologram(sensor, car=80.37, ref_phase=1.2)   # unknown piston, fractional carrier

carrier = W.estimate_carrier_2d(holo, dc_radius=30.0 / W.N)  # sub-bin, from the hologram
field0 = W.reconstruct_field(holo, carrier=carrier)
res = W.find_focus(field0, np.linspace(-10e-3, 0.0, 101))    # metric="realness"
print(res["z"] * 1e3, "mm")      # -4.998 (truth -5.000)
```

This is `examples/quickstart.py` (`python3 examples/quickstart.py` also writes a
figure). `examples/waveholo_demo.ipynb` is the same walk-through as an executed
notebook with plots, plus one phase-object scene. Both are synthetic.

## Run

```
python3 wave_hologram.py                       # demo: video + still + curve
python3 wave_hologram.py --ref-phase 1.2 --carrier 80.37   # realistic: unknown
                                               # piston, fractional carrier
python3 wave_hologram.py --metric ac_realness  # phase / mixed objects
python3 wave_hologram.py --no-video
python3 studies.py                             # S1–S10 research studies (~8 min)
python3 real_data/dhm_v013.py --data <dir>     # real DHM validation
python3 -m unittest discover -s tests -v       # core, game, NetCDF, and browser tests
```

Requires numpy and matplotlib. ffmpeg is needed for the MP4 and Pillow
for the DHM TIFFs.

For the complete development environment, install `pip install -e ".[test]"`.
The browser test requires Playwright's Chromium runtime (`playwright install
chromium`). CI provisions it in the Playwright container. The sandbox loads
Pyodide and Plotly from their public CDNs and performs all numerical work in
the browser; uploaded NetCDF-3 focus curves are decoded locally. The
sandbox's default cell plots a hand-written analytic curve (minimum at
−5.0 mm by construction). It is a plotting demo, not output of
`wave_hologram.py`, and it says nothing about accuracy.

## What v0.1.3 measures

- **Real microscope (Cuenat et al. DHM, 10×, 530 holograms, 6 stage
  positions):** focus tracks the stage with slope −96.4 ± 1.3 µm/µm
  (standard error of a per-hologram fit that treats the 530 holograms as
  independent; they sit at only 6 distinct stage positions). A fit to the
  6 class means gives −96.3 µm/µm with a 95% interval of about ±5 (t,
  4 degrees of freedom). That is consistent with M² ≈ 100 for a nominal
  10× objective, but the other two metrics give −85 and −108, and the
  objective's calibrated magnification is not available. Class means lie
  0.19 µm (in stage units) from the line. The focus curves are shallow
  (median/min contrast about 1.4×).
- **Synthetic demo** (fractional carrier 80.37, reference piston
  1.2 rad, both unknown to the pipeline): truth −5.000 mm, measured
  −4.998 mm (2 µm error), recovery correlation 0.956. This is one clean
  scene; sharp objects carry a bias of −50 to +64 µm (see "Honest limits").
- **Noise** (synthetic ring/cross/dots scene, random piston in every
  trial, 12 trials per level, realness metric): 12/12 within 0.1 mm up to
  0.35× gaussian noise and down to 100 *peak* photons per pixel (about 17
  on average for this hologram); 9/12 at 30 peak photons. 12/12 is a count,
  not a guarantee (Wilson 95% lower bound for 12/12 is 76%). The v0.1.2
  metric (the project's own earlier one, not prior art) manages 0–2/12 at
  every level.

The full write-up, with all tables and figures, is
**[RESULTS_v0.1.3.md](RESULTS_v0.1.3.md)**. Changes are listed in
[CHANGELOG.md](CHANGELOG.md).

## Related work

Focus criteria for numerical refocusing in digital holography are an
established field, and WaveHolo's metrics sit inside it. They were not
benchmarked against the following here:

- Dubois et al., Opt. Express 14, 5895 (2006): integrated amplitude, with a
  minimum for amplitude objects and a maximum for phase objects; and
  Dubois et al., Opt. Lett. 39, 4286 (2014): a high-pass filtered variant
  that gives a minimum for both. `ac_realness` (mean-subtracted field) is
  conceptually close to the second.
- Langehanenberg et al., Appl. Opt. 47, D176 (2008): Fourier-weighted and
  gradient-based criteria for pure phase objects.
- Memmolo et al., Opt. Lett. 36, 1945 (2011) and 39, 4719 (2014): Tamura and
  sparsity criteria. Zhang et al., Opt. Lett. 42, 3824 (2017): sparsity of the
  gradient, which removes the polarity problem. Tamamitsu et al. (2020):
  Gini index versus Tamura coefficient.

`realness` is also a criterion that is minimised over distance, not a direct
measurement of position. Whether it is new has not been verified. It has
been compared only with the project's own earlier metric, a plain Tamura
coefficient and intensity variance.

## Default pipeline

1. Sub-bin carrier from the hologram (`estimate_carrier_2d`, error
   ≤ 0.03 bins).
2. Sideband filter, optionally cropped losslessly (`reconstruct_field`).
3. For objective-based DHM, remove the lens curvature first
   (`remove_quadratic_phase`).
4. `find_focus`: coarse sweep, then golden-section refinement of
   **realness** = |Σ U²| / Σ |U|². This metric is invariant to the
   unknown reference phase. Use `ac_realness` for phase objects and
   `local_focus_map` for multi-plane samples.

## Find the focus (game)

[`docs/game/index.html`](docs/game/index.html) is a small offline game: five
simulated holograms, a depth slider, and a score against the true synthetic
distance and against WaveHolo's estimated one. Open the file in a browser
(`file://` works), or play it on GitHub Pages:
<https://sparkainlp-x.github.io/waveholo/game/>.

![Find the focus, level 4 after locking in](docs/game/preview.png)

- **Banner in the page:** "Synthetic hologram game. Not a hardware or biology
  validation."
- **Where the frames come from:** `scripts/make_game_frames.py` runs
  `record_hologram`, `estimate_carrier_2d`, `reconstruct_field` and
  `find_focus` on seeded objects with an unknown reference piston and a
  fractional carrier. It stores 41 frames per level (0.1 mm apart,
  24 gray levels, amplitude and wrapped phase) plus the true focus and
  WaveHolo's estimated focus. `python3 scripts/make_game_frames.py` rewrites the
  data inside the page; `--check` compares it with a regeneration.
- **Offline and private:** one HTML file, no network requests, no trackers,
  best score kept only in `localStorage`.

Recorded results (this repository, seed-fixed; `error` = |measured − true|):

| Level | Scene | Noise | WaveHolo metric | WaveHolo error | Plain intensity-variance maximum |
|---|---|---|---|---:|---:|
| 1 | ring, cross, dots | Gaussian 0.05× | realness | 17.2 µm | 63.0 µm |
| 2 | random discs | Gaussian 0.2× | realness | 1.3 µm | 113.0 µm |
| 3 | random diamonds | Poisson, 100 photons | realness | 39.5 µm | 2582.0 µm |
| 4 | binary phase π/2 object | Gaussian 0.1× | ac_realness | 4.4 µm | 2680.0 µm |
| 5 | random bars | Gaussian 0.8× | realness | 315.9 µm | 150.0 µm |

On level 4, plain `realness` is 1320 µm off, so that level uses `ac_realness`
(as the "Default pipeline" below recommends for phase objects). On level 5 the
plain intensity-variance maximum is closer than WaveHolo. Neither result is
tuned.

Game limits: five fixed scenes prove nothing about general performance.
The scenes were chosen after a short exploratory pass over about eleven
candidates, so the table is a recorded demonstration, not an unbiased
benchmark. The "plain intensity-variance maximum" baseline has the wrong
polarity for phase objects by design (Dubois et al. 2006), so level 4 is a
weak comparison.
Frames are 0.1 mm apart, so a human pick can be up to 50 µm from the truth,
while WaveHolo's `find_focus` refines continuously (levels 1, 2 and 4 cannot be
beaten by any frame pick). The answers are in the page source. It is not a
hardware, microscope or biology validation.

## Honest limits

- The DHM test covers 6 of 100 stage positions (Z −50…−37 µm) and
  checks the slope, not absolute focus. Metric slopes differ (85 / 96
  / 108) because curvature removal rescales coordinates.
- The residual synthetic bias (−50 to +64 µm for sharp objects, ≤ 1.6 µm
  for band-limited ones) is sideband cross-talk, set by the recording
  geometry.
- No metric wins on every object. Smooth phase objects and two
  materials on separate supports remain unsolved.
- `imag_energy` and its sub-micron "comb" (v0.1.2) are superseded: the
  comb was the λ/2 piston oscillation, aliased by 1 µm sampling. The
  metric is kept only for comparison.
- The MP4 stores rendered views, not the complex field.
- This is OPTICAL holography (wavefronts), not associative holographic
  memory and not the LED fan display.
- What the tests check against a reference, and what has no such check, is
  listed in [docs/validation.md](docs/validation.md).

## Data credits
- pyDHM fly-head hologram `real_data/fly-20x.bmp` (MIT),
  https://github.com/catrujilla/pyDHM.
- Cuenat et al. DHM dataset, Zenodo 6337535 (CC-BY-4.0),
  https://doi.org/10.5281/zenodo.6337535.

Details and licence texts are in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). Project direction is in
[WAVEHOLO_DIRECTION.md](WAVEHOLO_DIRECTION.md).

## Cite

Archived on Zenodo. Concept DOI (always the latest version):
[10.5281/zenodo.23644416](https://doi.org/10.5281/zenodo.23644416). DOI of
v0.1.5: [10.5281/zenodo.23653721](https://doi.org/10.5281/zenodo.23653721);
v0.1.4: [10.5281/zenodo.23644417](https://doi.org/10.5281/zenodo.23644417).
`CITATION.cff` cites v0.1.5 (its date is the Zenodo publication date, UTC).

License: AGPL-3.0-only (see `LICENSE`).

## Publishing to PyPI (maintainers)

The package is **not** on PyPI, and `.github/workflows/publish.yml` will not
publish it unless you opt in. A `guard` job must pass before the
`pypi-publish` job can run:

1. The repository variable `PYPI_PUBLISH` must be exactly `true`
   (Settings, Secrets and variables, Actions, Variables). It is unset by default.
2. The run must be on a tag `vX.Y.Z` whose version equals `version` in
   `pyproject.toml` (checked by `scripts/check_release_tag.py`). A mismatch
   fails the run when the variable is set.

The tests in `test-gate` still run on every version tag.

To enable publishing:

1. On PyPI, add a trusted publisher for project `wave-hologram`: owner
   `sparkainlp-x`, repository `waveholo`, workflow `publish.yml`, environment `pypi`.
2. In the repository settings, add **required reviewers** to the `pypi`
   environment, so a person approves each upload.
3. Set the `PYPI_PUBLISH` variable to `true`, then push a matching tag.

Set the variable back to anything else (or delete it) to turn publishing off again.
