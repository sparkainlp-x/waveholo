# Changelog

## Unreleased

Fixes after an independent review. No change to the focus pipeline's results.

### Fixed
- `pair_spectrum_realness` built one phase kernel for the whole sweep
  (`len(sweep) * ny * nx` complex numbers, about 1.6 GB for 101 frames of
  1024^2). It now works in chunks of at most about 64 MB. Results are
  unchanged (tests compare with propagation to 1e-12).
- The WebAssembly sandbox printed "WaveHolo Measured Focus Plane" for a
  hand-written analytic curve whose minimum is at -5.0 mm by construction. It is
  now labelled as an illustrative placeholder, and the browser test follows.

### Changed (wording and metadata)
- README: "estimate" instead of "measure"; the real-data result is described as
  a slope check on 6 of 100 stage positions, with the class-mean 95% interval;
  "100 photons per pixel" now says peak (about 17 mean); the synthetic bias
  range is -50 to +64 µm; new Related work section; PyPI statement; game scene
  selection and baseline caveats; WaveHolo / `wave-hologram` naming note.
- RESULTS_v0.1.3.md: an addendum (originals unchanged). REVIEW.md and
  WAVEHOLO_DIRECTION.md: historical notes.
- CITATION.cff: title "WaveHolo (wave-hologram): ...", ORCID. The Zenodo record
  title keeps its wording until it is edited or a new version is archived.
- THIRD_PARTY_NOTICES.md: the arXiv title is "Tiny Transformer Networks".
- `docs/index.html` so https://sparkainlp-x.github.io/waveholo/ no longer returns
  404; the game page footer states the scene-selection caveat.
- `tests/test_review_fixes.py`.

## 0.1.4 (2026-10-25)

Packaging, CI, examples and a game. **No change to the focus pipeline, the
metrics or any measured result**; only the version string was bumped.

### Added
- `pyproject.toml` for the flat layout (single module `wave_hologram`,
  numpy and matplotlib; extras `dhm` and `notebook`). `pip install -e .`
  works on Python 3.9 to 3.14. No PyPI release.
- GitHub Actions workflow `.github/workflows/ci.yml`: unittest on
  Python 3.9 to 3.14, plus a check that the game data matches a regeneration.
- `examples/quickstart.py` and an executed `examples/waveholo_demo.ipynb`
  (synthetic hologram, focus curve, reconstructed intensity and phase, one
  phase-object scene).
- "Find the focus" game: `docs/game/index.html` (offline, single file) and
  `scripts/make_game_frames.py`, which builds its frames and recorded
  results with `wave_hologram.py`. Preview in `docs/game/preview.png`.
- `tests/test_game.py`: 13 tests (frame generator determinism and shape, PNG
  round trip, recorded truths and measured foci match a regeneration, page is
  self-contained, quickstart). The suite is now 56 tests.
- README sections: Install, Quickstart, Find the focus (game).

### Unchanged
- `wave_hologram.py` apart from `__version__`, the 43 existing tests, all
  results and limits.

### Notes
- Zenodo archive of v0.1.4: https://doi.org/10.5281/zenodo.23644417 (concept DOI
  10.5281/zenodo.23644416), added to README and CITATION.cff after the release.
- The game is served by GitHub Pages (from `main` `/docs`) at
  https://sparkainlp-x.github.io/waveholo/game/ . It is a synthetic game, not a
  hardware or biology validation.

## 0.1.3 (2026-10-25)

### Fixed
- **Focus metric depended on the reference phase.** `imag_energy`
  follows sin² of the unknown object/reference piston, so its focus
  error reaches 0.79 mm. The new default `realness` = |Σ U²| / Σ |U|²
  is exactly piston-invariant.
- **The "comb" of sub-micron focus teeth (DISCOVERY finding 5) was an
  artefact.** It is the λ/2 = 0.266 µm piston oscillation aliased by
  1 µm sampling into 4.16 µm teeth. That finding is retracted.
- **Integer-bin carrier.** It is replaced by a two-pass local-DTFT
  sub-bin estimate (`estimate_carrier_2d`, `estimate_carrier`), with
  1-D and 2-D versions and an optional hint, so lattice peaks in real
  data can be avoided.
- Three comb-dependent tests were replaced. The suite is now 43 tests,
  all passing.
- ffmpeg is called with `-nostdin`, because it stopped when run from a
  background shell.

### Added
- `ac_realness` (for phase and mixed objects), `amplitude_variance`,
  `tamura`, `phase_variance` and `best_piston`.
- `find_focus` (coarse sweep plus golden-section refinement, off-grid)
  and `golden_min`.
- `local_focus_map`: per-tile depth maps for multi-plane samples.
- `remove_quadratic_phase`: removes microscope-objective curvature
  without unwrapping (needed for realness on DHM data).
- `reconstruct_field(..., crop=True)` with `cropped_pitch`: a lossless,
  faster sideband crop.
- `record_hologram(..., ref_phase=, carrier_y=)` with a fractional
  carrier.
- `studies.py`: ten reproducible studies (S1–S10) with figures in
  `results/`.
- `real_data/dhm_v013.py`: multi-class validation on the Cuenat et al.
  DHM dataset (530 holograms) in `real_data/dhm_v013/`.
- `RESULTS_v0.1.3.md`: the write-up.

### Changed
- Default metric: `imag_energy` is replaced by `realness`. The demo
  estimates the carrier by default (`--assume-carrier` restores the old
  behaviour).
