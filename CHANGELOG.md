# Changelog

## 0.1.5 (2026-10-26)

Release preparation; no GitHub release or new Zenodo archive is published
by these repository changes.

### Added
- `pair_spectrum_focus_curve`: exact realness and AC-realness curves from
  the S9 pair-spectrum identity, evaluated in bounded distance blocks.
- Optional CuPy backend and `find_focus_batch`, processing one field at a
  time; NumPy remains the default and GPU dependencies are not required.
- Explicit `find_focus(..., dhm_mode=True)` quadratic-phase correction,
  disabled by default. Corrected coordinates are shifted/scaled, not an
  absolute focus calibration.
- Acceleration regression tests for direct-propagation equivalence,
  bounded blocks, custom metrics, DHM, batches, invalid input, and extreme
  finite amplitudes. CUDA tests skip without CuPy or a usable device.
- `.zenodo.json` with creator ORCID, affiliation, and AGPL-3.0-only license;
  matching creator details in `CITATION.cff`.

### Changed
- `find_focus` uses the pair spectrum for realness/AC-realness, including
  refinement, preserving existing positional arguments and result keys.
  Other metrics and custom callables retain direct propagation.
- `studies.py` delegates its S9 curve to the library implementation.
- Local coarse-curve CPU speedups measured 1.33× for realness and 1.49× for
  AC-realness, with maximum differences of 1.7 × 10⁻¹⁶; these are not
  universal speed guarantees. Benchmark details are in the README.

### Unchanged
- Geometry-cache capacity, existing game data, and historical DHM evidence.
- The v0.1.4 Zenodo archive and its DOI. Correcting that record's metadata
  requires the owner's Zenodo login; no v0.1.5 DOI is claimed here.

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
