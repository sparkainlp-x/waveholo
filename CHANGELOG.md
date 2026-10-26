# Changelog

## Unreleased

### Changed
- `CITATION.cff` now cites v0.1.5 with its Zenodo version DOI
  10.5281/zenodo.23653721 (concept DOI 10.5281/zenodo.23644416); it still said 0.1.4.
- `publish.yml` has a `guard` job: PyPI publishing needs the repository variable
  `PYPI_PUBLISH` set to `true` and a `vX.Y.Z` tag equal to the `pyproject.toml`
  version. Off by default. See the README.
- `scripts/check_release_tag.py` and `tests/test_release_metadata.py`: the guard's
  check, and a test that `pyproject.toml`, `__version__` and `CITATION.cff` agree.
- `CITATION.cff` author entry gains the ORCID (0009-0000-9778-5374), as in membrane-arena-3x3.

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
