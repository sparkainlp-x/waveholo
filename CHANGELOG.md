# Changelog

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
