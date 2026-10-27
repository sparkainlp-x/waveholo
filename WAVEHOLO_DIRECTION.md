# WAVEHOLO — DIRECTION

## WHAT IT IS

Off-axis digital holography with focus estimation. Angular-spectrum
propagation, sideband reconstruction, and focus metrics — every claim
measured, every failure documented. Classical computation only. No
hardware, no speedup, no product claim.

## WHAT IT IS NOT

Optical holography (wavefronts) is not membrane holographic memory
(W = sum |p><p|) and is not an LED fan display. Same word, three
things. Never merge the claims.

## WHERE WE ARE — v0.1.3, 43/43 tests

- **Metric fixed.** `realness` = |Σ U²| / Σ |U|² replaces `imag_energy`.
  It is exactly invariant to the unknown object/reference piston;
  `imag_energy` followed sin² of it and erred up to 0.79 mm. Demo with
  a 1.2 rad piston and a fractional carrier (80.37), both unknown to
  the pipeline: truth −5.000 mm, measured −4.998 mm, recovery corr
  0.956.
- **Carrier.** Sub-bin estimate (two-pass local DTFT), error within
  0.03 bins. It tracks the ±0.1–0.2 bin shot-to-shot carrier jitter in
  the real DHM series.
- **Noise.** With a random piston in every trial, realness lands within
  0.1 mm in 12/12 trials up to 0.35× Gaussian noise and down to 100
  photons per pixel; `imag_energy` manages 0–2/12 at every level.
  v0.1.2's "holds to 2x Gaussian" was true only at piston = 0 —
  withdrawn.
- **Off-grid focus.** `find_focus` (coarse sweep + golden section):
  2–49 µm error for off-grid truths. The remaining bias is sideband
  cross-talk set by the recording geometry, not the metric.
- **Comb retracted.** The v0.1.2 sub-micron "focus comb" was the
  λ/2 = 0.266 µm piston oscillation aliased by 1 µm sampling into
  4.16 µm teeth. Not optics.
- **3D sectioning.** `local_focus_map` (per-tile focus) recovers both
  planes of the two-plane sample (5.0–5.1 mm and 7.5–7.7 mm), where a
  single global focus reports one number.
- **Phase objects.** `ac_realness` finds 9/12 synthetic phase and mixed
  objects (binary phase, single-material amplitude+phase). Smooth phase
  objects and two materials on separate supports are still open.
- **Real data 1:** pyDHM fly-head hologram. Three metrics
  (ac_realness, amplitude_variance, tamura) agree within 0.3 mm near
  +1.1 mm. Agreement between metrics only — no ground-truth distance,
  no error claimed.
- **Real data 2:** Cuenat et al. DHM series (Zenodo 6337535), 530
  holograms at 6 logged stage positions (Z −50…−37 µm). Realness
  (curvature removed) gives slope −96.4 ± 1.3 µm focus per µm stage,
  consistent with M² ≈ 100 for a 10x objective; per-position means sit
  0.19 µm from the fit. Three metrics give slopes 85 / 96 / 108, so the
  true magnification is unresolved. Slope only — no absolute
  per-hologram error quotable from it.
- **Open problems, stated not hidden:** phase-object autofocus (smooth
  phase), absolute focus never tested against independent truth, only
  6 of 100 DHM stage positions used, Talbot replica planes on the
  periodic target not investigated, sideband cross-talk bias up to
  ±50 µm on sharp objects.

## RULES

- License AGPL-3.0-only. Sole author: Jean-François Brisson,
  Spark AI NLP, Fredericton, ORCID 0009-0000-9778-5374.
- Third-party data credited: pyDHM (MIT), Cuenat et al. Zenodo
  6337535 (CC-BY-4.0).
- Never report a focus plane without its metric and contrast.
- Never claim an absolute error without independent ground truth.
- A single hologram never separates offset from error — only a
  logged series makes a test.

## NEXT, IN ORDER

Done: v0.1.3 pushed to github.com/sparkainlp-x/waveholo (branch main).
(Update: v0.1.3 was never tagged. v0.1.4 was tagged and released, and its Zenodo
record is 10.5281/zenodo.23644417. The three steps below are historical.)

1. Switch waveholo ON in Zenodo GitHub settings BEFORE tagging —
   otherwise no DOI mints.
2. Tag v0.1.3.
3. Absolute test: a new lensless USAF amplitude recording, distance
   logged, per real_data/RECORDING_PROTOCOL.md.
4. Phase-object autofocus: still the open research problem.
   `ac_realness` solves binary phase objects; smooth phase objects
   defeat every metric, and an amplitude-variance minimum is a hint,
   not a solution.
