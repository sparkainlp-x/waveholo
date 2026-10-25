# RECORDING_PROTOCOL.md — a calibration-target hologram that closes the gap for good

Goal: one (or a small series of) real off-axis hologram(s) for which the
focus distance is known INDEPENDENTLY of any focusing algorithm, so the
wave-hologram pipeline can finally be scored as measured-vs-truth
instead of measured-vs-nothing. This protocol is written from the
failures of the two previous real-data attempts:

- The pyDHM fly-head hologram reconstructed crisply but had NO
  documented distance, so no error could be claimed (see REAL_DATA.md).
- The JoVE DLHM "database" turned out to host only an index CSV, no
  images, and is in-line anyway.
- The Cuenat et al. DHM dataset (Zenodo 6337535, 40,400 holograms)
  DOES carry a per-class ground truth (hexapod Z, −50…+42 µm), but
  only as a 48.6 GB tarball, of a PHASE object, through a 10×
  microscope objective — so its truth lives in sample space and the
  hologram-plane scaling (×1 vs ×M²) has to be established from the
  data, not assumed. See `dhm_tiny/` for that run.

## 1. Target choice

- **First choice: a chrome-on-glass USAF 1951 resolution target, used
  as an AMPLITUDE object** (opaque bars on clear glass, transmission).
  Amplitude matters: the package's imag-energy metric assumes a real
  amplitude object and bottoms out only for one. A phase object (like
  the Cuenat pattern or a fly) inverts the metrics and leaves only a
  shallow amplitude-variance dip — usable as a hint, not a test.
- **Second choice (same session, cheap): a pinhole (1–5 µm) and a
  sparse microsphere slide.** A pinhole is the closest thing to a
  point object: its reconstruction is a direct check of the
  propagator, and its focus curve is the cleanest possible.
- **Do NOT use:** a biological sample, a thick sample, or a phase-only
  object as the primary target. Those are interesting later targets,
  not calibration targets. If you record one anyway, label it
  phase/thick and do not score imag-energy on it.
- Record the target's identity (manufacturer, part number, group /
  element used) so the bar widths are known: they give an independent
  lateral-scale check (reconstructed bar width in pixels × pitch
  must equal the known width, after magnification if any).

## 2. Geometry: the distance IS the experiment

The single most important number is the object–sensor distance, and it
must be measured by a ruler, not by focusing:

- **Preferred geometry: lensless, plane-wave illumination, no
  microscope objective between object and sensor.** Then the
  hologram-plane pitch IS the sensor pitch, the propagation distance
  IS the object–sensor distance, and there is no magnification
  ambiguity at all. Put the target on a translation stage, measure
  the target-to-sensor distance with calipers / a micrometer stage
  readout, and log it to ±0.1 mm or better.
- **If a microscope objective is unavoidable:** log the objective
  magnification M (measured, not nominal — see §4), and record the
  sample displacement from a known reference position with the
  stage encoder. Remember the scaling law: a sample displacement δ
  corresponds to M²·δ of propagation in hologram-plane coordinates
  (longitudinal magnification), and a hologram-plane pitch dx
  corresponds to dx/M in sample space. The pipeline's measured
  distance must be converted with the SAME convention as the logged
  truth — state the convention in the log file, and test both
  scalings if unsure (the Cuenat run does exactly this: class steps
  of ~1 µm in sample space appear as ~1 µm or ~100 µm steps in
  hologram-plane coordinates depending on the convention, and the
  data decides).
- **Record a SERIES, not a single hologram:** the same target at
  5–10 stage positions spanning at least a few mm (lensless) or
  ~100 µm (microscope), with the stage position logged per file.
  A series tests the slope (measured vs logged step size), which
  survives an unknown absolute offset; a single hologram tests
  only one number and cannot separate offset from error.
- **Log the sign convention:** which stage direction moves the
  target AWAY from the sensor, and whether the pipeline's
  back-propagation distance is negative (as in this package, where
  focus of an object at +z appears at −z). A sign flip is the
  cheapest possible false "failure".

## 3. What to log per hologram (a checklist, in the file header)

Every hologram file ships with a sidecar (same basename, `.json` or
`.txt`) containing ALL of these — a hologram without its sidecar is
not a test, it is a picture:

- [ ] wavelength λ (nm) + source type + how verified (datasheet /
      wavemeter / filter spec). If a filtered SLD/LED, log the
      filter bandwidth: it sets the coherence length, which caps
      the usable distance range.
- [ ] sensor pixel pitch dx, dy (µm) + camera model + how verified
      (datasheet AND a pitch check: image a known grating and
      confirm the fringe period in pixels matches λ/angle
      predictions, or image the USAF target bars).
- [ ] object–sensor distance z (mm) + how measured (calipers /
      stage encoder / micrometer) + uncertainty.
- [ ] stage position per file (if a series) + step size + direction
      convention (see §2).
- [ ] objective magnification M, if any, + how measured (see §4).
      If none, write "lensless, M = 1".
- [ ] reference-beam angle, or equivalently the carrier: log the
      tilt if known, but ALWAYS re-derive it from the spectrum
      (the package's `estimate_carrier` approach) — the pipeline
      must never trust an assumed carrier (a 2-bin error collapsed
      focus contrast to 1.2× in the sensor test; see REVIEW.md).
- [ ] exposure time, gain, bit depth. Saturation clips fringes and
      halves metric contrast (sensor test: 27.8× vs 50.8× clean).
- [ ] dark frame and flat field, if taken (subtract/divide before
      reconstruction; vignetting and dust cost contrast, not focus,
      but a clean hologram is a cleaner test).
- [ ] file format: lossless only (TIFF/PNG, 8- or preferably
      16-bit). NEVER JPEG for a hologram: lossy compression eats
      exactly the high-frequency fringe information the
      reconstruction lives on.
- [ ] target identity + whether amplitude or phase (see §1).
- [ ] ambient notes: vibration isolation (tilt jitter is the one
      vibration that kills focus — see DISCOVERY.md finding 3),
      temperature drift if the session is long.

## 4. Verifying pixel pitch and magnification (do not trust labels)

- **Pitch:** the fly hologram's repo quoted 2.4 µm but its own
  propagation example used 5.4 µm unexplained — and every z in the
  analysis scales as dx², so a pitch error of ×2.25 is a distance
  error of ×5.06. Verify pitch by imaging a known-period grating
  or the USAF bars and comparing measured vs nominal feature size
  in pixels.
- **Magnification:** image the USAF target, measure a known bar
  group in sensor pixels, and compute M = (pixels × dx) /
  (known width). Log the measured M, not the engraved M.

## 5. How to run it through this package (the check procedure)

From `~/workspace/wave-hologram/real_data/`, adapting
`real_fly_focus.py` (2-D sideband generalisation) for the new file:

1. Load the hologram (lossless, float). Subtract the mean.
2. **Estimate the carrier from the spectrum** (strongest off-DC
   peak, DC blanked, centroid-refined — the `estimate_carrier`
   idea, generalised to 2-D). Record the estimated carrier; if it
   disagrees with the logged tilt by more than ~1 bin, say so —
   that disagreement is itself a finding about the setup.
3. Crop the sideband around the estimated carrier, shift to
   baseband. Window size: as large as possible without touching
   DC or the twin sideband (note the truncation cap: a window
   narrower than the object spectrum caps recovery correlation,
   exactly as the synthetic 0.951 ceiling).
4. **Coarse sweep first** (e.g. 0.1 mm steps, or 1 µm steps in
   sample-space units) across a range bracketing the logged
   distance, tracking imag-energy (amplitude target) AND
   amplitude variance (phase-target diagnostic). Find the basin.
5. **Sub-micron scan around each candidate dip.** The imag-energy
   curve is a COMB (teeth every ~0.004 mm at 532 nm / 6.5 µm;
   tooth spacing scales with λ and pitch — recompute for your
   parameters, do not copy 0.004 mm). A coarse grid that misses
   the central tooth lands on the broad dip ~0.3 mm off (see
   DISCOVERY.md finding 5). Procedure: coarse sweep → find basin
   → fine scan with steps ≪ tooth spacing around the dip →
   report the tooth. Never report the coarse-grid argmin as the
   measurement, and never use parabolic interpolation on the
   comb (it fits the broad shoulder, not the tooth).
6. Report, per hologram: logged z, measured z (fine-scan tooth),
   error, metric contrast (median/min), and — for the USAF
   target — the reconstructed bar widths vs known widths
   (lateral check). For a series, ALSO report the fitted slope
   of measured-vs-logged (slope ≈ 1 is the strong claim; a
   constant offset is a calibration constant, a wrong slope is
   a physics error) and the residual scatter around the fit.
7. State the scaling convention used (lensless: direct;
   microscope: which factor of M was applied) and, if the
   convention was uncertain, show both scalings and which one
   the data supports.

## 6. What "passing" looks like

- Amplitude target, lensless: measured-vs-logged slope 1.00 ± 0.05
  across the series, per-point error within the fine-scan tooth
  width, imag-energy contrast ≥ 10×.
- Microscope geometry: same, AFTER the stated M²/pitch
  conversion — and the conversion itself documented as part of
  the result, not buried.
- Phase target (secondary): amplitude-variance minimum within
  the series' step size of the logged plane, with the shallow-dip
  warning attached. A phase target that "passes" imag-energy
  should be treated as suspicious, not celebrated — the metric
  is not supposed to work there.

## 7. Minimal session (if time/equipment is short)

One afternoon, one target, one series: USAF amplitude target,
lensless, 532 nm (or any logged λ), sensor pitch logged and
verified, 6 stage positions 5 mm apart starting ~20 mm from the
sensor, lossless 16-bit TIFFs, one JSON sidecar per file with
§3's checklist. That session alone would give this package its
first measured-vs-truth score — the thing two archive hunts
have not yet produced.
