# FINDINGS.md — Cuenat et al. DHM holograms through the pipeline (2026-10-25)

## What was obtained

- **Source:** Cuenat et al., "Fast Autofocusing using Tiny Networks
  for Digital Holographic Microscopy" (arXiv:2203.07772), dataset
  Zenodo 6337535 (CC-BY-4.0).
- **What it is:** 40,400 EXPERIMENTAL holograms (1024×1024, 8-bit)
  of a pseudo-periodic PHASE pattern (9 µm period, 12-bit encoding,
  100 nm Al layer), recorded on a Lyncee Tec DHM with a Leica 10×
  objective (NA 0.32), Basler camera (5.86 µm pitch), λ = 674.99 nm
  (SLD + 5 nm filter, coherence ≈ 100 µm). 100 classes × 400
  holograms; the mapping file gives the hexapod Z per class,
  −50.000 … +42.015 µm in ~1 µm steps.
- **What was actually downloaded here:** the mapping file
  (`mapping_z.txt`, complete) plus the FIRST 2 GB of the 48.6 GB
  `dataset.tar.gz` (HTTP range request), stream-extracted to
  `extracted/`: 1,983 complete holograms from directories
  class_0, class_1, class_10, class_11, class_12 (tar order is
  alphabetical, hence the jump). The full 48.6 GB was NOT
  downloaded — unnecessary for this test, and the chunk already
  spans 5 classes.
- **Scripts:** `run_dhm.py` (single hologram: carrier, coarse
  sweep, sub-micron fine scan), `multi_dhm.py` (multi-class
  basin test). The core package is untouched.

## The headline: a real measured-vs-logged slope

One hologram per available class, carrier (281, −268) bins,
window 90, sweep −10…+10 mm in 0.05 mm steps, basin centre from a
parabolic fit on the smoothed amplitude-variance curve:

| dir class | logged Z (µm) | file used | basin centre (mm) |
|---:|---:|---|---:|
| 0 | −50.000 | class_1_holo_0.tif | −8.113 |
| 1 | −49.000 | class_2_holo_399.tif | −8.189 |
| 10 | −40.065 | class_11_holo_3999.tif | −9.105 |
| 11 | −39.079 | class_12_holo_4399.tif | −9.163 |
| 12 | −38.095 | class_13_holo_4799.tif | −9.148 |

Linear fit of basin centre vs logged Z:

- **Slope: −93.7 µm of basin per µm of logged Z**
  (−0.0937 mm/µm), residual RMS 0.047 mm.
- Reading: the M² hypothesis predicts +100 µm/µm in magnitude
  (10× objective → longitudinal magnification 100); the direct
  hypothesis predicts +1. The measured |slope| = 93.7 is the M²
  scaling, 6% low — consistent with a true magnification of
  √(93.7) ≈ 9.68 rather than the nominal 10× (objectives' nominal
  magnification is approximate; the tube-lens factor also enters).
  The SIGN is negative because the pipeline reports
  back-propagation distances (focus of an object at +z appears at
  −z) and the stage direction runs opposite to it; the sign
  convention is a bookkeeping fact, not a failure.
- The intercept (−12.794 mm) is the DHM image-plane calibration
  constant: where the microscope's image plane sits in
  hologram-plane coordinates. It is not an error term.
- **This is the first measured-vs-truth score this package has
  ever had on real data** — and it is a SLOPE (differential),
  which survives the two absolute ambiguities below. A single
  hologram could never have shown this; the series is what makes
  it a test.

## Ambiguities, stated not hidden

1. **Class index ±1.** Tar directories are class_0… but files
   inside are named class_1… (1-indexed originals repacked
   0-indexed, or vice versa). This report uses the DIRECTORY
   index as the mapping class. If the filename index is the
   truth, add exactly +1 µm to every logged Z above. The slope
   is identical either way; only the intercept shifts by 0.094 mm.
2. **Absolute offset.** The intercept absorbs where the DHM
   image plane sits; without an independent measurement of that
   plane, no absolute per-hologram error (measured − logged)
   can be quoted. The slope + residual RMS is the honest score.
3. **Basin, not tooth.** The amplitude-variance dips are 9.0–10.4%
   (shallow) and ~1.5 mm wide in hologram-plane coordinates
   (15 µm sample DoF × M²). There is NO sub-micron tooth here:
   fine scans (±20 µm around the basin) are flat (0.0% local
   dip). That is expected, not a failure — the comb of
   DISCOVERY.md finding 5 belongs to imag-energy on amplitude
   objects, not to amplitude variance on phase objects. Do NOT
   report a "focus error" finer than the basin width for this
   object type; the honest resolution here is the residual RMS
   (47 µm hologram-plane ≈ 0.5 µm sample-space at M² = 100).
4. **Phase object.** Imag-energy is flat on every run here
   (1.02–1.48× contrast), exactly as predicted for a phase
   object (REAL_DATA.md finding B). Nothing about this dataset
   tests imag-energy; it tests amplitude variance + the M²
   scaling, and both behave as predicted.

## Wrong turns on the way (documented, not hidden)

1. **Auto carrier (70, −70), window 200:** the strongest
   spectral peak is a PATTERN LATTICE point, not the carrier
   (the 9 µm periodic pattern puts diffraction orders on a
   70-bin lattice centred on DC). Window 200 overlapped DC and
   the twin — contaminated reconstruction, meaningless metrics
   (amp-variance "min" +8.25 mm, Laplacian max −0.25 mm).
2. **Same peak, window 45:** clean of DC but cuts the pattern
   fundamental (66.7 bins from sideband centre) out — mush
   reconstruction, amp-variance min at the sweep edge.
3. **Off-lattice peak search** (masking the 70-bin lattice)
   found the true carrier candidate: (−281, +268) /
   (+281, −268), 388 bins, magnitude 1.43e6 vs lattice 3.24e6.
   With window 90 the reconstruction SHOWS the pattern
   (see dhm_c0_car281_best.png) — qualitative validation of
   the carrier. Lesson for future data: on a periodic object,
   never trust the strongest spectral peak as the carrier;
   mask the object's own lattice first, and always check the
   reconstruction by eye.
4. **Smoothed-dip print bug (multi_dhm.py, first version):**
   mode='same' convolution edges fake a ~49% dip. Corrected
   to interior-only: true dips ~9%. Basin centres were
   unaffected (edges were already masked in the fit) and were
   independently re-verified for class 0 (raw min −8.10 mm,
   raw dip 8.9%).

## Why this dataset is NOT the ideal calibration target

It works as a differential test, but it is a phase object
through a microscope (M² ambiguity, 15 µm DoF, lattice/carrier
collision, 48.6 GB single tarball). The ideal remains
RECORDING_PROTOCOL.md §7: a USAF AMPLITUDE target, lensless,
a logged stage series, lossless files, one JSON sidecar per
hologram. That session would give an absolute per-hologram
error, which this dataset structurally cannot.

## Other routes checked this round (all blocked, briefly)

- **arXiv 2308.13115 (USAF, 532 nm, 3.45 µm, 86 mm):** parameters
  fully documented in text, but NO data-availability statement
  and no downloadable hologram. Blocked at "no file".
- **P-TDHM (Jin et al., PMC10847153; Zenodo 8206220):** the
  Zenodo record holds only the software/STL zip, no holograms;
  the paper's 150 µm figure is an IN-LINE particle-tracking
  result. Blocked at "no hologram files".
- **Zenodo generic "holography" search:** dominated by analogue
  holography, electron holography, and software; no off-axis
  optical hologram dataset with distances surfaced. The
  "off-axis digital holography" phrase query returns a single
  unrelated theory paper.
- **pyDHM samples:** re-checked; star-target and USAF
  phase-shifting samples carry no documented distance (the
  "+3cm.tif" numerical-propagation sample is an object image,
  not a hologram). Blocked at "no distance".
- **JoVE DLHM:** already closed in the previous round (OSF
  hosts only the index CSV; in-line anyway).

## Files in dhm_tiny/

- `mapping_z.txt` — class→Z mapping (complete, 100 classes)
- `chunk.bin` — first 2 GB of dataset.tar.gz (range download)
- `extracted/` — 1,983 holograms, classes 0, 1, 10, 11, 12
- `run_dhm.py`, `multi_dhm.py` — run scripts (package untouched)
- `dhm_c0_*`, `multi_basin.png` — outputs for class 0 + slope plot
- `holo_c0.png`, `holo_c0_crop.png` — raw hologram views
- `FINDINGS.md` (this file)
