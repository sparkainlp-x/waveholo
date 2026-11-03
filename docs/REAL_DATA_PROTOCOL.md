# Real-data protocol: limited microscope slope check

**Evidence label:** REPORTED. **Reproducibility:** PARTIALLY REPRODUCIBLE.

The historical result is retained; full independent rerun requires additional
protocol metadata and access to the external dataset. The project analysis reports
530 holograms drawn from six of 100 stage positions in the Cuenat et al. DHM dataset
(Zenodo 6337535, CC-BY-4.0). The complete archive is not stored in Git. The
historical analysis identifies wavelength 674.99 nm, 5.86 µm camera pixels, a
10× objective, and phase-target imaging, but does not establish every effective
plane pitch and magnification convention needed to independently reconstruct the
focus axis.

**Data source/access:** Zenodo 6337535; download/use under CC-BY-4.0; analysis entry
point `real_data/dhm_v013.py`. **Count:** 530 analyzed frames. **Stage positions:**
six reported positions from the first partial archive subset, not the full 100.
The historical data acquisition downloaded the first approximately 2.2 GB of the
48.6 GB archive, covering six available positions. This is an availability-limited
subset, not a randomized sample of all positions; the exact frame-selection log is
insufficiently frozen. Wavelength is reported as 674.99 nm; camera pixel pitch as
5.86 µm; objective as 10×. Effective hologram-plane pitch and full geometry
calibration remain insufficiently documented.
Do not treat 530 frames as 530 independent stage positions.

**Preprocessing/reconstruction:** the historical script preprocesses the
off-axis holograms, estimates the carrier, removes an objective-curvature term,
then reconstructs using angular-spectrum propagation. Exact versioned dependency
locks, the full frame manifest/checksums, all exclusions, and some calibration
metadata are absent from this historical record. The reported focus measure is
`realness`; slope is fit between metric-derived focus distance and logged stage
position. Historical reported slopes vary by metric and coordinate convention;
the class-mean `realness` fit is about −96.3 µm/µm (reported per-hologram fit
−96.4 ± 1.3, where uncertainty treated frames as independent).

Uncertainty must be recomputed by resampling stage positions or a hierarchical
model, not treating repeated frames at one position as independent positions.
Exclusions and preprocessing deviations must be recorded before examining slopes.
Expected interpretation is consistency of a relative trend under the locked
coordinate mapping; no absolute intercept or clinical/biological implication is
claimed. The historical result is retained and cannot be overwritten by a rerun.

> The real-data check tests consistency of a focus trend under a locked protocol.
> It does not establish absolute focus accuracy, generalization, or production readiness.

Until the missing source-frame identifiers, processing environment, and coordinate
mapping are frozen and rerun, label this check:

> PARTIALLY REPRODUCIBLE — historical result retained; full independent rerun requires additional protocol metadata.
