# Data and external inputs

Do not commit large raw datasets or generated experiment outputs. Synthetic fixtures
used by tests are small and deterministic.

The historical microscope result used a subset of the Cuenat et al. DHM dataset
(Zenodo 6337535, CC-BY-4.0); the full archive is external and its per-frame checksum
manifest is not part of this repository. Independent evaluation must acquire the
licensed source directly, retain its original metadata, and publish a file manifest
and SHA-256 checksums outside this Git repository. See
[`REAL_DATA_PROTOCOL.md`](../docs/REAL_DATA_PROTOCOL.md).

`fly-20x.bmp` is a small credited external example whose source and licence are
documented in `THIRD_PARTY_NOTICES.md`. It has no ground-truth focus distance.
