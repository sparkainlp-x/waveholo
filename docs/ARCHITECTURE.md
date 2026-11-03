# Architecture

The maintained numerical implementation is currently the flat module
`wave_hologram.py`; compatibility and scientific regression coverage are retained
there rather than forcing an unreviewed package migration. It contains simulation,
angular-spectrum propagation, hologram formation, sideband reconstruction, focus
metrics/search, and plotting. Historical/real-data scripts live under `real_data/`.

The benchmark layer belongs under `benchmarks/`: versioned YAML configuration,
deterministic scene generation, isolated metric evaluation, result serialization,
manifest/checksum creation, and figure generation. It must not alter historical
reports or raw data. A future `src/wave_hologram/` migration is a separately
reviewable API-compatibility change.
