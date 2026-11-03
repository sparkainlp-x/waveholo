# Reproducibility

Each run is identified by its protocol version, exact config, commit SHA, UTC
timestamp, Python/dependency versions, machine details, random seeds, input/output
checksums, and runtime. Outputs use stable JSON key ordering and UTF-8.

Run the smoke case with:

```bash
make benchmark CONFIG=benchmarks/configs/waveholo_smoke.yaml
```

Run the nominal case with:

```bash
make benchmark CONFIG=benchmarks/configs/waveholo_nominal.yaml
```

The full benchmark is manual/scheduled, not per-commit CI. Every output directory
must contain `config.yaml`, `manifest.json`, `environment.json`, `code_commit.txt`,
`metrics.csv`, `per_seed_results.csv`, `summary.json`, `figures/`, and
`checksums.sha256`. A run is invalid when a required seed/metric is absent, any
checksum fails, config differs from the locked protocol, a claim gate fails, or the
code commit is unknown. Historical real-data results are not silently recreated
from current code.

Large external data and raw acquisitions stay out of Git. Consult
[`data/README.md`](../data/README.md) for source, licence, access, and checksum
instructions.
