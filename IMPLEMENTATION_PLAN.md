# Implementation plan

The changes are intended as small, independently reviewable pull requests; each
pull request must retain the evidence labels and claim policy.

1. `docs: add scientific protocol and claims policy` — protocol, metrics, real-data
   boundary, validation matrix, limitations, and claims.
2. `feat(metrics): add reconstruction and focus metrics` — separately testable
   focus criteria and image metric formulae.
3. `feat(simulation): add deterministic synthetic benchmark fixtures` — seeded
   point, extended, and textured objects and small smoke config.
4. `feat(experiments): add factorial benchmark runner` — versioned configs,
   reconstruction-distance sweeps, seeds, and condition outputs.
5. `test: add unit, property, regression, and integration tests` — preserve legacy
   physics tests and extend required layers.
6. `feat(reporting): add manifests, checksums, figures, and claim gates` — complete
   metadata, CSV summaries, generated figures, and fail-closed eligibility.
7. `docs: add limitations, validation matrix, and reproducibility guide` — finalize
   operational documentation after checking actual runner behavior.
8. `benchmark: run smoke benchmark and archive results` — run only after commit
   identity is known; mark failed/partial conditions explicitly.
9. `benchmark: run full protocol benchmark and publish summary` — manual/scheduled
   execution; never fabricate or replace historical outputs.
10. `docs: update README with measured results and failure envelope` — only after
    the full run and independent review.
