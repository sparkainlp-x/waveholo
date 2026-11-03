# Validation matrix

| Capability | Evidence type | Status | Metric | Protocol | Limitation |
|---|---|---|---|---|---|
| Synthetic hologram generation | Unit/integration test | Demonstrated | Shape, dtype, deterministic fixture | Existing tests; PROTOCOL-v1.0.0 extension UNRUN | Synthetic only |
| Reconstruction | Benchmark | Demonstrated for existing fixtures; factorial target UNRUN | RMSE, PSNR, SSIM | PROTOCOL-v1.0.0 | Algorithm-dependent |
| Relative focus estimation | Benchmark | Demonstrated for selected prior synthetic scenes; factorial target UNRUN | Focus error, trend agreement | PROTOCOL-v1.0.0 | Not absolute focus |
| Real-data trend check | REPORTED | Existing; PARTIALLY REPRODUCIBLE | Slope sign and uncertainty | REAL_DATA_PROTOCOL.md | Six stage positions; not absolute focus |
| Zernike recovery | Benchmark | Not in scope | Not applicable | Not applicable | See companion telescope repository |
| Telescope focus recovery | Benchmark | Not in scope | Not applicable | Not applicable | See companion telescope repository |
| Operational telescope validation | External measurement | UNRUN | Not defined | Not defined | No real data |
