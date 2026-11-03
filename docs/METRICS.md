# Metrics

Distances are metres internally; reports may additionally show micrometres.
Metrics describe numerical similarity, not physical correctness.

| Metric | Definition | Units / interpretation | Implementation |
|---|---|---|---|
| Reconstruction RMSE | `sqrt(mean((A-B)^2))` | Input amplitude units; lower is closer to reference | `benchmarks/run_benchmark.py` |
| PSNR | `20 log10(MAX/RMSE)` | dB; `MAX` is the declared amplitude range | `benchmarks/run_benchmark.py` |
| SSIM | Mean local luminance/contrast/structure score with declared window and constants | Dimensionless; 1 is identical | `benchmarks/run_benchmark.py` |
| Normalized cross-correlation | `sum((A-mean A)(B-mean B))/sqrt(sum((A-mean A)^2)sum((B-mean B)^2))` | Dimensionless; 1 is linear agreement | `benchmarks/run_benchmark.py` |
| Focus-curve contrast | For a minimized positive criterion, `median(m(z))/min(m(z))`; for a maximized non-negative criterion, `max(m(z))/median(m(z))` when denominator is nonzero | Dimensionless; 1 is flat, larger is more contrast; report the exact polarity | `wave_hologram.py:find_focus` (`contrast` for realness) |
| Focus-metric peak location | `z[argmin(m(z))]` for minimized metrics, `z[argmax(m(z))]` for maximized metrics | m; selected plane is the reconstruction distance where the chosen criterion reaches its declared min/max | `wave_hologram.py:find_focus` |
| Absolute focus error | `abs(z_selected-z_true)` | m | Benchmark runner |
| Relative focus error | `abs(z_selected-z_true)/abs(z_true)` | Dimensionless; undefined at zero truth | Benchmark runner |
| Real-data trend agreement | Sign of fitted focus-versus-stage slope, reported with slope and interval | Sign and µm/µm; not absolute focus accuracy | `real_data/dhm_v013.py` |
| Runtime / memory | Wall time and process peak RSS per experiment | seconds / bytes; machine dependent | Benchmark runner |

Focus metric polarity is not interchangeable: `realness` and `ac_realness` are
minimized, while some listed criteria are maximized. PSNR and SSIM are widely used
reconstruction-quality metrics; report them alongside RMSE, without treating them
as proof of physical correctness. [16][22][27] Focus selection is the
reconstruction distance at which the selected criterion reaches its maximum or
minimum, depending on its definition. [20]

No single scalar is reported without its condition, N, and uncertainty. Existing
metric implementations are in `wave_hologram.py`; the complete cross-condition
metric export is **UNRUN** pending the benchmark runner.

References numbered [16], [20], [22], [27] must be mapped to full bibliographic
entries before publication. This repository does not currently identify those
numbered citations; they are source placeholders, not verified references.
