"""Image reconstruction metrics used by the benchmark result schema."""
from __future__ import annotations

import numpy as np


def _pair(reference, estimate):
    a = np.asarray(reference, dtype=float)
    b = np.asarray(estimate, dtype=float)
    if a.shape != b.shape or a.ndim != 2 or a.size == 0:
        raise ValueError("images must be matching non-empty 2D arrays")
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("images must be finite")
    return a, b


def rmse(reference, estimate):
    a, b = _pair(reference, estimate)
    return float(np.sqrt(np.mean((a - b) ** 2)))


def psnr(reference, estimate, data_range=1.0):
    error = rmse(reference, estimate)
    if not np.isfinite(data_range) or data_range <= 0:
        raise ValueError("data_range must be finite and positive")
    if error == 0:
        return float("inf")
    return float(20.0 * np.log10(data_range / error))


def global_ssim(reference, estimate, data_range=1.0):
    """SSIM equation applied to global image moments (not a local-window SSIM)."""
    a, b = _pair(reference, estimate)
    if not np.isfinite(data_range) or data_range <= 0:
        raise ValueError("data_range must be finite and positive")
    c1, c2 = (0.01 * data_range) ** 2, (0.03 * data_range) ** 2
    mean_a, mean_b = float(a.mean()), float(b.mean())
    var_a, var_b = float(a.var()), float(b.var())
    cov = float(np.mean((a - mean_a) * (b - mean_b)))
    return float(
        ((2 * mean_a * mean_b + c1) * (2 * cov + c2))
        / ((mean_a**2 + mean_b**2 + c1) * (var_a + var_b + c2))
    )


def normalized_cross_correlation(reference, estimate):
    a, b = _pair(reference, estimate)
    da, db = a - a.mean(), b - b.mean()
    denom = float(np.sqrt(np.sum(da * da) * np.sum(db * db)))
    if denom == 0:
        return 1.0 if np.array_equal(a, b) else 0.0
    return float(np.sum(da * db) / denom)
