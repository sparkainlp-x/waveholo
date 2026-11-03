"""Focus criteria, intentionally independent from hologram reconstruction."""
from __future__ import annotations

import numpy as np


def _amplitude(field):
    value = np.asarray(field)
    if value.ndim != 2 or value.size == 0 or not np.isfinite(value).all():
        raise ValueError("field must be a non-empty finite 2D array")
    return np.abs(value)


def variance_focus(field):
    """Intensity variance; intended as a maximized amplitude-object criterion."""
    intensity = _amplitude(field) ** 2
    return float(np.var(intensity))


def tenengrad_focus(field):
    """Mean squared amplitude gradient; larger means more high-frequency detail."""
    amp = _amplitude(field)
    gy, gx = np.gradient(amp)
    return float(np.mean(gx * gx + gy * gy))


def entropy_focus(field, bins=64):
    """Shannon entropy (bits) of a fixed-range normalized amplitude histogram."""
    if not isinstance(bins, int) or bins < 2:
        raise ValueError("bins must be an integer >= 2")
    amp = _amplitude(field)
    peak = float(amp.max())
    if peak == 0:
        return 0.0
    counts, _ = np.histogram(amp / peak, bins=bins, range=(0.0, 1.0))
    probability = counts[counts > 0] / counts.sum()
    return float(-np.sum(probability * np.log2(probability)))


FOCUS_CRITERIA = {
    "realness": ("min", None),
    "variance": ("max", variance_focus),
    "tenengrad": ("max", tenengrad_focus),
    "entropy": ("max", entropy_focus),
}


def select_focus(distances, scores, polarity):
    """Select a global score extremum and return its index and location."""
    z = np.asarray(distances, dtype=float)
    score = np.asarray(scores, dtype=float)
    if z.ndim != 1 or score.shape != z.shape or z.size == 0:
        raise ValueError("distances and scores must be matching non-empty 1D arrays")
    if not np.isfinite(z).all() or not np.isfinite(score).all():
        raise ValueError("distances and scores must be finite")
    if polarity not in {"min", "max"}:
        raise ValueError("polarity must be 'min' or 'max'")
    index = int(np.argmin(score) if polarity == "min" else np.argmax(score))
    return index, float(z[index])
