#!/usr/bin/env python3
"""Noise robustness of the wave-hologram focus measurement.

Question: how much sensor noise can the focus measurement take before it
breaks? Two noise models, both realistic for a sensor:

- gaussian: additive read noise, sigma as a fraction of the hologram mean.
- poisson: shot noise. The hologram is scaled to `photons` counts, sampled
  from Poisson, and scaled back; fewer photons = noisier. Intensity is
  clipped at zero in both models, as a real sensor cannot record a
  negative intensity.

Method: for each level, add noise, re-run the off-axis reconstruction +
focus sweep, and report the measured focus against the truth. The breaking
point is the result, not a failure to hide.

Run:  python3 noise_robustness.py [--model gaussian|poisson] [--seed N]
"""

import argparse
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
from matplotlib.figure import Figure

from wave_hologram import (
    angular_spectrum, make_object, record_hologram, reconstruct_field,
    imag_energy, correlation, Z_OBJ, OUTDIR_DEFAULT,
)

GAUSS_LEVELS = [0.0, 0.02, 0.05, 0.10, 0.20, 0.40, 0.80, 1.20, 2.00, 4.00, 8.00]
POISSON_LEVELS = [1e7, 1e6, 1e5, 1e4, 1e3, 1e2]  # peak photons per pixel
SWEEP = np.linspace(-10e-3, 0.0, 101)


def add_noise(hologram, model, level, rng):
    """Return a noisy hologram. Levels: gaussian -> sigma fraction of the
    hologram mean; poisson -> peak photons per pixel."""
    if model == "gaussian":
        if level == 0.0:
            return hologram
        return np.maximum(
            hologram + rng.normal(0, level * float(np.mean(hologram)),
                                  hologram.shape), 0)
    if model == "poisson":
        scale = level / float(np.max(hologram))
        return rng.poisson(hologram * scale) / scale
    raise ValueError(f"unknown model {model!r}")


def measure(field0, obj, sweep=SWEEP):
    """Focus sweep on a baseband field; return (best_mm, error_mm, contrast, corr)."""
    fields = [angular_spectrum(field0, d) for d in sweep]
    metrics = np.array([imag_energy(f) for f in fields])
    best_i = int(np.argmin(metrics))
    best_d = float(sweep[best_i])
    return (
        best_d * 1e3,
        abs(best_d + Z_OBJ) * 1e3,
        float(np.median(metrics) / metrics[best_i]),
        correlation(np.abs(fields[best_i]) ** 2, obj),
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--model", choices=["gaussian", "poisson"],
                        default="gaussian")
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--outdir", default=OUTDIR_DEFAULT)
    args = parser.parse_args(argv)

    rng = np.random.default_rng(args.seed)
    obj = make_object()
    sensor = angular_spectrum(obj.astype(complex), Z_OBJ)
    clean = record_hologram(sensor)
    levels = GAUSS_LEVELS if args.model == "gaussian" else POISSON_LEVELS

    os.makedirs(args.outdir, exist_ok=True)
    print(f"model {args.model}, seed {args.seed}, truth -5.00 mm")
    print(f"{'level':>10} {'best focus':>11} {'error':>8} {'contrast':>9} {'corr':>6}")
    rows = []
    for level in levels:
        noisy = add_noise(clean, args.model, level, rng)
        field0 = reconstruct_field(noisy)
        best_mm, err, contrast, corr = measure(field0, obj)
        rows.append((level, err, contrast, corr))
        print(f"{level:>10.2g} {best_mm:>+10.2f} {err:>7.2f} "
              f"{contrast:>8.1f}x {corr:>6.3f}")

    fig = Figure(figsize=(8, 4), dpi=100)
    ax = fig.subplots(1, 1)
    ax.plot([r[0] for r in rows], [r[1] for r in rows], marker="o")
    ax.set_xscale("log")
    ax.set_xlabel("noise level (sigma fraction, or peak photons)")
    ax.set_ylabel("focus error [mm]")
    ax.set_title(f"Focus error vs sensor noise ({args.model}, truth -5.00 mm)")
    fig.tight_layout()
    path = os.path.join(args.outdir, f"noise_curve_{args.model}.png")
    fig.savefig(path)
    print(f"curve: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
