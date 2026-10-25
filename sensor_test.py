#!/usr/bin/env python3
"""Closest-to-physical test: a virtual sensor.

The demo so far used a perfect sensor. A real camera is not perfect, so
this script degrades the hologram the way hardware would, one artifact
at a time and then all together, and re-measures focus:

- vignetting: Gaussian illumination falloff on the hologram intensity.
- dust: random opaque disks on the sensor (intensity x 0.1 inside).
- gain: per-pixel fixed-pattern gain, sigma 0.02.
- quantize: 12-bit ADC rounding (peak mapped to 4095, no clipping).
- saturate: 12-bit ADC with the top 1% of pixels clipped at 4095.
- carrier: reconstruction assumes a reference tilt 2 bins off the true
  one (miscalibrated reference angle, a classic lab error).
- combined: vignetting + dust + gain + quantize + small read noise.

Each case reports the measured focus against the -5.00 mm truth and the
recovery correlation. The question is not whether artifacts hurt (they
do) but which ones move the focus plane and which ones only cost image
quality -- the same distinction the noise study found.

Run:  python3 sensor_test.py [--seed N]
"""

import argparse

import numpy as np

from wave_hologram import (
    angular_spectrum, make_object, record_hologram, reconstruct_field,
    imag_energy, correlation, Z_OBJ, N,
)
from noise_robustness import measure

SWEEP = np.linspace(-10e-3, 0.0, 101)


def vignette(shape, sigma_px=90.0):
    y, x = np.mgrid[0:shape[0], 0:shape[1]]
    cy, cx = shape[0] // 2, shape[1] // 2
    r2 = (x - cx) ** 2 + (y - cy) ** 2
    return np.exp(-r2 / (2 * sigma_px ** 2))


def dust_mask(shape, rng, count=25, rmin=2, rmax=4):
    mask = np.ones(shape)
    y, x = np.mgrid[0:shape[0], 0:shape[1]]
    for _ in range(count):
        cx, cy = rng.integers(0, shape[1]), rng.integers(0, shape[0])
        r = rng.integers(rmin, rmax + 1)
        mask[(x - cx) ** 2 + (y - cy) ** 2 < r ** 2] = 0.1
    return mask


def quantize(hologram, bits=12, saturate_pct=None):
    levels = 2 ** bits - 1
    if saturate_pct is None:
        scale = levels / float(np.max(hologram))
    else:
        scale = levels / float(np.percentile(hologram, saturate_pct))
    return np.round(np.clip(hologram * scale, 0, levels)) / scale


def run_case(name, hologram, obj, car=80, win=48):
    field0 = reconstruct_field(hologram, car=car, win=win)
    best_mm, err, contrast, corr = measure(field0, obj)
    print(f"{name:>10} {best_mm:>+9.2f} {err:>6.2f} {contrast:>8.1f}x {corr:>6.3f}")
    return best_mm, err, contrast, corr


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--seed", type=int, default=11)
    args = parser.parse_args(argv)

    rng = np.random.default_rng(args.seed)
    obj = make_object()
    sensor = angular_spectrum(obj.astype(complex), Z_OBJ)
    clean = record_hologram(sensor)

    print("case            best     err  contrast   corr   (truth -5.00 mm)")
    run_case("clean", clean, obj)
    run_case("vignette", clean * vignette(clean.shape), obj)
    run_case("dust", clean * dust_mask(clean.shape, rng), obj)
    run_case("gain", clean * rng.normal(1.0, 0.02, clean.shape), obj)
    run_case("quantize", quantize(clean), obj)
    run_case("saturate", quantize(clean, saturate_pct=99.0), obj)
    run_case("carrier+2", clean, obj, car=82)
    combined = clean * vignette(clean.shape) * dust_mask(clean.shape, rng)
    combined = combined * rng.normal(1.0, 0.02, clean.shape)
    combined = quantize(combined)
    combined = np.maximum(
        combined + rng.normal(0, 0.05 * float(np.mean(combined)),
                              combined.shape), 0)
    run_case("combined", combined, obj)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
