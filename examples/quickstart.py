#!/usr/bin/env python3
# wave-hologram: off-axis digital holography with focus estimation.
# Copyright (C) 2026 Jean-François Brisson
# SPDX-License-Identifier: AGPL-3.0-only
"""Quickstart: record a synthetic off-axis hologram, then measure its focus.

Everything here is SIMULATED with wave_hologram.py itself. The reference
beam has an unknown piston (1.2 rad) and a fractional carrier (80.37 bins),
and the pipeline is told neither.

Run from the repository root, or after `pip install -e .`:

    python3 examples/quickstart.py [--outdir DIR]

Prints the measured focus and writes quickstart.png (focus curve,
reconstructed intensity and phase) into DIR (default: examples/output).
"""

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))  # lets this run without installing

import wave_hologram as W  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--outdir", default=os.path.join(HERE, "output"))
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args(argv)

    z_obj = 5.0e-3                         # true object distance (m)
    obj = W.make_object()                  # ring, cross and three dots
    sensor = W.angular_spectrum(obj.astype(complex), z_obj)
    holo = W.record_hologram(sensor, car=80.37, ref_phase=1.2)

    # The carrier is estimated from the hologram itself (sub-bin).
    carrier = W.estimate_carrier_2d(holo, dc_radius=30.0 / W.N)
    field0 = W.reconstruct_field(holo, carrier=carrier)

    # Focus: coarse sweep, then golden-section refinement of "realness".
    sweep = np.linspace(-10e-3, 0.0, 101)
    res = W.find_focus(field0, sweep, metric="realness")
    best = W.angular_spectrum(field0, res["z"])

    print(f"estimated carrier (ky, kx): ({carrier[0]:+.4f}, {carrier[1]:+.4f}) bins"
          "  (true: (0, -80.37))")
    print(f"truth   : {-z_obj * 1e3:+.4f} mm")
    print(f"measured: {res['z'] * 1e3:+.4f} mm   "
          f"error {abs(res['z'] + z_obj) * 1e6:.1f} um   "
          f"contrast {res['contrast']:.1f}x")
    print(f"correlation of |U|^2 with the object: "
          f"{W.correlation(np.abs(best) ** 2, obj):.3f}")

    if not args.no_plot:
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib.figure import Figure

        os.makedirs(args.outdir, exist_ok=True)
        fig = Figure(figsize=(12, 3.6), dpi=100)
        a0, a1, a2 = fig.subplots(1, 3)
        a0.plot(sweep * 1e3, res["curve"])
        a0.axvline(res["z"] * 1e3, color="r", ls="--", label="measured")
        a0.axvline(-z_obj * 1e3, color="g", ls=":", label="truth")
        a0.set(xlabel="reconstruction distance [mm]", ylabel="realness",
               title="Focus curve")
        a0.legend()
        a1.imshow(np.abs(best) ** 2, cmap="gray")
        a1.set_title("Intensity at the measured focus")
        a2.imshow(np.angle(best), cmap="twilight", vmin=-np.pi, vmax=np.pi)
        a2.set_title("Wrapped phase (offset = reference piston)")
        for a in (a1, a2):
            a.axis("off")
        fig.tight_layout()
        path = os.path.join(args.outdir, "quickstart.png")
        fig.savefig(path)
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
