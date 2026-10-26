#!/usr/bin/env python3
"""Optional live DHM benchmark harness for camera and stage adapters.

The hardware dependencies are deliberately optional.  Supply frame and stage
adapters from an experiment, or use ``--synthetic`` to validate the logging
and focus-search path without hardware.
"""
from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

import numpy as np

# Support ``python real_data/live_dhm_benchmark.py`` from the repository root.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import wave_hologram as W


def run(frames, stage_positions, output):
    """Measure focus for each 2-D hologram and write a CSV result table."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    sweep = np.linspace(-10e-3, 0.0, 101)
    rows = []
    for i, (hologram, stage_z) in enumerate(zip(frames, stage_positions)):
        carrier = W.estimate_carrier_2d(hologram, dc_radius=30.0 / hologram.shape[0])
        field = W.reconstruct_field(hologram, carrier=carrier)
        result = W.find_focus(field, sweep)
        rows.append({"frame": i, "stage_z_m": float(stage_z),
                     "focus_z_m": result["z"], "contrast": result["contrast"]})
    with output.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys() if rows else
                                ["frame", "stage_z_m", "focus_z_m", "contrast"])
        writer.writeheader(); writer.writerows(rows)
    return rows


def synthetic(count):
    obj = W.make_object().astype(complex)
    frames, stages = [], []
    for stage in np.linspace(4.8e-3, 5.2e-3, count):
        sensor = W.angular_spectrum(obj, stage)
        frames.append(W.record_hologram(sensor, car=80.37, ref_phase=1.2))
        stages.append(stage)
    return frames, stages


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--synthetic", action="store_true", help="run without camera hardware")
    parser.add_argument("--count", type=int, default=3)
    parser.add_argument("--output", default="live_dhm_benchmark.csv")
    args = parser.parse_args(argv)
    if not args.synthetic:
        parser.error("provide hardware adapters in a calling script, or use --synthetic")
    started = time.time()
    frames, stages = synthetic(args.count)
    rows = run(frames, stages, args.output)
    print(f"wrote {len(rows)} measurements to {args.output} in {time.time()-started:.2f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
