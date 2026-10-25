#!/usr/bin/env python3
"""Multi-class basin test on Cuenat et al. DHM holograms.

Why basins, not teeth: the pattern is a PHASE object seen through a
10x objective with a 15 um depth of field (paper). In hologram-plane
coordinates that DoF is M^2 x 15 um = 1.5 mm wide, so the
amplitude-variance focus feature is a BROAD BASIN, not a comb
tooth — fine scans (+/-20 um) correctly show nothing. The testable
claim is therefore the SLOPE: as the logged class Z steps by ~1 um
(sample space), does the basin centre move by ~1 um (direct
scaling) or ~100 um (M^2 scaling), or not at all (no signal)? A
constant offset between basin centre and logged Z is the DHM
image-plane calibration, not an error; a wrong slope is.

Carrier (281, -268) bins, window 90: validated qualitatively —
the reconstruction shows the pseudo-periodic pattern (see
dhm_c0_car281_best.png). The auto-found (70, -70) peak is a
pattern lattice point, NOT the carrier; the win=200 run
overlapped DC and the win=45 run cut the pattern fundamental
out. Both are kept in the record as wrong turns, not hidden.

Class ambiguity note: tar directories are class_0..class_12 but
files inside are named class_1..class_13 (1-indexed originals
repacked 0-indexed, or vice versa). The mapping file is
0-indexed (classes 0..99). This script uses the DIRECTORY index
as the mapping class; if the filename index is the truth, every
logged Z below is off by exactly +1 um — the SLOPE is unaffected
either way, and the report states both readings.

Usage: python3 multi_dhm.py  (uses first hologram of each class
directory present under extracted/)
"""
import glob
import os

import numpy as np
from PIL import Image

import run_dhm as R

HERE = os.path.dirname(os.path.abspath(__file__))
MAPPING = [
    -50.000000, -49.000000, -48.001440, -47.004320, -46.008640,
    -45.014400, -44.021600, -43.030240, -42.040320, -41.051840,
    -40.064800, -39.079200, -38.095040, -37.112320, -36.131040,
]


def smooth(y, w=11):
    return np.convolve(y, np.ones(w) / w, mode="same")


def basin_center(sweep, curve):
    """Parabolic fit around the smoothed argmin; returns centre in
    the sweep's units. Falls back to the argmin if degenerate."""
    sm = smooth(curve)
    # ignore convolution edge artefacts
    sm[:6] = np.inf
    sm[-6:] = np.inf
    i = int(np.argmin(sm))
    if i == 0 or i == len(sm) - 1:
        return float(sweep[i]), float(sm[i])
    y0, y1, y2 = sm[i - 1], sm[i], sm[i + 1]
    denom = y0 - 2 * y1 + y2
    if denom <= 0:
        return float(sweep[i]), float(y1)
    step = float(sweep[i + 1] - sweep[i])
    off = 0.5 * (y0 - y2) / denom
    if abs(off) > 1:
        return float(sweep[i]), float(y1)
    return float(sweep[i] + off * step), float(y1)


def main():
    sweep = np.arange(-10.0, 10.0 + 1e-9, 0.05) * 1e-3  # 401 points
    rows = []
    for d in sorted(glob.glob(os.path.join(HERE, "extracted", "class_*"))):
        cls = int(os.path.basename(d).split("_")[1])
        if cls >= len(MAPPING):
            continue
        files = sorted(glob.glob(os.path.join(d, "*.tif")))
        if not files:
            continue
        f = files[0]
        holo = np.asarray(Image.open(f)).astype(float)
        field0, pr, pc, _ = R.isolate_sideband(holo, 90,
                                               carrier=(281, -268))
        av = np.array([R.amplitude_variance(R.propagate(field0, z))
                       for z in sweep])
        centre, _ = basin_center(sweep, av)
        # NOTE (corrected 2026-10-25): an earlier version took the min
        # of the mode='same' smoothed curve INCLUDING its edges, where
        # partial sums fake a ~49% dip. Interior only: true dip ~9%.
        sm_interior = smooth(av)[11:-11]
        dip = (1 - float(np.min(sm_interior))
               / float(np.median(sm_interior))) * 100
        rows.append((cls, MAPPING[cls], os.path.basename(f), centre, dip))
        print(f"class {cls:>2} logged {MAPPING[cls]:+8.3f} um  "
              f"file {os.path.basename(f)}  basin {centre*1e3:+8.3f} mm  "
              f"smoothed dip {dip:4.1f}%", flush=True)

    print("\n--- slope fit: basin centre (mm) vs logged Z (um) ---")
    z = np.array([r[1] for r in rows])
    c = np.array([r[3] * 1e3 for r in rows])
    if len(rows) >= 2:
        A = np.polyfit(z, c, 1)
        pred = np.polyval(A, z)
        print(f"slope: {A[0]*1e3:+.1f} um of basin per um of logged Z "
              f"(direct scaling predicts +1, M^2 predicts +100)")
        print(f"  (in mm/um: {A[0]:+.4f}; M^2 scaling predicts "
              f"+0.100 mm/um, direct predicts +0.001 mm/um)")
        print(f"intercept: {A[1]:+.3f} mm (image-plane calibration)")
        print(f"residuals (mm): "
              + " ".join(f"{v:+.3f}" for v in (c - pred)))
        print(f"residual RMS: {float(np.sqrt(np.mean((c-pred)**2))):.3f} mm")
    print("\nIf filename index (not directory) is the true class, add "
          "+1 um to every logged Z above; slope unchanged.")

    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.figure import Figure
    fig = Figure(figsize=(8, 4), dpi=100)
    ax = fig.subplots(1, 1)
    ax.plot(z, c, "o-")
    ax.set_xlabel("logged class Z [um]")
    ax.set_ylabel("amp-variance basin centre [mm, hologram plane]")
    ax.set_title("DHM multi-class basin test (carrier 281,-268, win 90)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "multi_basin.png"))
    print("plot: multi_basin.png")


if __name__ == "__main__":
    main()
