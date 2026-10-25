#!/usr/bin/env python3
"""Real-data validation on the Cuenat et al. DHM holograms (v0.1.3).

Data: Zenodo 6337535 (CC-BY-4.0), Cuenat et al., arXiv:2203.07772.
Lyncee Tec DHM, 10x NA 0.32, lambda 674.99 nm, 5.86 um pixels,
1024x1024 8-bit holograms of a 9 um pseudo-periodic PHASE pattern.
Each class is one hexapod Z position (mapping_z.txt); 400 holograms
per class. Only the classes contained in the first ~2 GB of the 48.6 GB
archive are used (tar order: class_0, class_1, class_10, class_11, ...).

For every hologram: sub-bin carrier near (281, -268), sideband win 90
cropped (lossless), then focus by four metrics:
  amplitude_variance, tamura      on the raw field
  realness, ac_realness           after removing the objective's
                                  quadratic phase (remove_quadratic_phase)
Scores: linear fit of focus vs logged stage Z (slope, residual RMS, R^2)
and WITHIN-class scatter (repeatability over many holograms at one Z).

Run: python3 dhm_v013.py --data /path/with/class_dirs [--per-class 40]
"""

import argparse
import csv
import glob
import os
import re
import sys
import time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import wave_hologram as W  # noqa: E402

WL = 674.99e-9
PITCH = 5.86e-6
HINT = (281, -268)
WIN = 90
SWEEP = np.arange(-20.0, 4.0001, 0.1) * 1e-3
METRICS = (("amplitude_variance", "raw"), ("tamura", "raw"),
           ("realness", "flat"), ("ac_realness", "flat"))


def load_mapping(path):
    z = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 2 and parts[0].isdigit():
                z[int(parts[0])] = float(parts[1])
    return z


def process(path):
    from PIL import Image
    try:
        h = np.asarray(Image.open(path)).astype(float)
    except Exception:  # truncated last member of a partial tar
        return None
    if h.shape != (1024, 1024):
        return None
    car = W.estimate_carrier_2d(h, dc_radius=0.1, half="any", hint=HINT,
                                hint_radius=12)
    raw = W.reconstruct_field(h, carrier=car, win=WIN, crop=True)
    pitch = W.cropped_pitch(h.shape, WIN, PITCH)
    flat, _ = W.remove_quadratic_phase(raw)
    out = {"file": os.path.basename(path), "ky": car[0], "kx": car[1]}
    curves = {}
    for m, which in METRICS:
        f = raw if which == "raw" else flat
        r = W.find_focus(f, SWEEP, metric=m, dx=pitch, wl=WL, tol=1e-7)
        out[m] = r["z"] * 1e3
        out[m + "_contrast"] = r["contrast"]
        curves[m] = r["curve"]
    return out, curves


def fit(x, y):
    A = np.vstack([x, np.ones_like(x)]).T
    (a, b), res, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = a * x + b
    rms = float(np.sqrt(np.mean((y - pred) ** 2)))
    n = len(x)
    se = float(np.sqrt(np.sum((y - pred) ** 2) / max(n - 2, 1)
                       / np.sum((x - x.mean()) ** 2)))
    r2 = 1 - np.sum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2)
    return float(a), float(b), rms, se, float(r2)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="/workspace/dhm_data/ex")
    ap.add_argument("--per-class", type=int, default=40)
    ap.add_argument("--outdir", default=os.path.join(HERE, "dhm_v013"))
    ap.add_argument("--mapping", default=os.path.join(HERE, "dhm_tiny",
                                                      "mapping_z.txt"))
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    zmap = load_mapping(a.mapping)
    jobs = []
    for d in sorted(glob.glob(os.path.join(a.data, "class_*"))):
        ci = int(re.findall(r"class_(\d+)$", d)[0])
        files = sorted(glob.glob(os.path.join(d, "*.tif")),
                       key=lambda p: int(re.findall(r"(\d+)\.tif$", p)[0]))
        if len(files) < 5:
            print(f"skip {d}: only {len(files)} files")
            continue
        pick = np.unique(np.linspace(0, len(files) - 1,
                                     min(a.per_class, len(files))).astype(int))
        jobs += [(ci, files[i]) for i in pick]
    print(f"{len(jobs)} holograms from classes "
          f"{sorted(set(c for c, _ in jobs))}", flush=True)
    t0 = time.time()
    with Pool(os.cpu_count()) as pool:
        results = pool.map(process, [p for _, p in jobs])
    print(f"processed in {time.time() - t0:.0f} s", flush=True)

    rows = []
    mean_curves = {}
    skipped = sum(res is None for res in results)
    if skipped:
        print(f"skipped {skipped} unreadable files")
    for (ci, _), res in zip(jobs, results):
        if res is None:
            continue
        r, curves = res
        r = dict(r, dir_class=ci, logged_z_um=zmap[ci])
        rows.append(r)
        for m, c in curves.items():
            mean_curves.setdefault((ci, m), []).append(c)
    keys = ["dir_class", "logged_z_um", "file", "ky", "kx"] + \
        [k for m, _ in METRICS for k in (m, m + "_contrast")]
    with open(os.path.join(a.outdir, "dhm_per_hologram.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in keys})

    lines = [f"DHM validation (v{W.__version__}): {len(rows)} holograms, "
             f"classes {sorted(set(r['dir_class'] for r in rows))}",
             f"carrier ky {np.mean([r['ky'] for r in rows]):+.3f} +- "
             f"{np.std([r['ky'] for r in rows]):.3f}, kx "
             f"{np.mean([r['kx'] for r in rows]):+.3f} +- "
             f"{np.std([r['kx'] for r in rows]):.3f} bins", ""]
    lines.append(f"{'metric':>20} {'slope mm/um':>12} {'+-SE':>7} "
                 f"{'|slope|':>7} {'R^2':>6} {'resid RMS':>9} "
                 f"{'within-class SD':>15} {'contrast':>8}")
    x = np.array([r["logged_z_um"] for r in rows])
    summary = {}
    for m, _ in METRICS:
        y = np.array([r[m] for r in rows])
        a_, b_, rms, se, r2 = fit(x, y)
        sds = [np.std([r[m] for r in rows if r["dir_class"] == c], ddof=1)
               for c in sorted(set(r["dir_class"] for r in rows))]
        con = np.median([r[m + "_contrast"] for r in rows])
        summary[m] = (a_, b_, rms, se, r2, float(np.median(sds)))
        lines.append(f"{m:>20} {a_:>+12.4f} {se:>7.4f} {abs(a_) * 1e3:>7.1f} "
                     f"{r2:>6.3f} {rms:>9.3f} {np.median(sds):>15.4f} {con:>8.3f}")
    lines.append("")
    lines.append("|slope| in um of hologram-plane focus per um of stage Z;"
                 " M^2 = 100 for a nominal 10x objective.")
    lines.append("")
    lines.append("per-class mean focus (mm) / SD:")
    for c in sorted(set(r["dir_class"] for r in rows)):
        sel = [r for r in rows if r["dir_class"] == c]
        lines.append(f"  class {c:>3} Z {zmap[c]:+8.3f} um n={len(sel):>3}: " +
                     "  ".join(f"{m} {np.mean([r[m] for r in sel]):+.3f}/"
                               f"{np.std([r[m] for r in sel], ddof=1):.3f}"
                               for m, _ in METRICS))
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(a.outdir, "dhm_summary.txt"), "w") as fh:
        fh.write(text + "\n")

    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.figure import Figure
    fig = Figure(figsize=(13, 9), dpi=100)
    axes = fig.subplots(2, 2).ravel()
    for ax, (m, _) in zip(axes, METRICS):
        y = np.array([r[m] for r in rows])
        ax.plot(x, y, ".", ms=4, alpha=0.5)
        a_, b_, rms, se, r2, sd = summary[m]
        xs = np.linspace(x.min(), x.max(), 10)
        ax.plot(xs, a_ * xs + b_, "r-",
                label=f"slope {a_ * 1e3:+.1f} um/um, R^2 {r2:.3f}, "
                      f"within-class SD {sd * 1e3:.0f} um")
        ax.set_title(m)
        ax.set_xlabel("logged stage Z [um]")
        ax.set_ylabel("focus [mm, hologram plane]")
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "dhm_focus_vs_stage.png"))

    fig = Figure(figsize=(13, 9), dpi=100)
    axes = fig.subplots(2, 2).ravel()
    for ax, (m, _) in zip(axes, METRICS):
        for c in sorted(set(r["dir_class"] for r in rows)):
            cur = np.mean(mean_curves[(c, m)], axis=0)
            ax.plot(SWEEP * 1e3, cur / np.median(cur), label=f"Z {zmap[c]:+.1f} um")
        ax.set_title(f"{m}: class-mean curve / median")
        ax.set_xlabel("z [mm]")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(a.outdir, "dhm_mean_curves.png"))
    print("outputs in", a.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
