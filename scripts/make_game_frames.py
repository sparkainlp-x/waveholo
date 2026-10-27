#!/usr/bin/env python3
# wave-hologram: off-axis digital holography with focus estimation.
# Copyright (C) 2026 Jean-François Brisson
# SPDX-License-Identifier: AGPL-3.0-only
"""Build the data for the "Find the focus" game (docs/game/index.html).

Every level is a SYNTHETIC scene made with the repo's own pipeline
(`wave_hologram.py`):

  object (seeded)  ->  angular-spectrum propagation to the sensor
  ->  off-axis hologram with an UNKNOWN reference piston and a FRACTIONAL
      carrier (record_hologram)  ->  sensor noise (noise_robustness.add_noise)
  ->  carrier estimated from the hologram (estimate_carrier_2d)
  ->  sideband reconstruction (reconstruct_field)
  ->  focus sweep + golden-section refinement (find_focus).

For each level this script stores
  * a z-stack of reconstructed frames (amplitude | wrapped phase), quantised
    to 24 levels and packed in one 8-bit grayscale PNG (base64);
  * the TRUE synthetic focus distance and WaveHolo's measured one;
  * the focus-metric curve and a plain intensity-variance curve, which the
    game reveals after the player locks in an answer.

The frames come from the band-limited (cropped) field, `crop=True`, which
`reconstruct_field` documents as lossless for the selected sideband.

Usage:
    python3 scripts/make_game_frames.py                   # rewrite the data
                                                          # block of index.html
    python3 scripts/make_game_frames.py --check           # verify index.html
    python3 scripts/make_game_frames.py --json out.json   # only dump JSON

Standard library + numpy only (the PNG is written with zlib).
"""

import argparse
import base64
import json
import os
import re
import struct
import sys
import zlib

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import wave_hologram as W  # noqa: E402
from noise_robustness import add_noise  # noqa: E402

HTML_DEFAULT = os.path.join(ROOT, "docs", "game", "index.html")
DATA_BEGIN = '<script id="game-data" type="application/json">'
DATA_END = "</script><!-- /game-data -->"

QLEVELS = 24          # frames are stored with 24 gray levels
N_FRAMES = 41         # frames per level
STEP_MM = 0.1         # depth step between frames
WIN = W.WIN_DEFAULT   # sideband half-width (bins), as in the default pipeline
N = W.N

#: The five levels. All numbers here are inputs; every result is computed.
LEVELS = [
    dict(id=1, title="Warm-up", family="ring, cross and dots (the repo's own object)",
         blurb="Amplitude object, almost clean sensor.",
         obj=("repo", 0), z_obj_mm=5.337, carrier=80.37, carrier_y=0.0,
         ref_phase=1.2, noise=("gaussian", 0.05), noise_seed=101,
         metric="realness", start_mm=-6.9),
    dict(id=2, title="Scattered discs", family="seven random discs",
         blurb="A different amplitude object and 0.2x Gaussian sensor noise.",
         obj=("disc", 11), z_obj_mm=6.213, carrier=78.61, carrier_y=0.4,
         ref_phase=2.4, noise=("gaussian", 0.2), noise_seed=102,
         metric="realness", start_mm=-7.6),
    dict(id=3, title="Photon-starved", family="seven random diamonds",
         blurb="Amplitude object recorded with only 100 photons per pixel (peak).",
         obj=("diamond", 13), z_obj_mm=4.482, carrier=81.25, carrier_y=-0.3,
         ref_phase=4.0, noise=("poisson", 100), noise_seed=103,
         metric="realness", start_mm=-5.9),
    dict(id=4, title="Phase object", family="binary phase pi/2 (ring, cross, dots)",
         blurb=("A transparent object: it changes the phase of the light, not "
                "the brightness. WaveHolo switches to ac_realness."),
         obj=("phase_repo", 0), z_obj_mm=5.78, carrier=79.4, carrier_y=0.2,
         ref_phase=0.7, noise=("gaussian", 0.1), noise_seed=104,
         metric="ac_realness", start_mm=-7.1),
    dict(id=5, title="Boss: heavy noise", family="seven random bars",
         blurb="Amplitude object under 0.8x Gaussian sensor noise.",
         obj=("bar", 14), z_obj_mm=5.05, carrier=80.8, carrier_y=0.1,
         ref_phase=5.5, noise=("gaussian", 0.8), noise_seed=105,
         metric="realness", start_mm=-6.3),
]


# ------------------------------------------------------------------ objects
def make_scene_object(kind, seed, n=N):
    """Seeded synthetic object (float array, or complex for phase objects)."""
    if kind in ("repo", "phase_repo"):
        amp = W.make_object(n)
        return np.exp(1j * (np.pi / 2) * amp) if kind == "phase_repo" else amp
    y, x = np.mgrid[0:n, 0:n]
    rng = np.random.default_rng(seed)
    obj = np.zeros((n, n))
    for _ in range(7):
        cx, cy = rng.integers(50, n - 50, 2)
        if kind == "disc":
            obj[np.hypot(x - cx, y - cy) < rng.integers(8, 20)] = 1.0
        elif kind == "bar":
            w, h = rng.integers(4, 40, 2)
            obj[cy - h // 2:cy + h // 2, cx - w // 2:cx + w // 2] = 1.0
        elif kind == "diamond":
            obj[(np.abs(x - cx) + np.abs(y - cy)) < rng.integers(10, 24)] = 1.0
        else:
            raise ValueError(f"unknown object kind {kind!r}")
    return obj


# --------------------------------------------------------------------- PNG
def _chunk(tag, data):
    body = tag + data
    return struct.pack(">I", len(data)) + body + struct.pack(">I", zlib.crc32(body))


def encode_png_gray8(img):
    """Encode a 2-D uint8 array as an 8-bit grayscale PNG (stdlib only).

    Each row uses PNG filter 0 or 1, whichever has the smaller sum of
    absolute (signed) bytes, which usually compresses best.
    """
    img = np.asarray(img, dtype=np.uint8)
    h, w = img.shape
    rows = []
    for r in img:
        sub = (r.astype(np.int16) - np.concatenate(([0], r[:-1]))).astype(np.uint8)
        cost = lambda a: int(np.abs(a.astype(np.int8).astype(np.int16)).sum())  # noqa: E731
        if cost(sub) < cost(r):
            rows.append(b"\x01" + sub.tobytes())
        else:
            rows.append(b"\x00" + r.tobytes())
    ihdr = struct.pack(">IIBBBBB", w, h, 8, 0, 0, 0, 0)
    return (b"\x89PNG\r\n\x1a\n" + _chunk(b"IHDR", ihdr)
            + _chunk(b"IDAT", zlib.compress(b"".join(rows), 9))
            + _chunk(b"IEND", b""))


def decode_png_gray8(data):
    """Decode a PNG written by `encode_png_gray8` (filters 0-4, 8-bit gray)."""
    assert data[:8] == b"\x89PNG\r\n\x1a\n", "not a PNG"
    pos, idat, w, h = 8, b"", 0, 0
    while pos < len(data):
        (ln,) = struct.unpack(">I", data[pos:pos + 4])
        tag, body = data[pos + 4:pos + 8], data[pos + 8:pos + 8 + ln]
        if tag == b"IHDR":
            w, h, depth, ctype = struct.unpack(">IIBB", body[:10])
            assert depth == 8 and ctype == 0, "only 8-bit grayscale supported"
        elif tag == b"IDAT":
            idat += body
        pos += 12 + ln
    raw = zlib.decompress(idat)
    out = np.zeros((h, w), dtype=np.uint8)
    prev = np.zeros(w, dtype=np.int32)
    for i in range(h):
        ft = raw[i * (w + 1)]
        cur = np.frombuffer(raw[i * (w + 1) + 1:(i + 1) * (w + 1)],
                            dtype=np.uint8).astype(np.int32)
        if ft == 1:
            cur = np.cumsum(cur) & 255
        elif ft == 2:
            cur = (cur + prev) & 255
        elif ft in (3, 4):
            left = 0
            res = np.zeros(w, dtype=np.int32)
            for j in range(w):
                up = prev[j]
                ul = prev[j - 1] if j else 0
                if ft == 3:
                    pred = (left + up) // 2
                else:
                    p = left + up - ul
                    pa, pb, pc = abs(p - left), abs(p - up), abs(p - ul)
                    pred = left if (pa <= pb and pa <= pc) else (up if pb <= pc else ul)
                left = (cur[j] + pred) & 255
                res[j] = left
            cur = res
        out[i] = cur
        prev = cur
    return out


# -------------------------------------------------------------------- build
def sweep_for(spec):
    """Depth of each frame in metres (reconstruction distance, negative)."""
    return (spec["start_mm"] + STEP_MM * np.arange(N_FRAMES)) * 1e-3


def simulate_level(spec):
    """Run the real pipeline for one level; return the numbers and frames."""
    kind, oseed = spec["obj"]
    obj = make_scene_object(kind, oseed)
    z = spec["z_obj_mm"] * 1e-3
    sensor = W.angular_spectrum(obj.astype(complex), z)
    holo = W.record_hologram(sensor, spec["carrier"], carrier_y=spec["carrier_y"],
                             ref_phase=spec["ref_phase"])
    model, level = spec["noise"]
    holo = add_noise(holo, model, level, np.random.default_rng(spec["noise_seed"]))
    carrier = W.estimate_carrier_2d(holo, dc_radius=30.0 / N)
    field0 = W.reconstruct_field(holo, carrier=carrier, win=WIN)
    sweep = sweep_for(spec)

    res = W.find_focus(field0, sweep, metric=spec["metric"])
    res_real = W.find_focus(field0, sweep, metric="realness")
    prop = W.Propagator(field0)
    intensity_var = np.array([np.var(np.abs(prop(d)) ** 2) for d in sweep])
    naive_i = int(np.argmax(intensity_var))

    # Frames: the same field on its lossless band-limited (cropped) grid.
    crop0 = W.reconstruct_field(holo, carrier=carrier, win=WIN, crop=True)
    pitch = W.cropped_pitch(holo.shape, WIN, W.DX)
    cprop = W.Propagator(crop0, dx=pitch)
    stack = [cprop(d) for d in sweep]
    truth = -z
    return dict(spec=spec, truth_m=truth, carrier=carrier, sweep=sweep, res=res,
                res_real=res_real, intensity_var=intensity_var, naive_i=naive_i,
                stack=stack)


def quantise_stack(stack):
    """Amplitude and wrapped phase of the stack -> uint8 level indices."""
    amp = np.stack([np.abs(f) for f in stack])
    scale = float(np.percentile(amp, 99.5))
    qa = np.rint(np.clip(amp / scale, 0, 1) * (QLEVELS - 1)).astype(np.uint8)
    ph = np.stack([np.angle(f) for f in stack])
    qp = np.floor((ph + np.pi) / (2 * np.pi) * QLEVELS).astype(np.int64) % QLEVELS
    return qa, qp.astype(np.uint8), scale


def _round(a, nd=6):
    return [round(float(v), nd) for v in np.asarray(a).ravel()]


def build_level(spec):
    sim = simulate_level(spec)
    qa, qp, scale = quantise_stack(sim["stack"])
    nf, fh, fw = qa.shape
    # atlas: one row-band per frame, amplitude on the left, phase on the right
    atlas = np.concatenate([qa, qp], axis=2).reshape(nf * fh, 2 * fw)
    png = encode_png_gray8(atlas)
    res, truth = sim["res"], sim["truth_m"]
    sweep_mm = sim["sweep"] * 1e3
    nearest = int(np.argmin(np.abs(sim["sweep"] - truth)))
    naive_z = float(sim["sweep"][sim["naive_i"]])
    return {
        "id": spec["id"], "title": spec["title"], "family": spec["family"],
        "blurb": spec["blurb"],
        "noise": {"model": spec["noise"][0], "level": spec["noise"][1]},
        "reference_piston_rad": spec["ref_phase"],
        "true_carrier_bins": [-spec["carrier_y"], -spec["carrier"]],
        "estimated_carrier_bins": [round(sim["carrier"][0], 4),
                                   round(sim["carrier"][1], 4)],
        "metric": spec["metric"],
        "start_mm": spec["start_mm"], "step_mm": STEP_MM, "n_frames": nf,
        "frame_w": fw, "frame_h": fh, "levels": QLEVELS,
        "amp_scale": round(scale, 6),
        "truth_mm": round(truth * 1e3, 6),
        "wavehologram": {
            "z_mm": round(res["z"] * 1e3, 6),
            "error_um": round(abs(res["z"] - truth) * 1e6, 3),
            "contrast": round(min(float(res["contrast"]), 1e6), 3),
            "realness_z_mm": round(sim["res_real"]["z"] * 1e3, 6),
            "realness_error_um": round(abs(sim["res_real"]["z"] - truth) * 1e6, 3),
        },
        "naive_intensity_variance": {
            "z_mm": round(naive_z * 1e3, 6),
            "error_um": round(abs(naive_z - truth) * 1e6, 3),
        },
        "nearest_frame_error_um": round(abs(float(sim["sweep"][nearest]) - truth) * 1e6, 3),
        "curve_wavehologram": _round(res["curve"]),
        "curve_intensity_variance": _round(sim["intensity_var"], 8),
        "atlas_png_b64": base64.b64encode(png).decode("ascii"),
        "sweep_mm": _round(sweep_mm, 4),
    }


def build_game_data(levels=None):
    return {
        "game": "Find the focus",
        "version": W.__version__,
        "banner": "Synthetic hologram game. Not a hardware or biology validation.",
        "pipeline": ("wave_hologram.py: record_hologram, estimate_carrier_2d, "
                     "reconstruct_field, find_focus"),
        "levels": [build_level(s) for s in (levels or LEVELS)],
    }


# --------------------------------------------------------------------- HTML
def dumps(data):
    # "</" must not appear inside a <script> block.
    return json.dumps(data, separators=(",", ":"), ensure_ascii=False).replace("</", "<\\/")


def extract_data(html):
    i = html.index(DATA_BEGIN) + len(DATA_BEGIN)
    j = html.index(DATA_END, i)
    return json.loads(html[i:j])


def inject(html, data):
    i = html.index(DATA_BEGIN) + len(DATA_BEGIN)
    j = html.index(DATA_END, i)
    return html[:i] + dumps(data) + html[j:]


def summary(data):
    lines = [f"{'lvl':>3} {'title':<20} {'truth mm':>10} {'WaveHolo mm':>12} "
             f"{'err um':>8} {'intensity-var err um':>21} {'metric':>12}"]
    for lv in data["levels"]:
        lines.append(f"{lv['id']:>3} {lv['title']:<20} {lv['truth_mm']:>10.3f} "
                     f"{lv['wavehologram']['z_mm']:>12.3f} "
                     f"{lv['wavehologram']['error_um']:>8.1f} "
                     f"{lv['naive_intensity_variance']['error_um']:>21.1f} "
                     f"{lv['metric']:>12}")
    return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--html", default=HTML_DEFAULT,
                    help="game HTML whose data block is rewritten")
    ap.add_argument("--json", help="write the data to this JSON file instead")
    ap.add_argument("--check", action="store_true",
                    help="regenerate and compare with the data in --html")
    args = ap.parse_args(argv)

    data = build_game_data()
    print(summary(data))
    if args.json:
        with open(args.json, "w", encoding="utf-8") as fh:
            fh.write(dumps(data))
        print(f"wrote {args.json}")
        return 0
    with open(args.html, encoding="utf-8") as fh:
        html = fh.read()
    if args.check:
        old = extract_data(html)
        ok = all(abs(a["truth_mm"] - b["truth_mm"]) < 1e-9
                 and abs(a["wavehologram"]["z_mm"] - b["wavehologram"]["z_mm"]) < 1e-3
                 for a, b in zip(old["levels"], data["levels"]))
        print("index.html data matches regeneration" if ok else "MISMATCH")
        return 0 if ok else 1
    with open(args.html, "w", encoding="utf-8") as fh:
        fh.write(inject(html, data))
    print(f"updated {args.html} ({os.path.getsize(args.html) / 1e6:.2f} MB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
