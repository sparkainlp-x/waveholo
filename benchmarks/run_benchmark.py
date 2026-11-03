#!/usr/bin/env python3
"""Run deterministic synthetic focus/reconstruction conditions."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import resource
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import wave_hologram as W

from benchmarks.focus_metrics import (
    FOCUS_CRITERIA,
    entropy_focus,
    select_focus,
    tenengrad_focus,
    variance_focus,
)
from benchmarks.metrics import global_ssim, normalized_cross_correlation, psnr, rmse

REPOSITORY = "https://github.com/sparkainlp-x/waveholo"


def _object(kind, n, rng):
    y, x = np.mgrid[-1:1:complex(n), -1:1:complex(n)]
    if kind == "point":
        return np.exp(-(x * x + y * y) / (2 * 0.025**2))
    if kind == "extended":
        return (((np.abs(x) < 0.35) & (np.abs(y) < 0.12))
                | ((x * x + y * y > 0.20**2) & (x * x + y * y < 0.27**2))).astype(float)
    if kind == "textured":
        image = np.zeros((n, n), dtype=float)
        yy, xx = rng.integers(0, n, size=(2, 80))
        radii = rng.integers(max(2, n // 100), max(3, n // 35), size=80)
        for cy, cx, radius in zip(yy, xx, radii):
            image[(x * n / 2 - (cx - n / 2)) ** 2
                  + (y * n / 2 - (cy - n / 2)) ** 2 < radius**2] = 1.0
        return image
    raise ValueError(f"unknown object class: {kind}")


def _carrier(value, n, fractions):
    if isinstance(value, str):
        return int(round(float(fractions[value]) * n))
    return int(value)


def _run_one(config, *, n, carrier, kind, noise_sigma, seed):
    started = time.perf_counter()
    rng = np.random.default_rng(seed)
    obj = _object(kind, n, rng)
    distance = float(config["true_distance_m"])
    sensor = W.angular_spectrum(obj.astype(complex), distance)
    piston = float(rng.uniform(-np.pi, np.pi))
    hologram = W.record_hologram(sensor, car=carrier, ref_phase=piston)
    sigma = float(noise_sigma) * float(np.sqrt(np.mean(hologram**2)))
    noisy = hologram + rng.normal(0.0, sigma, hologram.shape)
    if config.get("clip_negative_intensity", False):
        noisy = np.maximum(0.0, noisy)
    input_checksum = hashlib.sha256(np.ascontiguousarray(noisy).tobytes()).hexdigest()
    measured_carrier = W.estimate_carrier_2d(noisy, dc_radius=0.08, half="left", refine=False)
    baseband = W.reconstruct_field(noisy, carrier=measured_carrier, win=int(config["sideband_half_width_bins"]))
    sweep = np.linspace(
        -distance - float(config["sweep_half_range_m"]),
        -distance + float(config["sweep_half_range_m"]),
        int(config["sweep_points"]),
    )
    fields = [W.angular_spectrum(baseband, float(z)) for z in sweep]
    definitions = {
        "realness": ("min", W.realness),
        "variance": ("max", variance_focus),
        "tenengrad": ("max", tenengrad_focus),
        "entropy": ("max", entropy_focus),
    }
    focus = {}
    curves = {}
    for name, (polarity, function) in definitions.items():
        scores = np.asarray([function(field) for field in fields], dtype=float)
        index, selected = select_focus(sweep, scores, polarity)
        focus[name] = {
            "polarity": polarity,
            "selected_distance_m": selected,
            "signed_error_m": selected + distance,
            "absolute_error_m": abs(selected + distance),
            "at_boundary": bool(index in (0, len(sweep) - 1)),
        }
        curves[name] = scores.tolist()

    best = fields[int(np.argmin(np.abs(sweep + distance)))]
    recovered = np.abs(best)
    rmse_value = rmse(obj, recovered)
    return {
        "seed": int(seed),
        "resolution": int(n),
        "carrier_bins": int(carrier),
        "object_class": kind,
        "noise_sigma_relative_rms": float(noise_sigma),
        "runtime_seconds": time.perf_counter() - started,
        "input_checksum": input_checksum,
        "true_distance_m": distance,
        "focus": focus,
        "metrics": {
            "rmse": rmse_value,
            "psnr_db": psnr(obj, recovered),
            "ssim_global": global_ssim(obj, recovered),
            "normalized_cross_correlation": normalized_cross_correlation(obj, recovered),
        },
        "sweep_m": sweep.tolist(),
        "curves": curves,
    }


def _write_csv(path, rows):
    columns = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def _claim_gate(rows, config, provenance_ok):
    groups = {}
    for row in rows:
        key = (
            row["resolution"], row["carrier_bins"], row["object_class"],
            row["noise_sigma_relative_rms"],
        )
        groups.setdefault(key, []).append(row)
    step = 2 * float(config["sweep_half_range_m"]) / (int(config["sweep_points"]) - 1)
    checks = []
    for condition, selected in sorted(groups.items()):
        successes = sum(
            not row["focus"]["realness"]["at_boundary"]
            and row["focus"]["realness"]["absolute_error_m"] <= step
            for row in selected
        )
        checks.append({
            "condition": condition,
            "n": len(selected),
            "successes": int(successes),
            "success_fraction": successes / len(selected),
            "minimum_seeds_met": len(selected) >= 10,
            "threshold_met": successes / len(selected) >= 0.9,
        })
    passed = provenance_ok and bool(checks) and all(
        check["minimum_seeds_met"] and check["threshold_met"] for check in checks
    )
    return {
        "scope": "synthetic_only",
        "claim_gate_status": "pass" if passed else "fail",
        "physical_claim": False,
        "provenance_valid": provenance_ok,
        "checks": checks,
    }


def _sha(path):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return digest


def _git_provenance():
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    dirty = bool(subprocess.check_output(
        ["git", "status", "--porcelain", "--untracked-files=all"], text=True
    ).strip())
    return commit, not dirty


def _write_figures(rows, output):
    figures = output / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    noise_levels = sorted({row["noise_sigma_relative_rms"] for row in rows})
    stats = []
    for level in noise_levels:
        values = np.asarray([
            row["focus"]["realness"]["absolute_error_m"] * 1e6 for row in rows
            if row["noise_sigma_relative_rms"] == level
        ])
        rng = np.random.default_rng(20261103)
        boot = values[rng.integers(0, values.size, (5000, values.size))].mean(axis=1)
        lo, hi = np.quantile(boot, [0.025, 0.975])
        stats.append((float(np.median(values)), float(lo), float(hi), values.size))
    width, height = 760, 440
    left, right, top, bottom = 85, 25, 45, 80
    xmax = max(noise_levels) or 1.0
    ymax = max(max(item[2] for item in stats) * 1.1, 1.0)
    points = []
    for noise, (median, lo, hi, _) in zip(noise_levels, stats):
        x = left + (width - left - right) * noise / xmax
        y = height - bottom - (height - top - bottom) * median / ymax
        points.append((x, y))
    polyline = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
    marks = []
    for noise, (median, lo, hi, n) in zip(noise_levels, stats):
        x = left + (width - left - right) * noise / xmax
        ylo = height - bottom - (height - top - bottom) * lo / ymax
        yhi = height - bottom - (height - top - bottom) * hi / ymax
        y = height - bottom - (height - top - bottom) * median / ymax
        marks.append(
            f'<path d="M{x:.1f},{ylo:.1f} V{yhi:.1f}" stroke="#225ea8"/>'
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5" fill="#225ea8"/>'
            f'<text x="{x:.1f}" y="{height-bottom+20}" text-anchor="middle" font-size="12">{noise:g} (n={n})</text>'
        )
    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="white"/>
<text x="{width/2}" y="22" text-anchor="middle" font-size="15">SYNTHETIC: focus error vs Gaussian noise</text>
<text x="{width/2}" y="{height-18}" text-anchor="middle" font-size="13">Noise sigma / hologram RMS</text>
<text transform="translate(20 {height/2}) rotate(-90)" text-anchor="middle" font-size="13">Median absolute focus error (µm)</text>
<path d="M{left},{top} V{height-bottom} H{width-right}" fill="none" stroke="#333"/>
<polyline points="{polyline}" fill="none" stroke="#225ea8" stroke-width="2"/>
{''.join(marks)}
<text x="{left}" y="{height-45}" font-size="11">Error bars: bootstrap 95% CI of the mean; aggregate configured scenes. Synthetic only, no physical accuracy claim.</text>
</svg>'''
    (figures / "focus_error_vs_noise.svg").write_text(svg, encoding="utf-8")
    if rows:
        row = rows[0]
        curve_points = []
        z, score = np.asarray(row["sweep_m"]), np.asarray(row["curves"]["realness"])
        zmin, zmax = float(z.min()), float(z.max())
        smin, smax = float(score.min()), float(score.max())
        for zi, si in zip(z, score):
            x = left + (width - left - right) * (zi - zmin) / max(zmax - zmin, 1e-30)
            y = height - bottom - (height - top - bottom) * (si - smin) / max(smax - smin, 1e-30)
            curve_points.append(f"{x:.1f},{y:.1f}")
        curve = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">
<rect width="100%" height="100%" fill="white"/>
<text x="{width/2}" y="22" text-anchor="middle" font-size="15">SYNTHETIC focus curve; seed {row["seed"]}</text>
<text x="{width/2}" y="{height-18}" text-anchor="middle" font-size="13">Reconstruction distance (m)</text>
<text transform="translate(20 {height/2}) rotate(-90)" text-anchor="middle" font-size="13">Realness (minimized)</text>
<path d="M{left},{top} V{height-bottom} H{width-right}" fill="none" stroke="#333"/>
<polyline points="{" ".join(curve_points)}" fill="none" stroke="#225ea8" stroke-width="2"/>
<text x="{left}" y="{height-45}" font-size="11">One deterministic synthetic realization; CI not estimable (N=1); resolution {row["resolution"]}, carrier {row["carrier_bins"]} bins, noise sigma {row["noise_sigma_relative_rms"]:g}; not physical validation.</text>
</svg>'''
        (figures / "focus_metric_curve.svg").write_text(curve, encoding="utf-8")


def run(config_path, output_path):
    config_path, output = Path(config_path), Path(output_path)
    config = json.loads(config_path.read_text(encoding="utf-8"))
    started = time.perf_counter()
    rows = []
    for n in config["resolutions"]:
        carriers = config["carriers"]
        for carrier_value in carriers:
            carrier = _carrier(carrier_value, int(n), config.get(
                "carrier_fractions", {"low": 0.24, "nominal": 0.30, "high": 0.36}
            ))
            for kind in config["objects"]:
                for noise in config["noise_sigmas"]:
                    for seed in config["seeds"]:
                        rows.append(_run_one(
                            config, n=int(n), carrier=carrier, kind=kind,
                            noise_sigma=float(noise), seed=int(seed),
                        ))
    output.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).isoformat()
    commit, tree_clean = _git_provenance()
    config_sha = _sha(config_path)
    result_rows = []
    metric_rows = []
    for row in rows:
        for method, result in row["focus"].items():
            result_rows.append({
                "protocol_version": config["protocol_version"],
                "code_commit_sha": commit,
                "timestamp_utc": timestamp,
                "seed": row["seed"],
                "software_versions": json.dumps({"python": sys.version.split()[0], "numpy": np.__version__}, sort_keys=True),
                "config_json": json.dumps(config, sort_keys=True),
                "input_checksum": row["input_checksum"],
                "resolution": row["resolution"],
                "carrier_bins": row["carrier_bins"],
                "object_class": row["object_class"],
                "noise_sigma_relative_rms": row["noise_sigma_relative_rms"],
                "metric": method,
                "focus_error_m": result["signed_error_m"],
                "focus_absolute_error_m": result["absolute_error_m"],
                "failure": result["at_boundary"],
                "runtime_seconds": row["runtime_seconds"],
            })
        for metric, value in row["metrics"].items():
            metric_rows.append({
                "protocol_version": config["protocol_version"],
                "code_commit_sha": commit,
                "timestamp_utc": timestamp,
                "seed": row["seed"],
                "software_versions": json.dumps({"python": sys.version.split()[0], "numpy": np.__version__}, sort_keys=True),
                "config_json": json.dumps(config, sort_keys=True),
                "input_checksum": row["input_checksum"],
                "resolution": row["resolution"],
                "carrier_bins": row["carrier_bins"],
                "object_class": row["object_class"],
                "noise_sigma_relative_rms": row["noise_sigma_relative_rms"],
                "metric": metric,
                "value": value,
                "runtime_seconds": row["runtime_seconds"],
            })
    _write_csv(output / "per_seed_results.csv", result_rows)
    _write_csv(output / "metrics.csv", metric_rows)
    claims = _claim_gate(rows, config, bool(commit) and tree_clean)

    summary = []
    conditions = sorted({
        (
            row["resolution"], row["carrier_bins"], row["object_class"],
            row["noise_sigma_relative_rms"],
        )
        for row in rows
    })
    for resolution, carrier, kind, noise in conditions:
        selected = [
            row for row in rows
            if (row["resolution"], row["carrier_bins"], row["object_class"],
                row["noise_sigma_relative_rms"]) == (resolution, carrier, kind, noise)
        ]
        for metric in ("rmse", "psnr_db", "ssim_global", "normalized_cross_correlation"):
            values = np.asarray([row["metrics"][metric] for row in selected], dtype=float)
            rng = np.random.default_rng(int(config["seeds"][0]))
            samples = values[rng.integers(0, values.size, size=(10000, values.size))]
            lower, upper = np.quantile(samples.mean(axis=1), [0.025, 0.975])
            summary.append({
                "metric": metric, "resolution": resolution, "carrier_bins": carrier,
                "object_class": kind, "noise_sigma_relative_rms": noise, "n": int(values.size),
                "protocol_version": config["protocol_version"], "code_commit_sha": commit,
                "timestamp_utc": timestamp, "config_json": json.dumps(config, sort_keys=True),
                "software_versions": json.dumps({"python": sys.version.split()[0], "numpy": np.__version__}, sort_keys=True),
                "input_checksum": config_sha,
                "mean": float(values.mean()), "median": float(np.median(values)),
                "standard_deviation": float(values.std(ddof=1)) if values.size > 1 else 0.0,
                "ci95_mean_lower": float(lower), "ci95_mean_upper": float(upper),
                "minimum": float(values.min()), "maximum": float(values.max()),
                "failure_rate": "",
                "runtime_seconds_mean": float(np.mean([row["runtime_seconds"] for row in selected])),
            })
        for method in FOCUS_CRITERIA:
            values = np.asarray([row["focus"][method]["signed_error_m"] for row in selected])
            failures = sum(row["focus"][method]["at_boundary"] for row in selected)
            rng = np.random.default_rng(int(config["seeds"][0]))
            samples = values[rng.integers(0, values.size, size=(10000, values.size))]
            lower, upper = np.quantile(samples.mean(axis=1), [0.025, 0.975])
            summary.append({
                "metric": f"focus_error_{method}", "resolution": resolution,
                "carrier_bins": carrier, "object_class": kind,
                "noise_sigma_relative_rms": noise, "n": int(values.size),
                "protocol_version": config["protocol_version"], "code_commit_sha": commit,
                "timestamp_utc": timestamp, "config_json": json.dumps(config, sort_keys=True),
                "software_versions": json.dumps({"python": sys.version.split()[0], "numpy": np.__version__}, sort_keys=True),
                "input_checksum": selected[0]["input_checksum"], "mean": float(values.mean()),
                "median": float(np.median(values)),
                "standard_deviation": float(values.std(ddof=1)) if values.size > 1 else 0.0,
                "ci95_mean_lower": float(lower), "ci95_mean_upper": float(upper),
                "minimum": float(values.min()), "maximum": float(values.max()),
                "failure_rate": failures / len(selected),
                "runtime_seconds_mean": float(np.mean([row["runtime_seconds"] for row in selected])),
            })
    _write_csv(output / "summary.csv", summary)
    _write_figures(rows, output)
    elapsed = time.perf_counter() - started
    environment = {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "cpu": platform.processor() or platform.machine(),
        "peak_rss_platform_units": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    (output / "config.yaml").write_text(config_path.read_text(encoding="utf-8"), encoding="utf-8")
    (output / "code_commit.txt").write_text(commit + "\n", encoding="utf-8")
    (output / "environment.json").write_text(json.dumps(environment, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {
        "protocol_version": config["protocol_version"],
        "repository_url": REPOSITORY,
        "commit_sha": commit,
        "timestamp_utc": timestamp,
        "python_version": sys.version.split()[0],
        "numpy_version": np.__version__,
        "cpu": environment["cpu"],
        "random_seeds": config["seeds"],
        "config_sha256": config_sha,
        "runtime_seconds": elapsed,
        "status": "COMPLETE" if tree_clean else "INVALID",
        "working_tree_clean": tree_clean,
        "output_checksums": {},
        "claim_gate_status": claims["claim_gate_status"],
        "software_versions": {"python": sys.version.split()[0], "numpy": np.__version__},
        "input_checksums": {
            "config": config_sha,
            "holograms": [row["input_checksum"] for row in rows],
        },
        "n_per_seed_rows": len(result_rows),
        "limitations": "SYNTHETIC only; no absolute-focus or product validation.",
        "source_hashes": {
            path.relative_to(Path(__file__).resolve().parents[1]).as_posix(): _sha(path)
            for path in [Path(__file__).resolve().parents[1] / "wave_hologram.py"]
            + sorted((Path(__file__).resolve().parent).glob("*.py"))
        },
        "peak_memory_bytes": int(
            resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            * (1024 if platform.system() == "Linux" else 1)
        ),
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output / "claims.json").write_text(json.dumps(claims, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    raw_files = [
        path for path in sorted(item for item in output.rglob("*") if item.is_file())
        if path.name not in {"manifest.json", "checksums.sha256", "raw_checksums.sha256"}
    ]
    manifest["output_checksums"] = {
        path.relative_to(output).as_posix(): _sha(path) for path in raw_files
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output / "raw_checksums.sha256").write_text(
        "\n".join(f"{_sha(path)}  {path.relative_to(output).as_posix()}" for path in raw_files) + "\n",
        encoding="utf-8",
    )
    checksums = []
    for path in sorted(item for item in output.rglob("*") if item.is_file()):
        if path.name != "checksums.sha256":
            checksums.append(f"{_sha(path)}  {path.relative_to(output).as_posix()}")
    (output / "checksums.sha256").write_text("\n".join(checksums) + "\n", encoding="utf-8")
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    rows = run(args.config, args.output)
    failures = sum(row["focus"]["realness"]["at_boundary"] for row in rows)
    claims = json.loads((Path(args.output) / "claims.json").read_text(encoding="utf-8"))
    print(
        f"Completed {len(rows)} synthetic trials; realness boundary failures={failures}; "
        f"claim gate={claims['claim_gate_status']}"
    )
    return 0 if claims["claim_gate_status"] == "pass" else 2


if __name__ == "__main__":
    raise SystemExit(main())
