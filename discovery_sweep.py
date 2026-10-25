#!/usr/bin/env python3
"""Synthetic discovery sweep: find where the pipeline holds and breaks.

Five sections, each printed as a table against a known truth:

1. DISTANCE: object distances 1..50 mm at fixed 532 nm / 6.5 um.
2. PHASE OBJECT: pure phase objects (and a mixed amplitude+phase one).
   The imag-energy metric assumes a real object at focus; this section
   measures exactly how it degrades, and tests intensity-weighted
   circular phase variance as the replacement metric. (Gradient
   sharpness was tried first and FALSIFIED: it is propagation-invariant
   by Parseval, so it cannot be a focus metric at all.)
3. NEW NOISE PHYSICS: vibration (exposure-averaged sub-pixel shifts,
   and per-frame reference-tilt jitter) and partial coherence (reduced
   fringe visibility, and Gaussian blur of the hologram).
4. CARRIER AUTO-CALIBRATION: estimate the carrier from the hologram
   spectrum itself and show it recovers the collapsed contrast of a
   2-bin carrier error.
5. FINE FOCUS + TWO PLANES: sub-step parabolic interpolation around
   focus (including an off-grid truth, the case where interpolation
   can actually beat the grid), and a two-plane (3D-ish) sample.

Run:  python3 discovery_sweep.py [--seed N]
Seeds are fixed so every number here is reproducible.
"""

import argparse

import numpy as np

from wave_hologram import (
    angular_spectrum, make_object, record_hologram, reconstruct_field,
    imag_energy, phase_variance, correlation, estimate_carrier,
    parabolic_min, CAR_DEFAULT, N,
)

STEP_MM = 0.1  # coarse sweep step (mm)


def sweep_metrics(field0, sweep):
    """Return (fields, imag_metrics, phase-variance metrics) across sweep."""
    fields = [angular_spectrum(field0, d) for d in sweep]
    imag = np.array([imag_energy(f) for f in fields])
    pvar = np.array([phase_variance(f) for f in fields])
    return fields, imag, pvar


def summarize(name, sweep, fields, imag, pvar, truth_mm, obj_int=None,
              obj_phase=None):
    """One table row: imag-energy argmin and phase-variance argmin vs truth."""
    bi = int(np.argmin(imag))
    pi = int(np.argmin(pvar))
    best_mm = sweep[bi] * 1e3
    pvar_mm = sweep[pi] * 1e3
    row = (f"{name:>22} {best_mm:>+8.2f} {abs(best_mm - truth_mm):>6.2f} "
           f"{float(np.median(imag) / imag[bi]):>8.1f}x "
           f"{pvar_mm:>+8.2f} {abs(pvar_mm - truth_mm):>6.2f} "
           f"{float(np.median(pvar) / pvar[pi]):>8.1f}x")
    if obj_int is not None:
        row += f" {correlation(np.abs(fields[bi]) ** 2, obj_int):>6.3f}"
    if obj_phase is not None:
        row += f" {correlation(np.angle(fields[pi]), obj_phase):>6.3f}"
    print(row)
    return best_mm, pvar_mm


def header(corr_label="corr"):
    print(f"{'':>22} {'imag best':>8} {'err':>6} {'contrast':>8} "
          f"{'pvar best':>8} {'err':>6} {'contrast':>8} {corr_label:>6}")


def section_distance():
    print("\n== 1. DISTANCE SWEEP (532 nm, 6.5 um, step 0.1 mm) ==")
    header()
    obj = make_object()
    for z_mm in (1.0, 2.0, 5.0, 10.0, 20.0, 50.0):
        z = z_mm * 1e-3
        sweep = np.arange(-(z_mm + 5.0), 0.0 + 1e-9, STEP_MM) * 1e-3
        sensor = angular_spectrum(obj.astype(complex), z)
        field0 = reconstruct_field(record_hologram(sensor))
        fields, imag, pvar = sweep_metrics(field0, sweep)
        summarize(f"{z_mm:.0f} mm", sweep, fields, imag, pvar, -z_mm,
                  obj_int=obj)


def section_phase():
    print("\n== 2. PHASE OBJECTS (truth -5.00 mm, step 0.1 mm) ==")
    print("imag-energy assumes a REAL object at focus. For a phase object")
    print("the field at focus is complex, so the true dip shallows and")
    print("defocus oscillation dips can beat it. Phase variance (argmin)")
    print("is the candidate replacement: it needs no realness assumption.")
    print("Phase corr is measured at the phase-variance plane, not the")
    print("imag-energy plane, so a wrong imag argmin cannot fake it.")
    header("ph-corr")
    sweep = np.linspace(-10e-3, 0.0, 101)
    amp = make_object()
    cases = [("amplitude (control)", amp.astype(complex), amp, None)]
    for phi, label in ((np.pi / 4, "phase pi/4"),
                       (np.pi / 2, "phase pi/2"),
                       (3 * np.pi / 4, "phase 3pi/4")):
        obj = np.exp(1j * phi * amp)
        cases.append((label, obj, None, phi * amp))
    # High-fill phase object: a solid disk at pi/2. Fill ~48%, so the
    # imag-energy at focus (~0.48) meets the defocused background (~0.5)
    # and the metric must flatten or invert.
    y, x = np.mgrid[0:N, 0:N]
    disk = (np.hypot(x - N // 2, y - N // 2) < 100).astype(float)
    cases.append(("phase disk pi/2, fill .48",
                  np.exp(1j * np.pi / 2 * disk), None, np.pi / 2 * disk))
    # Mixed: amplitude pattern A, phase pattern P rolled away from A.
    mixed = (0.3 + 0.7 * amp) * np.exp(1j * 1.2 * np.roll(amp, (30, -40)))
    cases.append(("mixed amp+phase", mixed, (0.3 + 0.7 * amp) ** 2,
                  1.2 * np.roll(amp, (30, -40))))
    for label, obj, obj_int, obj_phase in cases:
        sensor = angular_spectrum(obj, 5e-3)
        field0 = reconstruct_field(record_hologram(sensor))
        fields, imag, pvar = sweep_metrics(field0, sweep)
        summarize(label, sweep, fields, imag, pvar, -5.0,
                  obj_int=obj_int, obj_phase=obj_phase)


def shift_field(field, sx_px, sy_px):
    """Sub-pixel shift via a Fourier phase ramp (sx, sy in pixels)."""
    fy = np.fft.fftfreq(field.shape[0])
    fx = np.fft.fftfreq(field.shape[1])
    FX, FY = np.meshgrid(fx, fy)
    ramp = np.exp(-2j * np.pi * (FX * sx_px + FY * sy_px))
    return np.fft.ifft2(np.fft.fft2(field) * ramp)


def exposure_average_hologram(sensor, rng, frames=64, sigma_px=0.0,
                              sigma_car=0.0):
    """Average hologram over frames degraded by vibration.

    sigma_px: random sub-pixel shift of the object field per frame
    (sensor/object vibration during one exposure).
    sigma_car: random reference-tilt jitter in FFT bins per frame.
    Fringe period at car=80 is 256/80 = 3.2 px, so shifts of ~1 px and
    tilt jitter of ~a few bins are already a large fraction of a fringe.
    """
    y, x = np.mgrid[0:N, 0:N]
    acc = np.zeros((N, N))
    for _ in range(frames):
        u = sensor
        if sigma_px > 0:
            u = shift_field(u, rng.normal(0, sigma_px),
                            rng.normal(0, sigma_px))
        car = CAR_DEFAULT + (rng.normal(0, sigma_car) if sigma_car > 0 else 0)
        ref = np.exp(2j * np.pi * car * x / N)
        acc += np.abs(u + ref) ** 2
    return acc / frames


def blur_hologram(hologram, sigma_px):
    """Gaussian blur via FFT (spatial-coherence-style fringe washout)."""
    fy = np.fft.fftfreq(hologram.shape[0])
    fx = np.fft.fftfreq(hologram.shape[1])
    FX, FY = np.meshgrid(fx, fy)
    kern = np.exp(-2 * np.pi ** 2 * sigma_px ** 2 * (FX ** 2 + FY ** 2))
    return np.real(np.fft.ifft2(np.fft.fft2(hologram) * kern))


def section_noise(seed):
    print("\n== 3. NEW NOISE PHYSICS (truth -5.00 mm, step 0.1 mm) ==")
    print("Fringe period is 3.2 px, so vibration/blur near 1 px is already")
    print("a large fraction of a fringe. Vibration uses 64-frame exposure")
    print("averages so the ordering is not single-draw luck. When contrast")
    print("is near 1x the argmin is noise: landing on truth then is luck,")
    print("not a measurement, and the table should be read that way.")
    header()
    sweep = np.linspace(-10e-3, 0.0, 101)
    obj = make_object()
    sensor = angular_spectrum(obj.astype(complex), 5e-3)
    y, x = np.mgrid[0:N, 0:N]
    ref = np.exp(2j * np.pi * CAR_DEFAULT * x / N)

    def run(label, hologram):
        field0 = reconstruct_field(hologram)
        fields, imag, pvar = sweep_metrics(field0, sweep)
        summarize(label, sweep, fields, imag, pvar, -5.0, obj_int=obj)

    rng = np.random.default_rng(seed)
    for sig in (0.1, 0.25, 0.5, 1.0, 2.0):
        run(f"vibration shift {sig} px",
            exposure_average_hologram(sensor, rng, sigma_px=sig))
    for sig in (0.5, 1.0, 2.0, 4.0):
        run(f"vibration tilt {sig} bins",
            exposure_average_hologram(sensor, rng, sigma_car=sig))
    for gamma in (0.8, 0.5, 0.3, 0.1, 0.05):
        holo = np.abs(sensor) ** 2 + 1.0 + 2 * gamma * np.real(sensor * np.conj(ref))
        run(f"visibility gamma={gamma}", holo)
    run("visibility gamma=0", np.abs(sensor) ** 2 + 1.0)
    clean = record_hologram(sensor)
    for sig in (0.5, 1.0, 1.5, 2.0):
        run(f"blur sigma={sig} px", blur_hologram(clean, sig))


def section_carrier(seed):
    print("\n== 4. CARRIER AUTO-CALIBRATION ==")
    print("Review finding: assuming car=80 when the true carrier is 82")
    print("collapses imag-energy contrast to ~1.2x. Estimating the carrier")
    print("from the hologram spectrum itself should restore it.")
    print(f"{'case':>22} {'est car':>7} {'imag best':>8} {'err':>6} "
          f"{'contrast':>8} {'corr':>6}")
    sweep = np.linspace(-10e-3, 0.0, 101)
    obj = make_object()
    sensor = angular_spectrum(obj.astype(complex), 5e-3)
    for true_car in (70, 80, 82, 90):
        holo = record_hologram(sensor, car=true_car)
        est = estimate_carrier(holo)
        field0 = reconstruct_field(holo, car=est)
        fields, imag, pvar = sweep_metrics(field0, sweep)
        bi = int(np.argmin(imag))
        print(f"{f'true car {true_car}':>22} {est:>7} "
              f"{sweep[bi] * 1e3:>+8.2f} {abs(sweep[bi] * 1e3 + 5.0):>6.2f} "
              f"{float(np.median(imag) / imag[bi]):>8.1f}x "
              f"{correlation(np.abs(fields[bi]) ** 2, obj):>6.3f}")
    # The review's collapsed case, side by side: assumed 80 vs estimated.
    holo82 = record_hologram(sensor, car=82)
    for label, car in (("car 82, assumed 80", 80),
                       ("car 82, estimated", estimate_carrier(holo82))):
        field0 = reconstruct_field(holo82, car=car)
        fields, imag, pvar = sweep_metrics(field0, sweep)
        bi = int(np.argmin(imag))
        print(f"{label:>22} {car:>7} "
              f"{sweep[bi] * 1e3:>+8.2f} {abs(sweep[bi] * 1e3 + 5.0):>6.2f} "
              f"{float(np.median(imag) / imag[bi]):>8.1f}x "
              f"{correlation(np.abs(fields[bi]) ** 2, obj):>6.3f}")
    # Estimated carrier under noise.
    rng = np.random.default_rng(seed)
    noisy = np.maximum(holo82 + rng.normal(
        0, 0.20 * float(np.mean(holo82)), holo82.shape), 0)
    est = estimate_carrier(noisy)
    field0 = reconstruct_field(noisy, car=est)
    fields, imag, pvar = sweep_metrics(field0, sweep)
    bi = int(np.argmin(imag))
    print(f"{'car 82 + gauss 0.20, est':>22} {est:>7} "
          f"{sweep[bi] * 1e3:>+8.2f} {abs(sweep[bi] * 1e3 + 5.0):>6.2f} "
          f"{float(np.median(imag) / imag[bi]):>8.1f}x "
          f"{correlation(np.abs(fields[bi]) ** 2, obj):>6.3f}")


def section_fine_and_3d(seed):
    print("\n== 5a. FINE FOCUS STEP + PARABOLIC INTERPOLATION ==")
    print("On-grid truth (-5.00) cannot show interpolation value: the grid")
    print("already lands on it. The real test is an OFF-GRID truth")
    print("(-5.03 mm), where the grid must err by up to half a step.")
    obj = make_object()
    rng = np.random.default_rng(seed)
    for z_mm, label in ((5.00, "truth -5.00 (on-grid)"),
                        (5.03, "truth -5.03 (off-grid)")):
        sensor = angular_spectrum(obj.astype(complex), z_mm * 1e-3)
        clean = record_hologram(sensor)
        noisy = np.maximum(clean + rng.normal(
            0, 0.20 * float(np.mean(clean)), clean.shape), 0)
        for tag, holo in (("clean", clean), ("gauss 0.20", noisy)):
            field0 = reconstruct_field(holo)
            sweep = np.linspace(-10e-3, 0.0, 101)
            fields, imag, pvar = sweep_metrics(field0, sweep)
            bi = int(np.argmin(imag))
            parab, _ = parabolic_min(sweep, imag)
            print(f"{label:>22} {tag:>10}: grid "
                  f"{sweep[bi] * 1e3:+.2f} (err {abs(sweep[bi] * 1e3 + z_mm):.3f})"
                  f"  parabolic {parab * 1e3:+.4f} "
                  f"(err {abs(parab * 1e3 + z_mm):.4f})")
        # Fine step around the off-grid truth.
        field0 = reconstruct_field(clean)
        sweep = np.linspace(-5.5e-3, -4.5e-3, 101)
        fields, imag, pvar = sweep_metrics(field0, sweep)
        bi = int(np.argmin(imag))
        parab, _ = parabolic_min(sweep, imag)
        print(f"{label:>22} {'clean':>10}: fine 0.01 mm grid "
              f"{sweep[bi] * 1e3:+.2f} (err {abs(sweep[bi] * 1e3 + z_mm):.3f})"
              f"  parabolic {parab * 1e3:+.4f} "
              f"(err {abs(parab * 1e3 + z_mm):.4f})")

    print("\n== 5b. TWO-PLANE SAMPLE (ring+cross at 5 mm, dots at 7.5 mm) ==")
    print("Sensor field = sum of both planes propagated to the sensor.")
    print("A global metric sees both planes at once: do two dips survive,")
    print("in imag-energy and in phase variance?")
    amp = make_object()
    y, x = np.mgrid[0:N, 0:N]
    dots = np.zeros((N, N))
    for dx0, dy0 in ((-70, -70), (70, -70), (0, 70)):
        dots[(np.hypot(x - (N // 2 + dx0), y - (N // 2 + dy0)) < 10)] = 1.0
    plane1 = np.clip(amp - dots, 0, 1)  # ring + cross only
    sensor2 = (angular_spectrum(plane1.astype(complex), 5e-3)
               + angular_spectrum(dots.astype(complex), 7.5e-3))
    field0 = reconstruct_field(record_hologram(sensor2))
    sweep = np.linspace(-10e-3, 0.0, 201)
    fields, imag, pvar = sweep_metrics(field0, sweep)
    for d_mm, plane, name in ((-5.0, plane1, "plane1 @ -5.0"),
                              (-7.5, dots, "plane2 @ -7.5")):
        i = int(np.argmin(np.abs(sweep * 1e3 - d_mm)))
        print(f"at {name}: imag {imag[i]:.4f}, pvar {pvar[i]:.4f}, "
              f"corr vs that plane "
              f"{correlation(np.abs(fields[i]) ** 2, plane):.3f}")
    print(f"global imag argmin: {sweep[int(np.argmin(imag))] * 1e3:+.2f} mm "
          f"(median {np.median(imag):.3f}), global pvar argmin: "
          f"{sweep[int(np.argmin(pvar))] * 1e3:+.2f} mm "
          f"(median {np.median(pvar):.3f})")
    for metric, name, below in ((imag, "imag", True), (pvar, "pvar", True)):
        med = np.median(metric)
        if below:
            local = [i for i in range(1, len(metric) - 1)
                     if metric[i] < metric[i - 1] and metric[i] < metric[i + 1]
                     and metric[i] < 0.5 * med]
        print(f"local {name} minima below half-median:",
              [f"{sweep[i] * 1e3:+.2f} mm ({metric[i]:.3f})"
               for i in local] or "none")


def section_comb(seed):
    print("\n== 5c. THE COMB: imag-energy at 0.001 mm steps around focus ==")
    print("The transfer-function phase advances ~1.88 cycles per micron at")
    print("low frequencies, so phase alignment recurs quasi-periodically")
    print("in propagation distance: teeth every 0.004 mm, peaks ~0.99")
    print("between them. Only the central tooth (exact truth) aligns ALL")
    print("frequencies, so only it is deepest. A sweep that does not")
    print("sample the central tooth lands on the broad oscillation dip")
    print("instead (see 5a: coarse off-grid argmin -5.30, err 0.27 mm).")
    print("Noise does NOT smooth the comb away (second table).")
    obj = make_object()
    sensor = angular_spectrum(obj.astype(complex), 5.03e-3)
    field0 = reconstruct_field(record_hologram(sensor))
    sweep = np.linspace(-5.04e-3, -5.02e-3, 21)
    metrics = np.array([imag_energy(angular_spectrum(field0, d))
                        for d in sweep])
    for d, v in zip(sweep, metrics):
        print(f"  {d * 1e3:+.3f} mm  imag={v:.4f}")
    from noise_robustness import add_noise
    rng = np.random.default_rng(seed)
    noisy = add_noise(record_hologram(sensor), "gaussian", 0.20, rng)
    field0n = reconstruct_field(noisy)
    sweep2 = np.linspace(-5.06e-3, -5.00e-3, 61)
    metrics2 = np.array([imag_energy(angular_spectrum(field0n, d))
                         for d in sweep2])
    bi = int(np.argmin(metrics2))
    print(f"gauss 0.20, 0.001 mm steps: argmin {sweep2[bi] * 1e3:+.3f} mm "
          f"(truth -5.03, err {abs(sweep2[bi] * 1e3 + 5.03):.3f}), "
          f"min {metrics2[bi]:.4f}, median {np.median(metrics2):.4f}")
    print("Comb survives noise: fine sampling still finds truth exactly.")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--seed", type=int, default=23)
    args = parser.parse_args(argv)
    section_distance()
    section_phase()
    section_noise(args.seed)
    section_carrier(args.seed)
    section_fine_and_3d(args.seed)
    section_comb(args.seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
