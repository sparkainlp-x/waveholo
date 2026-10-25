#!/usr/bin/env python3
"""Run a Cuenat et al. DHM hologram (Zenodo 6337535) through the
wave-hologram pipeline and score it against its logged ground truth.

Dataset facts (paper arXiv:2203.07772 + Zenodo record + mapping file):
  * Experimental holograms of a pseudo-periodic PHASE pattern
    (9 um period, 100 nm Al layer), 1024x1024 px, 8-bit.
  * Lyncee Tec DHM (Fresnel off-axis family, cf. Cuche et al.),
    Leica 10x MO, NA 0.32, camera Basler acA1920-155um,
    pixel 5.86 um, wavelength 674.99 nm (SLD + 5 nm filter).
  * 100 classes x 400 holograms; mapping_z.txt gives the hexapod Z
    (um) per class, -50.000 ... +42.015 um in ~1 um steps.
  * The pattern is a PHASE object, so the package's imag-energy
    metric is NOT expected to bottom out (see REAL_DATA.md finding B);
    amplitude variance (minimum at focus for phase objects) and a
    Laplacian focus measure on intensity are tracked instead.

Scaling ambiguity (read this before believing any number):
  The logged Z is a SAMPLE-space displacement. The pipeline
  propagates the hologram-plane field with dx = 5.86 um. If the
  DHM software's reconstruction distance already lives in sample
  space, measured steps per class should be ~1 um; if the
  hologram-plane propagation must be converted by longitudinal
  magnification M^2 = 100, measured steps per class should be
  ~100 um. The script therefore reports the measured focus in
  hologram-plane coordinates AND the implied step-per-class under
  both scalings; the DATA decides which scaling holds. A constant
  offset between measured and logged planes is a calibration
  constant (where the DHM's image plane sits); a wrong STEP SIZE
  is a physics error. Score the slope, not just one point.

Comb lesson (DISCOVERY.md finding 5): the focus curve is a comb,
so a coarse grid argmin is NOT the measurement. Procedure here:
coarse sweep -> find candidate dips -> sub-micron fine scan
around each candidate dip -> report the tooth.

Usage:
  python3 run_dhm.py HOLO_FILE CLASS_Z_UM [--coarse-min MM ...]
Outputs: printed table + PNGs (spectrum, coarse curves, best focus)
in the same directory as this script.
"""
import argparse
import os
import sys

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))
from wave_hologram import imag_energy, amplitude_variance

HERE = os.path.dirname(os.path.abspath(__file__))
WL = 674.99e-9
DX = 5.86e-6


def propagate(field, z, dx=DX, wl=WL):
    """Angular spectrum, uncached (see real_fly_focus.py: the package
    cache would hold one 1024^2 transfer function per distance)."""
    n = field.shape[0]
    fx = np.fft.fftfreq(n, d=dx)
    FX, FY = np.meshgrid(fx, fx)
    k = 2 * np.pi / wl
    arg = 1.0 - (wl * FX) ** 2 - (wl * FY) ** 2
    H = np.where(arg > 0, np.exp(1j * k * z * np.sqrt(np.maximum(arg, 0))), 0)
    return np.fft.ifft2(np.fft.fft2(field) * H)


def isolate_sideband(holo, win, carrier=None):
    """2-D generalisation of the package's crop-and-shift: strongest
    off-DC spectral peak (or the forced carrier, if given), window
    cropped around it, shifted to baseband. Returns (field,
    peak_row_bins, peak_col_bins, mag)."""
    n = holo.shape[0]
    spec = np.fft.fftshift(np.fft.fft2(holo - holo.mean()))
    mag = np.abs(spec)
    c = n // 2
    yy, xx = np.mgrid[0:n, 0:n]
    r = np.sqrt((yy - c) ** 2 + (xx - c) ** 2)
    search = mag.copy()
    search[r < 60] = 0.0  # blank DC + skirt
    if carrier is not None:
        pr, pc = c + carrier[0], c + carrier[1]
    else:
        pr, pc = np.unravel_index(int(np.argmax(search)), search.shape)
    r0, r1 = max(pr - win, 0), min(pr + win, n - 1)
    c0, c1 = max(pc - win, 0), min(pc + win, n - 1)
    crop = spec[r0:r1 + 1, c0:c1 + 1]
    full = np.zeros((n, n), dtype=complex)
    full[c - (pr - r0):c + (r1 - pr) + 1,
         c - (pc - c0):c + (c1 - pc) + 1] = crop
    return np.fft.ifft2(np.fft.ifftshift(full)), int(pr - c), int(pc - c), mag


def laplacian_measure(field):
    """Mean squared Laplacian of intensity: the paper's LAP-style
    focus criterion family (their Fig. 4d uses the intensity-image
    Laplacian). For a phase object it should PEAK near focus, where
    defocus has not yet converted phase into amplitude texture...
    in practice it is tracked alongside amplitude variance and the
    data decides; do not assume its sign."""
    inten = np.abs(field) ** 2
    gy, gx = np.gradient(inten)
    gyy, _ = np.gradient(gy)
    _, gxx = np.gradient(gx)
    return float(np.mean((gxx + gyy) ** 2))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("holo", help="hologram image file")
    ap.add_argument("class_z_um", type=float,
                    help="logged hexapod Z for this hologram's class (um)")
    ap.add_argument("--win", type=int, default=200,
                    help="sideband half-width in bins (default 200)")
    ap.add_argument("--carrier", type=int, nargs=2, default=None,
                    metavar=("ROW", "COL"),
                    help="force carrier peak bins from DC, e.g. 281 -268")
    ap.add_argument("--coarse-min", type=float, default=-10.0,
                    help="coarse sweep start, mm (default -10)")
    ap.add_argument("--coarse-max", type=float, default=10.0,
                    help="coarse sweep end, mm (default +10)")
    ap.add_argument("--coarse-step", type=float, default=0.25,
                    help="coarse sweep step, mm (default 0.25)")
    ap.add_argument("--fine-halfwidth-um", type=float, default=20.0,
                    help="fine scan half-width around candidate, um")
    ap.add_argument("--fine-step-um", type=float, default=0.2,
                    help="fine scan step, um (must be << comb tooth)")
    ap.add_argument("--tag", default="dhm",
                    help="output filename tag")
    args = ap.parse_args()

    holo = np.asarray(Image.open(args.holo)).astype(float)
    if holo.ndim == 3:
        holo = holo[:, :, 0]
    print(f"hologram: {args.holo} shape {holo.shape} "
          f"range [{holo.min():.0f}, {holo.max():.0f}]")
    print(f"logged class Z: {args.class_z_um:+.3f} um "
          f"(lambda {WL*1e9:.2f} nm, dx {DX*1e6:.2f} um)")

    field0, pr, pc, mag = isolate_sideband(holo, args.win,
                                           carrier=args.carrier)
    dc_dist = float(np.hypot(pr, pc))
    print(f"carrier peak at ({pr}, {pc}) bins from DC, distance "
          f"{dc_dist:.1f} bins; window half-width {args.win} bins")
    if dc_dist < args.win + 60:
        print("WARNING: sideband window may overlap DC/twin; "
              "treat reconstruction with suspicion.")
    # Off-axis check: a clean separated sideband peak far from DC is
    # the signature of an off-axis recording; an in-line hologram
    # has no such peak.
    print(f"off-axis evidence: sideband/DC magnitude ratio "
          f"{mag[pr + holo.shape[0]//2, pc + holo.shape[0]//2] / mag.max():.3f}")

    sweep = np.arange(args.coarse_min, args.coarse_max + 1e-9,
                      args.coarse_step) * 1e-3
    ie, av, lap = [], [], []
    for d in sweep:
        f = propagate(field0, d)
        ie.append(imag_energy(f))
        av.append(amplitude_variance(f))
        lap.append(laplacian_measure(f))
    ie, av, lap = map(np.array, (ie, av, lap))
    i_ie, i_av, i_lap = map(int, (np.argmin(ie), np.argmin(av),
                                  np.argmax(lap)))
    print("\ncoarse sweep (hologram-plane mm):")
    print(f"  imag-energy min : {sweep[i_ie]*1e3:+.2f} mm "
          f"({ie[i_ie]:.4f} vs median {np.median(ie):.4f}, "
          f"contrast {np.median(ie)/ie[i_ie]:.2f}x)")
    print(f"  amp-variance min: {sweep[i_av]*1e3:+.2f} mm "
          f"({av[i_av]:.4e} vs median {np.median(av):.4e}, "
          f"dip {(1 - av[i_av]/np.median(av))*100:.1f}%)")
    print(f"  laplacian max   : {sweep[i_lap]*1e3:+.2f} mm "
          f"({lap[i_lap]:.4e} vs median {np.median(lap):.4e})")

    # Fine scan around the amp-variance dip (the phase-object
    # candidate), sub-micron steps per the comb lesson.
    z0 = sweep[i_av]
    fine = (z0 + np.arange(-args.fine_halfwidth_um,
                            args.fine_halfwidth_um + 1e-9,
                            args.fine_step_um) * 1e-6)
    av_f = np.array([amplitude_variance(propagate(field0, d))
                     for d in fine])
    i_f = int(np.argmin(av_f))
    print(f"\nfine amp-variance min around {z0*1e3:+.2f} mm: "
          f"{fine[i_f]*1e3:+.4f} mm "
          f"(dip {(1 - av_f[i_f]/np.median(av_f))*100:.1f}% "
          f"vs local median)")
    print(f"\nNOTE: logged Z is sample-space ({args.class_z_um:+.3f} um). "
          f"Compare fine focus in hologram-plane mm against it under "
          f"BOTH scalings: direct (1 um step/class) and M^2=100 "
          f"(100 um step/class). A constant offset = image-plane "
          f"calibration; report the step size across classes.")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.image
    from matplotlib.figure import Figure

    fig = Figure(figsize=(10, 4), dpi=100)
    ax = fig.subplots(1, 1)
    ax.plot(sweep * 1e3, av / np.median(av), label="amp variance / median")
    ax.plot(sweep * 1e3, ie / np.median(ie), label="imag-energy / median")
    ax.axvline(sweep[i_av] * 1e3, ls="--", label="coarse av min")
    ax.set_xlabel("propagation distance [mm, hologram plane]")
    ax.set_ylabel("metric / median")
    ax.set_title(f"DHM hologram, logged class Z {args.class_z_um:+.2f} um")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, f"{args.tag}_coarse.png"))

    f_best = propagate(field0, fine[i_f])
    inten = np.abs(f_best) ** 2
    matplotlib.image.imsave(os.path.join(HERE, f"{args.tag}_best.png"),
                            inten / inten.max())
    matplotlib.image.imsave(
        os.path.join(HERE, f"{args.tag}_spectrum.png"),
        np.log10(mag + 1) / np.log10(mag.max() + 1))
    print(f"\noutputs in {HERE}")


if __name__ == "__main__":
    main()
