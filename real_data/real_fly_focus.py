#!/usr/bin/env python3
"""First real-data run of the wave-hologram pipeline.

Hologram: 'fly 20x.bmp', an EXPERIMENTAL off-axis telecentric DHM
hologram shipped in the public pyDHM repository (MIT),
https://github.com/catrujilla/pyDHM, data/compensation samples/.
Parameters from pyDHM's own examples_compensation.py (FRS/ERS/CFS
examples): wavelength 0.532 um, pixel size 2.4 um x 2.4 um.

Adaptations vs the synthetic pipeline (wave_hologram.py), kept minimal
and documented here, not hidden:
  * Carrier is 2-D (diagonal tilt), not x-only: the sideband peak is
    located from the hologram's own Fourier spectrum and the window is
    cropped/shifted in both axes. The package's reconstruct_field only
    handles an x-only carrier, so this script implements the 2-D
    generalisation directly (same crop-and-shift idea).
  * Two focus metrics are tracked: the package's imag-energy fraction
    (assumes a real amplitude object) AND the package's gradient-energy
    sharpness (no such assumption). A fly is a thick biological sample
    (amplitude + phase), so imag-energy is NOT expected to bottom out
    at focus; sharpness is the honest metric for this object type.
  * NO ground-truth focus distance is documented for this hologram
    anywhere in pyDHM (it is an image-plane microscope recording, so
    focus is expected near z = 0, but that is an expectation, not a
    truth value). Results below are therefore MEASURED focus planes,
    not errors vs truth. Any claim of "error vs truth" on this file
    would be invented.
  * Amplitude variance is tracked as a third, diagnostic metric: for a
    PHASE object it is minimal (not maximal) at focus, because defocus
    converts phase structure into amplitude contrast. See REAL_DATA.md
    for the full reading, including why the package's gradient-energy
    metric is invariant under angular-spectrum propagation.

Outputs (this directory): fly_focus_curve.png, fly_amp_z0.png,
fly_amp_best.png, fly_spectrum_log.png, and a printed table.
"""
import os
import sys

import numpy as np
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from wave_hologram import imag_energy, field_sharpness

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = "/tmp/pydhm/data/compensation samples/fly 20x.bmp"
WL = 0.532e-6
DX = 2.4e-6
WIN = 200  # sideband half-width (FFT bins)


def propagate(field, z, dx=DX, wl=WL):
    """Angular spectrum, same formula as wave_hologram.angular_spectrum
    but WITHOUT its per-distance transfer cache: at 2048x2048 each cached
    transfer function is 67 MB, and an 81-step sweep OOM-killed the run
    (exit 137) through the cache. One distance at a time, no cache."""
    n = field.shape[0]
    fx = np.fft.fftfreq(n, d=dx)
    FX, FY = np.meshgrid(fx, fx)
    k = 2 * np.pi / wl
    arg = 1.0 - (wl * FX) ** 2 - (wl * FY) ** 2
    H = np.where(arg > 0, np.exp(1j * k * z * np.sqrt(np.maximum(arg, 0))), 0)
    return np.fft.ifft2(np.fft.fft2(field) * H)


def load_hologram(path):
    return np.asarray(Image.open(path)).astype(float)


def isolate_sideband(holo, win=WIN):
    """Find the strongest off-DC spectral peak, crop a window around it,
    shift it to baseband. Returns (field, peak_row, peak_col)."""
    n = holo.shape[0]
    spec = np.fft.fftshift(np.fft.fft2(holo - holo.mean()))
    mag = np.abs(spec)
    c = n // 2
    yy, xx = np.mgrid[0:n, 0:n]
    r = np.sqrt((yy - c) ** 2 + (xx - c) ** 2)
    search = mag.copy()
    search[r < 100] = 0.0  # blank DC + skirt
    pr, pc = np.unravel_index(int(np.argmax(search)), search.shape)
    r0, r1 = max(pr - win, 0), min(pr + win, n - 1)
    c0, c1 = max(pc - win, 0), min(pc + win, n - 1)
    crop = spec[r0:r1 + 1, c0:c1 + 1]
    full = np.zeros((n, n), dtype=complex)
    full[c - (pr - r0):c + (r1 - pr) + 1,
         c - (pc - c0):c + (c1 - pc) + 1] = crop
    field = np.fft.ifft2(np.fft.ifftshift(full))
    return field, int(pr - c), int(pc - c), mag


def main():
    holo = load_hologram(SRC)
    print(f"hologram: {SRC} shape {holo.shape} "
          f"range [{holo.min():.0f}, {holo.max():.0f}]")
    field0, pr, pc, mag = isolate_sideband(holo)
    print(f"carrier peak at ({pr}, {pc}) bins from DC, "
          f"distance {np.hypot(pr, pc):.1f} bins; window half-width {WIN}")

    sweep = np.arange(-10.0, 10.01, 0.25) * 1e-3  # metres
    ie, sh, av = [], [], []
    for d in sweep:
        f = propagate(field0, d, dx=DX, wl=WL)
        ie.append(imag_energy(f))
        sh.append(field_sharpness(f))
        av.append(float(np.var(np.abs(f) ** 2)))
    ie = np.array(ie)
    sh = np.array(sh)
    av = np.array(av)
    i_ie = int(np.argmin(ie))
    i_sh = int(np.argmax(sh))
    i_av = int(np.argmin(av))
    print(f"\nimag-energy minimum : z = {sweep[i_ie]*1e3:+.2f} mm "
          f"(value {ie[i_ie]:.4f}, median {np.median(ie):.4f})")
    print(f"sharpness maximum   : z = {sweep[i_sh]*1e3:+.2f} mm "
          f"(value {sh[i_sh]:.4e}, median {np.median(sh):.4e})")
    print(f"amp-variance minimum: z = {sweep[i_av]*1e3:+.2f} mm "
          f"(value {av[i_av]:.4e}, at z=0 {av[np.argmin(np.abs(sweep))]:.4e})")
    print(f"imag-energy at z=0  : {ie[np.argmin(np.abs(sweep))]:.4f}")
    print(f"sharpness at z=0    : {sh[np.argmin(np.abs(sweep))]:.4e}")

    # Refine around the amp-variance minimum with a 0.05 mm step
    # (sharpness is flat by construction here; refining it is noise).
    z0 = sweep[i_av]
    fine = np.arange(z0 - 0.25e-3, z0 + 0.2501e-3, 0.05e-3)
    av_f = [float(np.var(np.abs(propagate(field0, d, dx=DX, wl=WL)) ** 2))
            for d in fine]
    i_f = int(np.argmin(av_f))
    print(f"amp-variance refined: z = {fine[i_f]*1e3:+.2f} mm")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.image
    from matplotlib.figure import Figure

    fig = Figure(figsize=(10, 4), dpi=100)
    ax = fig.subplots(1, 1)
    ax.plot(sweep * 1e3, ie, label="imag-energy (min = focus if real object)")
    ax2 = ax.twinx()
    ax2.plot(sweep * 1e3, sh, color="r", label="gradient sharpness")
    ax.axvline(sweep[i_sh] * 1e3, color="r", ls="--",
               label=f"sharpness best {sweep[i_sh]*1e3:+.2f} mm")
    ax.set_xlabel("propagation distance [mm]")
    ax.set_ylabel("imag-energy fraction")
    ax2.set_ylabel("gradient sharpness")
    ax.set_title("Real fly hologram: focus curves (NO ground truth on file)")
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "fly_focus_curve.png"))

    def save_amp(z, name):
        f = propagate(field0, z, dx=DX, wl=WL)
        inten = np.abs(f) ** 2
        matplotlib.image.imsave(os.path.join(HERE, name),
                                inten / inten.max())

    save_amp(0.0, "fly_amp_z0.png")
    save_amp(float(fine[i_f]), "fly_amp_best.png")
    matplotlib.image.imsave(
        os.path.join(HERE, "fly_spectrum_log.png"),
        np.log10(mag + 1) / np.log10(mag.max() + 1))
    print(f"\noutputs in {HERE}")
    print("NOTE: no ground-truth distance exists for this hologram; "
          "the numbers above are measured planes, not errors.")


if __name__ == "__main__":
    main()
