#!/usr/bin/env python3
"""wave-hologram: off-axis digital holography with a measured focus.

Pipeline (v0.1.3)
-----------------
1. An object at z = 0 (amplitude, phase, or mixed).
2. Angular-spectrum propagation to a sensor plane; an off-axis hologram
   H = |U + R|^2 is recorded with a tilted reference R. R may carry an
   arbitrary constant phase (piston), as every real reference does.
3. Reconstruction: the carrier is estimated from the hologram itself to
   SUB-BIN precision, the object sideband is demodulated in real space
   by that fractional carrier, and then low-passed. Nothing assumes an
   integer carrier or an x-only tilt.
4. Focus metric: REALNESS, the global-phase-invariant imaginary-energy
   fraction  0.5 * (1 - |sum U^2| / sum |U|^2).  It is 0 for any real
   field times any constant phase factor. AC-REALNESS applies the same
   test to the mean-subtracted field, which also makes binary phase
   objects and single-material mixed objects "real up to a constant".
5. Focus search: a coarse sweep finds the basin, then golden-section
   refinement finds the continuous minimum. The metric is smooth, so
   no sub-micron scan is needed.

Why the change from v0.1.2
--------------------------
v0.1.2 used imag_energy = mean(Im^2)/mean(|U|^2). That quantity depends
on the field's global phase: at perfect focus it equals sin^2(theta) for
a field e^{i theta} * (real). Propagation multiplies the low frequencies
by e^{ikz}, so the metric oscillates with period lambda/2 = 0.266 um in
z. Sampled at 1 um steps, that oscillation aliases to the "teeth every
0.004 mm" reported as a physical comb in v0.1.2. The simulator happened
to hand the reconstruction a zero piston at true focus; a real reference
does not. imag_energy is kept for comparison, with that warning.

Why off-axis, not in-line
--------------------------
With an in-line hologram (R = 1) the recorded field is real, so the
reconstruction intensity at +d equals that at -d: the twin image makes
the focus distance ambiguous in sign. Off-axis recording separates the
terms in Fourier space (Leith-Upatnieks) and removes the twin.

Dependencies: numpy, matplotlib (plots), ffmpeg (MP4 only).
Run:   python3 wave_hologram.py [--z-obj MM] [--metric NAME] [--no-video]
Tests: python3 -m unittest discover -s tests -v
"""

import argparse
import os
import subprocess
import sys
from collections import OrderedDict

import numpy as np

__version__ = "0.1.4"

# ------------------------------------------------------------------ config
N = 256                # grid size (pixels, square) for the synthetic demo
DX = 6.5e-6            # pixel pitch (m)
WL = 532e-9            # wavelength (m)
Z_OBJ_DEFAULT = 5.0e-3  # object -> sensor distance (m)
Z_OBJ = Z_OBJ_DEFAULT  # alias kept for the tests and earlier scripts
CAR_DEFAULT = 80       # off-axis carrier (FFT bins, x tilt)
WIN_DEFAULT = 48       # sideband half-width (FFT bins)
FPS_DEFAULT = 20

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR_DEFAULT = os.path.join(HERE, "output")


# ------------------------------------------------------------- propagation
# v0.1.2 cached one full transfer function per distance, which OOM-killed
# a 2048^2 real-data sweep (81 x 67 MB). v0.1.3 caches only the kz grid
# per GEOMETRY (shape, pitch, wavelength) in a small LRU; the transfer
# function for a distance is then a single complex exponential.
_KZ_CACHE = OrderedDict()
_KZ_CACHE_MAX = 6


def _pitch(dx):
    """Return (dy, dx) from a scalar or a 2-tuple pitch."""
    if np.ndim(dx) == 0:
        return float(dx), float(dx)
    dy, dxx = dx
    return float(dy), float(dxx)


def _kz(shape, dx=DX, wl=WL):
    """Longitudinal wavenumber grid and propagating-wave mask (cached)."""
    if isinstance(shape, (int, np.integer)):
        shape = (int(shape), int(shape))
    dy, dxx = _pitch(dx)
    key = (tuple(shape), dy, dxx, float(wl))
    if key in _KZ_CACHE:
        _KZ_CACHE.move_to_end(key)
        return _KZ_CACHE[key]
    fy = np.fft.fftfreq(shape[0], d=dy)
    fx = np.fft.fftfreq(shape[1], d=dxx)
    FX, FY = np.meshgrid(fx, fy)
    arg = 1.0 - (wl * FX) ** 2 - (wl * FY) ** 2
    mask = arg > 0
    kz = (2 * np.pi / wl) * np.sqrt(np.maximum(arg, 0.0))
    _KZ_CACHE[key] = (kz, mask)
    while len(_KZ_CACHE) > _KZ_CACHE_MAX:
        _KZ_CACHE.popitem(last=False)
    return kz, mask


def _transfer(z, n=N, dx=DX, wl=WL):
    """Angular-spectrum transfer function for distance z.

    Evanescent components (argument of the square root <= 0) are set to
    zero: they carry no energy to the sensor.
    """
    kz, mask = _kz(n, dx, wl)
    return np.where(mask, np.exp(1j * z * kz), 0)


def angular_spectrum(field, z, dx=DX, wl=WL):
    """Propagate a complex 2-D field by distance z (negative = back).

    Rectangular fields are allowed (v0.1.3); dx may be a scalar or a
    (dy, dx) pair. Raises ValueError on non-finite input or non-2-D input.
    """
    field = np.asarray(field, dtype=complex)
    if field.ndim != 2:
        raise ValueError("field must be a 2D array")
    if not np.isfinite(field).all() or not np.isfinite(z):
        raise ValueError("field and z must be finite")
    return np.fft.ifft2(np.fft.fft2(field) * _transfer(z, field.shape, dx, wl))


class Propagator:
    """Repeated propagation of ONE field to many distances.

    The forward FFT is done once; each distance costs one exponential
    and one inverse FFT. Use this for focus sweeps.
    """

    def __init__(self, field0, dx=DX, wl=WL):
        field0 = np.asarray(field0, dtype=complex)
        if field0.ndim != 2 or not np.isfinite(field0).all():
            raise ValueError("field0 must be a finite 2D array")
        self.shape = field0.shape
        self.kz, mask = _kz(field0.shape, dx, wl)
        self.spec = np.fft.fft2(field0) * mask
        self.dx, self.wl = dx, wl

    def __call__(self, z):
        if not np.isfinite(z):
            raise ValueError("z must be finite")
        return np.fft.ifft2(self.spec * np.exp(1j * z * self.kz))


# ------------------------------------------------------------------ objects
def make_object(n=N):
    """Amplitude object: a ring, a cross, and three dots (one plane)."""
    y, x = np.mgrid[0:n, 0:n]
    cx = cy = n // 2
    r = np.hypot(x - cx, y - cy)
    obj = np.zeros((n, n))
    obj[(r > 40) & (r < 48)] = 1.0
    obj[cy - 3:cy + 4, cx - 60:cx + 61] = 1.0
    obj[cy - 60:cy + 61, cx - 3:cx + 4] = 1.0
    for dx0, dy0 in ((-70, -70), (70, -70), (0, 70)):
        obj[(np.hypot(x - (cx + dx0), y - (cy + dy0)) < 10)] = 1.0
    return obj


# ------------------------------------------------------- record / reconstruct
def record_hologram(sensor_field, car=CAR_DEFAULT, n=None, dx=DX,
                    carrier_y=0.0, ref_phase=0.0, ref_amp=1.0):
    """Off-axis hologram intensity H = |U + R|^2.

    R = ref_amp * exp(i (2 pi (car x / nx + carrier_y y / ny) + ref_phase)).
    car and carrier_y may be fractional (real carriers always are).
    ref_phase is the reference piston: the reconstructed field carries
    exp(-i ref_phase), which is exactly what breaks imag_energy and what
    realness is invariant to. dx is unused (kept for the v0.1.2 API).
    """
    sensor_field = np.asarray(sensor_field, dtype=complex)
    ny, nx = sensor_field.shape
    n = nx if n is None else n
    if not (0 < car < n // 2):
        raise ValueError(f"carrier must satisfy 0 < car < {n // 2}")
    y, x = np.mgrid[0:ny, 0:nx]
    ref = ref_amp * np.exp(1j * (2 * np.pi * (car * x / nx + carrier_y * y / ny)
                                 + ref_phase))
    return np.abs(sensor_field + ref) ** 2


def _smooth3(a):
    pad = np.pad(a, 1, mode="edge")
    return (pad[:-2, :-2] + pad[:-2, 1:-1] + pad[:-2, 2:]
            + pad[1:-1, :-2] + pad[1:-1, 1:-1] + pad[1:-1, 2:]
            + pad[2:, :-2] + pad[2:, 1:-1] + pad[2:, 2:]) / 9.0


def _local_dtft(h, ky0, kx0, half, step):
    """|DTFT| of h on a local (ky, kx) grid around (ky0, kx0), in bins."""
    ny, nx = h.shape
    v = ky0 + np.arange(-half, half + step / 2, step)
    u = kx0 + np.arange(-half, half + step / 2, step)
    Ey = np.exp(-2j * np.pi * np.outer(v, np.arange(ny)) / ny)
    Ex = np.exp(-2j * np.pi * np.outer(u, np.arange(nx)) / nx)
    F = np.abs(Ey @ h @ Ex.T)
    i, j = np.unravel_index(int(np.argmax(F)), F.shape)
    return float(v[i]), float(u[j]), F


def estimate_carrier_2d(hologram, dc_radius=0.12, half="left", hint=None,
                        hint_radius=6, refine=True, exclude=None):
    """Estimate the sideband position (ky, kx) in FFT bins, sub-bin precise.

    Convention: the returned (ky, kx) is where the U R* sideband sits in
    the hologram spectrum (signed bins from DC). For the package's x-only
    reference with carrier `car`, this is (0, -car).

    Steps: mean-subtract, 3x3-smoothed |FFT|, blank DC within dc_radius
    (cycles/pixel), keep one half-plane (half='left': kx < 0, 'top':
    ky < 0, 'any'), optionally blank `exclude` = list of ((ky, kx), r)
    discs (e.g. a periodic object's own lattice), or search only within
    hint_radius of a `hint` (ky, kx). The integer argmax is then refined
    by a local upsampled DFT (two passes: 0.05 bin, then 0.0025 bin),
    which locates the continuous spectral peak, not the nearest bin.
    """
    h = np.asarray(hologram, dtype=float)
    if h.ndim != 2:
        raise ValueError("hologram must be a 2D array")
    h = h - h.mean()
    ny, nx = h.shape
    cy, cx = ny // 2, nx // 2
    mag = _smooth3(np.abs(np.fft.fftshift(np.fft.fft2(h))))
    V, U = np.mgrid[0:ny, 0:nx]
    V = V - cy
    U = U - cx
    rad = np.hypot(V / ny, U / nx)
    search = mag.copy()
    search[rad < dc_radius] = 0.0
    if half == "left":
        search[U >= 0] = 0.0
    elif half == "top":
        search[V >= 0] = 0.0
    elif half != "any":
        raise ValueError("half must be 'left', 'top' or 'any'")
    if exclude:
        for (ey, ex), er in exclude:
            search[np.hypot(V - ey, U - ex) < er] = 0.0
    if hint is not None:
        search[np.hypot(V - hint[0], U - hint[1]) > hint_radius] = 0.0
    if not np.any(search > 0):
        raise ValueError("no sideband candidate left after masking")
    pr, pc = np.unravel_index(int(np.argmax(search)), search.shape)
    ky, kx = float(pr - cy), float(pc - cx)
    if refine:
        ky, kx, _ = _local_dtft(h, ky, kx, 1.0, 0.05)
        ky, kx, _ = _local_dtft(h, ky, kx, 0.05, 0.0025)
    return ky, kx


def estimate_carrier(hologram, subbin=True):
    """x-carrier of the package's x-only reference, in FFT bins.

    v0.1.3: returns the SUB-BIN estimate (float) by default. v0.1.2
    rounded to an integer, which leaves a residual phase ramp of up to
    half a bin; that ramp moved imag-energy focus by up to 0.8 mm.
    subbin=False restores the old integer behaviour.
    """
    hologram = np.asarray(hologram, dtype=float)
    if hologram.ndim != 2 or hologram.shape[0] != hologram.shape[1]:
        raise ValueError("hologram must be a square 2D array")
    n = hologram.shape[0]
    # Same search band as v0.1.2: 30..120 bins left of DC.
    excl = None
    ky, kx = estimate_carrier_2d(hologram, dc_radius=30.0 / n, half="left",
                                 refine=subbin, exclude=excl)
    if not subbin:
        return int(round(-kx))
    return float(-kx)


def reconstruct_field(hologram, car=CAR_DEFAULT, win=WIN_DEFAULT, n=None,
                      carrier=None, window="square", crop=False):
    """Isolate the object sideband and bring it to baseband.

    carrier : (ky, kx) sideband position in bins (from estimate_carrier_2d).
              If None, the x-only legacy carrier (0, -car) is used.
              Fractional values are handled EXACTLY by demodulating in
              real space (multiplying by exp(-2 pi i (ky y/ny + kx x/nx)))
              before the low-pass; there is no rounding to a bin.
    win     : half-width of the low-pass window in bins (int, or (wy, wx)).
    window  : 'square' (v0.1.2 behaviour) or 'circle'.
    crop    : if True, return the band-limited field on a reduced grid of
              (2*wy+2, 2*wx+2) samples. This is lossless for the selected
              band and makes real-data sweeps 10-100x cheaper. The pixel
              pitch of the cropped field is dx * n / m (see `cropped_pitch`).

    Frequencies whose ORIGINAL position (before demodulation) lies outside
    the measured spectrum stay zero instead of wrapping in from the other
    side (the v0.1.2 edge-window rule, kept).
    """
    h = np.asarray(hologram, dtype=float)
    ny, nx = h.shape
    if carrier is None:
        nn = nx if n is None else n
        if not (0 < car < nn) or not (0 <= np.max(win) < nn // 2):
            raise ValueError("need 0 < carrier < n and 0 <= window < n/2")
        carrier = (0.0, -float(car))
    ky, kx = float(carrier[0]), float(carrier[1])
    wy, wx = (win, win) if np.ndim(win) == 0 else win
    if wy < 0 or wx < 0 or wy >= ny // 2 or wx >= nx // 2:
        raise ValueError("window must satisfy 0 <= win < n/2")
    if ky == round(ky) and kx == round(kx):
        # Integer carrier: an exact circular shift of the spectrum,
        # bit-identical to v0.1.2's crop-and-shift.
        spec = np.roll(np.fft.fftshift(np.fft.fft2(h)),
                       (-int(round(ky)), -int(round(kx))), axis=(0, 1))
    else:
        y = np.arange(ny)[:, None]
        x = np.arange(nx)[None, :]
        demod = h * np.exp(-2j * np.pi * (ky * y / ny + kx * x / nx))
        spec = np.fft.fftshift(np.fft.fft2(demod))
    cy, cx = ny // 2, nx // 2
    V, U = np.mgrid[0:ny, 0:nx]
    V = V - cy
    U = U - cx
    if window == "square":
        keep = (np.abs(V) <= wy) & (np.abs(U) <= wx)
    elif window == "circle":
        keep = (V / max(wy, 1)) ** 2 + (U / max(wx, 1)) ** 2 <= 1.0
    else:
        raise ValueError("window must be 'square' or 'circle'")
    oy, ox = V + ky, U + kx
    keep &= (oy >= -ny / 2) & (oy < ny / 2) & (ox >= -nx / 2) & (ox < nx / 2)
    if not keep.any():
        raise ValueError("sideband window does not overlap the spectrum")
    spec = np.where(keep, spec, 0)
    if not crop:
        return np.fft.ifft2(np.fft.ifftshift(spec))
    my, mx = 2 * wy + 2, 2 * wx + 2
    sub = spec[cy - my // 2:cy + my // 2, cx - mx // 2:cx + mx // 2]
    return np.fft.ifft2(np.fft.ifftshift(sub)) * (my * mx) / (ny * nx)


def cropped_pitch(shape, win, dx):
    """Pixel pitch (dy, dx) of a crop=True reconstruction."""
    ny, nx = shape
    wy, wx = (win, win) if np.ndim(win) == 0 else win
    dy0, dx0 = _pitch(dx)
    return dy0 * ny / (2 * wy + 2), dx0 * nx / (2 * wx + 2)


def offaxis_field(sensor_field, car=CAR_DEFAULT, win=WIN_DEFAULT,
                  ref_phase=0.0):
    """Record an off-axis hologram and return the baseband object field."""
    return reconstruct_field(record_hologram(sensor_field, car,
                                             ref_phase=ref_phase), car, win)


def remove_quadratic_phase(field, iters=3):
    """Fit and remove a 2-D quadratic + tilt phase (no unwrapping needed).

    Microscope-objective DHMs (e.g. Lyncee Tec) add a large spherical
    phase curvature to the reconstructed field; with it no field is real
    at any plane, so realness / ac_realness are flat (measured 0.43-0.45
    on the Cuenat data). The fit uses wrapped phase DIFFERENCES between
    neighbours, weighted by amplitude, so it needs no phase unwrapping:
    a quadratic phase has a linear gradient, fitted by least squares,
    iterated on the residual. Returns (corrected_field, removed_phase).
    A removed quadratic acts as a thin lens, so focus distances measured
    after removal are in a shifted/scaled coordinate: compare them
    DIFFERENTIALLY (slopes), as the DHM study does.
    """
    g = np.asarray(field, dtype=complex).copy()
    ny, nx = g.shape
    y, x = np.mgrid[0:ny, 0:nx]
    y = (y - ny / 2) / ny
    x = (x - nx / 2) / nx
    total = np.zeros(g.shape)
    xx = (x[:, 1:] + x[:, :-1]) / 2
    yx = y[:, 1:]
    xy = x[1:, :]
    yy = (y[1:, :] + y[:-1, :]) / 2
    one_x, one_y = np.ones_like(xx), np.ones_like(yy)
    A1 = np.stack([2 * xx / nx, 0 * xx, yx / nx, one_x / nx, 0 * xx], -1).reshape(-1, 5)
    A2 = np.stack([0 * yy, 2 * yy / ny, xy / ny, 0 * yy, one_y / ny], -1).reshape(-1, 5)
    for _ in range(iters):
        gx = np.angle(g[:, 1:] * np.conj(g[:, :-1]))
        gy = np.angle(g[1:, :] * np.conj(g[:-1, :]))
        wx = np.sqrt(np.abs(g[:, 1:] * g[:, :-1])).reshape(-1, 1)
        wy = np.sqrt(np.abs(g[1:, :] * g[:-1, :])).reshape(-1, 1)
        A = np.vstack([A1 * wx, A2 * wy])
        b = np.concatenate([gx.reshape(-1, 1) * wx, gy.reshape(-1, 1) * wy])[:, 0]
        p = np.linalg.lstsq(A, b, rcond=None)[0]
        ph = p[0] * x ** 2 + p[1] * y ** 2 + p[2] * x * y + p[3] * x + p[4] * y
        total += ph
        g *= np.exp(-1j * ph)
    return g, total


# ------------------------------------------------------------------ metrics
def imag_energy(field):
    """Fraction of energy in the imaginary part (v0.1.2 metric).

    WARNING: depends on the field's global phase. For a real field times
    exp(i theta) it equals sin^2(theta), so on real data (arbitrary
    reference piston) its value at true focus is arbitrary, and in z it
    oscillates with period lambda/2. Use `realness`. Kept so the v0.1.2
    results can be reproduced and compared.
    """
    return float(np.mean(np.imag(field) ** 2) / np.mean(np.abs(field) ** 2))


def realness(field):
    """Global-phase-invariant imaginary-energy fraction, in [0, 0.5].

    realness(U) = min over theta of mean(Im(e^{i theta} U)^2) / mean|U|^2
                = 0.5 * (1 - |sum U^2| / sum |U|^2).
    It is 0 exactly when U is a real field times a constant phase factor,
    and 0.5 when the phases are circularly spread. A real amplitude object
    is real at focus whatever the reference piston, so this keeps the
    v0.1.2 design idea while dropping its accidental phase dependence.
    """
    f = np.asarray(field, dtype=complex)
    e = float(np.sum(np.abs(f) ** 2))
    if e <= 0 or not np.isfinite(e):
        return 0.5
    return float(0.5 * (1.0 - abs(np.sum(f * f)) / e))


def ac_realness(field):
    """Realness of the mean-subtracted (AC) field.

    For a thin object whose transmission is t = A + c * r(x), with r real
    and A, c complex constants, the AC part c * (r - mean r) is real up to
    a constant phase at focus. That covers amplitude masks on dark OR
    bright backgrounds, binary phase objects of ANY depth (t = 1 +
    (e^{i phi} - 1) m(x)), weak phase objects (t ~ 1 + i phi(x)), and
    single-material mixed objects. Propagation keeps the mean (DC) as a
    plane wave, so subtracting it commutes with propagation.
    """
    f = np.asarray(field, dtype=complex)
    return realness(f - f.mean())


def best_piston(field):
    """Constant phase theta that makes exp(i theta) * field most real."""
    f = np.asarray(field, dtype=complex)
    return float(-0.5 * np.angle(np.sum(f * f)))


def amplitude_variance(field):
    """Variance of the reconstructed intensity |field|^2.

    PHASE-OBJECT diagnostic: defocus converts phase into amplitude
    contrast, so for a phase object this is MINIMAL near focus. The dip
    is shallow (about 10% on the real fly), so never use it alone.
    Do not use it for amplitude objects.
    """
    inten = np.abs(field) ** 2
    return float(np.var(inten))


def tamura(field):
    """Tamura coefficient sqrt(std|U| / mean|U|) of the amplitude.

    A standard DHM sparsity criterion (Memmolo et al., 2011): maximal at
    focus for amplitude objects, minimal for phase objects. Included as a
    comparison metric for the phase-object study.
    """
    a = np.abs(np.asarray(field))
    m = float(a.mean())
    return float(np.sqrt(a.std() / m)) if m > 0 else 0.0


def field_sharpness(field):
    """Mean squared gradient of the complex field.

    FALSIFIED as a focus metric: angular-spectrum propagation multiplies
    each spatial frequency by a unit-modulus phase, so by Parseval this
    quantity is invariant under propagation. More generally, any
    Hermitian form <U, K U> with a shift-invariant K is invariant.
    (v0.1.2 said "any quadratic functional", which is too strong:
    sum Im(U)^2 and sum U^2 are quadratic but not Hermitian, and
    they do change with z.) Kept as the documented dead end.
    """
    gy, gx = np.gradient(field)
    return float(np.mean(np.abs(gx) ** 2 + np.abs(gy) ** 2))


def phase_variance(field):
    """Intensity-weighted circular variance of the phase.

    Measured to FAIL for phase objects (v0.1.2 DISCOVERY finding 2): the
    uniform background phasor dominates the mean at every plane, so the
    curve is flat. Works for amplitude objects, degrades earlier than
    realness under every noise type tested. Kept for comparison.
    """
    field = np.asarray(field, dtype=complex)
    w = np.abs(field) ** 2
    total = float(w.sum())
    if total <= 0 or not np.isfinite(total):
        return 1.0
    mean_phasor = np.sum(w * np.exp(1j * np.angle(field))) / total
    return float(1.0 - abs(mean_phasor))


METRICS = {
    "realness": realness,
    "ac_realness": ac_realness,
    "imag_energy": imag_energy,
    "phase_variance": phase_variance,
    "amplitude_variance": amplitude_variance,
    "neg_tamura": lambda f: -tamura(f),   # amplitude objects: max -> min
    "tamura": tamura,                     # phase objects: min at focus
}


# ------------------------------------------------------------- focus search
def parabolic_min(sweep, metrics):
    """Parabolic interpolation of the metric minimum around the argmin.

    Valid for SMOOTH metrics (realness, ac_realness). On imag_energy it
    fits the lambda/2 piston oscillation, not the focus basin. Returns
    (distance, metric); falls back to the grid minimum at an edge or for
    a degenerate fit; the offset is clamped to one step.
    """
    sweep = np.asarray(sweep, dtype=float)
    metrics = np.asarray(metrics, dtype=float)
    i = int(np.argmin(metrics))
    if i == 0 or i == len(metrics) - 1:
        return float(sweep[i]), float(metrics[i])
    y0, y1, y2 = metrics[i - 1], metrics[i], metrics[i + 1]
    denom = y0 - 2 * y1 + y2
    if denom <= 0:
        return float(sweep[i]), float(y1)
    step = float(sweep[i + 1] - sweep[i])
    offset = 0.5 * (y0 - y2) / denom
    if abs(offset) > 1:
        return float(sweep[i]), float(y1)
    return (float(sweep[i] + offset * step),
            float(y1 - 0.25 * (y0 - y2) * offset))


def golden_min(fun, a, b, tol=1e-8, max_iter=200):
    """Golden-section minimisation of a unimodal function on [a, b]."""
    g = (np.sqrt(5) - 1) / 2
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = fun(c), fun(d)
    for _ in range(max_iter):
        if abs(b - a) < tol:
            break
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = fun(c)
        else:
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = fun(d)
    x = 0.5 * (a + b)
    return float(x), float(fun(x))


def _metric(metric):
    if callable(metric):
        return metric
    try:
        return METRICS[metric]
    except KeyError as exc:
        raise ValueError(f"unknown metric {metric!r}; "
                         f"choose from {sorted(METRICS)}") from exc


def focus_sweep(field0, sweep, metric=realness, dx=DX, wl=WL,
                keep_fields=True):
    """Propagate across sweep; return (fields, metrics, best_index).

    metric: callable or a METRICS key. v0.1.3 default is `realness`
    (v0.1.2 used imag_energy; pass metric=imag_energy to reproduce it).
    keep_fields=False returns fields=None (saves memory on big arrays).
    """
    sweep = np.asarray(sweep, dtype=float)
    if sweep.size < 2 or not np.isfinite(sweep).all():
        raise ValueError("sweep needs at least two finite distances")
    fn = _metric(metric)
    prop = Propagator(field0, dx, wl)
    fields, vals = [], []
    for d in sweep:
        f = prop(d)
        vals.append(fn(f))
        if keep_fields:
            fields.append(f)
    vals = np.array(vals)
    return (fields if keep_fields else None), vals, int(np.argmin(vals))


def find_focus(field0, sweep, metric=realness, dx=DX, wl=WL, refine=True,
               tol=1e-8):
    """Coarse sweep + golden-section refinement of the focus distance.

    Returns a dict: z (refined, m), coarse_z, value, contrast
    (median / min of the coarse curve), curve (coarse metric values).
    The refinement brackets the coarse argmin by its two neighbours, so
    it can never jump to another basin.
    """
    sweep = np.asarray(sweep, dtype=float)
    fn = _metric(metric)
    prop = Propagator(field0, dx, wl)
    curve = np.array([fn(prop(d)) for d in sweep])
    i = int(np.argmin(curve))
    z, val = float(sweep[i]), float(curve[i])
    if refine and len(sweep) >= 3:
        a = sweep[max(i - 1, 0)]
        b = sweep[min(i + 1, len(sweep) - 1)]
        z, val = golden_min(lambda d: fn(prop(d)), min(a, b), max(a, b), tol)
        if val > curve[i]:
            z, val = float(sweep[i]), float(curve[i])
    med = float(np.median(curve))
    return {"z": z, "coarse_z": float(sweep[i]), "value": val,
            "contrast": med / val if val > 0 else np.inf, "curve": curve}


def local_focus_map(field0, sweep, tile=32, metric="realness", dx=DX, wl=WL,
                    min_energy_frac=0.05):
    """Per-tile focus distance: a depth map for thick / multi-plane samples.

    The field is cut into tile x tile blocks; for each z the per-block
    realness (or ac_realness) is computed from block sums, so one sweep
    gives every block's focus curve. Blocks whose energy at their best
    plane is below min_energy_frac of the brightest block are NaN.
    Returns (zmap [ty, tx], curves [len(sweep), ty, tx]).
    """
    if metric not in ("realness", "ac_realness"):
        raise ValueError("local_focus_map supports realness/ac_realness")
    sweep = np.asarray(sweep, dtype=float)
    prop = Propagator(field0, dx, wl)
    ny, nx = prop.shape
    ty, tx = ny // tile, nx // tile
    curves = np.empty((len(sweep), ty, tx))
    energy = np.empty((len(sweep), ty, tx))
    npx = tile * tile
    for k, d in enumerate(sweep):
        f = prop(d)[:ty * tile, :tx * tile].reshape(ty, tile, tx, tile)
        e = np.sum(np.abs(f) ** 2, axis=(1, 3))
        s2 = np.sum(f * f, axis=(1, 3))
        if metric == "ac_realness":
            s1 = np.sum(f, axis=(1, 3))
            s2 = s2 - s1 ** 2 / npx
            e = e - np.abs(s1) ** 2 / npx
        with np.errstate(invalid="ignore", divide="ignore"):
            curves[k] = 0.5 * (1 - np.abs(s2) / e)
        energy[k] = e
    best = np.argmin(np.nan_to_num(curves, nan=1.0), axis=0)
    zmap = sweep[best]
    e_best = np.take_along_axis(energy, best[None], axis=0)[0]
    zmap = np.where(e_best >= min_energy_frac * e_best.max(), zmap, np.nan)
    return zmap, curves


def correlation(a, b):
    """Pearson correlation between two images."""
    return float(np.corrcoef(np.asarray(a).ravel(), np.asarray(b).ravel())[0, 1])


# ------------------------------------------------------------------ outputs
def plot_focus_curve(sweep, metrics, best_d, truth_d, path,
                     ylabel="realness (phase-invariant)"):
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib.figure import Figure

    fig = Figure(figsize=(8, 4), dpi=100)
    ax = fig.subplots(1, 1)
    ax.plot(sweep * 1e3, metrics)
    ax.axvline(best_d * 1e3, color="r", ls="--",
               label=f"best {best_d * 1e3:+.3f} mm")
    ax.axvline(truth_d * 1e3, color="g", ls=":",
               label=f"truth {truth_d * 1e3:+.3f} mm")
    ax.set_xlabel("reconstruction distance [mm]")
    ax.set_ylabel(ylabel)
    ax.set_title("Focus curve: minimum = field real up to a constant phase")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path)
    return path


def render_video(fields, sweep, fps, out_path, frames_dir):
    """Render intensity + phase frames and encode them with ffmpeg."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.image
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    os.makedirs(frames_dir, exist_ok=True)
    peak = max(float(np.max(np.abs(f) ** 2)) for f in fields)
    thresh = 0.0025 * peak
    fig = Figure(figsize=(12, 6), dpi=100)
    canvas = FigureCanvasAgg(fig)
    axes = fig.subplots(1, 2)
    first = fields[0]
    inten = np.abs(first) ** 2
    im0 = axes[0].imshow(inten, origin="lower", cmap="inferno", vmin=0, vmax=peak)
    im1 = axes[1].imshow(
        np.ma.masked_where(inten <= thresh, np.angle(first)),
        origin="lower", cmap="twilight", vmin=-np.pi, vmax=np.pi,
    )
    axes[0].set_title("Recovered intensity")
    axes[1].set_title("Wrapped phase (piston removed)")
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
    fig.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04,
                 label="Relative intensity")
    fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04, label="Phase [rad]")
    title = fig.suptitle("")
    fig.subplots_adjust(left=0.04, right=0.94, bottom=0.10, top=0.88, wspace=0.35)

    for i, (d, rec) in enumerate(zip(sweep, fields)):
        rec = rec * np.exp(1j * best_piston(rec))
        inten = np.abs(rec) ** 2
        im0.set_data(inten)
        im1.set_data(np.ma.masked_where(inten <= thresh, np.angle(rec)))
        title.set_text(
            f"Numerical reconstruction - z = {d * 1e3:+.4f} mm "
            f"| frame {i + 1}/{len(sweep)}"
        )
        canvas.draw()
        rgb = np.asarray(canvas.buffer_rgba())[:, :, :3]
        matplotlib.image.imsave(os.path.join(frames_dir, f"f{i:04d}.png"), rgb)

    try:
        subprocess.run(
            ["ffmpeg", "-nostdin", "-y", "-v", "error", "-framerate", str(fps),
             "-i", os.path.join(frames_dir, "f%04d.png"),
             "-c:v", "libx264", "-pix_fmt", "yuv420p", out_path],
            check=True, stdin=subprocess.DEVNULL,
        )
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        raise RuntimeError(
            "ffmpeg failed or is missing; re-run with --no-video to keep "
            "the frames, still, and curve"
        ) from exc
    return out_path


# ---------------------------------------------------------------------- run
def run(z_obj=Z_OBJ_DEFAULT, sweep_min=-10e-3, sweep_max=0.0, frames=101,
        fps=FPS_DEFAULT, car=CAR_DEFAULT, win=WIN_DEFAULT,
        outdir=OUTDIR_DEFAULT, make_video=True, metric="realness",
        ref_phase=0.0, estimate=True):
    """Run the full demo; return a dict of measured results.

    estimate=True: the carrier is estimated from the hologram (sub-bin)
    rather than assumed; the true carrier `car` may be fractional.
    """
    if z_obj <= 0:
        raise ValueError("z_obj must be positive")
    if not (sweep_min < sweep_max) or frames < 2:
        raise ValueError("need sweep_min < sweep_max and frames >= 2")
    os.makedirs(outdir, exist_ok=True)
    sweep = np.linspace(sweep_min, sweep_max, frames)

    obj = make_object()
    sensor_field = angular_spectrum(obj.astype(complex), z_obj)
    holo = record_hologram(sensor_field, car, ref_phase=ref_phase)
    if estimate:
        carrier = estimate_carrier_2d(holo, dc_radius=30.0 / N)
    else:
        carrier = (0.0, -float(car))
    field0 = reconstruct_field(holo, win=win, carrier=carrier)
    res = find_focus(field0, sweep, metric=metric)
    fields, metrics, best_i = focus_sweep(field0, sweep, metric=metric,
                                          keep_fields=make_video)
    best_field = angular_spectrum(field0, res["z"])
    corr = correlation(np.abs(best_field) ** 2, obj)

    results = {
        "truth_mm": -z_obj * 1e3,
        "carrier_est": carrier,
        "best_mm": res["z"] * 1e3,
        "coarse_mm": res["coarse_z"] * 1e3,
        "error_um": abs(res["z"] + z_obj) * 1e6,
        "metric": metric,
        "metric_best": res["value"],
        "metric_median": float(np.median(res["curve"])),
        "contrast": res["contrast"],
        "corr": corr,
    }

    plot_focus_curve(sweep, metrics, res["z"], -z_obj,
                     os.path.join(outdir, "focus_curve.png"), ylabel=metric)
    if make_video:
        render_video(fields, sweep, fps,
                     os.path.join(outdir, "focus_sweep.mp4"),
                     os.path.join(outdir, "frames"))
    import matplotlib.image

    inten = np.abs(best_field) ** 2
    matplotlib.image.imsave(os.path.join(outdir, "best_focus.png"),
                            inten / inten.max())
    return results


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--z-obj", type=float, default=Z_OBJ_DEFAULT * 1e3,
                        help="object -> sensor distance in mm (default 5.0)")
    parser.add_argument("--sweep-min", type=float, default=-10.0,
                        help="sweep start in mm (default -10.0)")
    parser.add_argument("--sweep-max", type=float, default=0.0,
                        help="sweep end in mm (default 0.0)")
    parser.add_argument("--frames", type=int, default=101)
    parser.add_argument("--fps", type=int, default=FPS_DEFAULT)
    parser.add_argument("--carrier", type=float, default=CAR_DEFAULT,
                        help="TRUE carrier in bins (may be fractional)")
    parser.add_argument("--window", type=int, default=WIN_DEFAULT)
    parser.add_argument("--ref-phase", type=float, default=0.0,
                        help="reference piston in rad (real data: arbitrary)")
    parser.add_argument("--metric", default="realness", choices=sorted(METRICS))
    parser.add_argument("--assume-carrier", action="store_true",
                        help="do not estimate the carrier (v0.1.2 behaviour)")
    parser.add_argument("--outdir", default=OUTDIR_DEFAULT)
    parser.add_argument("--no-video", action="store_true")
    args = parser.parse_args(argv)

    r = run(
        z_obj=args.z_obj * 1e-3,
        sweep_min=args.sweep_min * 1e-3,
        sweep_max=args.sweep_max * 1e-3,
        frames=args.frames, fps=args.fps,
        car=args.carrier, win=args.window,
        outdir=args.outdir, make_video=not args.no_video,
        metric=args.metric, ref_phase=args.ref_phase,
        estimate=not args.assume_carrier,
    )
    ky, kx = r["carrier_est"]
    print(f"wave-hologram {__version__}  metric: {r['metric']}")
    print(f"carrier (ky, kx)     : ({ky:+.4f}, {kx:+.4f}) bins")
    print(f"true object distance : {r['truth_mm']:+.4f} mm (back-propagation)")
    print(f"coarse best focus    : {r['coarse_mm']:+.4f} mm")
    print(f"refined best focus   : {r['best_mm']:+.4f} mm")
    print(f"focus error          : {r['error_um']:.2f} um")
    print(f"metric at best       : {r['metric_best']:.5f} "
          f"(median {r['metric_median']:.5f}, contrast {r['contrast']:.1f}x)")
    print(f"corr with object     : {r['corr']:.3f}")
    print(f"outputs in {args.outdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
