"""Sideband separation rule for off-axis recording.

For an object field whose spectrum is confined to |kx|, |ky| <= B bins on an
n x n grid and an x-carrier of `car` bins:

* the zero-order terms (|R|^2 + |U|^2) occupy |k| <= 2B, because |U|^2 is the
  autocorrelation of the spectrum of U;
* the wanted sideband U R* occupies car - B .. car + B around the carrier;
* so a window of half-width B recovers U exactly when car - B > 2B, that is
  car >= 3B + 1 (bins are integers: at car = 3B the first sideband bin lies
  on the last zero-order bin);
* and only while the sideband stays below Nyquist, car + B < n / 2, that is
  car <= n/2 - B - 1.

Fixed seeds, exact band limits, integer carriers. The 1e-10 bound for "exact"
is a test bound; the measured error there is about 1e-15 (see
docs/validation.md). Below the threshold the cross-talk is deterministic for a
given seed, and the asserted lower bounds are well under the measured values.
"""

import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import wave_hologram as W  # noqa: E402

N = 128
PISTON = 0.9


def bandlimited_field(n, b, seed):
    """Complex field with spectrum exactly on |ky|, |kx| <= b (peak amplitude 0.5)."""
    rng = np.random.default_rng(seed)
    spec = np.zeros((n, n), complex)
    k = np.arange(-b, b + 1)
    coeff = rng.normal(size=(k.size, k.size)) + 1j * rng.normal(size=(k.size, k.size))
    spec[np.ix_(k % n, k % n)] = coeff
    u = np.fft.ifft2(spec)
    return 0.5 * u / np.abs(u).max()


def bandlimited_amplitude(n, b, seed):
    """Real, positive amplitude object, band-limited to b bins (no DC term drawn)."""
    rng = np.random.default_rng(seed)
    spec = np.zeros((n, n), complex)
    k = np.arange(-b, b + 1)
    coeff = rng.normal(size=(k.size, k.size)) + 1j * rng.normal(size=(k.size, k.size))
    coeff[b, b] = 0.0                                   # DC set separately below
    spec[np.ix_(k % n, k % n)] = coeff
    a = np.fft.ifft2(spec).real
    return 0.5 + 0.4 * a / np.abs(a).max()


def rel_error(u, car, win, piston=PISTON):
    holo = W.record_hologram(u, car=car, ref_phase=piston)
    rec = W.reconstruct_field(holo, car=car, win=win)
    want = u * np.exp(-1j * piston)
    return np.linalg.norm(rec - want) / np.linalg.norm(want)


class TestSeparationThreshold(unittest.TestCase):
    def test_threshold_is_three_b_plus_one_for_several_bandwidths(self):
        for b in (4, 6, 8, 10):
            u = bandlimited_field(N, b, seed=11 + b)
            lo, hi = 3 * b + 1, N // 2 - b - 1
            for car in (lo, (lo + hi) // 2, hi):
                with self.subTest(b=b, car=car, side="separated"):
                    self.assertLess(rel_error(u, car, b), 1e-10)
            with self.subTest(b=b, car=3 * b, side="one bin short"):
                err = rel_error(u, 3 * b, b)
                self.assertGreater(err, 1e-3)       # a single shared bin already leaks
                self.assertLess(err, 0.1)           # but only a little
            with self.subTest(b=b, car=N // 2 - b, side="past Nyquist"):
                self.assertGreater(rel_error(u, N // 2 - b, b), 0.05)

    def test_error_grows_as_the_carrier_falls_below_the_threshold(self):
        b = 8
        cars = (12, 16, 20, 23, 24, 25)               # 25 = 3b + 1
        for seed in range(5):
            u = bandlimited_field(N, b, seed)
            errs = [rel_error(u, c, b) for c in cars]
            with self.subTest(seed=seed):
                self.assertTrue(all(x > y for x, y in zip(errs[:-1], errs[1:])),
                                msg=f"not decreasing: {errs}")
                self.assertGreater(errs[0], 0.03)     # cross-talk of several percent
                self.assertLess(errs[-1], 1e-10)

    def test_carrier_at_or_below_b_overlaps_the_twin_image(self):
        b = 8
        u = bandlimited_field(N, b, 2026)
        for car in (2, 4, 8):
            with self.subTest(car=car):
                self.assertGreater(rel_error(u, car, b), 1.0)   # not a reconstruction

    def test_the_rule_applies_to_the_window_not_only_the_object(self):
        # At car = 25 (3b + 1 for b = 8) a window wider than b reaches the
        # zero-order terms; a narrower one cuts the object itself.
        u = bandlimited_field(N, 8, 2026)
        self.assertLess(rel_error(u, 25, 8), 1e-10)
        self.assertGreater(rel_error(u, 25, 9), 1e-3)
        self.assertGreater(rel_error(u, 25, 10), 1e-3)
        self.assertGreater(rel_error(u, 25, 7), 0.1)
        self.assertGreater(rel_error(u, 25, 6), 0.1)

    def test_the_piston_does_not_change_the_threshold(self):
        u = bandlimited_field(N, 8, 3)
        for piston in (0.0, 2.1, -4.0):
            self.assertLess(rel_error(u, 25, 8, piston), 1e-10)
            self.assertGreater(rel_error(u, 24, 8, piston), 1e-3)


class TestFocusBiasFollowsTheSeparationRule(unittest.TestCase):
    """The focus estimate is exact (to 1e-9 um) above the threshold and biased
    by cross-talk below it. Real, positive, band-limited amplitude objects,
    propagated 5 mm with the package defaults for pitch and wavelength."""

    B = 8
    SWEEP = np.linspace(-10e-3, 0.0, 101)

    def focus_error_um(self, seed, car):
        amp = bandlimited_amplitude(N, self.B, seed)
        sensor = W.angular_spectrum(amp.astype(complex), 5e-3)
        holo = W.record_hologram(sensor, car=car, ref_phase=1.1)
        f0 = W.reconstruct_field(holo, car=car, win=self.B)
        return abs(W.find_focus(f0, self.SWEEP)["z"] + 5e-3) * 1e6

    def test_above_the_threshold_the_focus_is_exact_in_the_noiseless_case(self):
        for seed in range(4):
            for car in (25, 40, 55):
                with self.subTest(seed=seed, car=car):
                    self.assertLess(self.focus_error_um(seed, car), 1e-3)

    def test_below_the_threshold_cross_talk_biases_the_focus(self):
        # Measured 1.4 to 187 um for these carriers and seeds; asserted > 0.5 um.
        # The size is not monotonic in the carrier (the leak is a random
        # interference term), so only the lower bound is asserted.
        for seed in range(4):
            for car in (12, 16, 20, 24):
                with self.subTest(seed=seed, car=car):
                    self.assertGreater(self.focus_error_um(seed, car), 0.5)


class TestDefaultGeometry(unittest.TestCase):
    def test_default_carrier_and_window_do_not_satisfy_the_rule_for_a_filled_window(self):
        # CAR_DEFAULT = 80, WIN_DEFAULT = 48, N = 256. For an object that fills
        # the window the rule needs 3 * 48 + 1 = 145 <= car <= 256 / 2 - 48 - 1 = 79,
        # which no carrier satisfies: the window is too wide for this grid.
        # The demo object is not band-limited to 48 bins, so this is a statement
        # about the geometry, not a measured error; the recorded residual bias is
        # in the README and RESULTS_v0.1.3.md.
        lo = 3 * W.WIN_DEFAULT + 1
        hi = W.N // 2 - W.WIN_DEFAULT - 1
        self.assertGreater(lo, hi)
        self.assertLess(W.CAR_DEFAULT, lo)
        self.assertGreater(W.CAR_DEFAULT, hi)


if __name__ == "__main__":
    unittest.main()
