"""Analytic and physical validation of the propagation and recording code.

Each test compares the code with a closed-form result, an exact symmetry, or a
conservation law. The tolerances are test bounds, not accuracy claims; the
measured agreement for each case is listed in docs/validation.md.

Conventions under test: W.angular_spectrum(field, z) multiplies the spectrum
by exp(+i kz z), kz = (2 pi / wl) sqrt(1 - (wl fx)^2 - (wl fy)^2), so a
positive z moves a wave with positive fx toward +x.
"""

import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import wave_hologram as W  # noqa: E402

DX, WL = W.DX, W.WL


def kz_of(fx):
    return 2 * np.pi / WL * np.sqrt(1 - (WL * fx) ** 2)


class TestPlaneWaveAndSymmetries(unittest.TestCase):
    def test_on_grid_plane_wave_gains_exactly_exp_i_kz_z(self):
        n, m, z = 64, 5, 3e-3
        x = np.arange(n)[None, :] * np.ones((n, 1))
        u = np.exp(2j * np.pi * m * x / n)
        fx = m / (n * DX)
        expect = u * np.exp(1j * kz_of(fx) * z)
        self.assertLess(np.abs(W.angular_spectrum(u, z) - expect).max(), 1e-12)

    def test_uniform_field_gains_the_free_space_phase_exp_i_k_z(self):
        z = 5e-3
        out = W.angular_spectrum(np.ones((16, 16), complex), z)
        self.assertLess(np.abs(out - np.exp(2j * np.pi / WL * z)).max(), 1e-12)

    def test_propagation_composes_additively_in_distance(self):
        rng = np.random.default_rng(0)
        f = rng.normal(size=(64, 64)) + 1j * rng.normal(size=(64, 64))
        two_steps = W.angular_spectrum(W.angular_spectrum(f, 1.3e-3), 2.1e-3)
        one_step = W.angular_spectrum(f, 3.4e-3)
        self.assertLess(np.abs(two_steps - one_step).max(), 1e-9)

    def test_energy_is_conserved_when_every_wave_propagates(self):
        # On this grid the highest frequency is 1/(2 dx) << 1/wl, so no
        # component is evanescent and sum |U|^2 must not change.
        rng = np.random.default_rng(1)
        f = rng.normal(size=(64, 64)) + 1j * rng.normal(size=(64, 64))
        g = W.angular_spectrum(f, 7e-3)
        ratio = np.sum(np.abs(g) ** 2) / np.sum(np.abs(f) ** 2)
        self.assertAlmostEqual(ratio, 1.0, places=12)

    def test_realness_of_a_real_object_is_even_in_z(self):
        # For real U0, F(-k) = conj F(k) and kz is even, so U_{-z} = conj U_z,
        # and realness is invariant under conjugation.
        obj = W.make_object().astype(complex)
        for z in (1e-3, 3e-3, 7e-3):
            a = W.realness(W.angular_spectrum(obj, z))
            b = W.realness(W.angular_spectrum(obj, -z))
            self.assertAlmostEqual(a, b, places=12)
            self.assertGreater(a, 0.01)  # not trivially zero away from focus


class TestGaussianBeam(unittest.TestCase):
    """Free-space Gaussian beam: w(z) = w0 sqrt(1 + (z / zR)^2), zR = pi w0^2 / wl.

    The Gaussian is a paraxial solution; with w0 = 20 pixels (130 um) at
    532 nm the paraxial error is far below the tolerance used here.
    """

    N = 512

    def setUp(self):
        n = self.N
        y, x = np.mgrid[0:n, 0:n]
        self.xs = (x - n / 2) * DX
        r2 = self.xs ** 2 + ((y - n / 2) * DX) ** 2
        self.w0 = 20 * DX
        self.zR = np.pi * self.w0 ** 2 / WL
        self.field = np.exp(-r2 / self.w0 ** 2).astype(complex)

    def test_width_and_peak_intensity_follow_the_analytic_beam(self):
        for frac in (0.5, 1.0, 2.0):
            z = frac * self.zR
            inten = np.abs(W.angular_spectrum(self.field, z)) ** 2
            c = (inten * self.xs).sum() / inten.sum()
            var = (inten * (self.xs - c) ** 2).sum() / inten.sum()
            w = 2 * np.sqrt(var)          # I ~ exp(-2 r^2 / w^2): var_x = w^2 / 4
            w_true = self.w0 * np.sqrt(1 + frac ** 2)
            self.assertLess(abs(w / w_true - 1), 1e-4, msg=f"z = {frac} zR")
            peak = inten.max() / (self.w0 / w_true) ** 2   # on-axis I = (w0/w)^2
            self.assertLess(abs(peak - 1), 1e-4, msg=f"z = {frac} zR")


class TestTalbotSelfImaging(unittest.TestCase):
    """A cosine amplitude grating of period p reproduces itself at the Talbot
    distance. With the exact (non-paraxial) dispersion the first-order phase
    difference is 2 pi at zT = wl / (1 - sqrt(1 - (wl / p)^2)); at zT / 2 the
    pattern is shifted by half a period.
    """

    def setUp(self):
        self.n, self.p, self.m = 256, 16, 0.3      # period of 16 px divides 256
        x = np.arange(self.n)[None, :] * np.ones((self.n, 1))
        self.x = x
        self.grating = (1 + self.m * np.cos(2 * np.pi * x / self.p)).astype(complex)
        self.zT = WL / (1 - np.sqrt(1 - (WL / (self.p * DX)) ** 2))

    def test_talbot_distance_is_close_to_the_paraxial_value(self):
        paraxial = 2 * (self.p * DX) ** 2 / WL
        self.assertLess(abs(self.zT / paraxial - 1), 1e-4)

    def test_pattern_repeats_at_the_talbot_distance(self):
        out = W.angular_spectrum(self.grating, self.zT)
        i0 = np.abs(self.grating) ** 2
        self.assertLess(np.abs(np.abs(out) ** 2 - i0).max(), 1e-9)

    def test_pattern_is_half_a_period_shifted_at_half_the_talbot_distance(self):
        out = W.angular_spectrum(self.grating, 0.5 * self.zT)
        shifted = (1 - self.m * np.cos(2 * np.pi * self.x / self.p)) ** 2
        self.assertLess(np.abs(np.abs(out) ** 2 - shifted).max(), 1e-9)


class TestCircularApertureFresnel(unittest.TestCase):
    """On-axis intensity behind a circular aperture of radius a (Fresnel
    diffraction): I / I0 = 4 sin^2(pi N_F / 2), N_F = a^2 / (wl z). It is
    a paraxial result and the aperture is a pixelated disc, so the comparison
    is made to a coarse absolute tolerance, using the radius implied by the
    pixel count.
    """

    def test_axial_intensity_follows_the_fresnel_oscillation(self):
        n = 512
        y, x = np.mgrid[0:n, 0:n]
        ap = (np.hypot(x - n / 2, y - n / 2) < 20).astype(float)
        a_eff = np.sqrt(ap.sum() / np.pi) * DX
        for nf in (0.5, 1.0, 1.5, 2.5, 3.0):
            z = a_eff ** 2 / (WL * nf)
            c = W.angular_spectrum(ap.astype(complex), z)[n // 2, n // 2]
            expect = 4 * np.sin(np.pi * nf / 2) ** 2
            self.assertAlmostEqual(abs(c) ** 2, expect, delta=0.05,
                                   msg=f"N_F = {nf}")

    def test_axial_intensity_vanishes_at_an_even_fresnel_number(self):
        n = 512
        y, x = np.mgrid[0:n, 0:n]
        ap = (np.hypot(x - n / 2, y - n / 2) < 20).astype(float)
        a_eff = np.sqrt(ap.sum() / np.pi) * DX
        z = a_eff ** 2 / (WL * 2.0)
        c = W.angular_spectrum(ap.astype(complex), z)[n // 2, n // 2]
        self.assertLess(abs(c) ** 2, 1e-3)


class TestSamplingLimits(unittest.TestCase):
    """The grid is periodic. A wave packet with spatial frequency fx moves
    laterally by z tan(theta), sin(theta) = wl fx; once that exceeds the
    field of view it wraps around instead of leaving. The same condition
    bounds the usable distance: z < N dx^2 / wl for fx up to 1/(2 dx).
    """

    N = 256

    def packet(self, m):
        y, x = np.mgrid[0:self.N, 0:self.N]
        env = np.exp(-((x - self.N / 2) ** 2 + (y - self.N / 2) ** 2) / (2 * 12 ** 2))
        return env * np.exp(2j * np.pi * m * x / self.N), x

    def centroid_x(self, inten, x):
        theta = 2 * np.pi * x / self.N              # circular mean on a periodic grid
        return (np.angle((inten * np.exp(1j * theta)).sum()) / (2 * np.pi) * self.N) % self.N

    def predicted_shift_px(self, m, z):
        fx = m / (self.N * DX)
        tan = WL * fx / np.sqrt(1 - (WL * fx) ** 2)
        return z * tan / DX

    def test_lateral_shift_matches_z_tan_theta(self):
        for m, z in ((20, 3e-3), (60, 8e-3), (60, 20e-3)):
            g, x = self.packet(m)
            c = self.centroid_x(np.abs(W.angular_spectrum(g, z)) ** 2, x)
            expect = (self.N / 2 + self.predicted_shift_px(m, z)) % self.N
            self.assertAlmostEqual(c, expect, delta=0.05, msg=f"m={m}, z={z}")

    def test_shift_beyond_the_field_of_view_wraps_around(self):
        m, z = 60, 60e-3
        shift = self.predicted_shift_px(m, z)
        self.assertGreater(shift, self.N / 2)       # leaves the half-field...
        g, x = self.packet(m)
        c = self.centroid_x(np.abs(W.angular_spectrum(g, z)) ** 2, x)
        self.assertAlmostEqual(c, (self.N / 2 + shift) % self.N, delta=0.05)
        # ...and comes back in from the other side: it is NOT at N/2 + shift
        self.assertLess(c, self.N / 2)

    def test_back_propagation_of_a_wrapped_packet_still_recovers_it(self):
        # Wrapping does not destroy information (propagation is unitary on the
        # grid); it only makes the apparent position ambiguous.
        g, _ = self.packet(60)
        back = W.angular_spectrum(W.angular_spectrum(g, 60e-3), -60e-3)
        self.assertLess(np.abs(back - g).max(), 1e-10)


class TestRecordingGeometry(unittest.TestCase):
    """Off-axis recording of a single tone: H = |U + R|^2 has a DC level of
    1 + a^2 and two conjugate fringe components of amplitude a at +-(car - m)
    bins (Fourier coefficient a each), so the sideband is exactly recoverable.
    """

    n = 128

    def setUp(self):
        y, x = np.mgrid[0:self.n, 0:self.n]
        self.x = x
        self.a, self.m, self.car, self.phi = 0.4, 3, 40, 1.1
        self.u = self.a * np.exp(2j * np.pi * self.m * x / self.n)

    def test_fourier_coefficients_of_the_hologram(self):
        h = W.record_hologram(self.u, car=self.car, ref_phase=self.phi)
        f = np.fft.fft2(h) / self.n ** 2
        self.assertAlmostEqual(f[0, 0].real, 1 + self.a ** 2, places=12)
        self.assertAlmostEqual(abs(f[0, (self.m - self.car) % self.n]), self.a, places=12)
        self.assertAlmostEqual(abs(f[0, (self.car - self.m) % self.n]), self.a, places=12)
        f[0, 0] = 0
        f[0, (self.m - self.car) % self.n] = 0
        f[0, (self.car - self.m) % self.n] = 0
        self.assertLess(np.abs(f).max(), 1e-12)     # nothing else is present

    def test_reconstruction_returns_the_object_times_the_reference_piston(self):
        h = W.record_hologram(self.u, car=self.car, ref_phase=self.phi)
        rec = W.reconstruct_field(h, car=self.car, win=8)
        self.assertLess(np.abs(rec - self.u * np.exp(-1j * self.phi)).max(), 1e-10)

    def test_fringes_above_nyquist_alias_to_a_lower_frequency(self):
        # Fringe frequency car - m = 70 > n/2 = 64 is indistinguishable, on this
        # grid, from n - 70 = 58 with the same (zero) piston: the recorded
        # hologram is identical, so the object frequency would be misread.
        a = 0.4
        above = W.record_hologram(a * np.exp(2j * np.pi * (-30) * self.x / self.n),
                                  car=40, ref_phase=0.0)
        below = W.record_hologram(a * np.exp(2j * np.pi * (-18) * self.x / self.n),
                                  car=40, ref_phase=0.0)
        self.assertLess(np.abs(above - below).max(), 1e-12)

    def test_carrier_at_or_above_nyquist_is_rejected(self):
        for car in (0, self.n // 2, self.n):
            with self.assertRaises(ValueError):
                W.record_hologram(self.u, car=car)


class TestDeterminismAndReadmeNumbers(unittest.TestCase):
    """Quickstart numbers quoted in the README, and run-to-run determinism."""

    @staticmethod
    def quickstart():
        obj = W.make_object()
        sensor = W.angular_spectrum(obj.astype(complex), 5.0e-3)
        holo = W.record_hologram(sensor, car=80.37, ref_phase=1.2)
        carrier = W.estimate_carrier_2d(holo, dc_radius=30.0 / W.N)
        field0 = W.reconstruct_field(holo, carrier=carrier)
        res = W.find_focus(field0, np.linspace(-10e-3, 0.0, 101))
        return obj, carrier, field0, res

    def test_readme_quickstart_focus_and_carrier(self):
        obj, carrier, field0, res = self.quickstart()
        self.assertEqual(round(res["z"] * 1e3, 3), -4.998)   # README: -4.998 mm
        self.assertLess(abs(carrier[1] + 80.37), 0.03)       # README: error <= 0.03 bins
        self.assertLess(abs(carrier[0]), 0.03)
        best = W.angular_spectrum(field0, res["z"])
        self.assertEqual(round(W.correlation(np.abs(best) ** 2, obj), 3), 0.956)

    def test_pipeline_is_bitwise_repeatable_in_one_process(self):
        _, c1, f1, r1 = self.quickstart()
        _, c2, f2, r2 = self.quickstart()
        self.assertEqual(c1, c2)
        self.assertTrue(np.array_equal(f1, f2))
        self.assertEqual(r1["z"], r2["z"])


if __name__ == "__main__":
    unittest.main()
