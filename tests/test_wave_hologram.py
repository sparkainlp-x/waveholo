"""Tests for wave-hologram v0.1.3. Run from the package directory:
    python3 -m unittest discover -s tests -v
Every test checks a claim the README or RESULTS make. Thresholds come
from measured behaviour (studies.py output, 2026-10-25) with margin.
Tests that encode a v0.1.2 claim later shown wrong are kept, rewritten
to assert the corrected statement, and say so in their docstring.
"""
import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import wave_hologram as W  # noqa: E402
from wave_hologram import (  # noqa: E402
    angular_spectrum, make_object, record_hologram, reconstruct_field,
    offaxis_field, imag_energy, realness, ac_realness, correlation,
    focus_sweep, find_focus, estimate_carrier, estimate_carrier_2d,
    local_focus_map, Propagator, Z_OBJ,
)

N = 256
SWEEP = np.linspace(-10e-3, 0.0, 101)
STEP_MM = 0.1


def pipeline(obj, z_obj=Z_OBJ, ref_phase=0.0, car=W.CAR_DEFAULT,
             estimate=True, noise=None, seed=7):
    sensor = angular_spectrum(np.asarray(obj, complex), z_obj)
    h = record_hologram(sensor, car, ref_phase=ref_phase)
    if noise is not None:
        from noise_robustness import add_noise
        h = add_noise(h, noise[0], noise[1], np.random.default_rng(seed))
    carrier = (estimate_carrier_2d(h, dc_radius=30.0 / N) if estimate
               else (0.0, -float(car)))
    return reconstruct_field(h, carrier=carrier)


def two_plane_field(ref_phase=1.3):
    amp = make_object()
    y, x = np.mgrid[0:N, 0:N]
    dots = np.zeros((N, N))
    for dx0, dy0 in ((-70, -70), (70, -70), (0, 70)):
        dots[np.hypot(x - (N // 2 + dx0), y - (N // 2 + dy0)) < 10] = 1.0
    plane1 = np.clip(amp - dots, 0, 1)
    sensor = (angular_spectrum(plane1.astype(complex), 5e-3)
              + angular_spectrum(dots.astype(complex), 7.5e-3))
    return pipeline_sensor(sensor, ref_phase), plane1, dots


def pipeline_sensor(sensor, ref_phase):
    h = record_hologram(sensor, ref_phase=ref_phase)
    return reconstruct_field(h, carrier=estimate_carrier_2d(h, dc_radius=30.0 / N))


class TestPropagation(unittest.TestCase):
    def test_round_trip_recovers_object(self):
        obj = make_object()
        back = angular_spectrum(angular_spectrum(obj.astype(complex), Z_OBJ),
                                -Z_OBJ)
        self.assertGreater(correlation(np.abs(back), obj), 0.999)

    def test_rectangular_fields_allowed(self):
        """v0.1.2 rejected non-square fields; real sensors are not square."""
        f = np.zeros((40, 64), complex)
        f[20, 32] = 1
        out = angular_spectrum(f, 1e-3, dx=(5e-6, 6e-6))
        self.assertEqual(out.shape, (40, 64))
        self.assertAlmostEqual(float(np.sum(np.abs(out) ** 2)), 1.0, places=6)

    def test_rejects_1d(self):
        with self.assertRaises(ValueError):
            angular_spectrum(np.zeros(10, complex), 1e-3)

    def test_rejects_non_finite(self):
        bad = np.zeros((N, N), dtype=complex)
        bad[0, 0] = np.inf
        with self.assertRaises(ValueError):
            angular_spectrum(bad, 1e-3)

    def test_rejects_bad_carrier(self):
        with self.assertRaises(ValueError):
            record_hologram(np.zeros((N, N), dtype=complex), car=0)
        with self.assertRaises(ValueError):
            record_hologram(np.zeros((N, N), dtype=complex), car=N)

    def test_propagator_matches_angular_spectrum(self):
        f0 = pipeline(make_object())
        p = Propagator(f0)
        for z in (-5e-3, -1.234e-3, 2e-3):
            np.testing.assert_allclose(p(z), angular_spectrum(f0, z),
                                       rtol=1e-10, atol=1e-12)

    def test_kz_cache_is_bounded(self):
        """v0.1.2 cached one transfer per distance (OOM on 2048^2)."""
        f = np.ones((32, 32), complex)
        for z in np.linspace(0, 1e-3, 50):
            angular_spectrum(f, z)
        self.assertLessEqual(len(W._KZ_CACHE), W._KZ_CACHE_MAX)


class TestPiston(unittest.TestCase):
    """The reference piston is arbitrary on real data (S1)."""

    def test_imag_energy_depends_on_piston(self):
        obj = make_object()
        vals = [imag_energy(angular_spectrum(pipeline(obj, ref_phase=p), -Z_OBJ))
                for p in (0.0, np.pi / 2)]
        self.assertLess(vals[0], 0.05)
        self.assertGreater(vals[1], 0.9)

    def test_realness_is_piston_invariant(self):
        obj = make_object()
        f = angular_spectrum(pipeline(obj), -Z_OBJ)
        for p in np.linspace(0, 2 * np.pi, 7):
            self.assertAlmostEqual(realness(f * np.exp(1j * p)), realness(f),
                                   places=12)

    def test_realness_definition(self):
        """realness = min over theta of imag-energy fraction of e^{i theta} U."""
        rng = np.random.default_rng(0)
        f = rng.normal(size=(16, 16)) + 0.3j * rng.normal(size=(16, 16))
        thetas = np.linspace(0, np.pi, 20001)
        brute = min(imag_energy(f * np.exp(1j * t)) for t in thetas)
        self.assertAlmostEqual(realness(f), brute, places=6)
        g = f * np.exp(1j * W.best_piston(f))
        self.assertAlmostEqual(imag_energy(g), realness(f), places=10)

    def test_focus_found_for_any_piston(self):
        obj = make_object()
        for p in (0.0, 0.8, np.pi / 2, 2.5):
            r = find_focus(pipeline(obj, ref_phase=p), SWEEP)
            self.assertAlmostEqual(r["z"] * 1e3, -5.0, delta=0.05)

    def test_imag_energy_fails_with_piston(self):
        """With a quarter-wave piston the v0.1.2 metric misses focus."""
        r = find_focus(pipeline(make_object(), ref_phase=np.pi / 2), SWEEP,
                       metric="imag_energy", refine=False)
        self.assertGreater(abs(r["z"] * 1e3 + 5.0), 0.1)


class TestFocus(unittest.TestCase):
    def test_focus_when_clean(self):
        """Refined realness focus is within 30 um of truth; the residual
        bias (+27 um at win 48) is sideband cross-talk (S4)."""
        obj = make_object()
        f0 = pipeline(obj)
        r = find_focus(f0, SWEEP)
        self.assertAlmostEqual(r["z"] * 1e3, -5.00, delta=0.03)
        best = angular_spectrum(f0, r["z"])
        self.assertGreater(correlation(np.abs(best) ** 2, obj), 0.9)

    def test_no_bias_without_crosstalk(self):
        """On the true sensor field (no twin, no |U|^2 term) the realness
        minimum is at the truth to < 0.1 um."""
        sensor = angular_spectrum(make_object().astype(complex), Z_OBJ)
        r = find_focus(sensor, SWEEP, tol=1e-10)
        self.assertLess(abs(r["z"] + Z_OBJ), 1e-7)

    def test_metric_bottoms_out_at_truth(self):
        f0 = pipeline(make_object(), ref_phase=1.0)
        at_truth = realness(angular_spectrum(f0, -Z_OBJ))
        off = realness(angular_spectrum(f0, -Z_OBJ + 1e-3))
        self.assertLess(3 * at_truth, off)

    def test_reconstruct_handles_edge_window(self):
        obj = make_object()
        hologram = record_hologram(angular_spectrum(obj.astype(complex), Z_OBJ))
        field = reconstruct_field(hologram, car=82)
        self.assertEqual(field.shape, (N, N))
        self.assertTrue(np.isfinite(field).all())

    def test_integer_carrier_reproduces_v012(self):
        """Integer carriers take the exact-shift path (v0.1.2 behaviour)."""
        h = record_hologram(angular_spectrum(make_object().astype(complex), Z_OBJ))
        a = reconstruct_field(h)
        b = reconstruct_field(h, carrier=(0.0, -80.0))
        y, x = np.mgrid[0:N, 0:N]
        demod = h * np.exp(2j * np.pi * 80 * x / N)
        spec = np.fft.fftshift(np.fft.fft2(demod))
        c = N // 2
        m = np.zeros_like(spec)
        m[c - 48:c + 49, c - 48:c + 49] = 1
        ref = np.fft.ifft2(np.fft.ifftshift(spec * m))
        np.testing.assert_allclose(a, b)
        np.testing.assert_allclose(a, ref, atol=1e-9)

    def test_crop_is_lossless_for_focus(self):
        h = record_hologram(angular_spectrum(make_object().astype(complex), Z_OBJ))
        full = find_focus(reconstruct_field(h), SWEEP)
        fc = reconstruct_field(h, crop=True)
        crop = find_focus(fc, SWEEP, dx=W.cropped_pitch(h.shape, 48, W.DX))
        self.assertEqual(fc.shape, (98, 98))
        self.assertAlmostEqual(full["z"], crop["z"], delta=1e-8)

    def test_inline_twin_ambiguity(self):
        obj = make_object()
        sensor = angular_spectrum(obj.astype(complex), Z_OBJ)
        hologram = np.abs(sensor + 1.0) ** 2
        plus = np.abs(angular_spectrum(hologram.astype(complex), 5e-3)) ** 2
        minus = np.abs(angular_spectrum(hologram.astype(complex), -5e-3)) ** 2
        np.testing.assert_allclose(plus, minus, rtol=1e-6, atol=1e-9)

    def test_golden_min(self):
        x, fx = W.golden_min(lambda t: (t - 0.3) ** 2 + 1, -1, 2, tol=1e-10)
        self.assertAlmostEqual(x, 0.3, places=7)


class TestNoComb(unittest.TestCase):
    """v0.1.2 reported a physical 'comb' (teeth every 4 um). It is the
    lambda/2 = 0.266 um piston oscillation of imag_energy, aliased by
    1 um sampling (S2). These tests replace v0.1.2's TestComb."""

    def test_imag_energy_oscillates_with_lambda_half(self):
        sensor = angular_spectrum(make_object().astype(complex), Z_OBJ)
        zs = -Z_OBJ + np.linspace(-0.6e-6, 0.6e-6, 241)
        p = Propagator(sensor)
        ie = np.array([imag_energy(p(z)) for z in zs])
        mins = [i for i in range(1, len(ie) - 1)
                if ie[i] < ie[i - 1] and ie[i] < ie[i + 1]]
        period = float(np.median(np.diff(zs[mins])))
        self.assertAlmostEqual(period, W.WL / 2, delta=0.01e-6)

    def test_realness_is_smooth_near_focus(self):
        f0 = pipeline(make_object(), z_obj=5.03e-3, ref_phase=0.3)
        p = Propagator(f0)
        r = np.array([realness(p(z)) for z in np.linspace(-5.04e-3, -5.02e-3, 81)])
        self.assertLess((r.max() - r.min()) / r.min(), 0.05)

    def test_coarse_plus_refine_finds_off_grid_truth(self):
        """v0.1.2 needed a 0.01 mm grid; a 0.1 mm grid + golden section
        now lands within 40 um of an off-grid truth."""
        f0 = pipeline(make_object(), z_obj=5.03e-3, ref_phase=2.0)
        r = find_focus(f0, SWEEP)
        self.assertAlmostEqual(r["z"] * 1e3, -5.03, delta=0.04)


class TestNoise(unittest.TestCase):
    def test_focus_holds_at_moderate_noise(self):
        obj = make_object()
        f0 = pipeline(obj, ref_phase=1.1, noise=("gaussian", 0.20))
        r = find_focus(f0, SWEEP)
        self.assertAlmostEqual(r["z"] * 1e3, -5.00, delta=STEP_MM)
        self.assertGreater(correlation(np.abs(angular_spectrum(f0, r["z"])) ** 2,
                                       obj), 0.9)

    def test_focus_breaks_at_extreme_noise(self):
        f0 = pipeline(make_object(), ref_phase=1.1, noise=("gaussian", 4.0))
        r = find_focus(f0, SWEEP)
        self.assertLess(r["contrast"], 1.5)

    def test_poisson_shot_noise_holds(self):
        f0 = pipeline(make_object(), ref_phase=0.4, noise=("poisson", 1e5))
        self.assertAlmostEqual(find_focus(f0, SWEEP)["z"] * 1e3, -5.0, delta=0.05)


class TestCarrierEstimation(unittest.TestCase):
    def setUp(self):
        self.sensor = angular_spectrum(make_object().astype(complex), Z_OBJ)

    def test_estimates_integer_carrier(self):
        for true_car in (70, 80, 82, 90):
            h = record_hologram(self.sensor, car=true_car)
            self.assertEqual(estimate_carrier(h, subbin=False), true_car)
            self.assertAlmostEqual(estimate_carrier(h), true_car, delta=0.05)

    def test_estimates_fractional_carrier(self):
        for true_car in (80.1, 80.25, 80.5, 83.77, 89.9):
            h = record_hologram(self.sensor, car=true_car, ref_phase=1.0)
            self.assertAlmostEqual(estimate_carrier(h), true_car, delta=0.05)

    def test_estimates_2d_carrier(self):
        """Sideband U R* sits at (-carrier_y, -car); half='any' may return
        the equally bright twin, so search the left half."""
        h = record_hologram(self.sensor, car=60.3, carrier_y=-40.7)
        ky, kx = estimate_carrier_2d(h, dc_radius=30.0 / N, half="left")
        self.assertAlmostEqual(ky, 40.7, delta=0.05)
        self.assertAlmostEqual(kx, -60.3, delta=0.05)

    def test_estimates_carrier_under_noise(self):
        from noise_robustness import add_noise
        h = add_noise(record_hologram(self.sensor, car=82.4), "gaussian", 0.2,
                      np.random.default_rng(23))
        self.assertAlmostEqual(estimate_carrier(h), 82.4, delta=0.1)

    def test_subbin_demodulation_fixes_half_bin_carrier(self):
        """True carrier 80.5: rounding leaves a half-bin ramp."""
        h = record_hologram(self.sensor, car=80.5)
        sub = find_focus(reconstruct_field(h, carrier=estimate_carrier_2d(
            h, dc_radius=30.0 / N)), SWEEP)
        self.assertAlmostEqual(sub["z"] * 1e3, -5.0, delta=0.05)

    def test_assumed_wrong_carrier_collapses_contrast(self):
        h = record_hologram(self.sensor, car=82)
        bad = find_focus(reconstruct_field(h, car=80), SWEEP)
        good = find_focus(reconstruct_field(h, carrier=estimate_carrier_2d(
            h, dc_radius=30.0 / N)), SWEEP)
        self.assertGreater(good["contrast"], 2 * bad["contrast"])


class TestPhaseMetrics(unittest.TestCase):
    def _err(self, obj, metric):
        r = find_focus(pipeline(obj, ref_phase=0.7), SWEEP, metric=metric)
        return abs(r["z"] * 1e3 + 5.0), r["contrast"]

    def test_ac_realness_finds_binary_phase_objects(self):
        """New in v0.1.3: binary phase objects of any depth are 'real up
        to a constant' after removing the mean, so ac_realness finds them."""
        for depth in (np.pi / 8, np.pi / 2, 3 * np.pi / 4, np.pi):
            err, _ = self._err(np.exp(1j * depth * make_object()), "ac_realness")
            self.assertLess(err, 0.05, msg=f"depth {depth:.3f}")

    def test_ac_realness_finds_single_material_object(self):
        amp = make_object()
        err, _ = self._err((1 - 0.4 * amp) * np.exp(1j * 0.8 * amp), "ac_realness")
        self.assertLess(err, 0.05)

    def test_ac_realness_fails_for_two_materials(self):
        """Documented limit: two materials on different supports are not
        real up to one constant."""
        amp = make_object()
        obj = (0.3 + 0.7 * amp) * np.exp(1j * 1.2 * np.roll(amp, (30, -40), (0, 1)))
        err, _ = self._err(obj, "ac_realness")
        self.assertGreater(err, 0.5)

    def test_realness_fails_for_weak_phase_object(self):
        err, contrast = self._err(np.exp(1j * np.pi / 8 * make_object()), "realness")
        self.assertGreater(err, 1.0)

    def test_imag_energy_misled_for_phase_object(self):
        err, _ = self._err(np.exp(1j * np.pi / 2 * make_object()), "imag_energy")
        self.assertGreater(err, 1.0)

    def test_phase_variance_flat_for_phase_object(self):
        _, contrast = self._err(np.exp(1j * np.pi / 2 * make_object()),
                                "phase_variance")
        self.assertLess(contrast, 2.0)

    def test_amplitude_variance_near_focus_for_phase_object(self):
        err, _ = self._err(np.exp(1j * np.pi / 2 * make_object()),
                           "amplitude_variance")
        self.assertLess(err, 0.2)


class TestSpectralIdentity(unittest.TestCase):
    def test_realness_from_pair_spectrum(self):
        """sum U_z^2 = sum_k F(k)F(-k) exp(2 i z kz): realness depends only
        on the pair spectrum (S9)."""
        f0 = pipeline(make_object(), ref_phase=0.4)
        F = np.fft.fft2(f0)
        Fm = np.roll(np.flip(F, (0, 1)), 1, (0, 1))
        kz, mask = W._kz(F.shape)
        E = np.sum(np.abs(F) ** 2)
        for z in (-3e-3, -5e-3, -7e-3):
            pred = 0.5 * (1 - abs(np.sum(F * Fm * mask * np.exp(2j * z * kz))) / E)
            self.assertAlmostEqual(pred, realness(angular_spectrum(f0, z)), places=12)


class TestFalsifiedSharpness(unittest.TestCase):
    def test_sharpness_invariant_under_propagation(self):
        from wave_hologram import field_sharpness
        obj = make_object().astype(complex)
        s0 = field_sharpness(obj)
        s5 = field_sharpness(angular_spectrum(obj, -5e-3))
        self.assertLess(abs(s5 - s0) / s0, 0.01)


class TestTwoPlane(unittest.TestCase):
    def test_global_metric_sees_only_bright_plane(self):
        f0, _, _ = two_plane_field()
        r = find_focus(f0, SWEEP)
        self.assertAlmostEqual(r["z"] * 1e3, -5.0, delta=0.1)

    def test_local_focus_map_recovers_both_planes(self):
        """New in v0.1.3: per-tile realness separates the two planes."""
        f0, plane1, dots = two_plane_field()
        zmap, _ = local_focus_map(f0, SWEEP, tile=32)
        # tiles that contain a dot centre
        for cx, cy in ((58, 58), (198, 58), (128, 198)):
            z = zmap[cy // 32, cx // 32]
            self.assertAlmostEqual(z * 1e3, -7.5, delta=0.3, msg=(cx, cy))
        # central tile (cross) is on plane 1
        self.assertAlmostEqual(zmap[4, 4] * 1e3, -5.0, delta=0.2)


if __name__ == "__main__":
    unittest.main()
