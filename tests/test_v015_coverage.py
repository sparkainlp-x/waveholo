"""Behavioural tests for code paths the original suite left unexercised.

Every test asserts a documented behaviour: a known-answer synthetic
hologram, a focus recovered within a stated tolerance, an exact identity,
or a specific error. No library code is changed by this module.
Run: python3 -m unittest discover -s tests -v
"""
import contextlib
import io
import os
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import wave_hologram as W  # noqa: E402

PNG_MAGIC = b"\x89PNG\r\n\x1a\n"


def sensor_field(z_obj=W.Z_OBJ_DEFAULT):
    return W.angular_spectrum(W.make_object().astype(complex), z_obj)


def baseband(z_obj=W.Z_OBJ_DEFAULT, **kw):
    """Object field at the sensor plane, recorded and reconstructed."""
    holo = W.record_hologram(sensor_field(z_obj), **kw)
    carrier = W.estimate_carrier_2d(holo, dc_radius=30.0 / W.N)
    return holo, carrier, W.reconstruct_field(holo, carrier=carrier)


class TestGeometryHelpers(unittest.TestCase):
    def test_pitch_scalar_and_pair(self):
        self.assertEqual(W._pitch(6.5e-6), (6.5e-6, 6.5e-6))
        self.assertEqual(W._pitch((1e-6, 2e-6)), (1e-6, 2e-6))

    def test_kz_accepts_int_shape_and_tuple_shape_identically(self):
        kz_a, mask_a = W._kz(64)
        kz_b, mask_b = W._kz((64, 64))
        np.testing.assert_array_equal(kz_a, kz_b)
        np.testing.assert_array_equal(mask_a, mask_b)

    def test_kz_cache_is_least_recently_used(self):
        W._KZ_CACHE.clear()
        first = (8, 8)
        W._kz(first)
        for n in range(9, 9 + W._KZ_CACHE_MAX - 1):  # fill to capacity
            W._kz((n, n))
        self.assertEqual(len(W._KZ_CACHE), W._KZ_CACHE_MAX)
        W._kz(first)                  # touch the oldest entry
        W._kz((100, 100))             # forces one eviction
        keys = [k[0] for k in W._KZ_CACHE]
        self.assertIn(first, keys)    # the touched entry survived
        self.assertNotIn((9, 9), keys)  # the next oldest was evicted
        self.assertEqual(len(W._KZ_CACHE), W._KZ_CACHE_MAX)

    def test_evanescent_waves_are_zeroed_in_the_transfer_function(self):
        # pitch below lambda/2: the highest frequencies cannot propagate
        h = W._transfer(1e-3, 32, dx=0.2 * W.WL)
        kz, mask = W._kz(32, 0.2 * W.WL)
        self.assertTrue((~mask).any())
        self.assertTrue(np.all(h[~mask] == 0))
        np.testing.assert_allclose(np.abs(h[mask]), 1.0)

    def test_anisotropic_pitch_propagation_round_trips(self):
        rng = np.random.default_rng(3)
        f = rng.normal(size=(32, 48)) * np.hanning(32)[:, None] * np.hanning(48)
        pitch = (6.5e-6, 8.0e-6)
        there = W.angular_spectrum(f, 2e-3, dx=pitch)
        back = W.angular_spectrum(there, -2e-3, dx=pitch)
        self.assertLess(np.max(np.abs(back - f)), 1e-9)
        # tuple pitch must matter: the scalar-pitch result differs
        self.assertGreater(np.max(np.abs(there - W.angular_spectrum(f, 2e-3))), 1e-4)

    def test_propagator_validates_input(self):
        with self.assertRaises(ValueError):
            W.Propagator(np.ones(8))
        bad = np.ones((8, 8), dtype=complex)
        bad[0, 0] = np.nan
        with self.assertRaises(ValueError):
            W.Propagator(bad)
        p = W.Propagator(np.ones((8, 8)))
        with self.assertRaises(ValueError):
            p(np.inf)

    def test_propagator_conserves_energy_of_propagating_waves(self):
        f = sensor_field()
        p = W.Propagator(f)
        e0 = np.sum(np.abs(p(0.0)) ** 2)
        self.assertAlmostEqual(np.sum(np.abs(p(-3e-3)) ** 2) / e0, 1.0, places=9)


class TestObjectAndRecording(unittest.TestCase):
    def test_make_object_known_geometry(self):
        obj = W.make_object()
        self.assertEqual(obj.shape, (W.N, W.N))
        self.assertTrue(set(np.unique(obj)) <= {0.0, 1.0})
        c = W.N // 2
        self.assertEqual(obj[c, c], 1.0)            # the cross centre
        self.assertEqual(obj[c + 20, c + 20], 0.0)  # empty between ring and cross
        self.assertEqual(obj[c + 70, c + 0], 1.0)   # dot at (0, +70)
        with self.assertWarns(UserWarning):          # dots are clipped below n = 160
            small = W.make_object(128)
        self.assertEqual(small.shape, (128, 128))

    def test_hologram_is_real_nonnegative_with_known_mean(self):
        obj = np.zeros((64, 64), dtype=complex)
        obj[20:28, 20:28] = 0.5  # 8 px wide = exactly one carrier period (64 / 8)
        holo = W.record_hologram(obj, car=8)
        self.assertEqual(holo.dtype.kind, "f")
        self.assertGreaterEqual(holo.min(), 0.0)
        # |U+R|^2 = |U|^2 + 1 + 2 Re(U R*); the cross term has zero mean here
        np.testing.assert_allclose(holo.mean(), np.mean(np.abs(obj) ** 2) + 1.0, atol=1e-9)

    def test_record_hologram_carrier_limits(self):
        f = np.ones((64, 64), dtype=complex)
        for bad in (0, -3, 32, 40):
            with self.assertRaises(ValueError):
                W.record_hologram(f, car=bad)
        W.record_hologram(f, car=31.5)  # fractional carriers just below Nyquist are fine

    def test_reference_amplitude_scales_the_dc_term(self):
        f = np.zeros((32, 32), dtype=complex)
        holo = W.record_hologram(f, car=5, ref_amp=3.0)
        np.testing.assert_allclose(holo, 9.0)

    def test_offaxis_field_recovers_the_sensor_intensity(self):
        sf = sensor_field()
        rec = W.offaxis_field(sf, car=W.CAR_DEFAULT, win=W.WIN_DEFAULT)
        corr = W.correlation(np.abs(rec) ** 2, np.abs(sf) ** 2)
        self.assertGreater(corr, 0.95)


class TestCarrierEstimationPaths(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.holo = W.record_hologram(sensor_field(), car=80)

    def test_estimate_carrier_2d_rejects_bad_arguments(self):
        with self.assertRaises(ValueError):
            W.estimate_carrier_2d(np.ones(16))
        with self.assertRaises(ValueError):
            W.estimate_carrier_2d(self.holo, half="diagonal")

    def test_default_search_sits_on_the_left_half_plane(self):
        ky, kx = W.estimate_carrier_2d(self.holo, dc_radius=30.0 / W.N)
        self.assertAlmostEqual(ky, 0.0, delta=0.05)
        self.assertAlmostEqual(kx, -80.0, delta=0.05)

    def test_half_any_finds_a_sideband_and_top_finds_a_vertical_one(self):
        ky, kx = W.estimate_carrier_2d(self.holo, dc_radius=30.0 / W.N, half="any")
        self.assertAlmostEqual(abs(kx), 80.0, delta=0.05)
        vert = W.record_hologram(sensor_field(), car=40, carrier_y=-60.0)
        ky, kx = W.estimate_carrier_2d(vert, dc_radius=0.1, half="top", hint=(-60, 40), hint_radius=8)
        self.assertAlmostEqual(ky, -60.0, delta=0.05)
        self.assertAlmostEqual(kx, 40.0, delta=0.05)

    def test_exclude_masks_out_the_true_peak(self):
        # Blank a disc around the real sideband: the estimate must move away.
        ky, kx = W.estimate_carrier_2d(self.holo, dc_radius=30.0 / W.N,
                                       exclude=[((0, -80), 6)], refine=False)
        self.assertGreaterEqual(np.hypot(ky - 0, kx + 80), 6.0)  # outside the blanked disc

    def test_hint_restricts_search_and_empty_search_raises(self):
        ky, kx = W.estimate_carrier_2d(self.holo, dc_radius=30.0 / W.N,
                                       hint=(0, -80), hint_radius=4)
        self.assertAlmostEqual(kx, -80.0, delta=0.05)
        with self.assertRaises(ValueError):
            # the hint is in the masked-out right half-plane: nothing is left
            W.estimate_carrier_2d(self.holo, dc_radius=30.0 / W.N,
                                  hint=(0, 80), hint_radius=3)

    def test_integer_estimate_without_refinement(self):
        k = W.estimate_carrier(self.holo, subbin=False)
        self.assertIsInstance(k, int)
        self.assertEqual(k, 80)
        self.assertAlmostEqual(W.estimate_carrier(self.holo), 80.0, delta=0.05)

    def test_estimate_carrier_requires_square_2d(self):
        with self.assertRaises(ValueError):
            W.estimate_carrier(np.ones((16, 32)))
        with self.assertRaises(ValueError):
            W.estimate_carrier(np.ones(16))

    def test_local_dtft_locates_a_fractional_tone(self):
        n = 64
        x = np.arange(n)
        h = np.cos(2 * np.pi * 9.37 * x / n)[None, :] * np.ones((n, 1))
        ky, kx, F = W._local_dtft(h, 0.0, 9.0, 1.0, 0.01)
        self.assertAlmostEqual(kx, 9.37, delta=0.02)
        self.assertAlmostEqual(ky, 0.0, delta=0.02)
        self.assertEqual(F.shape[0] * F.shape[1], F.size)


class TestReconstructionPaths(unittest.TestCase):
    def test_argument_validation(self):
        h = W.record_hologram(sensor_field(), car=80)
        with self.assertRaises(ValueError):
            W.reconstruct_field(h, car=0)             # carrier out of range
        with self.assertRaises(ValueError):
            W.reconstruct_field(h, car=300)
        with self.assertRaises(ValueError):
            W.reconstruct_field(h, win=W.N // 2)      # window too large
        with self.assertRaises(ValueError):
            W.reconstruct_field(h, win=(-1, 10))
        with self.assertRaises(ValueError):
            W.reconstruct_field(h, window="triangle")
        with self.assertRaises(ValueError):
            # a carrier far outside the measured spectrum leaves nothing to keep
            W.reconstruct_field(h, carrier=(0.0, 300.0))

    def test_circle_window_finds_the_same_focus_as_the_square_window(self):
        holo, carrier, _ = baseband()
        sweep = np.linspace(-10e-3, 0.0, 101)
        zs = []
        for window in ("square", "circle"):
            f0 = W.reconstruct_field(holo, carrier=carrier, window=window)
            zs.append(W.find_focus(f0, sweep)["z"])
        self.assertAlmostEqual(zs[0], -W.Z_OBJ_DEFAULT, delta=30e-6)  # documented <= 30 um bias
        # the circular window trims more sideband: allow its (measured) larger bias
        self.assertAlmostEqual(zs[1], -W.Z_OBJ_DEFAULT, delta=50e-6)

    def test_rectangular_window_tuple(self):
        holo, carrier, _ = baseband()
        f = W.reconstruct_field(holo, carrier=carrier, win=(30, 48))
        self.assertEqual(f.shape, holo.shape)
        crop = W.reconstruct_field(holo, carrier=carrier, win=(30, 48), crop=True)
        self.assertEqual(crop.shape, (62, 98))

    def test_cropped_pitch_scales_with_the_reduced_grid(self):
        dy, dx = W.cropped_pitch((256, 256), 48, 6.5e-6)
        self.assertAlmostEqual(dy, 6.5e-6 * 256 / 98, places=15)
        self.assertEqual(dy, dx)
        dy, dx = W.cropped_pitch((256, 512), (30, 48), (6.5e-6, 8e-6))
        self.assertAlmostEqual(dy, 6.5e-6 * 256 / 62, places=15)
        self.assertAlmostEqual(dx, 8e-6 * 512 / 98, places=15)

    def test_fractional_and_integer_carrier_agree_for_an_integer_value(self):
        holo = W.record_hologram(sensor_field(), car=80)
        a = W.reconstruct_field(holo, carrier=(0.0, -80.0))
        b = W.reconstruct_field(holo, carrier=(0.0, -80.0 + 1e-9))  # fractional branch
        self.assertLess(np.max(np.abs(a - b)), 1e-6)


class TestQuadraticPhaseRemoval(unittest.TestCase):
    def test_removes_a_known_curvature_and_restores_realness(self):
        n = 96
        y, x = np.mgrid[0:n, 0:n]
        amp = np.exp(-(((x - n / 2) ** 2 + (y - n / 2) ** 2) / (2 * 22.0 ** 2)))
        u = (x - n / 2) / n
        v = (y - n / 2) / n
        curvature = 25.0 * (u ** 2 + v ** 2) + 3.0 * u - 2.0 * v
        field = amp * np.exp(1j * curvature)
        self.assertGreater(W.realness(field), 0.05)   # the curvature ruins realness
        fixed, removed = W.remove_quadratic_phase(field)
        self.assertLess(W.realness(fixed), 1e-3)
        self.assertEqual(removed.shape, field.shape)
        # the removed phase equals the imposed phase up to a constant
        diff = removed - curvature
        self.assertLess(np.std(diff[amp > 0.3]), 0.05)

    def test_input_is_not_modified(self):
        f = np.exp(1j * np.linspace(0, 1, 64 * 64).reshape(64, 64) ** 2)
        before = f.copy()
        W.remove_quadratic_phase(f, iters=1)
        np.testing.assert_array_equal(f, before)


class TestMetricEdgeCases(unittest.TestCase):
    def test_realness_extremes(self):
        self.assertEqual(W.realness(np.zeros((8, 8))), 0.5)
        real_with_piston = np.random.default_rng(0).normal(size=(16, 16)) * np.exp(0.7j)
        self.assertLess(W.realness(real_with_piston), 1e-12)
        spread = np.exp(2j * np.pi * np.arange(256) / 256)  # phases evenly spread
        self.assertAlmostEqual(W.realness(spread), 0.5, places=12)

    def test_best_piston_makes_a_field_real(self):
        base = np.random.default_rng(1).normal(size=(16, 16))
        for theta in (0.3, 1.1, -0.8):
            f = base * np.exp(1j * theta)
            p = W.best_piston(f)
            self.assertLess(W.imag_energy(f * np.exp(1j * p)), 1e-12)

    def test_imag_energy_known_values(self):
        self.assertAlmostEqual(W.imag_energy(np.ones((4, 4))), 0.0)
        self.assertAlmostEqual(W.imag_energy(1j * np.ones((4, 4))), 1.0)
        theta = 0.4
        self.assertAlmostEqual(W.imag_energy(np.exp(1j * theta) * np.ones((4, 4))),
                               np.sin(theta) ** 2)

    def test_ac_realness_ignores_a_constant_offset(self):
        r = np.random.default_rng(2).normal(size=(32, 32))
        f = (3.0 + 2.0j) + (0.5 - 0.5j) * r
        self.assertLess(W.ac_realness(f), 1e-12)
        self.assertGreater(W.realness(f), 1e-3)

    def test_amplitude_variance_and_tamura_of_flat_fields(self):
        flat = np.full((8, 8), 2.0 + 1.0j)
        self.assertAlmostEqual(W.amplitude_variance(flat), 0.0, places=12)
        self.assertAlmostEqual(W.tamura(flat), 0.0, places=6)
        self.assertEqual(W.tamura(np.zeros((8, 8))), 0.0)
        varied = np.where(np.arange(64).reshape(8, 8) % 2, 1.0, 3.0)
        self.assertGreater(W.amplitude_variance(varied), 1.0)
        self.assertGreater(W.tamura(varied), 0.0)

    def test_phase_variance_limits(self):
        self.assertEqual(W.phase_variance(np.zeros((4, 4))), 1.0)
        self.assertAlmostEqual(W.phase_variance(np.exp(0.9j) * np.ones((4, 4))), 0.0, places=12)
        spread = np.exp(2j * np.pi * np.arange(64) / 64).reshape(8, 8)
        self.assertAlmostEqual(W.phase_variance(spread), 1.0, places=12)

    def test_field_sharpness_is_positive_for_structured_fields(self):
        self.assertAlmostEqual(W.field_sharpness(np.ones((8, 8), dtype=complex)), 0.0)
        self.assertGreater(W.field_sharpness(sensor_field()), 0.0)

    def test_registry_is_complete_and_finite(self):
        f = sensor_field()
        self.assertEqual(set(W.METRICS), {"realness", "ac_realness", "imag_energy",
                                          "phase_variance", "amplitude_variance",
                                          "neg_tamura", "tamura"})
        for name, fn in W.METRICS.items():
            self.assertTrue(np.isfinite(fn(f)), name)
        self.assertEqual(W.METRICS["neg_tamura"](f), -W.METRICS["tamura"](f))

    def test_metric_lookup(self):
        self.assertIs(W._metric("realness"), W.realness)
        fn = lambda f: 1.0  # noqa: E731
        self.assertIs(W._metric(fn), fn)
        with self.assertRaises(ValueError) as cm:
            W._metric("nope")
        self.assertIn("realness", str(cm.exception))

    def test_correlation_known_values(self):
        a = np.random.default_rng(4).normal(size=(10, 10))
        self.assertAlmostEqual(W.correlation(a, a), 1.0)
        self.assertAlmostEqual(W.correlation(a, -a), -1.0)
        self.assertAlmostEqual(W.correlation(a, 5 * a + 2), 1.0)


class TestPairSpectrum(unittest.TestCase):
    def test_validation_and_degenerate_inputs(self):
        ok = np.ones((16, 16), dtype=complex)
        with self.assertRaises(ValueError):
            W.pair_spectrum_realness(np.ones(16), [0.0])
        bad = ok.copy()
        bad[1, 1] = np.inf
        with self.assertRaises(ValueError):
            W.pair_spectrum_realness(bad, [0.0])
        with self.assertRaises(ValueError):
            W.pair_spectrum_realness(ok, [[0.0, 1e-3]])
        with self.assertRaises(ValueError):
            W.pair_spectrum_realness(ok, [0.0, np.nan])
        out = W.pair_spectrum_realness(np.zeros((16, 16), dtype=complex), [0.0, 1e-3])
        np.testing.assert_array_equal(out, [0.5, 0.5])

    def test_matches_direct_propagation_on_a_rectangular_complex_field(self):
        rng = np.random.default_rng(11)
        f = (rng.normal(size=(40, 56)) + 1j * rng.normal(size=(40, 56)))
        f *= np.hanning(40)[:, None] * np.hanning(56)[None, :]
        sweep = np.linspace(-3e-3, 3e-3, 7)
        fast = W.pair_spectrum_realness(f, sweep)
        direct = np.array([W.realness(W.angular_spectrum(f, z)) for z in sweep])
        np.testing.assert_allclose(fast, direct, atol=1e-10)


class TestFocusSearch(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        _, _, cls.field0 = baseband()
        cls.sweep = np.linspace(-10e-3, 0.0, 101)
        cls.truth = -W.Z_OBJ_DEFAULT

    def test_parabolic_min_recovers_an_exact_parabola(self):
        s = np.linspace(0.0, 1.0, 11)
        m = (s - 0.437) ** 2 + 0.25
        z, v = W.parabolic_min(s, m)
        self.assertAlmostEqual(z, 0.437, places=12)
        self.assertAlmostEqual(v, 0.25, places=12)

    def test_parabolic_min_edges_and_flat_curves_fall_back_to_the_grid(self):
        s = np.array([0.0, 1.0, 2.0, 3.0])
        self.assertEqual(W.parabolic_min(s, [0.0, 1.0, 2.0, 3.0]), (0.0, 0.0))
        self.assertEqual(W.parabolic_min(s, [3.0, 2.0, 1.0, 0.0]), (3.0, 0.0))
        # a symmetric two-point plateau: the vertex is the midpoint
        z, v = W.parabolic_min(s, [2.0, 1.0, 1.0, 2.0])
        self.assertAlmostEqual(z, 1.5)
        self.assertAlmostEqual(v, 0.875)
        flat = W.parabolic_min(s, [1.0, 1.0, 1.0, 1.0])
        self.assertEqual(flat, (0.0, 1.0))
        z, _ = W.parabolic_min(s, [1.0, 1.0, 1.0, 2.0])
        self.assertEqual(z, 0.0)  # argmin is the first of a flat run, which is an edge

    @unittest.expectedFailure
    def test_parabolic_min_with_a_nan_metric_should_return_a_finite_distance(self):
        """KNOWN ISSUE (reported, not fixed here): np.argmin returns the index of
        the first NaN, so a single NaN in the metric curve makes parabolic_min
        return (nan, nan) instead of the best finite grid point."""
        z, v = W.parabolic_min(np.arange(5.0), [3.0, 2.0, np.nan, 2.0, 3.0])
        self.assertTrue(np.isfinite(z) and np.isfinite(v))

    def test_parabolic_min_refines_a_real_focus_curve(self):
        _, vals, i = W.focus_sweep(self.field0, self.sweep, keep_fields=False)
        z, _ = W.parabolic_min(self.sweep, vals)
        self.assertLess(abs(z - self.truth), 25e-6)

    def test_golden_min_hits_a_known_minimum(self):
        x, fx = W.golden_min(lambda t: (t - 0.3) ** 2 + 1.0, -1.0, 2.0, tol=1e-9)
        self.assertAlmostEqual(x, 0.3, places=6)
        self.assertAlmostEqual(fx, 1.0, places=10)
        # asymmetric bracket pushes the iteration down both branches
        x, _ = W.golden_min(lambda t: (t + 0.9) ** 2, -1.0, 5.0, tol=1e-9)
        self.assertAlmostEqual(x, -0.9, places=6)

    def test_focus_sweep_validation_and_options(self):
        with self.assertRaises(ValueError):
            W.focus_sweep(self.field0, [0.0])
        with self.assertRaises(ValueError):
            W.focus_sweep(self.field0, [0.0, np.nan])
        with self.assertRaises(ValueError):
            W.focus_sweep(self.field0, self.sweep, metric="bogus")
        fields, vals, i = W.focus_sweep(self.field0, self.sweep[::10], keep_fields=True)
        self.assertEqual(len(fields), len(self.sweep[::10]))
        self.assertEqual(len(vals), len(fields))
        none, vals2, i2 = W.focus_sweep(self.field0, self.sweep[::10], keep_fields=False)
        self.assertIsNone(none)
        np.testing.assert_array_equal(vals, vals2)
        self.assertEqual(i, i2)

    def test_focus_sweep_accepts_a_callable_metric(self):
        calls = []

        def probe(f):
            calls.append(1)
            return W.realness(f)

        _, vals, _ = W.focus_sweep(self.field0, self.sweep[::20], metric=probe,
                                   keep_fields=False)
        self.assertEqual(len(calls), len(vals))

    def test_find_focus_without_refinement_returns_the_grid_point(self):
        r = W.find_focus(self.field0, self.sweep, refine=False)
        self.assertEqual(r["z"], r["coarse_z"])
        self.assertIn(r["z"], list(self.sweep))
        self.assertLessEqual(abs(r["z"] - self.truth), 0.1e-3 / 2 + 1e-12)

    def test_find_focus_refinement_beats_the_grid(self):
        r = W.find_focus(self.field0, self.sweep)
        self.assertLess(abs(r["z"] - self.truth), 30e-6)  # documented <= 30 um bias
        self.assertGreater(r["contrast"], 5.0)
        self.assertEqual(len(r["curve"]), len(self.sweep))

    def test_refinement_is_sub_micron_without_sideband_crosstalk(self):
        """On the true sensor field the 0.1 mm grid + golden section lands on an
        off-grid truth (5.037 mm) within 0.5 um."""
        sf = sensor_field(5.037e-3)
        r = W.find_focus(sf, self.sweep, tol=1e-10)
        self.assertLess(abs(r["z"] + 5.037e-3), 0.5e-6)
        self.assertGreater(abs(r["coarse_z"] + 5.037e-3), 5e-6)  # the grid alone is worse

    def test_find_focus_with_other_metrics_and_a_two_point_sweep(self):
        r = W.find_focus(self.field0, self.sweep, metric="ac_realness")
        self.assertLess(abs(r["z"] - self.truth), 40e-6)
        two = W.find_focus(self.field0, self.truth + np.array([-1e-3, 1e-3]))
        self.assertEqual(two["z"], two["coarse_z"])  # no refinement with < 3 points

    def test_find_focus_with_callable_realness_uses_the_same_answer(self):
        a = W.find_focus(self.field0, self.sweep, metric="realness")
        b = W.find_focus(self.field0, self.sweep, metric=W.realness)
        self.assertEqual(a["z"], b["z"])

    def test_find_focus_of_a_dark_field_reports_infinite_or_unit_contrast(self):
        zero = np.zeros((32, 32), dtype=complex)
        r = W.find_focus(zero, np.linspace(-1e-3, 0.0, 5))
        self.assertTrue(np.all(r["curve"] == 0.5))
        self.assertEqual(r["contrast"], 1.0)


class TestLocalFocusMapPaths(unittest.TestCase):
    def test_rejects_unsupported_metrics(self):
        with self.assertRaises(ValueError):
            W.local_focus_map(np.ones((32, 32), dtype=complex), [0.0, 1e-3],
                              metric="tamura")

    def test_ac_realness_map_and_dark_tiles(self):
        # One bright amplitude block on a dark background; the other tiles stay dark.
        obj = np.zeros((128, 128), dtype=complex)
        obj[40:88, 40:88] = 1.0
        sf = W.angular_spectrum(obj, 4e-3)
        sweep = np.linspace(-6e-3, -2e-3, 41)
        zmap, curves = W.local_focus_map(sf, sweep, tile=32, metric="ac_realness")
        self.assertEqual(zmap.shape, (4, 4))
        self.assertEqual(curves.shape, (41, 4, 4))
        self.assertTrue(np.isnan(zmap[0, 0]))          # corner tile carries no energy
        self.assertTrue(np.isfinite(zmap[1:3, 1:3]).all())

    def test_realness_map_reports_the_object_plane_for_a_single_plane(self):
        obj = np.zeros((128, 128), dtype=complex)
        obj[32:96, 32:96] = np.random.default_rng(5).uniform(0.5, 1.0, (64, 64))
        sf = W.angular_spectrum(obj, 3e-3)
        sweep = np.linspace(-5e-3, -1e-3, 81)
        zmap, _ = W.local_focus_map(sf, sweep, tile=32, metric="realness")
        inner = zmap[1:3, 1:3]
        self.assertTrue(np.all(np.abs(inner + 3e-3) < 0.2e-3))


class TestOutputs(unittest.TestCase):
    def test_plot_focus_curve_writes_a_png(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "curve.png")
            sweep = np.linspace(-10e-3, 0, 11)
            out = W.plot_focus_curve(sweep, (sweep + 5e-3) ** 2, -5e-3, -5e-3, path)
            self.assertEqual(out, path)
            with open(path, "rb") as fh:
                self.assertEqual(fh.read(8), PNG_MAGIC)

    def test_render_video_writes_frames_and_invokes_ffmpeg(self):
        rng = np.random.default_rng(6)
        fields = [rng.normal(size=(16, 16)) + 1j * rng.normal(size=(16, 16)) for _ in range(3)]
        sweep = np.array([-2e-3, -1e-3, 0.0])
        with tempfile.TemporaryDirectory() as d:
            frames = os.path.join(d, "frames")
            with mock.patch.object(W.subprocess, "run") as run:
                out = W.render_video(fields, sweep, 12, os.path.join(d, "v.mp4"), frames)
            self.assertEqual(out, os.path.join(d, "v.mp4"))
            names = sorted(os.listdir(frames))
            self.assertEqual(names, ["f0000.png", "f0001.png", "f0002.png"])
            with open(os.path.join(frames, names[0]), "rb") as fh:
                self.assertEqual(fh.read(8), PNG_MAGIC)
            cmd = run.call_args[0][0]
            self.assertEqual(cmd[0], "ffmpeg")
            self.assertIn("12", cmd)
            self.assertTrue(cmd[-1].endswith("v.mp4"))

    def test_render_video_reports_a_missing_or_failing_ffmpeg(self):
        fields = [np.ones((8, 8), dtype=complex)] * 2
        for exc in (FileNotFoundError("ffmpeg"),
                    subprocess.CalledProcessError(1, "ffmpeg")):
            with tempfile.TemporaryDirectory() as d:
                with mock.patch.object(W.subprocess, "run", side_effect=exc):
                    with self.assertRaises(RuntimeError) as cm:
                        W.render_video(fields, [0.0, 1e-3], 10,
                                       os.path.join(d, "v.mp4"), os.path.join(d, "f"))
                self.assertIn("--no-video", str(cm.exception))


class TestRunAndCli(unittest.TestCase):
    def test_run_recovers_the_focus_and_writes_outputs(self):
        with tempfile.TemporaryDirectory() as d:
            r = W.run(frames=51, outdir=d, make_video=False)
            self.assertLess(r["error_um"], 30.0)  # documented sideband bias
            self.assertAlmostEqual(r["truth_mm"], -5.0)
            self.assertAlmostEqual(r["carrier_est"][1], -80.0, delta=0.05)
            self.assertGreater(r["corr"], 0.9)
            self.assertGreater(r["contrast"], 5.0)
            for name in ("focus_curve.png", "best_focus.png"):
                with open(os.path.join(d, name), "rb") as fh:
                    self.assertEqual(fh.read(8), PNG_MAGIC)
            self.assertFalse(os.path.exists(os.path.join(d, "focus_sweep.mp4")))

    def test_run_with_a_reference_piston_and_assumed_carrier(self):
        with tempfile.TemporaryDirectory() as d:
            r = W.run(frames=51, outdir=d, make_video=False, ref_phase=1.9,
                      estimate=False, z_obj=4.0e-3, sweep_min=-8e-3)
            self.assertLess(r["error_um"], 30.0)  # realness is piston-invariant
            self.assertEqual(r["carrier_est"], (0.0, -float(W.CAR_DEFAULT)))

    def test_run_with_video_calls_the_encoder(self):
        with tempfile.TemporaryDirectory() as d:
            with mock.patch.object(W.subprocess, "run") as enc:
                W.run(frames=5, outdir=d, make_video=True)
            enc.assert_called_once()
            self.assertEqual(len(os.listdir(os.path.join(d, "frames"))), 5)

    def test_run_validates_its_arguments(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(ValueError):
                W.run(z_obj=0.0, outdir=d, make_video=False)
            with self.assertRaises(ValueError):
                W.run(z_obj=-1e-3, outdir=d, make_video=False)
            with self.assertRaises(ValueError):
                W.run(sweep_min=0.0, sweep_max=-1e-3, outdir=d, make_video=False)
            with self.assertRaises(ValueError):
                W.run(frames=1, outdir=d, make_video=False)

    def test_main_prints_a_report_and_returns_zero(self):
        with tempfile.TemporaryDirectory() as d:
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rc = W.main(["--no-video", "--frames", "41", "--outdir", d,
                             "--z-obj", "6.0", "--sweep-min", "-12", "--metric",
                             "ac_realness"])
            text = buf.getvalue()
            self.assertEqual(rc, 0)
            self.assertIn(f"wave-hologram {W.__version__}  metric: ac_realness", text)
            self.assertIn("true object distance : -6.0000 mm", text)
            self.assertIn("focus error", text)
            self.assertIn(d, text)
            err = float(text.split("focus error")[1].split(":")[1].split("um")[0])
            self.assertLess(err, 40.0)

    def test_main_rejects_an_unknown_metric(self):
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as cm:
                W.main(["--metric", "nope"])
        self.assertEqual(cm.exception.code, 2)

    def test_module_runs_as_a_script(self):
        proc = subprocess.run([sys.executable, os.path.join(HERE, "wave_hologram.py"), "--help"],
                              capture_output=True, text=True, timeout=120)
        self.assertEqual(proc.returncode, 0)
        self.assertIn("--no-video", proc.stdout)
        self.assertIn("--assume-carrier", proc.stdout)


if __name__ == "__main__":
    unittest.main()
