"""Spectral autofocus equivalence, memory bounds, and optional GPU coverage."""

import unittest
from unittest.mock import patch

import numpy as np

import wave_hologram as W


class TestPairSpectrum(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(42)
        self.field = rng.normal(size=(15, 18)) + 1j * rng.normal(size=(15, 18))
        self.sweep = np.linspace(-0.005, 0.001, 11)

    def test_matches_direct_for_shapes_metrics_and_evanescent_mask(self):
        for shape in ((15, 18), (16, 17), (16, 16)):
            rng = np.random.default_rng(17)
            field = rng.normal(size=shape) + 1j * rng.normal(size=shape)
            for dx in (W.DX, (0.2 * W.WL, 0.4 * W.WL)):
                for metric in ("realness", "ac_realness"):
                    with self.subTest(shape=field.shape, dx=dx, metric=metric):
                        prop = W.Propagator(field, dx=dx)
                        expected = [W.METRICS[metric](prop(z)) for z in self.sweep]
                        for chunk in (1, 3, 100):
                            actual = W.pair_spectrum_focus_curve(
                                field, self.sweep, metric, dx=dx, chunk_size=chunk)
                            np.testing.assert_allclose(actual, expected, atol=2e-12,
                                                       rtol=0)

    def test_zero_and_constant(self):
        for metric in ("realness", "ac_realness"):
            result = W.find_focus(np.zeros((8, 16)), self.sweep, metric)
            np.testing.assert_array_equal(result["curve"], 0.5)
            self.assertEqual(result["value"], 0.5)
        for metric, value in (("realness", 0.0), ("ac_realness", 0.5)):
            curve = W.pair_spectrum_focus_curve(
                np.ones((8, 16)) * (1 + 2j), self.sweep, metric)
            np.testing.assert_allclose(curve, value, atol=3e-16)

    def test_extreme_finite_amplitudes_do_not_overflow_or_underflow(self):
        for metric in ("realness", "ac_realness"):
            expected = W.pair_spectrum_focus_curve(self.field, self.sweep, metric)
            for scale in (1e-200, 1e152, 1e200):
                with self.subTest(metric=metric, scale=scale), np.errstate(all="raise"):
                    actual = W.pair_spectrum_focus_curve(
                        self.field * scale, self.sweep, metric)
                    np.testing.assert_allclose(actual, expected, atol=2e-12)
        field = np.ones((16, 16)) * 1e152
        prop = W.Propagator(field)
        direct = [W.realness(prop(z)) for z in self.sweep]
        np.testing.assert_allclose(
            W.pair_spectrum_focus_curve(field, self.sweep), direct, atol=1e-15)
        for scale in (1e-320, 1.7e308):
            field = np.full((16, 16), complex(scale, scale))
            with self.subTest(scale=scale), np.errstate(all="raise"):
                actual = W.pair_spectrum_focus_curve(field, self.sweep)
            np.testing.assert_allclose(actual, 0, atol=1e-15)

    def test_refinement_compatibility_and_custom_metrics(self):
        field = W.angular_spectrum(W.make_object(), 0.00513)
        for name in ("realness", "ac_realness", "amplitude_variance"):
            fn = W.METRICS[name]
            expected = W.find_focus(field, self.sweep, metric=lambda f: fn(f))
            actual = W.find_focus(field, self.sweep, metric=name)
            self.assertEqual(set(actual), set(expected))
            np.testing.assert_allclose(actual["curve"], expected["curve"], atol=2e-12)
            self.assertAlmostEqual(actual["z"], expected["z"], delta=2e-8)
            self.assertAlmostEqual(actual["value"], expected["value"], delta=2e-12)

    def test_spectral_path_uses_no_inverse_ffts(self):
        with patch.object(np.fft, "ifft2", side_effect=AssertionError("inverse FFT")):
            W.find_focus(self.field, self.sweep)
            W.find_focus(self.field, self.sweep, metric=W.ac_realness)

    def test_memory_bound_and_partial_blocks(self):
        field = np.random.default_rng(4).normal(size=(512, 512))
        solver = W._PairSpectrum(field, W.DX, W.WL, W.realness, np, 10000)
        self.assertLessEqual(solver.chunk_size * solver.k.size * 16, 16 * 1024 ** 2)
        with patch.object(np, "exp", wraps=np.exp) as exponential:
            curve = solver.curve(self.sweep)
        self.assertEqual(curve.shape, self.sweep.shape)
        for call in exponential.call_args_list:
            self.assertLessEqual(call.args[0].shape[0], solver.chunk_size)
        self.assertGreater(exponential.call_count, 1)

    def test_invalid_inputs(self):
        for sweep in ([], [[0, 1]], [np.nan], [np.inf]):
            with self.subTest(sweep=sweep), self.assertRaises(ValueError):
                W.find_focus(self.field, sweep)
        for field in ([], np.ones(4), np.empty((0, 4)), [[np.nan]], [[np.inf]]):
            with self.subTest(field=field), self.assertRaises(ValueError):
                W.find_focus(field, self.sweep)
        for chunk in (0, -1, 1.5, True):
            with self.subTest(chunk=chunk), self.assertRaises(ValueError):
                W.find_focus(self.field, self.sweep, chunk_size=chunk)
        for dx, wl in ((0, W.WL), ((W.DX, -1), W.WL), (W.DX, np.nan)):
            with self.subTest(dx=dx, wl=wl), self.assertRaises(ValueError):
                W.find_focus(self.field, self.sweep, dx=dx, wl=wl)
        for tol in (0, -1, np.nan):
            with self.subTest(tol=tol), self.assertRaises(ValueError):
                W.find_focus(self.field, self.sweep, tol=tol)
        with self.assertRaises(ValueError):
            W.find_focus(self.field, self.sweep, metric="typo")
        with self.assertRaises(ValueError):
            W.pair_spectrum_focus_curve(self.field, self.sweep, metric="tamura")
        with self.assertRaises(ValueError):
            W.find_focus(self.field, self.sweep, backend="typo")

    def test_single_and_descending_sweeps(self):
        for sweep in ([0.0], self.sweep[::-1]):
            result = W.find_focus(self.field, sweep, refine=False)
            direct = W.Propagator(self.field)
            np.testing.assert_allclose(result["curve"],
                                       [W.realness(direct(z)) for z in sweep],
                                       atol=2e-12)

    def test_dhm_matches_explicit_correction_and_does_not_mutate(self):
        y, x = np.mgrid[:15, :18]
        field = self.field * np.exp(1j * (0.04 * x ** 2 + 0.02 * y ** 2))
        original = field.copy()
        corrected, _ = W.remove_quadratic_phase(field)
        for metric in ("realness", "ac_realness", "amplitude_variance"):
            actual = W.find_focus(field, self.sweep, metric, dhm_mode=True)
            expected = W.find_focus(corrected, self.sweep, metric)
            np.testing.assert_array_equal(actual["curve"], expected["curve"])
            self.assertEqual(actual["z"], expected["z"])
        np.testing.assert_array_equal(field, original)

    def test_batch_matches_individual_calls(self):
        fields = np.stack((self.field, self.field * np.exp(0.7j)))
        original = fields.copy()
        for dhm_mode in (False, True):
            for metric in ("realness", "ac_realness", "tamura"):
                results = W.find_focus_batch(fields, self.sweep, metric,
                                             dhm_mode=dhm_mode, backend="auto")
                for field, actual in zip(fields, results):
                    expected = W.find_focus(field, self.sweep, metric, dhm_mode=dhm_mode)
                    np.testing.assert_array_equal(actual["curve"], expected["curve"])
                    self.assertEqual(actual["z"], expected["z"])
        np.testing.assert_array_equal(fields, original)
        for invalid in (self.field, np.empty((0, 8, 8)), np.empty((2, 0, 8))):
            with self.assertRaises(ValueError):
                W.find_focus_batch(invalid, self.sweep)

    def test_missing_optional_backend_is_explicit(self):
        with patch.dict("sys.modules", {"cupy": None}):
            with self.assertRaisesRegex(ImportError, "separately installed CuPy"):
                W.find_focus(self.field, self.sweep, backend="cupy")
            W.find_focus(self.field, self.sweep, backend="auto")

    def test_host_batch_is_not_converted_by_gpu_backend_before_iteration(self):
        fields = np.stack((self.field, self.field))
        with patch.object(W, "_array_backend") as backend, \
                patch.object(W, "find_focus", return_value={}) as focus:
            W.find_focus_batch(fields, self.sweep, backend="cupy")
        backend.return_value.asarray.assert_not_called()
        self.assertEqual(focus.call_count, len(fields))
        for call in focus.call_args_list:
            self.assertEqual(call.args[0].shape, self.field.shape)
            self.assertTrue(np.shares_memory(call.args[0], fields))


class TestCuPyFocus(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            import cupy
        except ImportError:
            raise unittest.SkipTest("optional CuPy is not installed")
        try:
            if cupy.cuda.runtime.getDeviceCount() < 1:
                raise unittest.SkipTest("no CUDA device")
        except cupy.cuda.runtime.CUDARuntimeError as exc:
            raise unittest.SkipTest(f"CUDA unavailable: {exc}")
        cls.cp = cupy

    def test_device_batch_curves_and_dhm_match_numpy(self):
        cp = self.cp
        rng = np.random.default_rng(9)
        fields = rng.normal(size=(2, 15, 18)) + 1j * rng.normal(size=(2, 15, 18))
        device = cp.asarray(fields)
        sweep = np.linspace(-0.003, 0, 7)
        for metric in ("realness", "ac_realness"):
            for dhm in (False, True):
                expected = W.find_focus_batch(fields, sweep, metric, dhm_mode=dhm)
                actual = W.find_focus_batch(device, sweep, metric, backend="auto",
                                            dhm_mode=dhm, chunk_size=2)
                for cpu, gpu in zip(expected, actual):
                    self.assertIsInstance(gpu["curve"], np.ndarray)
                    np.testing.assert_allclose(cpu["curve"], gpu["curve"], atol=2e-12)
                    self.assertAlmostEqual(cpu["z"], gpu["z"], delta=2e-8)
        with self.assertRaisesRegex(ValueError, "only realness/ac_realness"):
            W.find_focus(device[0], sweep, metric="tamura", backend="cupy")

    def test_host_batch_uploads_only_individual_fields(self):
        cp = self.cp
        fields = np.ones((3, 8, 8), dtype=complex)
        with patch.object(cp, "asarray", wraps=cp.asarray) as upload:
            results = W.find_focus_batch(fields, [-0.001, 0], backend="cupy",
                                         refine=False)
        self.assertEqual(len(results), 3)
        for call in upload.call_args_list:
            self.assertLessEqual(np.ndim(call.args[0]), 2)


if __name__ == "__main__":
    unittest.main()
