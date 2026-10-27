"""Input validation: nonsensical physical parameters and arrays must raise a
ValueError or TypeError that names the argument, never a ZeroDivisionError,
a bare numpy message, or a silently wrong result.

Valid input is unchanged; that is covered by the existing golden tests and by
the last test class here.
"""

import os
import sys
import unittest
import warnings

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import wave_hologram as W  # noqa: E402

NAN, INF = float("nan"), float("inf")
RNG = np.random.default_rng(7)
FIELD = (RNG.normal(size=(32, 32)) + 1j * RNG.normal(size=(32, 32)))
SWEEP = np.linspace(-1e-3, 0.0, 5)
HOLO = W.record_hologram(W.angular_spectrum(FIELD, 1e-3), car=8)


class Base(unittest.TestCase):
    def raises(self, exc, fn, *args, match=None, **kw):
        with self.assertRaises(exc) as cm:
            fn(*args, **kw)
        if match is not None:
            self.assertIn(match, str(cm.exception))
        return cm.exception


class TestPitchAndWavelength(Base):
    """dx and wl reach every propagation path through _kz."""

    CALLS = {
        "angular_spectrum": lambda **kw: W.angular_spectrum(FIELD, 1e-3, **kw),
        "Propagator": lambda **kw: W.Propagator(FIELD, **kw),
        "pair_spectrum_realness": lambda **kw: W.pair_spectrum_realness(FIELD, SWEEP, **kw),
        "focus_sweep": lambda **kw: W.focus_sweep(FIELD, SWEEP, **kw),
        "find_focus": lambda **kw: W.find_focus(FIELD, SWEEP, **kw),
        "local_focus_map": lambda **kw: W.local_focus_map(FIELD, SWEEP, tile=8, **kw),
    }

    def test_non_positive_or_non_finite_values_raise_value_error(self):
        for name, call in self.CALLS.items():
            for arg in ("dx", "wl"):
                for bad in (0.0, -1e-6, NAN, INF, float("-inf")):
                    with self.subTest(fn=name, arg=arg, value=bad):
                        self.raises(ValueError, call, match=arg, **{arg: bad})

    def test_non_numbers_raise_type_error(self):
        for name, call in self.CALLS.items():
            for arg in ("dx", "wl"):
                for bad in ("6.5e-6", None, 1j, True):
                    with self.subTest(fn=name, arg=arg, value=bad):
                        self.raises(TypeError, call, match=arg, **{arg: bad})

    def test_pitch_pair_must_be_two_positive_numbers(self):
        for bad in ((1e-6, -1e-6), (1e-6, 0.0), (1e-6, NAN)):
            with self.subTest(dx=bad):
                self.raises(ValueError, W.angular_spectrum, FIELD, 1e-3, dx=bad)
        for bad in ((1e-6, 1e-6, 1e-6), (1e-6,), [[1e-6, 1e-6]], ()):
            with self.subTest(dx=bad):
                self.raises(ValueError, W.angular_spectrum, FIELD, 1e-3, dx=bad,
                            match="pair")

    def test_valid_scalar_pair_and_numpy_values_are_accepted(self):
        a = W.angular_spectrum(FIELD, 1e-3)
        for ok in (W.DX, np.float64(W.DX), (W.DX, W.DX), [W.DX, W.DX], np.array(W.DX)):
            with self.subTest(dx=ok):
                self.assertLess(np.abs(W.angular_spectrum(FIELD, 1e-3, dx=ok) - a).max(), 1e-12)
        # a float32 pitch is a different number (6.5e-6 is not representable), not an error
        b = W.angular_spectrum(FIELD, 1e-3, dx=np.float32(W.DX))
        self.assertLess(np.abs(b - a).max(), 1e-5)

    def test_record_hologram_checks_its_unused_dx_too(self):
        self.raises(ValueError, W.record_hologram, FIELD, car=8, dx=-1.0, match="dx")

    def test_cropped_pitch_checks_its_shape(self):
        for bad in ((32,), (0, 32), (32, 32, 32), 32):
            self.raises(ValueError, W.cropped_pitch, bad, 4, 1e-6, match="shape")

    def test_cropped_pitch_checks_dx(self):
        self.raises(ValueError, W.cropped_pitch, (32, 32), 4, -1e-6, match="dx")
        self.raises(ValueError, W.cropped_pitch, (32, 32), 4, 0.0, match="dx")


class TestDistances(Base):
    def test_z_must_be_one_finite_real_number(self):
        for bad in (NAN, INF, float("-inf")):
            self.raises(ValueError, W.angular_spectrum, FIELD, bad, match="z")
            self.raises(ValueError, W.Propagator(FIELD), bad, match="z")
        self.raises(ValueError, W.angular_spectrum, FIELD, np.array([1e-3, 2e-3]),
                    match="scalar")
        for bad in ("1e-3", None, 1j, True):
            self.raises(TypeError, W.angular_spectrum, FIELD, bad, match="z")

    def test_sweeps_must_be_one_dimensional_finite_and_long_enough(self):
        for bad in ([NAN, 0.0, 1e-3], [[0.0, 1e-3], [2e-3, 3e-3]], [], [0.0, INF]):
            with self.subTest(sweep=bad):
                self.raises(ValueError, W.find_focus, FIELD, bad, match="sweep")
                self.raises(ValueError, W.focus_sweep, FIELD, bad, match="sweep")
                self.raises(ValueError, W.local_focus_map, FIELD, bad, tile=8,
                            match="sweep")
        self.raises(ValueError, W.focus_sweep, FIELD, [0.0], match="at least 2")

    def test_non_realness_metrics_reject_a_two_dimensional_sweep(self):
        self.raises(ValueError, W.find_focus, FIELD, np.ones((2, 2)),
                    metric="tamura", match="sweep")

    def test_find_focus_tolerance_must_be_positive(self):
        for bad in (0.0, -1e-8, NAN):
            self.raises(ValueError, W.find_focus, FIELD, SWEEP, tol=bad, match="tol")

    def test_run_checks_its_numeric_arguments(self):
        kw = dict(make_video=False, outdir="/nonexistent-must-not-be-created")
        for name, bad in (("z_obj", NAN), ("z_obj", -1e-3), ("z_obj", 0.0),
                          ("sweep_min", NAN), ("sweep_max", INF), ("frames", 2.5),
                          ("frames", 1), ("fps", 0), ("fps", -5), ("fps", NAN)):
            with self.subTest(arg=name, value=bad):
                self.raises(ValueError, W.run, match=name, **{name: bad}, **kw)
        self.assertFalse(os.path.exists("/nonexistent-must-not-be-created"))


class TestFieldsAndImages(Base):
    def test_fields_must_be_finite_non_empty_two_dimensional(self):
        bad_fields = {
            "nan": FIELD * NAN,
            "inf": FIELD + INF,
            "1d": np.ones(8, complex),
            "3d": np.ones((2, 4, 4), complex),
            "empty": np.ones((0, 4), complex),
        }
        users = {
            "angular_spectrum": lambda f: W.angular_spectrum(f, 1e-3),
            "Propagator": lambda f: W.Propagator(f),
            "pair_spectrum_realness": lambda f: W.pair_spectrum_realness(f, SWEEP),
            "find_focus": lambda f: W.find_focus(f, SWEEP),
            "record_hologram": lambda f: W.record_hologram(f, car=2),
            "remove_quadratic_phase": lambda f: W.remove_quadratic_phase(f),
        }
        for uname, use in users.items():
            for fname, f in bad_fields.items():
                with self.subTest(fn=uname, field=fname):
                    self.raises(ValueError, use, f)

    def test_non_numeric_arrays_are_rejected(self):
        self.raises((TypeError, ValueError), W.angular_spectrum,
                    np.full((4, 4), "a"), 1e-3)

    def test_holograms_must_be_real_finite_two_dimensional(self):
        for bad in (HOLO * NAN, HOLO + INF, np.ones(8), np.ones((0, 8))):
            with self.subTest():
                self.raises(ValueError, W.reconstruct_field, bad, car=8, win=4)
                self.raises(ValueError, W.estimate_carrier_2d, bad)
        self.raises(TypeError, W.reconstruct_field, HOLO + 1j, car=8, win=4,
                    match="real")
        self.raises(TypeError, W.estimate_carrier_2d, HOLO + 1j, match="real")

    def test_a_complex_hologram_with_zero_imaginary_part_is_accepted(self):
        a = W.reconstruct_field(HOLO, car=8, win=4)
        b = W.reconstruct_field(HOLO + 0j, car=8, win=4)
        self.assertTrue(np.array_equal(a, b))

    def test_correlation_needs_matching_shapes(self):
        self.raises(ValueError, W.correlation, np.ones((2, 2)), np.ones((3, 3)))
        self.raises(ValueError, W.correlation, [1.0], [1.0])


class TestRecordingParameters(Base):
    def rec(self, **kw):
        return W.record_hologram(FIELD, **{"car": 8, **kw})

    def test_carrier_must_be_a_finite_real_in_range(self):
        for bad in (NAN, INF, 0, -3, 16, 99):
            with self.subTest(car=bad):
                self.raises(ValueError, self.rec, car=bad)
        self.raises(TypeError, self.rec, car=8 + 1j)
        self.raises(TypeError, self.rec, car="8")

    def test_vertical_carrier_must_be_below_nyquist(self):
        for bad in (NAN, INF, 16.0, -16.0, 100.0):
            with self.subTest(carrier_y=bad):
                self.raises(ValueError, self.rec, carrier_y=bad, match="carrier_y")
        self.rec(carrier_y=-15.5)        # inside the limit: accepted

    def test_reference_amplitude_must_be_positive_and_finite(self):
        for bad in (0.0, -1.0, NAN, INF):
            with self.subTest(ref_amp=bad):
                self.raises(ValueError, self.rec, ref_amp=bad, match="ref_amp")

    def test_reference_phase_must_be_finite_and_real(self):
        for bad in (NAN, INF):
            self.raises(ValueError, self.rec, ref_phase=bad, match="ref_phase")
        self.raises(TypeError, self.rec, ref_phase=1j, match="ref_phase")
        self.rec(ref_phase=-7.0)         # any finite phase is fine


class TestReconstructionParameters(Base):
    def test_carrier_must_be_a_finite_pair(self):
        for bad in ((NAN, 0.0), (0.0, INF), (1.0,), (1.0, 2.0, 3.0), "ab"):
            with self.subTest(carrier=bad):
                self.raises((ValueError, TypeError), W.reconstruct_field, HOLO,
                            carrier=bad, win=4)

    def test_window_must_be_a_finite_scalar_or_pair(self):
        for bad in (NAN, -1, (1, 2, 3), (4, NAN), "4"):
            with self.subTest(win=bad):
                self.raises((ValueError, TypeError), W.reconstruct_field, HOLO,
                            car=8, win=bad)
                self.raises((ValueError, TypeError), W.cropped_pitch, (32, 32),
                            bad, 1e-6)

    def test_carrier_search_parameters(self):
        for bad in (-0.1, NAN, INF):
            self.raises(ValueError, W.estimate_carrier_2d, HOLO, dc_radius=bad,
                        match="dc_radius")
        for bad in ((NAN, 0.0), (1.0,), (1.0, 2.0, 3.0)):
            self.raises(ValueError, W.estimate_carrier_2d, HOLO, hint=bad,
                        match="hint")
        for bad in (-1.0, 0.0, NAN):
            self.raises(ValueError, W.estimate_carrier_2d, HOLO, hint=(0.0, -8.0),
                        hint_radius=bad, match="hint_radius")


class TestOtherNumericArguments(Base):
    def test_local_focus_map_tile_and_energy_fraction(self):
        for bad in (0, -4, 2.5, 33, NAN):
            with self.subTest(tile=bad):
                self.raises(ValueError, W.local_focus_map, FIELD, SWEEP, tile=bad,
                            match="tile")
        for bad in (-0.1, 1.5, NAN):
            self.raises(ValueError, W.local_focus_map, FIELD, SWEEP, tile=8,
                        min_energy_frac=bad, match="min_energy_frac")

    def test_golden_min_arguments(self):
        f = lambda t: (t - 0.3) ** 2  # noqa: E731
        self.raises(ValueError, W.golden_min, f, 1.0, -1.0, match="a must not exceed b")
        for a, b in ((NAN, 1.0), (0.0, INF)):
            self.raises(ValueError, W.golden_min, f, a, b)
        for bad in (0.0, -1e-8, NAN):
            self.raises(ValueError, W.golden_min, f, -1.0, 1.0, tol=bad, match="tol")
        for bad in (0, 2.5, -1):
            self.raises(ValueError, W.golden_min, f, -1.0, 1.0, max_iter=bad,
                        match="max_iter")
        x, _ = W.golden_min(f, -1.0, 1.0, tol=1e-10)
        self.assertAlmostEqual(x, 0.3, places=7)

    def test_remove_quadratic_phase_arguments(self):
        for bad in (0, -1, 2.5, NAN):
            self.raises(ValueError, W.remove_quadratic_phase, FIELD, iters=bad,
                        match="iters")
        self.raises(ValueError, W.remove_quadratic_phase, np.ones((2, 8), complex),
                    match="3 x 3")

    def test_make_object_size(self):
        self.raises(TypeError, W.make_object, 256.0)
        self.raises(TypeError, W.make_object, "256")
        self.raises(TypeError, W.make_object, True)
        for bad in (0, -5, 100, 121):
            self.raises(ValueError, W.make_object, bad, match="122")
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            obj = W.make_object(128)
        self.assertEqual(obj.shape, (128, 128))
        self.assertEqual(len(w), 1)
        self.assertTrue(issubclass(w[0].category, UserWarning))
        self.assertIn("clipped", str(w[0].message))
        with warnings.catch_warnings():
            warnings.simplefilter("error")       # no warning at the default size
            W.make_object()
            W.make_object(160)

    def test_shape_arguments_of_the_geometry_cache(self):
        self.raises(ValueError, W._kz, (0, 4))
        self.raises(ValueError, W._kz, (4, 4, 4))


class TestValidInputIsUnchanged(unittest.TestCase):
    def test_the_documented_quickstart_result_is_unchanged(self):
        obj = W.make_object()
        sensor = W.angular_spectrum(obj.astype(complex), 5.0e-3)
        holo = W.record_hologram(sensor, car=80.37, ref_phase=1.2)
        carrier = W.estimate_carrier_2d(holo, dc_radius=30.0 / W.N)
        field0 = W.reconstruct_field(holo, carrier=carrier)
        res = W.find_focus(field0, np.linspace(-10e-3, 0.0, 101))
        self.assertEqual(round(res["z"] * 1e3, 3), -4.998)

    def test_integer_numpy_and_list_arguments_still_work(self):
        W.record_hologram(FIELD, car=np.int64(8), carrier_y=np.float32(1.5))
        W.reconstruct_field(HOLO, car=8, win=(4, 4))
        W.reconstruct_field(HOLO, carrier=[0, -8], win=np.int64(4))
        W.local_focus_map(FIELD, SWEEP, tile=np.int64(8))
        W.find_focus(FIELD, list(SWEEP), tol=1e-6)


if __name__ == "__main__":
    unittest.main()
