import unittest

import numpy as np

from benchmarks.focus_metrics import entropy_focus, select_focus, tenengrad_focus, variance_focus
from benchmarks.metrics import global_ssim, normalized_cross_correlation, psnr, rmse


class ImageMetricTests(unittest.TestCase):
    def test_hand_computed_rmse_and_psnr(self):
        a = np.array([[0.0, 1.0]])
        b = np.array([[0.0, 0.0]])
        self.assertAlmostEqual(rmse(a, b), np.sqrt(0.5))
        self.assertAlmostEqual(psnr(a, b), 20 * np.log10(1 / np.sqrt(0.5)))

    def test_identical_images_have_perfect_similarity(self):
        a = np.arange(16, dtype=float).reshape(4, 4) / 15
        self.assertEqual(global_ssim(a, a), 1.0)
        self.assertEqual(normalized_cross_correlation(a, a), 1.0)

    def test_constant_mismatch_has_defined_correlation(self):
        self.assertEqual(normalized_cross_correlation(np.zeros((2, 2)), np.ones((2, 2))), 0.0)


class FocusMetricTests(unittest.TestCase):
    def test_metrics_are_finite_for_valid_fields(self):
        field = np.ones((8, 8), dtype=complex)
        self.assertTrue(np.isfinite(variance_focus(field)))
        self.assertTrue(np.isfinite(tenengrad_focus(field)))
        self.assertTrue(np.isfinite(entropy_focus(field)))

    def test_focus_selection_honors_polarity(self):
        z = np.array([-1.0, 0.0, 1.0])
        self.assertEqual(select_focus(z, [2, 0, 1], "min"), (1, 0.0))
        self.assertEqual(select_focus(z, [2, 0, 1], "max"), (0, -1.0))

    def test_invalid_shapes_and_polarity_fail_closed(self):
        with self.assertRaises(ValueError):
            rmse(np.ones(3), np.ones(3))
        with self.assertRaises(ValueError):
            select_focus([0, 1], [1, 2], "middle")


if __name__ == "__main__":
    unittest.main()
