import unittest

import numpy as np

from benchmarks.metrics import global_ssim, normalized_cross_correlation, rmse


class MetricGoldenTests(unittest.TestCase):
    def test_fixed_small_image_metrics(self):
        reference = np.array([[0.0, 0.5], [0.5, 1.0]])
        estimate = np.array([[0.0, 0.5], [0.0, 1.0]])
        self.assertAlmostEqual(rmse(reference, estimate), 0.25, places=14)
        self.assertAlmostEqual(normalized_cross_correlation(reference, estimate), 0.8528028654224417, places=14)
        self.assertAlmostEqual(global_ssim(reference, estimate), 0.8088878131423021, places=10)


if __name__ == "__main__":
    unittest.main()
