import unittest

import numpy as np

from benchmarks.focus_metrics import entropy_focus, select_focus, tenengrad_focus, variance_focus


class FocusProperties(unittest.TestCase):
    def test_valid_random_fields_produce_finite_scores_and_in_range_extrema(self):
        for seed in range(10):
            rng = np.random.default_rng(seed)
            fields = [
                rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
                for _ in range(5)
            ]
            for metric in (variance_focus, tenengrad_focus, entropy_focus):
                values = np.asarray([metric(field) for field in fields])
                self.assertTrue(np.isfinite(values).all())
                index, _ = select_focus(np.arange(5), values, "max")
                self.assertIn(index, range(5))


if __name__ == "__main__":
    unittest.main()
