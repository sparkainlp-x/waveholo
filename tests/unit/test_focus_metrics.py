import unittest

import numpy as np

from benchmarks.focus_metrics import entropy_focus, select_focus, tenengrad_focus, variance_focus


class FocusMetricUnitTests(unittest.TestCase):
    def test_constant_field_has_zero_variance_and_gradient(self):
        field = np.ones((5, 5), dtype=complex)
        self.assertEqual(variance_focus(field), 0.0)
        self.assertEqual(tenengrad_focus(field), 0.0)

    def test_entropy_of_constant_nonzero_image_is_zero(self):
        self.assertEqual(entropy_focus(np.ones((4, 4))), 0.0)

    def test_selection_rejects_nonfinite_scores(self):
        with self.assertRaises(ValueError):
            select_focus([0, 1], [0, np.nan], "min")


if __name__ == "__main__":
    unittest.main()
