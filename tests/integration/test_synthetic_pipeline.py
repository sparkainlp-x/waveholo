import unittest

import numpy as np

import wave_hologram as W


class SyntheticPipelineIntegrationTests(unittest.TestCase):
    def test_small_seeded_record_reconstruct_path_is_finite(self):
        n = 160
        y, x = np.mgrid[-1:1:complex(n), -1:1:complex(n)]
        obj = np.exp(-(x * x + y * y) / 0.01).astype(complex)
        sensor = W.angular_spectrum(obj, 5e-3)
        holo = W.record_hologram(sensor, car=48, ref_phase=0.73)
        carrier = W.estimate_carrier_2d(holo, dc_radius=0.08, half="left")
        restored = W.reconstruct_field(holo, carrier=carrier, win=20)
        self.assertEqual(holo.shape, (n, n))
        self.assertEqual(holo.dtype.kind, "f")
        self.assertTrue(np.isfinite(restored).all())


if __name__ == "__main__":
    unittest.main()
