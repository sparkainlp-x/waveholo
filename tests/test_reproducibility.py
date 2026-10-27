"""Golden-output tests for the documented synthetic pipeline.

Two kinds of check, chosen from what was measured on the maintainers' machine
(Python 3.9 to 3.14, numpy 2.0.2 to 2.5.4, x86-64):

* Values (focus, carrier, curve points) agree between those environments to
  about 1e-16 relative, but NOT bit for bit: `value` and `contrast` differ in the
  last printed digit between numpy 2.0.2 and 2.5.4. They are therefore compared
  with a tolerance of 1e-9 (relative), which is seven orders of magnitude above
  the observed spread and far below any physical effect.
* Arrays are compared through a SHA-256 of the intensity rounded to 1e-6 and
  stored as little-endian int64. In the cases below the closest value to a
  rounding boundary is 3.4e-6 quanta away (clean) and 1.4e-5 (noisy), against
  a floating-point spread of about 1e-10 quanta, so the hash is stable in
  practice. If a platform ever changes it, the tolerance tests above it will
  still say whether the physics changed.

Exact bit-for-bit equality of the raw floating-point arrays is NOT claimed
across platforms. See docs/validation.md and requirements-ci.txt.

If the pipeline is changed on purpose, regenerate the literals with the
snippet in docs/validation.md (Reproducibility) and say why in the change.
"""

import hashlib
import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import wave_hologram as W  # noqa: E402

SWEEP = np.linspace(-10e-3, 0.0, 101)


def quantised_sha(a, scale=1e6):
    q = np.round(np.asarray(a) * scale).astype("<i8")
    return hashlib.sha256(q.tobytes()).hexdigest()


def clean_case():
    obj = W.make_object()
    sensor = W.angular_spectrum(obj.astype(complex), 5.0e-3)
    holo = W.record_hologram(sensor, car=80.37, ref_phase=1.2)
    carrier = W.estimate_carrier_2d(holo, dc_radius=30.0 / W.N)
    field0 = W.reconstruct_field(holo, carrier=carrier)
    return holo, carrier, W.find_focus(field0, SWEEP)


def noisy_case():
    holo = clean_case()[0]
    rng = np.random.default_rng(1000)
    noisy = holo * (1 + 0.2 * rng.normal(size=holo.shape))
    carrier = W.estimate_carrier_2d(noisy, dc_radius=30.0 / W.N)
    field0 = W.reconstruct_field(noisy, carrier=carrier)
    return noisy, carrier, W.find_focus(field0, SWEEP)


class TestCleanGolden(unittest.TestCase):
    def test_values(self):
        holo, carrier, res = clean_case()
        self.assertAlmostEqual(carrier[0], 0.002500000000000932, delta=1e-9)
        self.assertAlmostEqual(carrier[1], -80.385, delta=1e-9)
        self.assertAlmostEqual(res["z"], -0.004997998171773545, delta=1e-12)   # metres
        self.assertAlmostEqual(res["value"], 0.03110957852247126, delta=1e-9)
        self.assertAlmostEqual(res["contrast"], 2.9107318318512734, delta=1e-8)
        golden = {0: 0.1610302509436733, 25: 0.08782590551321995,
                  50: 0.031109793111289474, 75: 0.0927527723617359,
                  100: 0.17703290577987596}
        for i, v in golden.items():
            self.assertAlmostEqual(res["curve"][i], v, delta=1e-9, msg=f"sweep[{i}]")

    def test_hash(self):
        holo, _, _ = clean_case()
        self.assertEqual(quantised_sha(holo),
                         "9d1dfa7e602c16ca29327f0967093a0090d24002a52cad3bd1c6b19685be446b")


class TestSeededNoiseGolden(unittest.TestCase):
    """Noise from numpy.random.default_rng(1000): same draws on numpy 2.0 to 2.5."""

    def test_values(self):
        _, carrier, res = noisy_case()
        self.assertAlmostEqual(carrier[0], 0.030000000000000956, delta=1e-9)
        self.assertAlmostEqual(carrier[1], -80.39750000000001, delta=1e-9)
        self.assertAlmostEqual(res["z"], -0.005003743086557775, delta=1e-12)

    def test_hash(self):
        noisy, _, _ = noisy_case()
        self.assertEqual(quantised_sha(noisy),
                         "57a9bc26f1312c6a9b4d358fe4af73dd8689f7621d435ca5535b523a9694623b")

    def test_the_same_seed_gives_the_same_result_and_another_seed_does_not(self):
        a = noisy_case()[2]["z"]
        b = noisy_case()[2]["z"]
        self.assertEqual(a, b)
        holo = clean_case()[0]
        other = holo * (1 + 0.2 * np.random.default_rng(1001).normal(size=holo.shape))
        self.assertNotEqual(quantised_sha(other), quantised_sha(noisy_case()[0]))


class TestHashHelper(unittest.TestCase):
    def test_the_hash_ignores_sub_quantum_differences_and_sees_larger_ones(self):
        x = np.linspace(0.1, 3.0, 50) + 0.25e-6
        self.assertEqual(quantised_sha(x), quantised_sha(x + 1e-12))
        self.assertNotEqual(quantised_sha(x), quantised_sha(x + 1e-5))


if __name__ == "__main__":
    unittest.main()
