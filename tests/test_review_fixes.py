"""Tests added after the independent review: pair-spectrum realness equivalence
and bounded memory, plus guards against the wording problems that were fixed."""
import os
import re
import sys
import tracemalloc
import unittest
from unittest import mock

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import wave_hologram as W  # noqa: E402


def read(*parts):
    with open(os.path.join(ROOT, *parts), encoding="utf-8") as fh:
        return fh.read()


def demo_field():
    obj = W.make_object()
    sensor = W.angular_spectrum(obj.astype(complex), 5e-3)
    holo = W.record_hologram(sensor, 80.37, ref_phase=1.2)
    carrier = W.estimate_carrier_2d(holo, dc_radius=30.0 / W.N)
    return W.reconstruct_field(holo, carrier=carrier, crop=True), W.cropped_pitch(
        holo.shape, W.WIN_DEFAULT, W.DX)


class TestPairSpectrum(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.field, cls.pitch = demo_field()
        cls.sweep = np.linspace(-10e-3, 0.0, 101)

    def test_matches_propagation(self):
        prop = W.Propagator(self.field, dx=self.pitch)
        ref = np.array([W.realness(prop(d)) for d in self.sweep])
        got = W.pair_spectrum_realness(self.field, self.sweep, dx=self.pitch)
        np.testing.assert_allclose(got, ref, atol=1e-12)

    def test_chunking_does_not_change_the_result(self):
        full = W.pair_spectrum_realness(self.field, self.sweep, dx=self.pitch)
        with mock.patch.object(W, "_PAIR_CHUNK_ELEMENTS", 1):
            tiny = W.pair_spectrum_realness(self.field, self.sweep, dx=self.pitch)
        np.testing.assert_allclose(tiny, full, atol=1e-14)

    def test_memory_is_bounded_by_the_chunk_size(self):
        ny, nx = self.field.shape
        full_bytes = self.sweep.size * ny * nx * 16
        with mock.patch.object(W, "_PAIR_CHUNK_ELEMENTS", 2 * ny * nx):
            tracemalloc.start()
            W.pair_spectrum_realness(self.field, self.sweep, dx=self.pitch)
            _, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
        self.assertLess(peak, 0.25 * full_bytes)


class TestWording(unittest.TestCase):
    def test_readme_scopes_the_real_data_claim(self):
        r = read("README.md")
        self.assertNotIn("validated against logged stage positions", r)
        self.assertIn("slope check, not an absolute-focus test", r)
        self.assertIn("−50 to +64 µm", r)
        self.assertNotIn("up to ±50 µm", r)
        self.assertIn("## Related work", r)
        self.assertIn("100 *peak* photons", r)

    def test_name_and_orcid_in_citation(self):
        c = read("CITATION.cff")
        self.assertIn("WaveHolo (wave-hologram)", c)
        self.assertIn("https://orcid.org/0009-0000-9778-5374", c)

    def test_sandbox_does_not_call_a_placeholder_a_measurement(self):
        page = read("docs", "sandbox", "index.html")
        self.assertNotIn("WaveHolo Measured Focus Plane", page)
        self.assertIn("NOT a wave_hologram.py measurement", page)

    def test_pages_root_links_the_game(self):
        idx = read("docs", "index.html")
        self.assertIn('href="game/"', idx)
        self.assertIsNone(re.search(r"<script", idx))

    def test_third_party_title(self):
        self.assertIn("Tiny Transformer Networks", read("THIRD_PARTY_NOTICES.md"))


if __name__ == "__main__":
    unittest.main()
