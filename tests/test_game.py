"""Tests for the "Find the focus" game data and frame generator.

Run from the repository root:
    python3 -m unittest discover -s tests -v

The game data in docs/game/index.html is regenerated here from
wave_hologram.py and compared with what is recorded in the page.
"""
import base64
import contextlib
import io
import os
import sys
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "scripts"))
sys.path.insert(0, os.path.join(ROOT, "examples"))
import make_game_frames as G  # noqa: E402
import wave_hologram as W  # noqa: E402

HTML = os.path.join(ROOT, "docs", "game", "index.html")


def read_html():
    with open(HTML, encoding="utf-8") as fh:
        return fh.read()


def decode_atlas(level):
    return G.decode_png_gray8(base64.b64decode(level["atlas_png_b64"]))


class TestPng(unittest.TestCase):
    def test_roundtrip(self):
        rng = np.random.default_rng(0)
        img = rng.integers(0, G.QLEVELS, (37, 53)).astype(np.uint8)
        img[5] = 7  # a flat row exercises the other PNG filter
        out = G.decode_png_gray8(G.encode_png_gray8(img))
        np.testing.assert_array_equal(out, img)


class TestFrameGenerator(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.a = G.build_level(G.LEVELS[0])
        cls.b = G.build_level(G.LEVELS[0])

    def test_deterministic(self):
        self.assertEqual(G.dumps(self.a), G.dumps(self.b))

    def test_shapes(self):
        side = 2 * W.WIN_DEFAULT + 2
        lv = self.a
        self.assertEqual((lv["frame_w"], lv["frame_h"], lv["n_frames"]),
                         (side, side, G.N_FRAMES))
        atlas = decode_atlas(lv)
        self.assertEqual(atlas.shape, (G.N_FRAMES * side, 2 * side))
        self.assertLess(int(atlas.max()), G.QLEVELS)
        self.assertEqual(len(lv["curve_wavehologram"]), G.N_FRAMES)
        self.assertEqual(len(lv["curve_intensity_variance"]), G.N_FRAMES)
        self.assertEqual(len(lv["sweep_mm"]), G.N_FRAMES)

    def test_truth_is_inside_the_stack_and_off_centre(self):
        for spec in G.LEVELS:
            sweep = G.sweep_for(spec) * 1e3
            truth = -spec["z_obj_mm"]
            self.assertGreater(truth, sweep[0])
            self.assertLess(truth, sweep[-1])
            self.assertGreater(abs(truth - sweep[len(sweep) // 2]), 0.2,
                               "truth should not sit in the middle frame")

    def test_object_kinds_differ(self):
        objs = [G.make_scene_object(*s["obj"]) for s in G.LEVELS[:3]]
        for i in range(3):
            for j in range(i + 1, 3):
                self.assertFalse(np.array_equal(objs[i], objs[j]))
        phase = G.make_scene_object("phase_repo", 0)
        self.assertTrue(np.iscomplexobj(phase))
        np.testing.assert_allclose(np.abs(phase), 1.0)


class TestShippedGameData(unittest.TestCase):
    """The data inside docs/game/index.html matches a fresh regeneration."""

    @classmethod
    def setUpClass(cls):
        cls.html = read_html()
        cls.shipped = G.extract_data(cls.html)
        cls.fresh = G.build_game_data()

    def test_levels_and_truths(self):
        self.assertEqual(len(self.shipped["levels"]), len(G.LEVELS))
        for rec, spec in zip(self.shipped["levels"], G.LEVELS):
            self.assertEqual(rec["id"], spec["id"])
            self.assertAlmostEqual(rec["truth_mm"], -spec["z_obj_mm"], places=6)
            self.assertEqual(rec["metric"], spec["metric"])

    def test_measured_focus_matches_regeneration(self):
        for rec, new in zip(self.shipped["levels"], self.fresh["levels"]):
            # 1 um tolerance: other NumPy/FFT builds differ in the last digits
            self.assertAlmostEqual(rec["wavehologram"]["z_mm"],
                                   new["wavehologram"]["z_mm"], delta=1e-3)
            self.assertAlmostEqual(rec["wavehologram"]["error_um"],
                                   new["wavehologram"]["error_um"], delta=1.0)
            self.assertAlmostEqual(rec["naive_intensity_variance"]["z_mm"],
                                   new["naive_intensity_variance"]["z_mm"], delta=1e-6)

    def test_recorded_errors_are_consistent(self):
        for rec in self.shipped["levels"]:
            w = rec["wavehologram"]
            self.assertAlmostEqual(w["error_um"],
                                   abs(w["z_mm"] - rec["truth_mm"]) * 1e3, delta=0.01)
            n = rec["naive_intensity_variance"]
            self.assertAlmostEqual(n["error_um"],
                                   abs(n["z_mm"] - rec["truth_mm"]) * 1e3, delta=0.01)

    def test_frames_match_regeneration(self):
        for rec, new in zip(self.shipped["levels"], self.fresh["levels"]):
            a, b = decode_atlas(rec).astype(int), decode_atlas(new).astype(int)
            self.assertEqual(a.shape, b.shape)
            diff = np.abs(a - b)
            # quantisation edges may flip by one level across FFT builds (a phase pixel
            # on the wrap boundary may flip by a full turn, hence a fraction test)
            self.assertGreater(float(np.mean(diff <= 1)), 0.995)

    def test_levels_1_and_4_are_measured_within_100_um(self):
        """Sanity bound on the recorded numbers (levels 2, 3, 5 are not bounded)."""
        errs = [lv["wavehologram"]["error_um"] for lv in self.shipped["levels"]]
        self.assertTrue(all(np.isfinite(errs)))
        self.assertLess(errs[0], 100.0)          # level 1, nearly clean
        self.assertLess(errs[3], 100.0)          # phase object with ac_realness

    def test_page_is_self_contained(self):
        h = self.html
        self.assertIn("Synthetic hologram game. Not a hardware or biology validation.", h)
        self.assertIn("AGPL-3.0-only", h)
        for bad in ("<script src", "fetch(", "XMLHttpRequest", "WebSocket", "sendBeacon",
                    "<link ", "@import", "google-analytics"):
            self.assertNotIn(bad, h)
        self.assertLess(os.path.getsize(HTML), 3_000_000)

    def test_check_mode(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = G.main(["--check", "--html", HTML])
        self.assertEqual(rc, 0)


class TestQuickstart(unittest.TestCase):
    def test_quickstart_runs(self):
        import quickstart
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = quickstart.main(["--no-plot"])
        self.assertEqual(rc, 0)
        out = buf.getvalue()
        self.assertIn("measured: -4.99", out)


if __name__ == "__main__":
    unittest.main()
