import os
import re
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "scripts"))
sys.path.insert(0, ROOT)

import check_release_tag  # noqa: E402


def read(*p):
    with open(os.path.join(ROOT, *p), encoding="utf-8") as f:
        return f.read()


class ReleaseMetadata(unittest.TestCase):
    def test_versions_agree(self):
        import wave_hologram
        py = check_release_tag.pyproject_version(os.path.join(ROOT, "pyproject.toml"))
        cff = re.search(r"^version:\s*(\S+)", read("CITATION.cff"), re.M).group(1)
        self.assertEqual(py, wave_hologram.__version__)
        self.assertEqual(py, cff)

    def test_citation_keeps_author_and_orcid_free_of_invented_doi(self):
        cff = read("CITATION.cff")
        self.assertIn("Brisson", cff)
        self.assertRegex(cff, r"doi: 10\.5281/zenodo\.\d+")

    def test_tag_check(self):
        path = os.path.join(ROOT, "pyproject.toml")
        v = check_release_tag.pyproject_version(path)
        self.assertTrue(check_release_tag.check(f"refs/tags/v{v}", path)[0])
        self.assertFalse(check_release_tag.check("refs/tags/v99.0.0", path)[0])
        self.assertFalse(check_release_tag.check("refs/heads/main", path)[0])
        self.assertFalse(check_release_tag.check(f"refs/tags/{v}", path)[0])
        with tempfile.TemporaryDirectory() as d:
            bad = os.path.join(d, "pyproject.toml")
            open(bad, "w").write("[project]\nname='x'\n")
            with self.assertRaises(ValueError):
                check_release_tag.check(f"refs/tags/v{v}", bad)

    def test_publish_workflow_is_guarded(self):
        w = read(".github", "workflows", "publish.yml")
        self.assertIn("vars.PYPI_PUBLISH", w)
        self.assertIn("check_release_tag.py", w)
        self.assertIn("needs: [guard, test-gate]", w)
        self.assertIn("needs.guard.outputs.publish == 'true'", w)
        self.assertIn("test-gate:", w)
        # the only step that uploads sits behind the guarded job
        self.assertEqual(w.count("uv publish"), 1)
        self.assertLess(w.index("pypi-publish:"), w.index("uv publish"))


if __name__ == "__main__":
    unittest.main()
