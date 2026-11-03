import csv
import json
import tempfile
import unittest
from pathlib import Path

from benchmarks.run_benchmark import _run_one, run


class BenchmarkIntegrationTests(unittest.TestCase):
    def test_seeded_condition_is_deterministic(self):
        config = json.loads(Path("benchmarks/configs/waveholo_smoke.yaml").read_text())
        kwargs = dict(
            config=config, n=160, carrier=48, kind="point",
            noise_sigma=0.05, seed=22,
        )
        first, second = _run_one(**kwargs), _run_one(**kwargs)
        first.pop("runtime_seconds")
        second.pop("runtime_seconds")
        self.assertEqual(first, second)

    def test_smoke_writes_provenance_machine_results_and_gate(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "result"
            rows = run("benchmarks/configs/waveholo_smoke.yaml", output)
            self.assertEqual(len(rows), 4)
            for name in (
                "manifest.json", "environment.json", "code_commit.txt", "metrics.csv",
                "per_seed_results.csv", "summary.csv", "summary.json", "claims.json",
                "checksums.sha256", "raw_checksums.sha256",
                "figures/focus_error_vs_noise.svg", "figures/focus_metric_curve.svg",
            ):
                self.assertTrue((output / name).is_file(), name)
            manifest = json.loads((output / "manifest.json").read_text())
            self.assertIn(manifest["status"], {"COMPLETE", "INVALID"})
            self.assertEqual(manifest["protocol_version"], "PROTOCOL-v1.0.0")
            self.assertTrue(manifest["commit_sha"])
            claims = json.loads((output / "claims.json").read_text())
            self.assertEqual(claims["claim_gate_status"], "fail")
            self.assertFalse(claims["provenance_valid"] if manifest["status"] == "INVALID" else False)
            with (output / "per_seed_results.csv").open(newline="") as handle:
                self.assertGreater(len(list(csv.DictReader(handle))), 0)


if __name__ == "__main__":
    unittest.main()
