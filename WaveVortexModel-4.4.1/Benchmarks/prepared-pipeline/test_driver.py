#!/usr/bin/env python3

import array
import json
import pathlib
import subprocess
import sys
import tempfile
import unittest

import run_prepared_pipeline as driver


class PayloadComparisonTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = pathlib.Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def payload(self, name, values):
        path = self.root / name
        with path.open("wb") as stream:
            array.array("d", values).tofile(stream)
        return path

    def test_equal_payloads_pass(self):
        result = driver.compare_payload(
            self.payload("left", [0.0, -3.5, 8.0]),
            self.payload("right", [0.0, -3.5, 8.0]))
        self.assertTrue(result["passed"])
        self.assertEqual(result["maximumAbsoluteDifference"], 0.0)

    def test_tolerance_boundary_passes(self):
        result = driver.compare_payload(
            self.payload("left", [1.0, 0.0]),
            self.payload("right", [1.0 + 5e-11, 5e-13]))
        self.assertTrue(result["passed"])
        self.assertGreater(result["maximumAbsoluteDifference"], 0.0)

    def test_excessive_error_fails(self):
        with self.assertRaisesRegex(RuntimeError, "scientific|flux"):
            driver.compare_payload(self.payload("left", [1.0]),
                                   self.payload("right", [1.0 + 2e-10]))

    def test_nonfinite_payload_fails(self):
        with self.assertRaisesRegex(RuntimeError, "Nonfinite"):
            driver.compare_payload(self.payload("left", [float("inf")]),
                                   self.payload("right", [float("inf")]))

    def test_size_mismatch_fails(self):
        with self.assertRaisesRegex(RuntimeError, "lengths differ"):
            driver.compare_payload(self.payload("left", [1.0]),
                                   self.payload("right", [1.0, 2.0]))


class CampaignFailureTests(unittest.TestCase):
    def write_worker(self, root, name, ready, producers):
        path = root / name
        source = f"""\
#!/usr/bin/env python3
import array
import json
import sys

print(json.dumps({ready!r}), flush=True)
for line in sys.stdin:
    command = json.loads(line)
    if command["command"] == "quit":
        break
    payload = command.get("payload", "")
    if payload:
        with open(payload, "wb") as stream:
            array.array("d", [1.0, -2.0]).tofile(stream)
        with open(payload + ".fields.bin", "wb") as stream:
            array.array("d", [3.0, 4.0]).tofile(stream)
    print(json.dumps({{"event": "result", "producerMetrics": {producers!r}}}),
          flush=True)
"""
        path.write_text(source, encoding="utf-8")
        path.chmod(0o755)
        return path

    def ready(self, first_library, second_library, option=8):
        return {
            "event": "ready", "compiler": "fake-1", "family": "B",
            "grid": [2, 2, 2], "options": {"pointwiseWorkers": option},
            "matrixBackend": "fake-matrix", "horizontalSchedule": "fake-h",
            "source": {
                key: driver.sha256(pathlib.Path(driver.__file__).with_name(filename))
                for key, filename in (("workerSha256", "WVPreparedPipelineWorker.cpp"),
                                      ("cmakeSha256", "CMakeLists.txt"))},
            "provider": {"version": "fake", "baseLibrary": str(first_library),
                         "threadLibrary": str(second_library)},
        }

    def run_failure(self, mismatch):
        with tempfile.TemporaryDirectory() as temporary:
            root = pathlib.Path(temporary)
            input_path = root / "input.nc"
            input_path.write_bytes(b"fixture")
            libraries = []
            for name in ("base", "threads", "other"):
                library = root / name
                library.write_bytes(name.encode())
                libraries.append(library)
            baseline_ready = self.ready(libraries[0], libraries[1])
            candidate_ready = self.ready(libraries[0], libraries[1])
            baseline_producers = {"tiled": 2}
            candidate_producers = dict(baseline_producers)
            if mismatch == "options":
                candidate_ready["options"] = {"pointwiseWorkers": 7}
            elif mismatch == "provider":
                candidate_ready["provider"] = dict(candidate_ready["provider"])
                candidate_ready["provider"]["baseLibrary"] = str(libraries[2])
            elif mismatch == "stale":
                candidate_ready["source"]["workerSha256"] = "stale"
            else:
                candidate_producers["tiled"] = 3
            baseline = self.write_worker(root, "baseline", baseline_ready,
                                         baseline_producers)
            candidate = self.write_worker(root, "candidate", candidate_ready,
                                          candidate_producers)
            output = root / "results"
            command = [
                sys.executable, str(pathlib.Path(driver.__file__).resolve()),
                "--input", str(input_path), "--worker", f"base={baseline}",
                "--worker", f"candidate={candidate}",
                "--sequence", "base,candidate", "--warmups", "1",
                "--samples", "1", "--output", str(output),
            ]
            result = subprocess.run(command, text=True, capture_output=True,
                                    check=False)
            self.assertNotEqual(result.returncode, 0)
            receipt = json.loads((output / "receipt.json").read_text())
            self.assertFalse(receipt["completed"])
            self.assertIn("failure", receipt)
            self.assertIn("stale" if mismatch == "stale" else "differ",
                          receipt["failure"]["message"].lower())

    def test_campaign_records_configuration_and_producer_mismatches(self):
        for mismatch in ("options", "provider", "producers", "stale"):
            with self.subTest(mismatch=mismatch):
                self.run_failure(mismatch)


if __name__ == "__main__":
    unittest.main()
