#!/usr/bin/env python3
"""Focused summary regressions for the variable full-model driver."""
import importlib.util
from pathlib import Path
import unittest


DRIVER = Path(__file__).with_name("run_variable_model_adoption.py")
SPEC = importlib.util.spec_from_file_location("variable_model_adoption", DRIVER)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def contract(qualification=False):
    return {"profileSet": "boussinesq-256", "declaredProfileCount": 1,
            "profileCount": 1, "blocks": 10,
            "qualificationOnly": qualification}


def block(index, warmup, baseline=10.0):
    runs = {"control-frozen": {
        "report": {"timingSeconds": {"integrate": baseline}}}}
    for selection, ratio in (("frozen", 1.0), ("interleaved", 0.25),
                             ("compact", 0.125)):
        runs[f"candidate-{selection}"] = {
            "report": {"timingSeconds": {"integrate": baseline * ratio}}}
    return {"profile": "variable-boussinesq-composite-256",
            "block": index, "warmup": warmup, "passed": True, "runs": runs}


class SummaryTests(unittest.TestCase):
    def test_warmup_only_defers_timing_ratios(self):
        summary = MODULE.summarize([block(0, True)], contract())
        self.assertTrue(summary["passed"])
        self.assertFalse(summary["campaignComplete"])
        self.assertEqual(summary["decision"], "warmup progress; no timing claim")
        self.assertEqual(summary["timing"], [{
            "profile": "variable-boussinesq-composite-256",
            "measuredBlocks": 0, "ratios": {}}])

    def test_measured_progress_reports_available_samples(self):
        blocks = [block(0, True), block(1, True), block(2, False)]
        summary = MODULE.summarize(blocks, contract())
        self.assertFalse(summary["campaignComplete"])
        self.assertEqual(summary["decision"],
                         "partial measured progress; no adoption timing claim")
        row = summary["timing"][0]
        self.assertEqual(row["measuredBlocks"], 1)
        self.assertEqual(row["ratios"]["compact"]["samples"], [0.125])
        self.assertEqual(
            row["ratios"]["interleaved"]["candidateOverIndependentControlMedian"],
            0.25)

    def test_complete_campaign_reports_eight_measured_blocks(self):
        blocks = [block(index, index < 2, baseline=10.0 + index)
                  for index in range(10)]
        summary = MODULE.summarize(blocks, contract())
        self.assertTrue(summary["campaignComplete"])
        self.assertEqual(summary["decision"],
                         "complete-model measurements; adoption requires the other declared gates")
        self.assertEqual(summary["timing"][0]["measuredBlocks"], 8)

    def test_qualification_never_reports_timing(self):
        qualification = contract(qualification=True)
        qualification["blocks"] = 1
        summary = MODULE.summarize([block(0, False)], qualification)
        self.assertTrue(summary["campaignComplete"])
        self.assertEqual(summary["decision"],
                         "correctness qualification; no timing claim")
        self.assertNotIn("timing", summary)


if __name__ == "__main__":
    unittest.main()
