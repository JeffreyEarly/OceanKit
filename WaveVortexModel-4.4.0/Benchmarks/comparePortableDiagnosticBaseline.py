#!/usr/bin/env python3
"""Compare complete, identical MATLAB-authored runs before/after diagnostics.

Each fixture consists of <family>-source.nc and <family>-request.json authored
by WVModel.writePortableRunRequest. Only file/report destinations are changed.
Run on an otherwise idle host; the first pair is warmup, then ordering alternates.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import statistics
import subprocess
import time


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def compare(args):
    root = Path(args.fixtures).resolve()
    executables = {"baseline": str(Path(args.baseline).resolve()),
                   "candidate": str(Path(args.candidate).resolve())}
    result = {"schema": "portable-diagnostic-nonregression-v1", "issue": 305,
              "maximumRegressionFraction": .03, "runsPerExecutable": args.runs,
              "method": "One excluded warmup pair, alternating process order, ratio of medians; complete integration includes scheduled output.",
              "executableSHA256": {key: digest(value) for key, value in executables.items()},
              "cases": []}
    for family in ("constant", "hydrostatic", "boussinesq"):
        source = root / f"{family}-source.nc"
        request_source = root / f"{family}-request.json"
        template = json.loads(request_source.read_text())
        samples = {key: [] for key in executables}
        for repetition in range(-1, args.runs):
            order = list(executables) if repetition % 2 == 0 else list(reversed(executables))
            for label in order:
                model = root / f"{family}-{label}.nc"
                report = root / f"{family}-{label}-report.json"
                request = root / f"{family}-{label}-request.json"
                shutil.copyfile(source, model)
                payload = dict(template, modelFiles=[str(model)], report=str(report))
                request.write_text(json.dumps(payload, indent=2) + "\n")
                started = time.perf_counter()
                process = subprocess.run([executables[label], "--request", str(request)],
                                         capture_output=True, text=True, check=False)
                elapsed = time.perf_counter() - started
                if process.returncode:
                    raise RuntimeError(f"{family}/{label}: {process.stdout}\n{process.stderr}")
                value = json.loads(report.read_text())
                sample = {"integrationSeconds": value["timingSeconds"]["integrate"],
                          "processSeconds": elapsed,
                          "retainedBytes": value["livenessBytes"]["fullModelRetained"],
                          "maximumLiveBytes": value["livenessBytes"]["fullModelMaximumLive"],
                          "stepCount": value["state"]["stepCount"],
                          "rhsEvaluationCount": value["state"]["rhsEvaluationCount"]}
                if repetition >= 0:
                    samples[label].append(sample)
        summary = {label: {"integrationSeconds": statistics.median(s["integrationSeconds"] for s in values),
                           "retainedBytes": max(s["retainedBytes"] for s in values)}
                   for label, values in samples.items()}
        for key in ("stepCount", "rhsEvaluationCount"):
            counts = {sample[key] for values in samples.values() for sample in values}
            if len(counts) != 1:
                raise RuntimeError(f"{family}: diagnostics changed {key}")
        time_change = summary["candidate"]["integrationSeconds"] / summary["baseline"]["integrationSeconds"] - 1
        memory_change = summary["candidate"]["retainedBytes"] / summary["baseline"]["retainedBytes"] - 1
        row = {"family": family, "sourceSHA256": digest(source), "requestSHA256": digest(request_source),
               "samples": samples, "medians": summary, "runtimeChangeFraction": time_change,
               "retainedMemoryChangeFraction": memory_change,
               "passes": time_change <= .03 and memory_change <= .03}
        result["cases"].append(row)
        print(f"{family}: runtime {100*time_change:+.3f}%, retained memory {100*memory_change:+.3f}%", flush=True)
    result["passes"] = all(row["passes"] for row in result["cases"])
    Path(args.output).write_text(json.dumps(result, indent=2) + "\n")
    return result["passes"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--fixtures", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--runs", type=int, default=8)
    arguments = parser.parse_args()
    if arguments.runs < 3:
        parser.error("At least three measured runs per executable are required.")
    raise SystemExit(0 if compare(arguments) else 1)
