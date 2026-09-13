#!/usr/bin/env python3
"""Independent-source #463 screen: complete flux and horizontal scalar work.

This does not decide #455 default adoption. Successful duplicate payloads are
removed after comparison; the first payload per profile/selection and every
failure remain. Original hashes, reports, logs and sample journals are retained.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess

import numpy as np


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def compare(actual, expected, complex_values):
    dtype = np.dtype("<c16" if complex_values else "<f8")
    if actual.stat().st_size != expected.stat().st_size or not actual.stat().st_size:
        raise ValueError("Payload lengths differ or are empty")
    if actual.stat().st_size % dtype.itemsize:
        raise ValueError("Payload length is not an element multiple")
    a = np.memmap(actual, dtype=dtype, mode="r")
    b = np.memmap(expected, dtype=dtype, mode="r")
    maximum = reference_maximum = error_squared = reference_squared = 0.0
    for start in range(0, a.size, 262144):
        aa, bb = a[start:start + 262144], b[start:start + 262144]
        if not np.isfinite(aa).all() or not np.isfinite(bb).all():
            raise ValueError("Nonfinite output")
        error = np.abs(aa - bb)
        reference = np.abs(bb)
        maximum = max(maximum, float(error.max()))
        reference_maximum = max(reference_maximum, float(reference.max()))
        error_squared += float(np.sum(error * error))
        reference_squared += float(np.sum(reference * reference))
    scaled = maximum / max(1.0, reference_maximum)
    relative = (error_squared / max(reference_squared, np.finfo(float).tiny)) ** 0.5
    return {"maximumScaleNormalizedError": scaled, "relativeL2Error": relative,
            "passed": scaled <= 1e-10 and relative <= 1e-10}


def source_receipt(source):
    paths = subprocess.check_output(["git", "-C", str(source), "ls-files", "-z"]).decode().split("\0")
    return {"path": str(source),
            "commit": subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip(),
            "status": subprocess.check_output(["git", "-C", str(source), "status", "--short"], text=True).strip(),
            "files": {p: digest(source / p) for p in paths
                      if p.startswith(("CompiledKernel/", "PortableRuntime/")) and (source / p).is_file()}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("control_source", type=Path)
    parser.add_argument("candidate_source", type=Path)
    parser.add_argument("control_build", type=Path)
    parser.add_argument("candidate_build", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("Workers must be positive")
    args.output.mkdir(parents=True, exist_ok=False)
    manifest = json.loads(args.manifest.read_text())
    blocks, warmups, samples = (1, 0, 1) if args.smoke else (4, 2, 4)
    builds = {"control": args.control_build.resolve(), "candidate": args.candidate_build.resolve()}
    workers = {side: {work: build / name for work, name in
               (("flux", "wv-variable-adoption"), ("scalar-horizontal", "wv-variable-derivative"))}
               for side, build in builds.items()}
    provenance = {"schema": "wvm-derivative-qualification-v1", "boundary": "direct-kernel",
                  "host": platform.platform(), "blocks": blocks, "warmups": warmups, "samples": samples,
                  "horizontalWorkers": args.workers, "fixtureManifest": manifest,
                  "manifestSHA256": digest(args.manifest), "runnerSHA256": digest(__file__),
                  "harnessSHA256": {name: digest(Path(__file__).parent / name) for name in
                                    ("WVVariableAdoptionWorker.cpp", "WVVariableDerivativeWorker.cpp", "CMakeLists.txt")},
                  "sources": {"control": source_receipt(args.control_source), "candidate": source_receipt(args.candidate_source)},
                  "workers": {side: {work: {"path": str(path), "sha256": digest(path)} for work, path in entries.items()}
                              for side, entries in workers.items()},
                  "buildCaches": {side: (path / "CMakeCache.txt").read_text() for side, path in builds.items()},
                  "environment": {k: v for k, v in os.environ.items() if k.startswith(("VECLIB_", "OMP_", "OPENBLAS_", "MKL_", "DYLD_"))},
                  "payloadRetention": "first per profile/selection plus every failure; verified duplicates removed"}
    save(args.output / "provenance.json", provenance)
    runs, retained, summary = [], {}, {}
    selections = [("production", "control", "frozen"), ("prior-pruned", "control", "pruned-streamed"),
                  ("candidate", "candidate", "pruned-streamed")]
    for case in manifest["cases"]:
        for path, sha in (("sourcePath", "sourceSHA256"), ("matlabFluxPath", "matlabFluxSHA256")):
            if digest(case[path]) != case[sha]:
                raise ValueError("Fixture identity changed: " + case[path])
        for work in ("flux", "scalar-horizontal"):
            profile = case["id"] + "-" + work
            for block in range(blocks):
                for label, side, selection in selections if block % 2 == 0 else reversed(selections):
                    name = f"{profile}-{block}-{label}"
                    payload = args.output / (name + ".bin")
                    command = [str(workers[side][work]), case["sourcePath"], selection, str(args.workers),
                               str(warmups), str(samples), str(payload)]
                    row = {"id": name, "profile": profile, "selection": label, "block": block, "command": command}
                    runs.append(row)
                    with (args.output / (name + ".stdout")).open("w") as out, (args.output / (name + ".stderr")).open("w") as err:
                        row["returncode"] = subprocess.run(command, stdout=out, stderr=err, check=False).returncode
                    save(args.output / "runs.json", runs)
                    if row["returncode"]:
                        raise RuntimeError("Worker failed: " + name)
                    report = json.loads((args.output / (name + ".stdout")).read_text())
                    row["report"] = report
                    expected_bytes = (case["Nj"] * case["Nkl"] * (1 if case["family"] == "stratified-qg" else 3) * 16
                                      if work == "flux" else int(np.prod(case["grid"])) * 8)
                    if payload.stat().st_size != expected_bytes:
                        raise ValueError("Output extent differs from declared work: " + name)
                    expected_family = {"stratified-qg": "WVTransformStratifiedQG", "hydrostatic": "WVTransformHydrostatic", "boussinesq": "WVTransformBoussinesq"}[case["family"]]
                    expected_schedule = "full-fft-gather" if selection == "frozen" else "fftw-streaming-pruned-tile16"
                    times = report["samplesSeconds"]
                    expected_schema = "wvm-variable-screening-v1" if work == "flux" else "wvm-variable-derivative-v1"
                    if (report["schema"] != expected_schema or report["selection"] != selection or report["family"] != expected_family or report["grid"] != case["grid"]
                            or report["horizontalSchedule"] != expected_schedule or report["warmups"] != warmups
                            or report["horizontalWorkers"] != (1 if selection == "frozen" else args.workers)
                            or report["matrixBackend"] != "accelerate" or report["provider"]["fftThreads"] != 1
                            or len(times) != samples or not all(np.isfinite(x) and x > 0 for x in times)):
                        raise ValueError("Worker metadata differs: " + name)
                    for key in ("baseLibrary", "threadLibrary"):
                        report["provider"][key + "SHA256"] = digest(report["provider"][key])
                    row["payloadSHA256"] = digest(payload)
                    production = retained.get((profile, "production"), payload)
                    if production == payload and label != "production":
                        raise ValueError("Production must establish comparison first")
                    row["productionComparison"] = compare(payload, production, work == "flux")
                    if work == "flux":
                        row["matlabComparison"] = compare(payload, Path(case["matlabFluxPath"]), True)
                    else:
                        analytic = report["analyticComparison"]
                        if (report.get("workload") != "scalar-horizontal" or report.get("antialias") is not False
                                or report.get("inputPreserved") is not True or not analytic["passed"]
                                or not all(np.isfinite(analytic[k]) and 0 <= analytic[k] <= 1e-10
                                           for k in ("maximumScaleNormalizedError", "relativeL2Error"))):
                            raise ValueError("Analytical scalar check failed: " + name)
                    save(args.output / "runs.json", runs)
                    if not row["productionComparison"]["passed"] or not row.get("matlabComparison", {"passed": True})["passed"]:
                        raise ValueError("Scientific comparison failed: " + name)
                    key = (profile, label)
                    if key in retained:
                        payload.unlink()
                        row["payloadRetained"] = False
                    else:
                        retained[key] = payload
                        row["payloadRetained"] = True
                    save(args.output / "runs.json", runs)
                    print(name + " passed", flush=True)
            selected = {label: [r["report"] for r in runs if r["profile"] == profile and r["selection"] == label]
                        for label, _, _ in selections}
            metrics = {label: {"processMedianSeconds": [float(np.median(r["samplesSeconds"])) for r in reports],
                              "maximumPeakRSSBytes": max(r["peakRSSBytes"] for r in reports),
                              "persistentBytes": [r["persistentBytes"] for r in reports]}
                       for label, reports in selected.items()}
            ratio = lambda label: float(np.exp(np.mean(np.log(np.array(metrics["candidate"]["processMedianSeconds"]) /
                                                      np.array(metrics[label]["processMedianSeconds"])))))
            summary[profile] = {"metrics": metrics, "candidateToPriorTime": ratio("prior-pruned"),
                                "candidateToProductionTime": ratio("production"),
                                "candidateToProductionPeakRSS": metrics["candidate"]["maximumPeakRSSBytes"] / metrics["production"]["maximumPeakRSSBytes"]}
            save(args.output / "summary.json", {"decision": "correctness-smoke" if args.smoke else "resource-screen-only",
                                               "profiles": summary})
    unchanged = (provenance["sources"] == {"control": source_receipt(args.control_source), "candidate": source_receipt(args.candidate_source)}
                 and all(digest(path) == provenance["workers"][side][work]["sha256"]
                         for side, entries in workers.items() for work, path in entries.items()))
    save(args.output / "postflight.json", {"sourcesAndWorkersUnchanged": unchanged})
    if not unchanged:
        raise ValueError("Sources or executables changed during the campaign")


if __name__ == "__main__":
    main()
