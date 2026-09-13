#!/usr/bin/env python3
"""Four-block direct-kernel screen; never makes a default-adoption decision."""
import argparse
import hashlib
import json
import os
import platform
from pathlib import Path
import subprocess
import time

import numpy as np


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def compare(actual, reference):
    a, r = np.fromfile(actual, dtype="<c16"), np.fromfile(reference, dtype="<c16")
    if a.shape != r.shape or not a.size or not np.isfinite(a).all() or not np.isfinite(r).all():
        raise ValueError("Flux payload shape or finiteness differs")
    error = a - r
    maximum = float(np.max(np.abs(error)) / max(1.0, np.max(np.abs(r))))
    norm = float(np.linalg.norm(error) / max(np.linalg.norm(r), np.finfo(float).tiny))
    return {"maximumScaleNormalizedError": maximum, "relativeL2Error": norm,
            "passed": maximum <= 1e-10 and norm <= 1e-10}



def validate_report(report, case, selection, workers, vertical_workers, samples):
    families = {"stratified-qg": "WVTransformStratifiedQG", "hydrostatic": "WVTransformHydrostatic", "boussinesq": "WVTransformBoussinesq"}
    expected_schedule = "fftw-streaming-pruned-tile16" if "pruned" in selection else "full-fft-gather"
    expected_workers = workers if "pruned" in selection else 1
    if (report["schema"] != "wvm-variable-screening-v1" or report["family"] != families[case["family"]]
            or report["grid"] != case["grid"] or report["selection"] != selection
            or report["horizontalSchedule"] != expected_schedule or report["horizontalWorkers"] != expected_workers
            or report["verticalGroupWorkers"] != (1 if selection == "frozen" else vertical_workers)
            or report["matrixBackend"] != "accelerate" or report["provider"]["fftThreads"] != 1):
        raise ValueError("Worker metadata differs from the declared workload")
    if case["family"] == "boussinesq":
        if (report.get("verticalOperatorsPerRHS", 0) <= 0 or
                report.get("verticalMatrixGroupsPerRHS", 0) <= 0 or
                report["verticalOperatorExecutionCount"] != report["verticalOperatorsPerRHS"] * samples or
                report["verticalMatrixGroupExecutionCount"] != report["verticalMatrixGroupsPerRHS"] * samples):
            raise ValueError("Boussinesq vertical execution counts are incomplete")
    times = report["samplesSeconds"]
    if len(times) != samples or not all(np.isfinite(x) and x > 0 for x in times):
        raise ValueError("Invalid timing sample inventory")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("worker", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--vertical-workers", type=int, default=1,
                        help="Prepared vertical matrix-group workers for optimized selections")
    parser.add_argument("--smoke", action="store_true", help="One block, no warmup, one sample; correctness only")
    args = parser.parse_args()
    if args.workers <= 0 or args.vertical_workers <= 0:
        parser.error("horizontal and vertical workers must be positive")
    args.output.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[2]
    worker = args.worker.resolve()
    manifest = json.loads(args.manifest.read_text())
    if args.vertical_workers > 1 and any(
            case["family"] != "boussinesq" for case in manifest["cases"]):
        raise ValueError("Vertical matrix-group screening requires a Boussinesq-only fixture manifest")
    write(args.output / "fixture-manifest.json", manifest)
    source_paths = [p for directory in ("CompiledKernel", "PortableRuntime", "Benchmarks/variable-adoption")
                    for p in (root / directory).rglob("*") if p.is_file() and p.suffix in (".cpp", ".hpp", ".py", ".m", ".txt")]
    provenance = {"schema": "wvm-variable-screening-provenance-v1", "host": platform.platform(),
                  "uname": list(platform.uname()), "worker": str(worker), "workerSHA256": sha(worker),
                  "manifestSHA256": sha(args.manifest), "fixtureSourceCommit": manifest["sourceCommit"], "sourceCommit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
                  "sources": {str(p.relative_to(root)): sha(p) for p in sorted(source_paths)},
                  "environment": {key: value for key, value in os.environ.items() if key.startswith(("VECLIB_", "OMP_", "OPENBLAS_", "MKL_", "DYLD_"))},
                  "workers": args.workers, "verticalWorkers": args.vertical_workers,
                  "blocks": 1 if args.smoke else 4, "warmups": 0 if args.smoke else 2,
                  "samples": 1 if args.smoke else 4, "decision": "correctness-smoke" if args.smoke else "screening-only"}
    (args.output / "source.diff").write_bytes(subprocess.check_output(["git", "diff", "HEAD"], cwd=root))
    write(args.output / "provenance.json", provenance)
    runs = []
    selections = ["frozen", "pruned", "streamed", "pruned-streamed"]
    for case in manifest["cases"]:
        for key, expected in (("sourcePath", "sourceSHA256"), ("matlabFluxPath", "matlabFluxSHA256")):
            if sha(case[key]) != case[expected]:
                raise ValueError("Fixture hash mismatch: " + case[key])
        first_frozen = None
        for block in range(provenance["blocks"]):
            for selection in selections if block % 2 == 0 else reversed(selections):
                name = f'{case["id"]}-{block}-{selection}'
                payload = args.output / (name + ".bin")
                command = [str(worker), case["sourcePath"], selection, str(args.workers),
                           str(provenance["warmups"]), str(provenance["samples"]), str(payload), "1",
                           str(1 if selection == "frozen" else args.vertical_workers)]
                start = time.time()
                run = {"id": name, "case": case["id"], "block": block, "selection": selection, "command": command}
                with (args.output / (name + ".stdout")).open("w") as stdout, (args.output / (name + ".stderr")).open("w") as stderr:
                    process = subprocess.run(command, stdout=stdout, stderr=stderr, check=False)
                run.update(returncode=process.returncode, elapsedSeconds=time.time()-start)
                runs.append(run)
                write(args.output / "runs.json", runs)
                if process.returncode:
                    raise RuntimeError("Worker failed; evidence retained: " + name)
                report = json.loads((args.output / (name + ".stdout")).read_text())
                run["report"] = report
                validate_report(report, case, selection, args.workers, args.vertical_workers,
                                provenance["samples"])
                for key in ("baseLibrary", "threadLibrary"):
                    report["provider"][key+"SHA256"] = sha(report["provider"][key])
                if first_frozen is None:
                    if selection != "frozen":
                        raise ValueError("First payload must establish the frozen control")
                    first_frozen = payload
                run["payloadSHA256"] = sha(payload)
                run["matlab"] = compare(payload, case["matlabFluxPath"])
                run["frozen"] = compare(payload, first_frozen)
                write(args.output / "runs.json", runs)
                if not run["matlab"]["passed"] or not run["frozen"]["passed"]:
                    raise ValueError("Scientific comparison failed: " + name)
                print(name + " passed", flush=True)
    summary = {"decision": provenance["decision"], "allScientificChecksPassed": True, "profiles": {}}
    for case in manifest["cases"]:
        profile = {}
        for selection in selections:
            selected = [r for r in runs if r["case"] == case["id"] and r["selection"] == selection]
            medians = [float(np.median(r["report"]["samplesSeconds"])) for r in selected]
            profile[selection] = {"processMediansSeconds": medians, "geometricMeanSeconds": float(np.exp(np.mean(np.log(medians)))),
                                  "maximumPeakRSSBytes": max(r["report"]["peakRSSBytes"] for r in selected)}
        for selection in selections:
            profile[selection]["ratioToFrozen"] = profile[selection]["geometricMeanSeconds"] / profile["frozen"]["geometricMeanSeconds"]
        summary["profiles"][case["id"]] = profile
    write(args.output / "summary.json", summary)


if __name__ == "__main__":
    main()
