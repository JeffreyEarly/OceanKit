#!/usr/bin/env python3
"""Source-bound #358 qualification; process RSS is collected after child exit.

Author-only dependency: NumPy. No runtime or MATLAB dependencies are introduced.
Each invocation creates a new evidence folder and refuses to overwrite it.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def compare(actual_path, expected_path, max_tolerance, l2_tolerance):
    if actual_path.stat().st_size != expected_path.stat().st_size:
        raise ValueError("Flux payload lengths differ")
    actual = np.memmap(actual_path, dtype="<c16", mode="r")
    expected = np.memmap(expected_path, dtype="<c16", mode="r")
    max_error = max_scale = squared_error = squared_scale = 0.0
    for begin in range(0, actual.size, 262144):
        a, b = actual[begin:begin + 262144], expected[begin:begin + 262144]
        if not (np.isfinite(a).all() and np.isfinite(b).all()):
            raise ValueError("Nonfinite flux coefficient")
        error = np.abs(a - b)
        scale = np.abs(b)
        max_error = max(max_error, float(error.max(initial=0)))
        max_scale = max(max_scale, float(scale.max(initial=0)))
        squared_error += float(np.sum(error * error))
        squared_scale += float(np.sum(scale * scale))
    tiny = np.finfo(float).tiny
    maximum = max_error / max(max_scale, tiny)
    relative_l2 = math.sqrt(squared_error / max(squared_scale, tiny))
    return {"maximumScaleNormalizedError": maximum, "relativeL2Error": relative_l2,
            "maximumTolerance": max_tolerance, "relativeL2Tolerance": l2_tolerance,
            "passed": maximum <= max_tolerance and relative_l2 <= l2_tolerance}


def run_child(command, prefix, parse_worker=True):
    started = time.monotonic()
    with prefix.with_suffix(".stdout").open("wb") as stdout, prefix.with_suffix(".stderr").open("wb") as stderr:
        child = subprocess.Popen(command, stdout=stdout, stderr=stderr)
        # Do not call poll/wait/communicate before wait4: those can reap the
        # process and lose this child's complete-lifetime resource accounting.
        _, status, usage = os.wait4(child.pid, 0)
        child.returncode = os.waitstatus_to_exitcode(status)
    result = {"command": command, "exitCode": child.returncode,
              "completeLifetimeSeconds": time.monotonic() - started,
              "completeLifetimePeakRSSBytes": int(usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024)),
              "userSeconds": usage.ru_utime, "systemSeconds": usage.ru_stime}
    if child.returncode == 0 and parse_worker:
        result["worker"] = json.loads(prefix.with_suffix(".stdout").read_text())
    return result


def summarize(records):
    profiles = sorted({r["profile"] for r in records})
    ratios, memory_ratios, details = [], [], []
    for profile in profiles:
        current = [r for r in records if r["profile"] == profile]
        pairs = sorted({r["pair"] for r in current})
        timing, memory = [], []
        for pair in pairs:
            by_role = {r["role"]: r for r in current if r["pair"] == pair}
            if set(by_role) != {"control", "candidate"}:
                continue
            control, candidate = by_role["control"], by_role["candidate"]
            timing.append(candidate["worker"]["timing"]["medianSeconds"] / control["worker"]["timing"]["medianSeconds"])
            memory.append(candidate["completeLifetimePeakRSSBytes"] / control["completeLifetimePeakRSSBytes"])
        if timing:
            ratios.append(np.log(timing))
            memory_ratios.append(np.log(memory))
            details.append({"profile": profile, "pairs": len(timing),
                            "fluxRatio": float(np.exp(np.mean(np.log(timing)))),
                            "peakRSSRatio": float(np.exp(np.mean(np.log(memory))))})
    if not ratios:
        return {"profiles": [], "decision": "qualification-only"}
    rng = np.random.default_rng(358)
    boot = np.zeros(10000)
    memory_boot = np.zeros(10000)
    for timing, memory in zip(ratios, memory_ratios):
        indices = rng.integers(0, len(timing), size=(10000, len(timing)))
        boot += timing[indices].mean(axis=1) / len(ratios)
        memory_boot += memory[indices].mean(axis=1) / len(ratios)
    return {"profiles": details,
            "geometricFluxRatio": float(np.exp(np.mean([r.mean() for r in ratios]))),
            "fluxPairedBootstrap95": np.exp(np.quantile(boot, [.025, .975])).tolist(),
            "geometricPeakRSSRatio": float(np.exp(np.mean([r.mean() for r in memory_ratios]))),
            "peakRSSPairedBootstrap95": np.exp(np.quantile(memory_boot, [.025, .975])).tolist(),
            "decision": "requires owned-memory, full-model and MATLAB-loaded gates"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixtures", required=True, type=Path)
    parser.add_argument("--control", required=True, type=Path)
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mode", choices=["qualification", "screen", "final"], default="qualification")
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--provenance", type=Path, help="Frozen source/build/provider/host receipt; required for timing")
    parser.add_argument("--interface", choices=["native", "matlab-loaded"], default="native")
    parser.add_argument("--source-root", type=Path, help="MATLAB authoring checkout, required for MATLAB-loaded execution")
    parser.add_argument("--provider-root", type=Path, help="Expected pinned FFTW prefix for MATLAB-loaded execution")
    parser.add_argument("--profile", action="append", help="Exact profile ID; default all frozen cases")
    args = parser.parse_args()
    if args.mode != "qualification" and not args.candidate:
        parser.error("timing campaigns require a candidate")
    if args.mode != "qualification" and not args.provenance:
        parser.error("timing campaigns require a frozen source/build/provider/host receipt")
    if args.interface == "matlab-loaded" and (not args.source_root or not args.provider_root):
        parser.error("MATLAB-loaded execution requires source-root and provider-root")
    args.output.mkdir(parents=True, exist_ok=False)
    fixtures = json.loads(args.fixtures.read_text())
    cases = [case for case in fixtures["cases"] if not args.profile or case["id"] in args.profile]
    if not cases or (args.profile and set(args.profile) != {case["id"] for case in cases}):
        raise ValueError("Unknown or empty profile selection")
    for case in cases:
        for path_key, hash_key in [("sourcePath", "sourceSHA256"), ("matlabFluxPath", "matlabFluxSHA256")]:
            if sha256(case[path_key]) != case[hash_key]:
                raise ValueError(f"Frozen fixture hash mismatch: {case[path_key]}")
    executables = {"control": args.control.resolve()}
    if args.candidate:
        executables["candidate"] = args.candidate.resolve()
    warmups, samples, pairs = {"qualification": (0, 1, 1), "screen": (2, 3, 1), "final": (2, 7, 8)}[args.mode]
    save(args.output / "manifest.json", {"schema": "wvm-constant-adoption-campaign-v1",
         "mode": args.mode, "platform": platform.platform(), "threads": args.threads,
         "interface": args.interface,
         "warmups": warmups, "samples": samples, "pairs": pairs,
         "fixtures": fixtures, "fixturesManifestSHA256": sha256(args.fixtures),
         "harnessSHA256": sha256(__file__),
         "provenance": json.loads(args.provenance.read_text()) if args.provenance else None,
         "provenanceSHA256": sha256(args.provenance) if args.provenance else None,
         "executables": {role: {"path": str(path), "sha256": sha256(path)} for role, path in executables.items()}})
    records = []
    for case in cases:
        for pair in range(pairs):
            roles = list(executables)
            if pair % 2:
                roles.reverse()
            pair_records = {}
            for role in roles:
                prefix = args.output / f"{case['id']}-pair{pair:02}-{role}"
                binary = prefix.with_suffix(".bin")
                print(f"{args.mode}: {case['id']} pair {pair + 1}/{pairs} {role}", flush=True)
                if args.interface == "native":
                    result = run_child([str(executables[role]), case["sourcePath"], str(args.threads), str(warmups), str(samples), str(binary)], prefix)
                else:
                    module = executables[role]
                    config_path = prefix.with_suffix(".config.json")
                    result_path = prefix.with_suffix(".worker.json")
                    save(config_path, {"moduleFolder": str(module.parent), "moduleName": module.stem,
                         "modulePath": str(module), "moduleSHA256": sha256(module),
                         "providerRoot": str(args.provider_root.resolve()), "sourcePath": case["sourcePath"],
                         "threads": args.threads, "warmups": warmups, "samples": samples,
                         "fluxPath": str(binary.resolve()), "resultPath": str(result_path.resolve())})
                    root = args.source_root.resolve()
                    def literal(path):
                        return "'" + str(path).replace("'", "''") + "'"
                    statement = (f"cd({literal(root)});addpath({literal(root / 'tools')});"
                                 f"configureCIEnvironment({literal(root)},{literal(root.parent / 'OceanKit')});"
                                 f"addpath({literal(root / 'Benchmarks/constant-adoption')});"
                                 f"constantAdoptionMexWorker({literal(config_path.resolve())});")
                    result = run_child(["matlab", "-batch", statement], prefix, parse_worker=False)
                    if result["exitCode"] == 0:
                        result["worker"] = json.loads(result_path.read_text())
                result.update(profile=case["id"], pair=pair, role=role)
                if result["exitCode"] == 0:
                    result["matlabComparison"] = compare(binary, Path(case["matlabFluxPath"]), 1e-10, 1e-10)
                    result["fluxSHA256"] = sha256(binary)
                save(prefix.with_suffix(".json"), result)
                records.append(result)
                save(args.output / "runs.json", records)
                if result["exitCode"] or not result["matlabComparison"]["passed"]:
                    raise RuntimeError(f"Qualification failed; evidence retained in {prefix}")
                pair_records[role] = (binary, result, prefix)
            if len(pair_records) == 2:
                comparison = compare(pair_records["candidate"][0], pair_records["control"][0], 2e-12, 1e-12)
                pair_records["candidate"][1]["controlComparison"] = comparison
                save(pair_records["candidate"][2].with_suffix(".json"), pair_records["candidate"][1])
                save(args.output / "runs.json", records)
                if not comparison["passed"]:
                    raise RuntimeError("Candidate/control numerical gate failed; evidence retained")
            # Keep one complete pair for independent reproduction. Later payloads
            # have hashes and both numerical comparisons retained before removal.
            if pair:
                for binary, _, _ in pair_records.values():
                    binary.unlink()
            if args.mode == "qualification":
                summary = {"decision": "correctness-only; concurrent host work, no timing claim",
                           "completedRuns": len(records), "allMatlabComparisonsPassed": True,
                           "allAvailableControlComparisonsPassed": all(r.get("controlComparison", {}).get("passed", True) for r in records)}
            else:
                summary = summarize(records)
                if args.mode == "screen":
                    summary.pop("fluxPairedBootstrap95", None)
                    summary.pop("peakRSSPairedBootstrap95", None)
                    summary["decision"] = "worker screening only; no inferential adoption claim"
            save(args.output / "summary.json", summary)


if __name__ == "__main__":
    main()
