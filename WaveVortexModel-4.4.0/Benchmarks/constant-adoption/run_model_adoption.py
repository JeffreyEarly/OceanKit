#!/usr/bin/env python3
"""Compare frozen complete-model runs with retained outputs and exit-time RSS.

Author-only dependencies: NumPy and netCDF4. Qualification runs one fresh pair
per profile and makes no timing claim. Final mode requires frozen provenance,
uses two excluded warmup pairs and eight alternating measured pairs. Every run
gets a fresh source copy. Requests, reports, logs and comparisons remain; the
explicit large-only retention policy may remove verified extra NetCDF payloads
while retaining the first measured pair and all failures.
"""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

import netCDF4
import numpy as np


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def equal(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape:
        return False
    if a.dtype.kind in "fc" and b.dtype.kind in "fc":
        return bool(np.array_equal(a, b, equal_nan=True))
    return bool(np.array_equal(a, b))


def slices(shape, itemsize):
    if not shape:
        yield ()
        return
    if 0 in shape:
        return
    budget = max(1, (8 * 1024 * 1024) // max(1, itemsize))
    counts = [1] * len(shape)
    for axis in reversed(range(len(shape))):
        counts[axis] = min(shape[axis], budget)
        budget = max(1, budget // counts[axis])
    for starts in itertools.product(*(range(0, n, c) for n, c in zip(shape, counts))):
        yield tuple(slice(s, min(s + c, n)) for s, c, n in zip(starts, counts, shape))


def compare_graph(expected_path, actual_path, relative, absolute, history_policy="strict"):
    if history_policy not in ("strict", "matlab-writer-provenance"):
        raise ValueError("Unknown history comparison policy")
    result = {"passed": True, "relativeTolerance": relative, "absoluteTolerance": absolute,
              "toleranceRule": "abs(actual-expected) <= absolute + relative*abs(expected)",
              "metadataRule": "exact names, types, dimensions and values; attribute/group/variable declaration order is immaterial",
              "volatileAttributes": ["history (presence and nonempty value still required unless explicitly excluded below)"],
              "historyPolicy": history_policy, "excludedAttributes": [],
              "differences": [], "variables": [], "groups": []}

    def differ(path, reason):
        result["differences"].append({"path": path, "reason": reason})

    def attributes(expected, actual, path):
        left, right = set(expected.ncattrs()), set(actual.ncattrs())
        # Existing MATLAB output-close provenance has no corresponding C++
        # writer attribute. Exclude only that class's history at this boundary;
        # root/variable history and all C++-to-C++ attribute presence stay strict.
        output_group_history = (history_policy == "matlab-writer-provenance" and
            isinstance(expected, netCDF4.Group) and path != "/" and
            "AnnotatedClass" in left and
            expected.getncattr("AnnotatedClass") == "WVModelOutputGroupEvenlySpaced")
        if output_group_history and "history" in left | right:
            result["excludedAttributes"].append({"path": path + "/@history",
                "expectedPresent": "history" in left, "actualPresent": "history" in right,
                "reason": "preexisting MATLAB-only per-output-group completion provenance"})
            left.discard("history")
            right.discard("history")
        if left != right:
            differ(path, "attribute names differ")
        for name in sorted(left & right):
            a, b = expected.getncattr(name), actual.getncattr(name)
            if name == "history":
                if not str(a).strip() or not str(b).strip():
                    differ(path + "/@history", "volatile history must remain nonempty")
            elif np.asarray(a).dtype != np.asarray(b).dtype or not equal(a, b):
                differ(path + "/@" + name, "attribute value differs")

    def variable(expected, actual, path):
        attributes(expected, actual, path)
        if (expected.dimensions != actual.dimensions or expected.shape != actual.shape or
                str(expected.datatype) != str(actual.datatype)):
            differ(path, "variable type, shape or ordered dimensions differ")
            return
        row = {"path": path, "shape": list(expected.shape), "dimensions": list(expected.dimensions),
               "datatype": str(expected.datatype), "passed": True, "maximumAbsoluteError": 0.0,
               "referenceMaximum": 0.0, "maximumToleranceRatio": 0.0,
               "nonfiniteExpectedCount": 0, "nonfiniteActualCount": 0}
        exact = expected.name == "t" or np.dtype(expected.dtype).kind not in "fc"
        row["comparison"] = "exact" if exact else "pointwise-tolerance"
        for part in slices(expected.shape, np.dtype(expected.dtype).itemsize):
            a, b = np.asarray(expected[part]), np.asarray(actual[part])
            if exact:
                row["passed"] = row["passed"] and equal(a, b)
                continue
            finite_a, finite_b = np.isfinite(a), np.isfinite(b)
            row["nonfiniteExpectedCount"] += int(np.count_nonzero(~finite_a))
            row["nonfiniteActualCount"] += int(np.count_nonzero(~finite_b))
            matching_nonfinite = ((np.isnan(a) & np.isnan(b)) |
                                  (np.isinf(a) & np.isinf(b) & (a == b)))
            finite = finite_a & finite_b
            if not np.all(finite | matching_nonfinite):
                row["passed"] = False
            # Saved evolving diagnostics must be finite, even if both runs
            # contain the same invalid value. Static provenance may have NaNs.
            if "t" in expected.dimensions and not np.all(finite):
                row["passed"] = False
            av, bv = a[finite].astype(np.complex128 if a.dtype.kind == "c" else np.float64), b[finite]
            error, scale = np.abs(bv - av), np.abs(av)
            threshold = absolute + relative * scale
            row["maximumAbsoluteError"] = max(row["maximumAbsoluteError"], float(error.max(initial=0)))
            row["referenceMaximum"] = max(row["referenceMaximum"], float(scale.max(initial=0)))
            row["maximumToleranceRatio"] = max(row["maximumToleranceRatio"], float((error / threshold).max(initial=0)))
            row["passed"] = row["passed"] and bool(np.all(error <= threshold))
        normalized = row["maximumAbsoluteError"] / max(row["referenceMaximum"], np.finfo(float).tiny)
        row["maximumScaleNormalizedError"] = normalized if np.isfinite(normalized) else None
        if not row["passed"]:
            differ(path, "variable payload differs or evolving output is nonfinite")
        result["variables"].append(row)

    def group(expected, actual):
        path = expected.path
        attributes(expected, actual, path)
        dimensions = {k: {"length": len(v), "unlimited": v.isunlimited()} for k, v in expected.dimensions.items()}
        actual_dimensions = {k: {"length": len(v), "unlimited": v.isunlimited()} for k, v in actual.dimensions.items()}
        if dimensions != actual_dimensions:
            differ(path, "local dimension identities, lengths or unlimited flags differ")
        result["groups"].append({"path": path, "dimensions": dimensions})
        for collection, label in (("variables", "variable"), ("groups", "group")):
            left, right = getattr(expected, collection), getattr(actual, collection)
            if set(left) != set(right):
                differ(path, label + " names differ: expected=" + repr(sorted(left)) + ", actual=" + repr(sorted(right)))
            for name in sorted(set(left) & set(right)):
                if collection == "groups":
                    group(left[name], right[name])
                else:
                    variable(left[name], right[name], path.rstrip("/") + "/" + name)

    with netCDF4.Dataset(expected_path) as expected, netCDF4.Dataset(actual_path) as actual:
        for dataset in (expected, actual):
            dataset.set_auto_maskandscale(False)
            dataset.set_auto_chartostring(False)
        if expected.data_model != actual.data_model:
            differ("/", "NetCDF data model differs")
        group(expected, actual)
    result["passed"] = not result["differences"]
    result["variableCount"] = len(result["variables"])
    return result


def run_child(command, folder):
    started = time.monotonic()
    with (folder / "stdout.log").open("wb") as stdout, (folder / "stderr.log").open("wb") as stderr:
        child = subprocess.Popen(command, stdout=stdout, stderr=stderr)
        # wait4 owns reaping: poll/wait/communicate would lose this child's
        # complete-lifetime resource usage, including close and destruction.
        while True:
            try:
                _, status, usage = os.wait4(child.pid, 0)
                break
            except InterruptedError:
                continue
        child.returncode = os.waitstatus_to_exitcode(status)
    return {"command": command, "exitCode": child.returncode,
            "completeLifetimeSeconds": time.monotonic() - started,
            "completeLifetimePeakRSSBytes": int(usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024)),
            "userSeconds": usage.ru_utime, "systemSeconds": usage.ru_stime}


def summarize(pairs, mode, profile_count=4):
    measured = [p for p in pairs if not p["warmup"]]
    expected_pairs = profile_count * (1 if mode == "qualification" else 10)
    complete = len(pairs) == expected_pairs
    result = {"passed": complete and all(p["passed"] for p in pairs), "completedPairs": len(pairs),
              "expectedPairs": expected_pairs, "status": "complete" if complete else "incomplete",
              "decision": "correctness-only; concurrent host work, no timing claim"}
    if mode == "qualification" or not result["passed"]:
        return result
    metrics = {"integration": lambda r: r["timingSeconds"]["integrate"],
               "completeLifetime": lambda r: r["completeLifetimeSeconds"],
               "peakRSS": lambda r: r["completeLifetimePeakRSSBytes"],
               "retained": lambda r: r["livenessBytes"]["fullModelRetained"],
               "maximumLive": lambda r: r["livenessBytes"]["fullModelMaximumLive"]}
    result["profiles"] = []
    all_ratios = {name: [] for name in metrics}
    for profile in sorted({p["profile"] for p in measured}):
        current = [p for p in measured if p["profile"] == profile]
        row = {"profile": profile, "pairs": len(current)}
        for name, value in metrics.items():
            ratios = np.array([value(p["runs"]["candidate"]) / value(p["runs"]["control"]) for p in current])
            row[name + "GeometricRatio"] = float(np.exp(np.log(ratios).mean()))
            all_ratios[name].append(np.log(ratios))
        result["profiles"].append(row)
    result["pairedBootstrap"] = {}
    for name, profiles in all_ratios.items():
        random = np.random.default_rng(358)
        bootstrap = np.zeros(10000)
        for values in profiles:
            indices = random.integers(0, len(values), size=(10000, len(values)))
            bootstrap += values[indices].mean(axis=1) / len(profiles)
        result["pairedBootstrap"][name] = {"geometricRatio": float(np.exp(np.mean([r.mean() for r in profiles]))),
            "confidence95": np.exp(np.quantile(bootstrap, [.025, .975])).tolist(), "resamples": 10000, "seed": 358}
    result["decision"] = "measured complete-model evidence; adoption also requires standalone and MATLAB-loaded gates"
    return result


def fixture_contract(fixtures):
    """Named immutable workloads; a different run length requires a new profile."""
    schema = fixtures["schema"]
    base = ["constant-hydrostatic-coefficient-only", "constant-hydrostatic-composite",
            "constant-nonhydrostatic-coefficient-only", "constant-nonhydrostatic-composite"]
    if schema == "wvm-constant-model-fixtures-v1":
        expected = base
        steps, final = 64, 32
    elif schema == "wvm-constant-model-fixtures-large-v1":
        expected = [name + "-" + str(size) for size in (256, 512) for name in
                    (base[0], base[2], base[3])]
        steps, final = 4, 2
        parameters = fixtures["parameters"]
        if (parameters["Nxyz"] != [[256, 256, 129], [512, 512, 257]] or
                parameters["fields"] != ["u"] or parameters["denseOutputTimes"] != .75 or
                parameters["fixedStepCount"] != steps or parameters["rhsEvaluationCount"] != 4 * steps):
            raise ValueError("Large fixture axes, fields, dense time or work differ")
    else:
        raise ValueError("Unknown frozen complete-model profile schema")
    if [p["id"] for p in fixtures["profiles"]] != expected:
        raise ValueError("Frozen complete-model profile inventory/order differs")
    return {"method": "fixed-rk4", "finalTime": final, "initialStep": .5}, steps


def retain_verified_pair(pair, output, policy):
    """Delete only verified extra payloads after hashes/comparisons are durable."""
    if policy == "all" or not pair["passed"] or pair["pair"] == 0:
        return
    # The caller has persisted both run receipts, comparisons and pairs.json.
    # All failed pairs and the first measured pair stay independently inspectable.
    paths = []
    for role in ("control", "candidate"):
        run = pair["runs"][role]
        path = Path(run["directory"]) / "output.nc"
        expected = run["artifacts"]["output.nc"]
        if digest(path) != expected:
            raise ValueError("Verified output changed before retention cleanup")
        paths.append({"path": str(path), "sha256": expected, "bytes": path.stat().st_size})
    receipt = {"policy": policy, "profile": pair["profile"], "pair": pair["pair"],
               "reason": "Both complete graphs, numerical comparisons, work and source identities passed; hashes and comparisons retained",
               "status": "verified-pending-removal", "outputs": paths}
    receipt_path = output / (pair["profile"] + f"-pair{pair['pair']:02}-retention.json")
    save(receipt_path, receipt)
    for item in paths:
        Path(item["path"]).unlink()
    receipt["status"] = "removed-after-verification"
    save(receipt_path, receipt)
    pair["payloadRetention"] = {"path": str(receipt_path), "sha256": digest(receipt_path),
                                "status": receipt["status"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixtures", required=True, type=Path, help="Frozen manifest.json")
    parser.add_argument("--control", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--mode", choices=["qualification", "final"], default="qualification")
    parser.add_argument("--payload-retention", choices=["all", "first-measured-pair"], default="all",
                        help="Large profiles only: retain first measured pair and every failed payload; remove only verified extra outputs")
    parser.add_argument("--provenance", type=Path, help="Frozen source/build/provider/host receipt; required for final")
    args = parser.parse_args()
    if args.mode == "final" and not args.provenance:
        parser.error("Final measurement requires frozen source/build/provider/host provenance")
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    fixtures = json.loads(args.fixtures.read_text())
    expected_integration, expected_steps = fixture_contract(fixtures)
    if args.payload_retention != "all" and fixtures["schema"] != "wvm-constant-model-fixtures-large-v1":
        raise ValueError("Historical small profiles retain every output")
    executables = {"control": args.control.resolve(), "candidate": args.candidate.resolve()}
    executable_hashes = {role: digest(path) for role, path in executables.items()}
    for profile in fixtures["profiles"]:
        for kind in ("source", "reference", "request"):
            if digest(profile[kind + "Path"]) != profile[kind + "SHA256"]:
                raise ValueError("Frozen " + kind + " fixture hash differs: " + profile["id"])
        request = json.loads(Path(profile["requestPath"]).read_text())
        if request["integration"] != expected_integration:
            raise ValueError("Fixed integration request differs from the named frozen workload")
    save(args.output / "manifest.json", {"schema": "wvm-constant-model-adoption-v1", "mode": args.mode,
        "platform": platform.platform(), "harnessSHA256": digest(__file__), "fixtures": fixtures,
        "fixturesManifestSHA256": digest(args.fixtures), "warmupPairs": 0 if args.mode == "qualification" else 2,
        "measuredPairs": 1 if args.mode == "qualification" else 8,
        "payloadRetention": args.payload_retention, "expectedStepCount": expected_steps,
        "expectedRHSEvaluationCount": 4 * expected_steps,
        "executables": {role: {"path": str(path), "sha256": executable_hashes[role]} for role, path in executables.items()},
        "provenance": json.loads(args.provenance.read_text()) if args.provenance else None,
        "provenanceSHA256": digest(args.provenance) if args.provenance else None,
        "dependencies": {"numpy": np.__version__, "netCDF4": netCDF4.__version__}})
    pairs = []
    for profile in fixtures["profiles"]:
        repetitions = range(1) if args.mode == "qualification" else range(-2, 8)
        for pair_index in repetitions:
            pair = {"profile": profile["id"], "pair": pair_index, "warmup": pair_index < 0, "runs": {}, "passed": True}
            roles = ["control", "candidate"] if pair_index % 2 == 0 else ["candidate", "control"]
            for role in roles:
                folder = args.output / (profile["id"] + f"-pair{pair_index:02}-{role}")
                folder.mkdir()
                output, request_path, report_path = folder / "output.nc", folder / "request.json", folder / "report.json"
                if digest(profile["sourcePath"]) != profile["sourceSHA256"] or digest(executables[role]) != executable_hashes[role]:
                    raise ValueError("Source or executable changed during campaign")
                shutil.copyfile(profile["sourcePath"], output)
                request = json.loads(Path(profile["requestPath"]).read_text())
                request.update(modelFiles=[str(output)], report=str(report_path))
                save(request_path, request)
                print(profile["id"], "pair", pair_index, role, flush=True)
                run = run_child([str(executables[role]), "--request", str(request_path)], folder)
                run.update(profile=profile["id"], pair=pair_index, role=role, directory=str(folder), passed=False)
                run["executableUnchanged"] = digest(executables[role]) == executable_hashes[role]
                if run["exitCode"] == 0 and report_path.is_file():
                    report = json.loads(report_path.read_text())
                    run.update(timingSeconds=report["timingSeconds"], livenessBytes=report["livenessBytes"], state=report["state"])
                    for name in ("fullModelRetained", "fullModelMaximumLive"):
                        if name not in run["livenessBytes"]:
                            raise ValueError("Missing complete-model liveness metric " + name)
                    comparison = compare_graph(profile["referencePath"], output, 1e-10, 1e-12,
                                               history_policy="matlab-writer-provenance")
                    save(folder / "matlab-comparison.json", comparison)
                    run["matlabComparison"] = {"passed": comparison["passed"], "differences": comparison["differences"],
                        "excludedAttributes": comparison["excludedAttributes"], "variableCount": comparison["variableCount"]}
                    run["fixedWorkPassed"] = (run["state"]["stepCount"] == expected_steps and
                        run["state"]["rhsEvaluationCount"] == 4 * expected_steps and
                        run["state"]["initialTime"] == 0 and run["state"]["finalTime"] == expected_integration["finalTime"] and
                        run["state"]["deltaT"] == .5 and run["state"]["rejectedStepCount"] == 0)
                    run["passed"] = comparison["passed"] and run["fixedWorkPassed"] and run["executableUnchanged"]
                run["artifacts"] = {p.name: digest(p) for p in sorted(folder.iterdir()) if p.is_file()}
                save(folder / "run.json", run)
                pair["runs"][role] = run
                pair["passed"] = pair["passed"] and run["passed"]
                save(args.output / "in-progress.json", {"pairs": pairs, "currentPair": pair})
            control, candidate = pair["runs"]["control"], pair["runs"]["candidate"]
            if control["exitCode"] == candidate["exitCode"] == 0 and "state" in control and "state" in candidate:
                comparison = compare_graph(Path(control["directory"]) / "output.nc", Path(candidate["directory"]) / "output.nc", 2e-12, 1e-13)
                comparison_path = args.output / (profile["id"] + f"-pair{pair_index:02}-control-comparison.json")
                save(comparison_path, comparison)
                pair["controlComparison"] = {"path": str(comparison_path), "sha256": digest(comparison_path), "passed": comparison["passed"], "differences": comparison["differences"]}
                pair["matchedWork"] = all(control["state"][key] == candidate["state"][key] for key in ("stepCount", "rhsEvaluationCount"))
                pair["passed"] = pair["passed"] and comparison["passed"] and pair["matchedWork"]
            pairs.append(pair)
            save(args.output / "pairs.json", pairs)
            save(args.output / "summary.json", summarize(pairs, args.mode, len(fixtures["profiles"])))
            retain_verified_pair(pair, args.output, args.payload_retention)
            save(args.output / "pairs.json", pairs)
            if not pair["passed"] and args.mode == "final":
                raise RuntimeError("Qualification failed; negative pair and all outputs retained")
    for profile in fixtures["profiles"]:
        for kind in ("source", "reference", "request"):
            if digest(profile[kind + "Path"]) != profile[kind + "SHA256"]:
                raise ValueError("Frozen inputs changed during campaign")
    return 0 if all(pair["passed"] for pair in pairs) else 1


if __name__ == "__main__":
    raise SystemExit(main())
