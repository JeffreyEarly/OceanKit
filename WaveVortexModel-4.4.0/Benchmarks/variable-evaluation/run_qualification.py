#!/usr/bin/env python3
"""Qualify frozen model executables with paired, full-model continuations.

Run on an idle host. The manifest supplies immutable source files and v2
requests; every execution receives its own copy. No experiment file is opened
for writing. Timing is meaningful only after the numerical comparisons pass.
"""

import argparse
import copy
import datetime
import importlib.util
import json
import math
import os
import platform
import pathlib
import shutil
import statistics
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=pathlib.Path, required=True)
    parser.add_argument("--candidate", type=pathlib.Path, required=True)
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--smoke", action="store_true",
                        help="One correctness pair, explicitly not performance qualification")
    parser.add_argument("--maximum-low-memory-ratio", type=float, default=1.03,
                        help="Explicit owned-memory growth budget; recorded in the protocol")
    parser.add_argument("--policies", nargs="+", choices=("reuse", "low-memory"),
                        default=["reuse", "low-memory"],
                        help="Candidate policies to qualify; reuse is required")
    args = parser.parse_args()
    if "reuse" not in args.policies or len(set(args.policies)) != len(args.policies):
        parser.error("Candidate policies must be unique and include reuse")
    if not math.isfinite(args.maximum_low_memory_ratio) or args.maximum_low_memory_ratio < 1:
        parser.error("The owned-memory ratio budget must be finite and at least one")
    root = pathlib.Path(__file__).resolve().parents[2]
    helper_path = root / "Benchmarks/constant-adoption/run_model_adoption.py"
    spec = importlib.util.spec_from_file_location("model_adoption", helper_path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    manifest = json.loads(args.manifest.read_text())
    profiles = manifest["profiles"]
    source_commit = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    source_paths = subprocess.check_output(
        ["git", "ls-files", "CompiledKernel", "PortableRuntime"],
        cwd=root, text=True).splitlines()
    source_hashes = {path: helper.digest(root / path) for path in source_paths}
    source_dirty = subprocess.check_output(
        ["git", "status", "--porcelain", "--", "CompiledKernel", "PortableRuntime"],
        cwd=root, text=True).strip()
    if not args.smoke:
        assert not source_dirty, "Commit the candidate runtime source before qualification"

    binaries = {"baseline": args.baseline.resolve(),
                **{policy: args.candidate.resolve() for policy in args.policies}}
    frozen = {str(path): helper.digest(path) for path in
              [*binaries.values(), args.manifest, pathlib.Path(__file__), helper_path]}
    for profile in profiles:
        source = pathlib.Path(profile["sourcePath"])
        assert helper.digest(source) == profile["sourceSHA256"], source
        frozen[str(source)] = profile["sourceSHA256"]
    args.output.mkdir(parents=True, exist_ok=False)
    measured_count = 1 if args.smoke else 8
    warmups = 0 if args.smoke else 2
    helper.save(args.output / "protocol.json", {
        "kind": "correctness-smoke" if args.smoke else "performance-qualification",
        "startedAtUTC": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "warmupPairs": warmups, "measuredPairs": measured_count,
        "order": "/".join(binaries)+", reversed on alternate pairs",
        "candidatePolicies": args.policies,
        "numericalGate": "exact metadata and integration decisions; manifest scientific tolerances (zero by default)",
        "runtimeGate": "EddyTide improvement; investigate each other reuse ratio above 1.03",
        "storageGate": "low-memory maximum-live owned bytes no more than baseline times the explicit ratio budget; retained capacity reported separately",
        "maximumLowMemoryRatio": args.maximum_low_memory_ratio,
        "profiles": profiles, "frozenSHA256": frozen,
        "host": {"platform": platform.platform(), "machine": platform.machine(),
                 "logicalCPUCount": os.cpu_count(),
                 "threadEnvironment": {name: os.environ.get(name) for name in
                     ("OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OPENBLAS_NUM_THREADS")}},
        "candidateSourceCommit": source_commit,
        "candidateSourceDirty": bool(source_dirty),
        "candidateSourceSHA256": source_hashes,
    })
    pairs = []
    try:
        for profile in profiles:
            for index in range(-warmups, measured_count):
                pair = {"profile": profile["id"], "pair": index,
                        "warmup": index < 0, "runs": {}}
                order = list(binaries) if index % 2 == 0 else list(reversed(binaries))
                for role in order:
                    directory = args.output / f"{profile['id']}-{index}-{role}"
                    directory.mkdir()
                    output = directory / "output.nc"
                    shutil.copyfile(profile["sourcePath"], output)
                    request = copy.deepcopy(profile["request"])
                    request["modelFiles"] = [str(output.resolve())]
                    request["report"] = str((directory / "report.json").resolve())
                    request.setdefault("execution", {}).pop("variableEvaluationPolicy", None)
                    if role != "baseline":
                        request["execution"]["variableEvaluationPolicy"] = role
                    request_path = directory / "request.json"
                    helper.save(request_path, request)
                    assert helper.digest(binaries[role]) == frozen[str(binaries[role])]
                    run = helper.run_child([str(binaries[role]), "--request",
                                            str(request_path.resolve())], directory)
                    assert run["exitCode"] == 0, (directory, run)
                    report = json.loads((directory / "report.json").read_text())
                    assert report["status"] == "complete", directory
                    assert report["provider"]["id"] == "native-fftw", directory
                    if role != "baseline":
                        if not args.smoke:
                            assert report["source"]["commit"] == source_commit, (
                                "Reconfigure and rebuild the candidate at the frozen source commit", directory)
                        evaluation = report["variableEvaluation"]
                        assert evaluation["effectivePolicy"] == role
                        if role == "reuse":
                            for scope in ("rightHandSide", "output"):
                                assert evaluation[scope]["duplicateExecutions"] == 0
                            producers = evaluation["kernelProducers"]
                            # Per-context tests provide exact key attribution. This
                            # independent aggregate bound also catches production
                            # work escaping the ledger in representative runs.
                            assert producers["phasePreparations"] <= producers["stateValidations"]
                            for producer in producers["reconstructions"]:
                                assert producer["count"] <= producers["stateValidations"], (
                                    "Actual reconstruction count exceeds immutable-state evaluations",
                                    directory, producer, producers["stateValidations"])
                    run.update(directory=str(directory.resolve()),
                               timingSeconds=report["timingSeconds"],
                               livenessBytes=report["livenessBytes"],
                               storageBytes=report["storageBytes"],
                               state=report["state"],
                               integrationRequest=report["integrationRequest"],
                               acceptedSteps=report["integrator"].get("acceptedSteps", []),
                               variableEvaluation=report.get("variableEvaluation"),
                               source=report["source"],
                               outputSHA256=helper.digest(output))
                    pair["runs"][role] = run
                    helper.save(directory / "run.json", run)
                    print(profile["id"], index, role,
                          report["timingSeconds"]["integrate"], flush=True)
                baseline = pair["runs"]["baseline"]
                for role in args.policies:
                    candidate = pair["runs"][role]
                    assert candidate["state"] == baseline["state"], (profile["id"], role, "state")
                    # Exclude the wall-clock duration of selecting controls.
                    for value in (candidate, baseline):
                        value["integrationRequest"].pop("candidateEvaluationSeconds", None)
                    assert candidate["integrationRequest"] == baseline["integrationRequest"]
                    baseline_steps=copy.deepcopy(baseline["acceptedSteps"])
                    candidate_steps=copy.deepcopy(candidate["acceptedSteps"])
                    error_differences=[]
                    for expected,actual in zip(baseline_steps,candidate_steps):
                        expected_error=expected.pop("normalizedError",0.0)
                        actual_error=actual.pop("normalizedError",0.0)
                        assert math.isfinite(expected_error) and math.isfinite(actual_error)
                        error_differences.append(abs(actual_error-expected_error))
                    assert candidate_steps == baseline_steps
                    candidate["maximumNormalizedErrorDiagnosticDifference"]=max(error_differences,default=0.0)
                    tolerance=profile.get("comparison",{})
                    comparison = helper.compare_graph(
                        pathlib.Path(baseline["directory"]) / "output.nc",
                        pathlib.Path(candidate["directory"]) / "output.nc",
                        tolerance.get("relativeTolerance",0),tolerance.get("absoluteTolerance",0))
                    comparison["bitwiseEqual"]=comparison["passed"] and all(
                        variable.get("maximumAbsoluteError",0)==0 for variable in comparison["variables"])
                    helper.save(args.output / f"{profile['id']}-{index}-{role}-comparison.json", comparison)
                    assert comparison["passed"], comparison["differences"]
                pairs.append(pair)
                helper.save(args.output / "pairs.json", pairs)
                if index != 0:
                    for run in pair["runs"].values():
                        (pathlib.Path(run["directory"]) / "output.nc").unlink()
        summaries = []
        for profile in profiles:
            measured = [pair for pair in pairs if pair["profile"] == profile["id"] and not pair["warmup"]]
            row = {"profile": profile["id"], "policies": {}}
            for role in args.policies:
                ratios = [pair["runs"][role]["timingSeconds"]["integrate"] /
                          pair["runs"]["baseline"]["timingSeconds"]["integrate"] for pair in measured]
                logs = helper.np.log(ratios)
                rng = helper.np.random.default_rng(470)
                bootstrap = logs[rng.integers(0, len(logs), (10000, len(logs)))].mean(axis=1)
                retained = [pair["runs"][role]["livenessBytes"]["fullModelRetained"] /
                            pair["runs"]["baseline"]["livenessBytes"]["fullModelRetained"] for pair in measured]
                lifetime_ratios = [pair["runs"][role]["completeLifetimeSeconds"] /
                                   pair["runs"]["baseline"]["completeLifetimeSeconds"]
                                   for pair in measured]
                lifetime_logs = helper.np.log(lifetime_ratios)
                lifetime_bootstrap = lifetime_logs[
                    rng.integers(0, len(lifetime_logs), (10000, len(lifetime_logs)))].mean(axis=1)
                row["policies"][role] = {
                    "completeLifetimeRatio": math.exp(statistics.mean(map(math.log, lifetime_ratios))),
                    "completeLifetimeBootstrap95": helper.np.exp(helper.np.quantile(lifetime_bootstrap, [.025, .975])).tolist(),
                    "integrationRatio": math.exp(statistics.mean(map(math.log, ratios))),
                    "pairedBootstrap95": helper.np.exp(helper.np.quantile(bootstrap, [.025, .975])).tolist(),
                    "ratios": ratios, "maximumRetainedRatio": max(retained),
                    "maximumLiveRatio": max(pair["runs"][role]["livenessBytes"]["fullModelMaximumLive"] /
                        pair["runs"]["baseline"]["livenessBytes"]["fullModelMaximumLive"] for pair in measured),
                    "maximumPeakRSSBytes": max(pair["runs"][role]["completeLifetimePeakRSSBytes"] for pair in measured),
                }
            row["runtimePassed"] = row["policies"]["reuse"]["integrationRatio"] < (1 if profile.get("requiresImprovement") else 1.03)
            row["requiresRuntimeInvestigation"] = (
                row["policies"]["reuse"]["integrationRatio"] > 1.03 or
                row["policies"]["reuse"]["completeLifetimeRatio"] > 1.03)
            row["lowMemoryPassed"] = (row["policies"]["low-memory"]["maximumLiveRatio"] <= args.maximum_low_memory_ratio
                                      if "low-memory" in row["policies"] else None)
            summaries.append(row)
        assert all(helper.digest(path) == digest for path, digest in frozen.items())
        assert all(helper.digest(root / path) == digest
                   for path, digest in source_hashes.items())
        helper.save(args.output / "summary.json", {
            "qualified": not args.smoke and all(row["runtimePassed"] and row["lowMemoryPassed"] is not False and not row["requiresRuntimeInvestigation"] for row in summaries),
            "finishedAtUTC": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "allNumericalComparisonsPassed": True, "allIntegrationDecisionsIdentical": True,
            "postflightUnchanged": True,
            "profiles": summaries,
        })
    except Exception as error:
        helper.save(args.output / "failure.json", {"status": "incomplete", "cause": repr(error)})
        raise


if __name__ == "__main__":
    main()
