#!/usr/bin/env python3
"""Qualify frozen, optimized interleaved, and compact variable WVModel runs.

Each run starts in a fresh process from an immutable copied source file. The
old-source worker is an independent frozen control; the candidate-source worker
is exercised with frozen, optimized interleaved, and compact split policies.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import sys


HERE = Path(__file__).resolve().parent
EXPECTED_CONTROL_COMMIT = "b157988355d2fabb6cbc62d56de18b203456897f"
EXPECTED_CANDIDATE_COMMIT = "a3acbf4dc9305c4d656399e9f98bc8360ddfa66f"
CONSTANT_RUNNER = HERE.parent.parent / "constant-adoption" / "run_model_adoption.py"
SPEC = importlib.util.spec_from_file_location("wvm_constant_model_adoption", CONSTANT_RUNNER)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("Unable to load the constant-adoption comparison machinery")
CONSTANT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CONSTANT)


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def source_snapshot(path, expected_commit):
    root = Path(path).resolve()
    commit = subprocess.check_output(
        ["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    worktree = subprocess.run(
        ["git", "-C", str(root), "diff", "--quiet", "HEAD", "--",
         "CompiledKernel", "PortableRuntime"], check=False).returncode
    index = subprocess.run(
        ["git", "-C", str(root), "diff", "--cached", "--quiet", "HEAD", "--",
         "CompiledKernel", "PortableRuntime"], check=False).returncode
    result = {"root": str(root), "commit": commit,
              "compiledSourceDirty": worktree != 0 or index != 0}
    if commit != expected_commit or result["compiledSourceDirty"]:
        raise ValueError("Source boundary is not the expected clean commit: " +
                         json.dumps(result))
    return result


def provider_snapshot(base_library, thread_library):
    result = {}
    for key, path in (("baseLibrary", base_library),
                      ("threadLibrary", thread_library)):
        resolved = Path(path).resolve()
        if not resolved.is_file():
            raise FileNotFoundError(resolved)
        result[key] = str(resolved)
        result[key + "SHA256"] = digest(resolved)
    return result


def fixture_contract(manifest):
    schema = manifest.get("schema")
    if schema == "wvm-variable-model-fixtures-small-v1":
        profile_set, grid, final_time, step, blocks, warmups = (
            "small", [32, 24, 33], 2.0, 0.5, 1, 0)
    elif schema == "wvm-variable-model-fixtures-large-v1":
        profile_set, grid, final_time, step, blocks, warmups = (
            "large", [256, 256, 129], 0.5, 0.5, 10, 2)
    elif schema == "wvm-variable-boussinesq-model-fixture-256-v1":
        profile_set, grid, final_time, step, blocks, warmups = (
            "boussinesq-256", [256, 256, 129], 83.5, 0.5, 10, 2)
    else:
        raise ValueError("Unknown variable full-model fixture schema")
    parameters = manifest["parameters"]
    if manifest.get("sourceHEAD") != EXPECTED_CANDIDATE_COMMIT:
        raise ValueError("Fixture numerical source commit differs from the frozen candidate")
    if (parameters["grid"] != grid or parameters["finalTime"] != final_time or
            parameters["initialStep"] != step or
            parameters["forcing"] != "WVNonlinearAdvection"):
        raise ValueError("Fixture axes, time span, step, or forcing differ")
    if profile_set == "boussinesq-256":
        expected_ids = ["variable-boussinesq-composite-256"]
    else:
        expected_ids = [f"variable-{family}-composite-{profile_set}" for family in
                        ("stratified-qg", "hydrostatic")]
    if [profile["id"] for profile in manifest["profiles"]] != expected_ids:
        raise ValueError("Fixture profile inventory or order differs")
    for profile in manifest["profiles"]:
        expected_steps = round((profile.get("finalTime", final_time) -
                                profile.get("initialTime", 0.0)) / step)
        if (profile["grid"] != grid or
                profile["expectedStepCount"] != expected_steps or
                profile["expectedRHSEvaluationCount"] != 4 * expected_steps):
            raise ValueError("Fixture work declaration differs")
    return {"profileSet": profile_set, "grid": grid, "finalTime": final_time,
            "step": step, "blocks": blocks, "warmupBlocks": warmups,
            "profileCount": len(expected_ids),
            "declaredProfileCount": len(expected_ids)}


def verify_fixture(profile):
    for key, hash_key in (("sourcePath", "sourceSHA256"),
                          ("referencePath", "referenceSHA256")):
        path = Path(profile[key])
        if not path.is_file():
            raise FileNotFoundError(path)
        if digest(path) != profile[hash_key]:
            raise ValueError(f"Frozen fixture hash differs: {path}")


def fixture_snapshot(manifest_path, profiles):
    return {
        "manifestPath": str(Path(manifest_path).resolve()),
        "manifestSHA256": digest(manifest_path),
        "files": [{"id": profile["id"],
                   "sourceSHA256": digest(profile["sourcePath"]),
                   "sourceBytes": Path(profile["sourcePath"]).stat().st_size,
                   "referenceSHA256": digest(profile["referencePath"]),
                   "referenceBytes": Path(profile["referencePath"]).stat().st_size}
                  for profile in profiles]}


def validate_report(report, role, selection, profile, horizontal_workers,
                    pointwise_workers, contract, expected_source_root,
                    expected_provider):
    failures = []

    def require(condition, message):
        if not condition:
            failures.append(message)

    source = report.get("source", {})
    policy = report.get("selection", {})
    work = report.get("work", {})
    workload = report.get("workload", {})
    require(report.get("status") == "complete", "worker did not complete")
    require(report.get("interface") == "standalone-cpp-wvmodel",
            "worker interface differs")
    require(source.get("compiledSourceDirty") is False,
            "compiled kernel/runtime source boundary was dirty at configure time")
    require(source.get("buildType") == "Release", "worker is not a Release build")
    require(source.get("workerSourceSHA256") == digest(HERE / "WVVariableModelAdoptionWorker.cpp"),
            "worker was built from a different harness source")
    require(Path(source.get("root", "/nonexistent")).resolve() == expected_source_root.resolve(),
            "worker CMake source root differs from the declared source boundary")
    require(workload.get("transform") == profile["family"],
            "restored transform differs")
    require(workload.get("finalTime") ==
            profile.get("finalTime", contract["finalTime"]) and
            workload.get("step") == contract["step"], "time work differs")
    require(policy.get("id") == selection, "selection differs")
    require(policy.get("fftwInternalThreads") == 1 and
            policy.get("generalVerticalWorkers") == 1,
            "fixed FFTW/general-vertical topology differs")
    require(work.get("integratorSteps") == profile["expectedStepCount"] and
            work.get("acceptedSteps") == profile["expectedStepCount"] and
            work.get("rejectedSteps") == 0, "fixed accepted-step work differs")
    require(work.get("rhsEvaluations") == profile["expectedRHSEvaluationCount"],
            "fixed RK4 RHS work differs")
    require(work.get("forcingEvaluations") == profile["expectedRHSEvaluationCount"],
            "forcing work differs")
    require(work.get("scalarAdvections", 0) > 0 and
            work.get("tracerEvaluations", 0) > 0,
            "active scalar tracer work was not observed")
    if profile["expectedParticleRecords"] > 0:
        require(work.get("particleVelocityEvaluations", 0) > 0,
                "particle observer work was not observed")
    else:
        require(work.get("particleVelocityEvaluations", 0) == 0,
                "SQG unexpectedly performed unsupported particle work")
    require(work.get("outputInterpolations", 0) > 0 and
            work.get("committedRecords", 0) > 0,
            "dense/output work was not observed")
    require(report.get("provider", {}).get("noFallback") is True,
            "provider fallback was not excluded")
    provider = report.get("provider", {})
    for key in ("baseLibrary", "threadLibrary"):
        observed = Path(provider.get(key, "/nonexistent")).resolve()
        require(observed == Path(expected_provider[key]).resolve(),
                f"observed {key} differs from the pinned provider")
        require(provider.get(key + "SHA256") == expected_provider[key + "SHA256"],
                f"observed {key} hash differs from the pinned provider")
    backend = policy.get("matrixBackendIdentifier", "").lower()
    if role == "control":
        require(source.get("commit") == EXPECTED_CONTROL_COMMIT,
                "control source commit differs from the frozen boundary")
        require(source.get("variableServicesAvailable") is False,
                "control is not an independent pre-services source build")
        require(selection == "frozen" and "scalar" in backend,
                "control did not report the frozen scalar boundary")
    else:
        require(source.get("commit") == EXPECTED_CANDIDATE_COMMIT,
                "candidate source commit differs from the frozen boundary")
        require(source.get("variableServicesAvailable") is True,
                "candidate source lacks variable service injection")
        require(policy.get("matrixBackendInstanceCount", 0) > 0,
                "candidate factory did not observe backend construction")
        expected_backend = "scalar" if selection == "frozen" else "accelerate"
        require(expected_backend in backend, "candidate matrix backend differs")
        expected_streamed = selection != "frozen"
        require(policy.get("streamedNonlinear") is expected_streamed,
                "streamed policy differs")
        require(policy.get("compactSplitViews") is (selection == "compact"),
                "compact policy differs")
        require(policy.get("horizontalWorkers") ==
                (horizontal_workers if expected_streamed else 1),
                "horizontal worker topology differs")
        require(policy.get("pointwiseWorkers") ==
                (pointwise_workers if expected_streamed else 1),
                "pointwise worker topology differs")
    return failures


def run_one(worker, role, selection, profile, contract, output, block,
            horizontal_workers, pointwise_workers, expected_source_root,
            expected_provider):
    run_folder = output / profile["id"] / f"block-{block:02d}" / f"{role}-{selection}"
    run_folder.mkdir(parents=True)
    output_path = run_folder / "output.nc"
    report_path = run_folder / "worker-report.json"
    shutil.copy2(profile["sourcePath"], output_path)
    command = [str(worker), str(output_path), selection,
               str(horizontal_workers if selection != "frozen" else 1),
               str(profile.get("finalTime", contract["finalTime"])),
               str(contract["step"]),
               str(report_path),
               str(pointwise_workers if selection != "frozen" else 1)]
    process = CONSTANT.run_child(command, run_folder)
    receipt = {"role": role, "selection": selection,
               "workerPath": str(worker), "workerSHA256": digest(worker),
               "process": process, "outputPath": str(output_path),
               "reportPath": str(report_path)}
    if process["exitCode"] == 0 and report_path.is_file():
        receipt["report"] = json.loads(report_path.read_text())
        for key in ("baseLibrary", "threadLibrary"):
            provider_path = receipt["report"].get("provider", {}).get(key)
            if provider_path and Path(provider_path).is_file():
                receipt["report"]["provider"][key + "SHA256"] = digest(provider_path)
        receipt["outputSHA256"] = digest(output_path)
        receipt["outputBytes"] = output_path.stat().st_size
        receipt["validationFailures"] = validate_report(
            receipt["report"], role, selection, profile,
            horizontal_workers, pointwise_workers, contract,
            expected_source_root, expected_provider)
    else:
        receipt["validationFailures"] = ["worker process failed or omitted its report"]
    receipt["passed"] = not receipt["validationFailures"]
    save(run_folder / "receipt.json", receipt)
    return receipt


def compare_block(profile, runs, block_folder):
    comparisons = {}
    for key, run in runs.items():
        if not run["passed"]:
            comparisons[f"matlab-vs-{key}"] = {"passed": False,
                                                "differences": ["worker validation failed"]}
            continue
        comparison = CONSTANT.compare_graph(
            profile["referencePath"], run["outputPath"], 1e-10, 1e-12,
            "matlab-writer-provenance")
        comparisons[f"matlab-vs-{key}"] = comparison
        save(block_folder / f"comparison-matlab-vs-{key}.json", comparison)
    control = runs.get("control-frozen")
    if control and control["passed"]:
        for selection in ("frozen", "interleaved", "compact"):
            key = f"candidate-{selection}"
            candidate = runs.get(key)
            if candidate and candidate["passed"]:
                comparison = CONSTANT.compare_graph(
                    control["outputPath"], candidate["outputPath"],
                    2e-12, 1e-13, "strict")
            else:
                comparison = {"passed": False,
                              "differences": ["worker validation failed"]}
            comparisons[f"control-vs-{key}"] = comparison
            save(block_folder / f"comparison-control-vs-{key}.json", comparison)
    return comparisons


def compare_incremental(profile, run, retained_control_path, block_folder, key):
    comparisons = {}
    if not run["passed"]:
        return {f"matlab-vs-{key}": {
            "passed": False, "differences": ["worker validation failed"]}}
    if retained_control_path is None or not Path(retained_control_path).is_file():
        return {f"retained-control-vs-{key}": {
            "passed": False,
            "differences": ["first-block independent control is unavailable"]}}
    matlab = CONSTANT.compare_graph(
        profile["referencePath"], run["outputPath"], 1e-10, 1e-12,
        "matlab-writer-provenance")
    comparisons[f"matlab-vs-{key}"] = matlab
    save(block_folder / f"comparison-matlab-vs-{key}.json", matlab)
    retained = CONSTANT.compare_graph(
        retained_control_path, run["outputPath"], 2e-12, 1e-13, "strict")
    comparisons[f"retained-control-vs-{key}"] = retained
    save(block_folder / f"comparison-retained-control-vs-{key}.json", retained)
    return comparisons


def capacity_preflight(profiles, contract, output, retention,
                       host_physical_memory_bytes):
    free = shutil.disk_usage(output.parent).free
    reference_sizes = [profile["referenceBytes"] for profile in profiles]
    if retention == "all":
        simultaneous_payload = 4 * contract["blocks"] * sum(reference_sizes)
    elif retention == "none":
        simultaneous_payload = 2 * max(reference_sizes)
    else:
        retained_first_blocks = 4 * sum(reference_sizes)
        largest_later_payload = max(reference_sizes) if contract["blocks"] > 1 else 0
        simultaneous_payload = retained_first_blocks + largest_later_payload
    if contract["profileSet"] == "large":
        margin = 2 * 1024 ** 3
    elif contract["profileSet"] == "boussinesq-256":
        margin = 512 * 1024 ** 2
    else:
        margin = 256 * 1024 ** 2
    required = simultaneous_payload + margin
    result = {"freeBytes": free, "requiredWorkingBytes": required,
              "simultaneousPayloadEstimateBytes": simultaneous_payload,
              "referenceBytesByProfile": reference_sizes,
              "retention": retention, "marginBytes": margin,
              "hostPhysicalMemoryBytes": host_physical_memory_bytes,
              "policy": "refuse; never substitute a smaller grid"}
    if free < required:
        raise RuntimeError("Insufficient capacity for the exact declared workload: " +
                           json.dumps(result))
    return result


def summarize(blocks, contract):
    passed = all(block["passed"] for block in blocks)
    expected_blocks = contract["profileCount"] * contract["blocks"]
    result = {"passed": passed, "profileSet": contract["profileSet"],
              "declaredProfileCount": contract["declaredProfileCount"],
              "selectedProfileCount": contract["profileCount"],
              "completedBlocks": len(blocks),
              "expectedBlocks": expected_blocks,
              "campaignComplete": len(blocks) == expected_blocks,
              "decision": "correctness qualification; no timing claim"}
    if contract.get("qualificationOnly") or not passed:
        return result
    result["timing"] = []
    measured_count = 0
    for profile in sorted({block["profile"] for block in blocks}):
        measured = [block for block in blocks if block["profile"] == profile and
                    not block["warmup"]]
        measured_count += len(measured)
        row = {"profile": profile, "measuredBlocks": len(measured), "ratios": {}}
        if not measured:
            result["timing"].append(row)
            continue
        for selection in ("frozen", "interleaved", "compact"):
            ratios = []
            for block in measured:
                baseline = block["runs"]["control-frozen"]["report"]["timingSeconds"]["integrate"]
                candidate = block["runs"][f"candidate-{selection}"]["report"]["timingSeconds"]["integrate"]
                ratios.append(candidate / baseline)
            row["ratios"][selection] = {
                "candidateOverIndependentControlMedian": statistics.median(ratios),
                "samples": ratios}
        result["timing"].append(row)
    if measured_count == 0:
        result["decision"] = "warmup progress; no timing claim"
    elif result["campaignComplete"]:
        result["decision"] = "complete-model measurements; adoption requires the other declared gates"
    else:
        result["decision"] = "partial measured progress; no adoption timing claim"
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--control-worker", required=True, type=Path)
    parser.add_argument("--candidate-worker", required=True, type=Path)
    parser.add_argument("--control-source", required=True, type=Path)
    parser.add_argument("--candidate-source", required=True, type=Path)
    parser.add_argument("--fftw-base-library", required=True, type=Path)
    parser.add_argument("--fftw-thread-library", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--horizontal-workers", required=True, type=int)
    parser.add_argument("--pointwise-workers", required=True, type=int)
    parser.add_argument("--host-physical-memory-bytes", required=True, type=int)
    parser.add_argument("--mode", choices=("qualification", "final"),
                        default="qualification")
    parser.add_argument("--retention", choices=("all", "first-block", "none"),
                        default="first-block")
    parser.add_argument("--profile-id",
                        help="Run one exact profile from a multi-profile manifest")
    args = parser.parse_args()
    if (args.horizontal_workers < 1 or args.pointwise_workers < 1 or
            args.host_physical_memory_bytes < 1):
        parser.error("Worker counts must be positive")
    if args.output.exists():
        parser.error("Refusing to replace an existing output directory")
    for worker in (args.control_worker, args.candidate_worker):
        if not worker.is_file() or not os.access(worker, os.X_OK):
            parser.error(f"Worker is absent or not executable: {worker}")

    manifest = json.loads(args.manifest.read_text())
    contract = fixture_contract(manifest)
    generator_path = Path(manifest["generatorPath"])
    generator_receipt = {
        "path": str(generator_path),
        "declaredSHA256": manifest["generatorSHA256"],
        "currentSourcePresent": generator_path.is_file()}
    if generator_path.is_file():
        generator_receipt["currentSHA256"] = digest(generator_path)
        generator_receipt["currentMatchesDeclared"] = (
            generator_receipt["currentSHA256"] == manifest["generatorSHA256"])
    if args.mode == "qualification":
        contract["blocks"], contract["warmupBlocks"] = 1, 0
        contract["qualificationOnly"] = True
    elif contract["profileSet"] not in ("large", "boussinesq-256"):
        parser.error("Final timing mode requires the frozen large fixture schema")
    if args.mode == "final" and args.retention == "none":
        parser.error("Final timing evidence must retain its first block")
    for profile in manifest["profiles"]:
        verify_fixture(profile)
    selected_profiles = manifest["profiles"]
    if args.profile_id:
        selected_profiles = [profile for profile in selected_profiles
                             if profile["id"] == args.profile_id]
        if len(selected_profiles) != 1:
            parser.error("The requested profile id is absent or ambiguous")
    contract["profileCount"] = len(selected_profiles)
    source_preflight = {
        "control": source_snapshot(args.control_source, EXPECTED_CONTROL_COMMIT),
        "candidate": source_snapshot(args.candidate_source, EXPECTED_CANDIDATE_COMMIT)}
    provider_preflight = provider_snapshot(args.fftw_base_library,
                                           args.fftw_thread_library)
    fixture_preflight = fixture_snapshot(args.manifest, selected_profiles)
    harness_preflight = {
        "driverSHA256": digest(__file__),
        "comparisonMachinerySHA256": digest(CONSTANT_RUNNER)}
    worker_preflight = {
        "controlSHA256": digest(args.control_worker),
        "candidateSHA256": digest(args.candidate_worker),
        "workerSourceSHA256": digest(HERE / "WVVariableModelAdoptionWorker.cpp")}
    args.output.mkdir(parents=True)
    preflight = capacity_preflight(selected_profiles, contract, args.output,
                                   args.retention,
                                   args.host_physical_memory_bytes)
    provenance = {
        "schema": "wvm-variable-model-adoption-run-v1",
        "mode": args.mode, "manifestPath": str(args.manifest.resolve()),
        "retentionPolicy": args.retention,
        "manifestSHA256": digest(args.manifest),
        "fixtureGeneratorReceipt": generator_receipt,
        "selectedProfiles": [profile["id"] for profile in selected_profiles],
        "driverPath": str(Path(__file__).resolve()),
        "driverSHA256": digest(__file__),
        "comparisonMachineryPath": str(CONSTANT_RUNNER),
        "comparisonMachinerySHA256": digest(CONSTANT_RUNNER),
        "controlWorker": {"path": str(args.control_worker.resolve()),
                          "sha256": digest(args.control_worker)},
        "candidateWorker": {"path": str(args.candidate_worker.resolve()),
                            "sha256": digest(args.candidate_worker)},
        "topology": {"horizontalWorkers": args.horizontal_workers,
                     "pointwiseWorkers": args.pointwise_workers,
                     "generalVerticalWorkers": 1, "fftwInternalThreads": 1},
        "capacityPreflight": preflight,
        "limitations": [
            "Standalone C++ WVModel service injection only; public MATLAB-loaded variable service injection is absent.",
            "The Boussinesq-512 model run is excluded by its measured 20.29 GiB scientific-array requirement before model-output copies.",
            "Qualification mode is correctness-only; final timing requires the frozen large schema."],
        "expectedSourceCommits": {"control": EXPECTED_CONTROL_COMMIT,
                                  "candidate": EXPECTED_CANDIDATE_COMMIT},
        "sourcePreflight": source_preflight,
        "providerPreflight": provider_preflight,
        "fixturePreflight": fixture_preflight,
        "harnessPreflight": harness_preflight,
        "workerPreflight": worker_preflight,
        "blocks": []}
    save(args.output / "run.json", provenance)

    labels = [("control", "frozen"), ("candidate", "frozen"),
              ("candidate", "interleaved"), ("candidate", "compact")]
    global_block = 0
    for profile in selected_profiles:
        retained_control_path = None
        for block_index in range(contract["blocks"]):
            if args.retention == "none":
                order = labels
            else:
                order = labels[global_block % len(labels):] + labels[:global_block % len(labels)]
            runs = {}
            comparisons = {}
            transient_control_path = None
            for role, selection in order:
                worker = args.control_worker if role == "control" else args.candidate_worker
                receipt = run_one(worker, role, selection, profile, contract,
                                  args.output, block_index,
                                  args.horizontal_workers, args.pointwise_workers,
                                  args.control_source if role == "control"
                                  else args.candidate_source, provider_preflight)
                key = f"{role}-{selection}"
                runs[key] = receipt
                incremental_mode = (args.retention == "none" or
                                    (block_index > 0 and args.retention == "first-block"))
                if args.retention == "none" and key == "control-frozen":
                    transient_control_path = receipt["outputPath"]
                if incremental_mode:
                    incremental = compare_incremental(
                        profile, receipt,
                        transient_control_path if args.retention == "none"
                        else retained_control_path,
                        args.output / profile["id"] / f"block-{block_index:02d}", key)
                    comparisons.update(incremental)
                    successful = (receipt["passed"] and
                                  all(value["passed"] for value in incremental.values()))
                    keep_transient_control = (args.retention == "none" and
                                              key == "control-frozen")
                    receipt["outputRetained"] = not successful or keep_transient_control
                    if successful and not keep_transient_control:
                        Path(receipt["outputPath"]).unlink()
                        receipt["outputDeletedAfterVerification"] = True
                    save(Path(receipt["reportPath"]).parent / "receipt.json", receipt)
                provenance["blocks"].append({"progress": True, "profile": profile["id"],
                                             "block": block_index,
                                             "completedRun": f"{role}-{selection}"})
                save(args.output / "run.json", provenance)
            block_folder = args.output / profile["id"] / f"block-{block_index:02d}"
            if ((block_index == 0 and args.retention != "none") or
                    args.retention == "all"):
                comparisons = compare_block(profile, runs, block_folder)
                if block_index == 0 and runs["control-frozen"]["passed"]:
                    retained_control_path = runs["control-frozen"]["outputPath"]
            block = {"profile": profile["id"], "block": block_index,
                     "warmup": block_index < contract["warmupBlocks"],
                     "executionOrder": [f"{role}-{selection}" for role, selection in order],
                     "runs": runs, "comparisons": comparisons}
            block["passed"] = (all(run["passed"] for run in runs.values()) and
                               all(comparison["passed"] for comparison in comparisons.values()))
            if args.retention == "none":
                control = runs["control-frozen"]
                control["outputRetained"] = not block["passed"]
                if block["passed"]:
                    Path(control["outputPath"]).unlink()
                    control["outputDeletedAfterVerification"] = True
                save(Path(control["reportPath"]).parent / "receipt.json", control)
            for run in runs.values():
                if ((block_index == 0 and args.retention != "none") or
                        args.retention == "all"):
                    run["outputRetained"] = True
                    save(Path(run["reportPath"]).parent / "receipt.json", run)
            save(block_folder / "block.json", block)
            provenance["blocks"] = [entry for entry in provenance["blocks"]
                                    if not entry.get("progress")]
            provenance["blocks"].append(block)
            provenance["summary"] = summarize(provenance["blocks"], contract)
            save(args.output / "run.json", provenance)
            global_block += 1

    control_commits = {block["runs"]["control-frozen"]["report"]["source"]["commit"]
                       for block in provenance["blocks"] if block["passed"]}
    candidate_commits = {run["report"]["source"]["commit"]
                         for block in provenance["blocks"] if block["passed"]
                         for key, run in block["runs"].items() if key.startswith("candidate-")}
    identities_pass = (control_commits == {EXPECTED_CONTROL_COMMIT} and
                       candidate_commits == {EXPECTED_CANDIDATE_COMMIT})
    provenance["sourceIdentity"] = {"controlCommits": sorted(control_commits),
                                    "candidateCommits": sorted(candidate_commits),
                                    "independent": identities_pass}
    provenance["summary"] = summarize(provenance["blocks"], contract)
    provenance["summary"]["passed"] = provenance["summary"]["passed"] and identities_pass
    if not identities_pass:
        provenance["summary"]["sourceIdentityFailure"] = (
            "Control and candidate must each use one stable, distinct source commit")
    source_postflight = {
        "control": source_snapshot(args.control_source, EXPECTED_CONTROL_COMMIT),
        "candidate": source_snapshot(args.candidate_source, EXPECTED_CANDIDATE_COMMIT)}
    worker_postflight = {
        "controlSHA256": digest(args.control_worker),
        "candidateSHA256": digest(args.candidate_worker),
        "workerSourceSHA256": digest(HERE / "WVVariableModelAdoptionWorker.cpp")}
    provider_postflight = provider_snapshot(args.fftw_base_library,
                                            args.fftw_thread_library)
    fixture_postflight = fixture_snapshot(args.manifest, selected_profiles)
    harness_postflight = {
        "driverSHA256": digest(__file__),
        "comparisonMachinerySHA256": digest(CONSTANT_RUNNER)}
    provenance["sourcePostflight"] = source_postflight
    provenance["workerPostflight"] = worker_postflight
    provenance["providerPostflight"] = provider_postflight
    provenance["fixturePostflight"] = fixture_postflight
    provenance["harnessPostflight"] = harness_postflight
    immutable = (source_postflight == source_preflight and
                 worker_postflight == worker_preflight and
                 provider_postflight == provider_preflight and
                 fixture_postflight == fixture_preflight and
                 harness_postflight == harness_preflight)
    provenance["summary"]["passed"] = provenance["summary"]["passed"] and immutable
    provenance["summary"]["sourceAndBuildIdentityStable"] = immutable
    provenance["summary"]["sourceBuildAndProviderIdentityStable"] = immutable
    save(args.output / "run.json", provenance)
    print(json.dumps(provenance["summary"], indent=2))
    return 0 if provenance["summary"]["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
