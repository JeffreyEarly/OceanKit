#!/usr/bin/env python3
"""Apply the prospectively declared #358 gates to complete final campaigns."""
import argparse
import json
from pathlib import Path


def read(folder, name):
    return json.loads((folder / name).read_text())


def decide(native, matlab, model, selection):
    checks = []
    campaigns = []

    def require(name, passed, actual=None, bound=None):
        checks.append({"gate": name, "passed": bool(passed), "actual": actual, "bound": bound})

    for boundary, folder in [("native", native), ("matlab-loaded", matlab)]:
        manifest, runs, summary = [read(folder, name) for name in ("manifest.json", "runs.json", "summary.json")]
        campaigns.append(manifest)
        require(boundary + ": correct execution boundary", manifest["interface"] == boundary)
        executable_provenance = ("control", "candidate") if boundary == "native" else ("controlMex", "candidateMex")
        for role, provenance_key in zip(("control", "candidate"), executable_provenance):
            recorded = manifest["executables"][role]["sha256"]
            source = manifest["provenance"][provenance_key]
            expected = source["executable"]["sha256"] if boundary == "native" else source["moduleSHA256"]
            require(boundary + ": "+role+" binary provenance", recorded == expected)
        require(boundary + ": final protocol", manifest["mode"] == "final" and manifest["pairs"] == 8 and manifest["warmups"] == 2 and manifest["samples"] == 7)
        require(boundary + ": all paired profiles", len(runs) == 64 and len(summary["profiles"]) == 4 and all(row["pairs"] == 8 for row in summary["profiles"]))
        if "currentIntegrationAndProtocol" in manifest["provenance"]:
            expected_flux_profiles = {"constant-" + kind + "-" + size for kind in
                                     ("hydrostatic", "nonhydrostatic") for size in
                                     ("256x256x129", "512x512x257")}
            expected_runs = {(profile, pair, role) for profile in expected_flux_profiles
                             for pair in range(8) for role in ("control", "candidate")}
            actual_runs = {(r["profile"], r["pair"], r["role"]) for r in runs}
            require(boundary + ": exact large run inventory", actual_runs == expected_runs and len(actual_runs) == len(runs))
            require(boundary + ": exact large summary inventory", {r["profile"] for r in summary["profiles"]} == expected_flux_profiles)
        require(boundary + ": independent MATLAB agreement", all(r["exitCode"] == 0 and r["matlabComparison"]["passed"] for r in runs))
        candidates = [r for r in runs if r["role"] == "candidate"]
        require(boundary + ": frozen control agreement", len(candidates) == 32 and all(r["controlComparison"]["passed"] for r in candidates))
        require(boundary + ": geometric flux improvement", summary["geometricFluxRatio"] <= .90, summary["geometricFluxRatio"], .90)
        require(boundary + ": no profile regression", all(r["fluxRatio"] <= 1.03 for r in summary["profiles"]), {r["profile"]: r["fluxRatio"] for r in summary["profiles"]}, 1.03)
        require(boundary + ": paired interval excludes tie", summary["fluxPairedBootstrap95"][1] < 1, summary["fluxPairedBootstrap95"], 1)
        require(boundary + ": complete-lifetime RSS", summary["geometricPeakRSSRatio"] <= 1 and summary["peakRSSPairedBootstrap95"][1] <= 1.03, [summary["geometricPeakRSSRatio"], summary["peakRSSPairedBootstrap95"]], [1, 1.03])
        controls = {(r["profile"], r["pair"]): r for r in runs if r["role"] == "control"}
        ratios = []
        for candidate in candidates:
            control = controls[candidate["profile"], candidate["pair"]]
            def owned(run):
                worker = run["worker"]
                return worker["adoption"]["kernelPersistentBytes"] if boundary == "native" else worker["metrics"]["persistentBytes"]
            ratios.append(owned(candidate) / owned(control))
        require(boundary + ": no owned-memory growth", all(r <= 1 for r in ratios), max(ratios), 1)
        policy = manifest["provenance"]["effectivePolicy"]
        require(boundary + ": selected shared policy", policy["pointwiseWorkers"] == selection["selectedPointwiseWorkers"] and policy["candidateHorizontalOuter"] == 12 and policy["candidateTypeIInternal"] == 16 and policy["coefficientWorkers"] == 2)

    require("boundaries: identical scientific fixtures", campaigns[0]["fixturesManifestSHA256"] == campaigns[1]["fixturesManifestSHA256"])
    require("boundaries: identical frozen provenance", campaigns[0]["provenanceSHA256"] == campaigns[1]["provenanceSHA256"])

    manifest, pairs, summary = [read(model, name) for name in ("manifest.json", "pairs.json", "summary.json")]
    require("model: identical frozen provenance", manifest["provenanceSHA256"] == campaigns[0]["provenanceSHA256"])
    large = manifest.get("fixtures", {}).get("schema") == "wvm-constant-model-fixtures-large-v1"
    if large:
        expected_profiles = [
            "constant-hydrostatic-coefficient-only-256",
            "constant-nonhydrostatic-coefficient-only-256",
            "constant-nonhydrostatic-composite-256",
            "constant-hydrostatic-coefficient-only-512",
            "constant-nonhydrostatic-coefficient-only-512",
            "constant-nonhydrostatic-composite-512",
        ]
        fixture = manifest["fixtures"]
        expected_parameters = fixture["parameters"]
        current_protocol = manifest["provenance"].get("currentIntegrationAndProtocol", {})
        require("model: large protocol marker", current_protocol.get("schema") == "wvm-large-grid-constant-adoption-v1" and current_protocol.get("protocolPath") == ".github/planning/issue-455-large-grid-qualification.md")
        require("model: large fixture axes", expected_parameters["Nxyz"] == [[256, 256, 129], [512, 512, 257]])
        require("model: large fixed work", expected_parameters["fixedStepCount"] == 4 and expected_parameters["rhsEvaluationCount"] == 16 and manifest["expectedStepCount"] == 4 and manifest["expectedRHSEvaluationCount"] == 16)
        require("model: large output workload", expected_parameters["finalTime"] == 2 and expected_parameters["initialStep"] == .5 and expected_parameters["fields"] == ["u"] and expected_parameters["denseOutputTimes"] == .75)
        require("model: large protocol", manifest["mode"] == "final" and manifest["measuredPairs"] == 8 and manifest["warmupPairs"] == 2 and len(pairs) == 6 * (8 + manifest["warmupPairs"]))
        profile_ids = [p["id"] for p in fixture["profiles"]]
        require("model: exact large profile inventory", profile_ids == expected_profiles)
        expected_pairs = {(profile, pair) for profile in expected_profiles for pair in range(-2, 8)}
        actual_pairs = {(pair["profile"], pair["pair"]) for pair in pairs}
        require("model: exact large pair inventory", actual_pairs == expected_pairs and len(actual_pairs) == len(pairs))
        require("model: large warmup classification", all(p["warmup"] == (p["pair"] < 0) for p in pairs))
        for role, provenance_key in (("control", "modelControl"), ("candidate", "modelCandidateExecutable")):
            source = manifest["provenance"][provenance_key]
            expected_hash = source["executable"]["sha256"] if role == "control" else source["sha256"]
            require("model: "+role+" binary provenance", manifest["executables"][role]["sha256"] == expected_hash)
        require("model: current protocol receipt", bool(current_protocol))
    else:
        require("model: complete final campaign", manifest["mode"] == "final" and manifest["measuredPairs"] == 8 and len(pairs) == 4 * (8 + manifest["warmupPairs"]))
    require("model: scientific output and matched work", summary["passed"] and all(p["passed"] for p in pairs))
    profiles = summary.get("profiles", [])
    declared_profile_ids = {p["profile"] for p in profiles}
    expected_profile_ids = set(expected_profiles) if large else None
    require("model: all declared profiles", len(profiles) == (6 if large else 4) and all(p["pairs"] == 8 for p in profiles) and (not large or declared_profile_ids == expected_profile_ids))
    require("model: no profile regression", bool(profiles) and all(p["integrationGeometricRatio"] <= 1.03 for p in profiles), {p["profile"]: p["integrationGeometricRatio"] for p in profiles}, 1.03)
    owned_ratios = []
    for pair in pairs:
        if pair["warmup"]:
            continue
        control, candidate = pair["runs"]["control"], pair["runs"]["candidate"]
        for key in ["fullModelRetained", "fullModelMaximumLive"]:
            owned_ratios.append(candidate["livenessBytes"][key] / control["livenessBytes"][key])
    require("model: no owned-memory growth", bool(owned_ratios) and all(r <= 1 for r in owned_ratios), max(owned_ratios, default=None), 1)
    rss = summary.get("pairedBootstrap", {}).get("peakRSS", {})
    require("model: complete-lifetime RSS", bool(rss) and rss["geometricRatio"] <= 1 and rss["confidence95"][1] <= 1.03, rss, {"geometricRatio": 1, "upper95": 1.03})
    passed = all(check["passed"] for check in checks)
    return {"schema": "wvm-constant-adoption-decision-v1", "passed": passed,
            "decision": "performance gates passed; require final source/default/CI handoff" if passed else "retain frozen default",
            "checks": checks, "selectedPointwiseWorkers": selection["selectedPointwiseWorkers"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["native", "matlab", "model", "selection", "output"]:
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = decide(args.native, args.matlab, args.model, json.loads(args.selection.read_text()))
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    raise SystemExit(0 if result["passed"] else 1)
