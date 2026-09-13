#!/usr/bin/env python3
"""Alternate prepared benchmark blocks between already resident workers."""

import argparse
import array
import hashlib
import json
import math
import os
import pathlib
import subprocess
import sys
import time

try:
    import numpy as np
except ImportError:
    np = None


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_worker(value):
    if "=" not in value:
        raise argparse.ArgumentTypeError("worker must be LABEL=EXECUTABLE")
    label, executable = value.split("=", 1)
    if not label or not executable:
        raise argparse.ArgumentTypeError("worker must be LABEL=EXECUTABLE")
    return label, pathlib.Path(executable).resolve()


def read_record(process, label):
    line = process.stdout.readline()
    if not line:
        stderr = process.stderr.read()
        raise RuntimeError(f"{label} worker exited without a record: {stderr}")
    return json.loads(line)


def write_receipt(path, receipt):
    temporary = path.with_suffix(".tmp")
    with open(temporary, "w", encoding="utf-8") as stream:
        json.dump(receipt, stream, indent=2, sort_keys=True)
        stream.write("\n")
    os.replace(temporary, path)


def compare_payload(reference, candidate):
    if reference.stat().st_size != candidate.stat().st_size:
        raise RuntimeError("Flux payload lengths differ")
    if reference == candidate:
        return {"passed": True, "maximumAbsoluteDifference": 0.0,
                "relativeTolerance": 1e-10, "absoluteTolerance": 1e-12}
    maximum = 0.0
    with reference.open("rb") as left, candidate.open("rb") as right:
        while block := left.read(1024 * 1024):
            other = right.read(len(block))
            if len(block) % 8 or len(other) != len(block):
                raise RuntimeError("Malformed double-precision payload")
            if np is not None:
                x, y = np.frombuffer(block, dtype="=f8"), np.frombuffer(other, dtype="=f8")
                if not np.isfinite(x).all() or not np.isfinite(y).all():
                    raise RuntimeError("Nonfinite payload")
                difference = np.abs(x - y)
                maximum = max(maximum, float(difference.max()))
                if np.any(difference > 1e-12 + 1e-10 * np.abs(x)):
                    raise RuntimeError("Cross-worker scientific comparison failed")
                continue
            a, b = array.array("d"), array.array("d")
            a.frombytes(block)
            b.frombytes(other)
            for x, y in zip(a, b):
                if not math.isfinite(x) or not math.isfinite(y):
                    raise RuntimeError("Nonfinite flux payload")
                difference = abs(x - y)
                maximum = max(maximum, difference)
                if difference > 1e-12 + 1e-10 * abs(x):
                    raise RuntimeError("Cross-worker flux comparison failed")
    return {"passed": True, "maximumAbsoluteDifference": maximum,
            "relativeTolerance": 1e-10, "absoluteTolerance": 1e-12}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=pathlib.Path)
    parser.add_argument("--worker", required=True, action="append", type=parse_worker)
    parser.add_argument("--sequence", required=True,
                        help="Comma-separated resident worker labels, in block order")
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--samples", type=int, default=8)
    parser.add_argument("--output", required=True, type=pathlib.Path)
    parser.add_argument("--expected-provider-root", type=pathlib.Path)
    args = parser.parse_args()
    if args.warmups <= 0 or args.samples <= 0:
        parser.error("warmups and samples must be positive")
    workers = dict(args.worker)
    sequence = [item for item in args.sequence.split(",") if item]
    if not sequence or any(label not in workers for label in sequence):
        parser.error("sequence contains an unknown or empty worker label")
    if len(workers) != len(args.worker):
        parser.error("worker labels must be unique")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    input_path = args.input.resolve()

    processes = {}
    started = time.monotonic()
    receipt = {
        "schema": "wvm-prepared-pipeline-campaign-v1",
        "completed": False,
        "driver": {"path": str(pathlib.Path(__file__).resolve()),
                   "sha256": sha256(__file__), "argv": sys.argv,
                   "numpy": None if np is None else np.__version__},
        "input": {"path": str(input_path), "sha256": sha256(input_path)},
        "sequence": sequence,
        "warmups": args.warmups,
        "samples": args.samples,
        "workers": {},
        "blocks": [],
    }
    try:
        for label, executable in workers.items():
            process = subprocess.Popen(
                [str(executable), str(input_path)], stdin=subprocess.PIPE,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                bufsize=1)
            processes[label] = process
            ready = read_record(process, label)
            if ready.get("event") != "ready":
                raise RuntimeError(f"{label} did not emit a ready record")
            for key, filename in (("workerSha256", "WVPreparedPipelineWorker.cpp"),
                                  ("cmakeSha256", "CMakeLists.txt")):
                if ready["source"][key] != sha256(pathlib.Path(__file__).with_name(filename)):
                    raise RuntimeError(f"Stale benchmark build: {label} {filename}")
            provider_files = {}
            for key in ("baseLibrary", "threadLibrary"):
                library = pathlib.Path(ready["provider"][key]).resolve()
                if args.expected_provider_root:
                    expected = args.expected_provider_root.resolve()
                    if os.path.commonpath((library, expected)) != str(expected):
                        raise RuntimeError(
                            f"{label} {key} is outside expected provider root")
                provider_files[key] = {
                    "path": str(library), "sha256": sha256(library)}
            receipt["workers"][label] = {
                "executable": str(executable),
                "executableSha256": sha256(executable),
                "providerFiles": provider_files,
                "ready": ready,
            }
            first = next(iter(receipt["workers"].values()))
            for key in ("compiler", "provider", "family", "grid", "options",
                        "matrixBackend", "horizontalSchedule"):
                if ready[key] != first["ready"][key]:
                    raise RuntimeError(f"Worker configuration differs: {key}")
            if provider_files != first["providerFiles"]:
                raise RuntimeError("Worker provider library identities differ")
            write_receipt(output / "receipt.json", receipt)

        payload_written = set()
        for index, label in enumerate(sequence):
            payload = ""
            if label not in payload_written:
                payload = str(output / f"{label}-flux.bin")
                payload_written.add(label)
            command = {
                "command": "run", "id": f"{index}-{label}",
                "warmups": args.warmups, "samples": args.samples,
                "payload": payload,
            }
            process = processes[label]
            process.stdin.write(json.dumps(command) + "\n")
            process.stdin.flush()
            record = read_record(process, label)
            if record.get("event") != "result":
                raise RuntimeError(f"{label} did not emit a result record")
            record["worker"] = label
            if payload:
                record["payloadSha256"] = sha256(payload)
                record["fieldPayloadSha256"] = sha256(payload + ".fields.bin")
            receipt["blocks"].append(record)
            with open(output / "blocks.jsonl", "a", encoding="utf-8") as stream:
                stream.write(json.dumps(record, sort_keys=True) + "\n")
            write_receipt(output / "receipt.json", receipt)
        reference = next(iter(workers))
        receipt["fluxComparisons"] = {
            label: compare_payload(output / f"{reference}-flux.bin",
                                   output / f"{label}-flux.bin")
            for label in workers}
        receipt["fieldComparisons"] = {
            label: compare_payload(output / f"{reference}-flux.bin.fields.bin",
                                   output / f"{label}-flux.bin.fields.bin")
            for label in workers}
        first = receipt["blocks"][0]
        for block in receipt["blocks"]:
            if block["producerMetrics"] != first["producerMetrics"]:
                raise RuntimeError("Producer counts differ between blocks")
        if sha256(input_path) != receipt["input"]["sha256"]:
            raise RuntimeError("Input file changed during campaign")
        for item in receipt["workers"].values():
            if sha256(item["executable"]) != item["executableSha256"]:
                raise RuntimeError("Executable changed during campaign")
            for library in item["providerFiles"].values():
                if sha256(library["path"]) != library["sha256"]:
                    raise RuntimeError("Provider changed during campaign")
        if sha256(__file__) != receipt["driver"]["sha256"]:
            raise RuntimeError("Campaign driver changed during execution")
        # The large field payloads are disposable validation outputs. Retain
        # their hashes and comparisons; compact flux payloads remain available.
        for label in workers:
            (output / f"{label}-flux.bin.fields.bin").unlink()
        receipt["fieldPayloadDisposition"] = "removed after successful comparison; hashes retained"
        receipt["completed"] = True
    except Exception as error:
        receipt["failure"] = {"type": type(error).__name__, "message": str(error)}
        raise
    finally:
        for process in processes.values():
            if process.poll() is None:
                try:
                    process.stdin.write('{"command":"quit"}\n')
                    process.stdin.flush()
                except BrokenPipeError:
                    pass
        for label, process in processes.items():
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                process.terminate()
                process.wait(timeout=10)
            stderr = process.stderr.read()
            if stderr:
                (output / f"{label}.stderr.log").write_text(stderr,
                                                            encoding="utf-8")
            if process.returncode != 0:
                receipt["completed"] = False
                receipt.setdefault("workerFailures", {})[label] = process.returncode
        receipt["processLifetimeSeconds"] = time.monotonic() - started
        write_receipt(output / "receipt.json", receipt)
    if not receipt.get("completed"):
        raise RuntimeError("Prepared campaign did not complete successfully")
    print(json.dumps({"output": str(output), "blocks": len(receipt["blocks"])},
                     sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(error, file=sys.stderr)
        raise
