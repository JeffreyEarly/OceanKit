#!/usr/bin/env python3
"""Losslessly archive retained first-block model payloads after a campaign."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import sys


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def decompressed_digest(path):
    value = hashlib.sha256()
    size = 0
    with gzip.open(path, "rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            value.update(block)
            size += len(block)
    return value.hexdigest(), size


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign", type=Path)
    args = parser.parse_args()
    run_path = args.campaign / "run.json"
    archive_path = args.campaign / "payload-archive.json"
    if archive_path.exists():
        parser.error(f"Archive receipt already exists: {archive_path}")
    run = json.loads(run_path.read_text())
    if not run.get("summary", {}).get("passed"):
        parser.error("The campaign is not complete and passing")
    retained = []
    for block in run["blocks"]:
        if block.get("block") != 0:
            continue
        for key, receipt in block["runs"].items():
            if receipt.get("outputRetained"):
                retained.append((block["profile"], key, receipt))
    if len(retained) != 4:
        parser.error("Expected exactly four retained first-block payloads")
    receipt = {
        "schema": "wvm-variable-model-payload-archive-v1",
        "campaign": str(args.campaign.resolve()),
        "campaignRunSHA256BeforeArchive": digest(run_path),
        "policy": "gzip each retained payload, verify decompressed SHA-256 and size, then remove the original",
        "entries": [], "complete": False}
    save(archive_path, receipt)
    for profile, key, run_receipt in retained:
        source = Path(run_receipt["outputPath"])
        if not source.is_file():
            raise FileNotFoundError(source)
        source_size = source.stat().st_size
        source_sha = digest(source)
        if (source_size != run_receipt["outputBytes"] or
                source_sha != run_receipt["outputSHA256"]):
            raise ValueError(f"Retained payload differs from its worker receipt: {source}")
        free = shutil.disk_usage(source.parent).free
        required = source_size + 16 * 1024 ** 2
        if free < required:
            raise RuntimeError(
                f"Insufficient temporary capacity to gzip {source}: "
                f"need {required} bytes, have {free}")
        archive = source.with_suffix(source.suffix + ".gz")
        temporary = archive.with_suffix(archive.suffix + ".partial")
        if archive.exists() or temporary.exists():
            raise FileExistsError(archive)
        with source.open("rb") as input_stream, temporary.open("wb") as raw_output:
            with gzip.GzipFile(filename="output.nc", mode="wb", fileobj=raw_output,
                               compresslevel=6, mtime=0) as output_stream:
                shutil.copyfileobj(input_stream, output_stream, 8 * 1024 * 1024)
        restored_sha, restored_size = decompressed_digest(temporary)
        if restored_sha != source_sha or restored_size != source_size:
            raise ValueError(f"Gzip verification failed: {temporary}")
        temporary.rename(archive)
        entry = {"profile": profile, "run": key,
                 "originalPath": str(source), "originalBytes": source_size,
                 "originalSHA256": source_sha, "archivePath": str(archive),
                 "archiveBytes": archive.stat().st_size,
                 "archiveSHA256": digest(archive),
                 "decompressedBytes": restored_size,
                 "decompressedSHA256": restored_sha,
                 "originalRemoved": False}
        receipt["entries"].append(entry)
        save(archive_path, receipt)
        source.unlink()
        entry["originalRemoved"] = True
        save(archive_path, receipt)
    receipt["complete"] = True
    save(archive_path, receipt)
    print(json.dumps({"complete": True, "archivedPayloads": len(retained),
                      "archiveReceipt": str(archive_path)}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
