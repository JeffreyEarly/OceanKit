#!/usr/bin/env python3
"""Freeze an authoring checkout, compiler flags and provider alongside a binary."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def command(*args):
    return subprocess.check_output(args, text=True).strip()


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def freeze(source, build, executable, provider):
    paths = subprocess.check_output(["git", "-C", str(source), "ls-files", "--cached", "--others", "--exclude-standard", "-z"]).decode().split("\0")
    files = {p: digest(source / p) for p in paths if p and (source / p).is_file()
             and (p.startswith(("CompiledKernel/", "PortableRuntime/", "Benchmarks/constant-adoption/", "@WVCompiledBackend/")))}
    flags = {str(p.relative_to(build)): p.read_text() for p in build.rglob("flags.make")}
    return {"source": {"path": str(source), "commit": command("git", "-C", str(source), "rev-parse", "HEAD"),
             "tree": command("git", "-C", str(source), "rev-parse", "HEAD^{tree}"),
             "status": command("git", "-C", str(source), "status", "--short"), "sourceSHA256": files},
            "build": {"path": str(build), "cache": (build / "CMakeCache.txt").read_text(), "targetFlags": flags},
            "executable": {"path": str(executable), "sha256": digest(executable)},
            "provider": {"root": str(provider), "libraries": {name: {"path": str((provider / "lib" / name).resolve()), "sha256": digest(provider / "lib" / name)} for name in ["libfftw3.dylib", "libfftw3_threads.dylib"]}}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ["source", "build", "executable", "provider", "output"]:
        parser.add_argument("--" + arg, required=True, type=Path)
    args = parser.parse_args()
    result = freeze(args.source.resolve(), args.build.resolve(), args.executable.resolve(), args.provider.resolve())
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2)
        stream.write("\n")
