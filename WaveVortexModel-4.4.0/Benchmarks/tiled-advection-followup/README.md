# Native tiled-stage screen

This screen isolates the actual prepared native advection operator. It excludes vertical MM, field assembly, forcing and integration. Use the complete-model campaign before making an integration-speed claim. Follow [the optimization workflow](../CPP-OPTIMIZATION.md), including its idle-host and source-freeze requirements.

Build `NativeStage.cpp` against each frozen candidate's native provider library and the common, unchanged core library. The original compiler command is retained as `verification/wvm491-build-stage.py` in the artifact archive named in [the qualification report](../TILED-SCRATCH-AND-ARITHMETIC.md). It uses C++17, `-O3 -DNDEBUG`, the existing FFTW provider and no fast-math. Freeze executable/library hashes and the source differences between roles before running.

The screen expects an archive directory containing `baseline/NativeStage`, `memory-only/NativeStage`, `simd/NativeStage`, `modes-28.txt` and `modes-129.txt`. A mode file begins with `Nx Ny Nz M Lx Ly`, followed by `M` integer `k l` pairs. Use the retained mode sets from the representative fixtures; this screen's comparison partitions assume `Nx=Ny=256`.

Run with Python and NumPy:

```sh
python3 Benchmarks/tiled-advection-followup/run_stage_screen.py "$benchmark_archive" --blocks 4
```

The output directory must not already exist. The script rotates/reverses three fresh-process roles over four blocks, with two warmups and nine samples per process. Each process restores the target-z input outside the measured interval. Full physical fields and split target outputs are compared before the generated binary payloads are deleted; output hashes, timing samples, workspace bytes and producer counts remain in `stage-screen/`.

The direct executable interface is `NativeStage mode-file targets output-prefix`, where `targets` is 3 or 4. It writes a full-output `.bin` payload and emits one JSON timing/count record on stdout. Keep the frozen role provenance with the results; executable names alone do not establish which implementation was measured.
