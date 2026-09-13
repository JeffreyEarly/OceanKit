# Constant-stratification adoption qualification

This author-only harness qualifies issue #358 against a frozen v4 control. It is not part of the portable runtime or the MATLAB public API. The prospective workload, sampling and acceptance rules are recorded in [the original verification ledger](../../.github/planning/issue-358-verification.md) and [the owner-directed large-grid protocol](../../.github/planning/issue-455-large-grid-qualification.md).

The [large-grid campaign](../../.github/ci-evidence/issue-455-large-grid-adoption/README.md) passes the remaining adoption gates and the compact schedule is now the default for new builds. The owner accepted the 32×24×33 integration penalty and declined small-grid tuning. The [original campaign](../../.github/ci-evidence/issue-358-adoption/README.md) remains a historical rejected decision under its original timing gate; its interrupted samples were excluded from the new complete campaigns.

`writeConstantAdoptionFixtures` writes deterministic current-MATLAB states and independent complete-flux references at both calibration sizes. `writeConstantModelFixtures` writes four fixed-RK4 coefficient-only or particle/tracer/dense-output cases and independent MATLAB integration references. `writeConstantModelFixtures(folder,profileSet="large")` writes six four-step integration profiles at both calibration sizes: hydrostatic/nonhydrostatic coefficient-only plus nonhydrostatic tracer, particle and dense-output cases. Both writers preserve frozen payload hashes; existing outputs are not overwritten.

Build the same `WVConstantAdoptionWorker.cpp` against separate control and candidate checkouts:

```sh
cmake -S Benchmarks/constant-adoption -B /tmp/wvm-adoption-build \
  -DWVM_SOURCE=/absolute/path/to/checkout \
  -DCMAKE_BUILD_TYPE=Release \
  -DWV_RUNTIME_FFTW_ROOT=/absolute/path/to/validated/native-neon-pthreads \
  -DWV_KERNEL_COMPACT_CONSTANT_CANDIDATE=ON \
  -DWV_KERNEL_COMPACT_HORIZONTAL_WORKERS=12 \
  -DWV_KERNEL_COMPACT_POINTWISE_WORKERS=12
cmake --build /tmp/wvm-adoption-build --target wv-constant-adoption wave-vortex-run
```

The explicit settings reproduce the selected policy. Keep the frozen control on its original source/build flags; the option now defaults to `ON` for new candidate builds. The original compiled artifacts may be reused only with verified source equivalence and their original build provenance.

`freeze_build.py` captures source files, checkout state, compiler flags, executable and provider identities. Combine those receipts with host topology and the declared worker policy before measurement. `run_adoption.py` supports correctness-only qualification, worker screening and eight-pair final campaigns. Its MATLAB-loaded mode uses isolated modules built by `buildConstantAdoptionMex` and executes one fresh `matlab -batch` process per run. Run those child processes outside the local macOS sandbox as required by the workspace MATLAB instructions.

`run_model_adoption.py` preserves complete NetCDF outputs and compares every group, dimension, variable and nonvolatile attribute. Record times are exact; numerical tolerances are separate for the C++ control and independent MATLAB reference. It checks fixed-step work as well as output. `decide_adoption.py` applies the declared timing and memory gates to complete final campaigns; a screening summary cannot establish adoption.

The Python tools need NumPy; the complete-model comparator also needs netCDF4. These are authoring dependencies only. Process peak RSS comes from `wait4` after each child closes output and destroys its state. C++ capacity metrics and process RSS are distinct: provider internals and native thread storage remain outside the known C++ capacity lower bound. Preserve both measurements.

During authoritative measurement, stop competing builds, tests and MATLAB jobs. Keep failed/interrupted evidence. Each campaign creates a new directory and refuses to overwrite an earlier result. Later successful flux payloads may be removed only after their hashes and both numerical comparisons are retained; the first complete pair remains available for independent inspection. Historical small-model campaigns retain every output. Large-model campaigns support `--payload-retention first-measured-pair`: preserve the first measured pair and every failed output, then remove only successful additional/warmup payloads after both scientific comparisons and hashes are persisted. The decision checker validates exact profile/pair inventories and executable provenance.
