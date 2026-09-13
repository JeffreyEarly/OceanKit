# Variable-transform optimization screen

This author-only harness compares explicit C++ kernel execution options for issue #455. Production constructors retain the frozen path. The [prospective protocol](../../.github/planning/issue-455-variable-matrix-qualification.md) separates this first screen from complete-model adoption.

`writeVariableAdoptionFixtures(folder)` creates independent MATLAB fluxes and model checkpoints for Stratified QG and Hydrostatic at 256×256×129. Optional `grid` and `families` select correctness fixtures or later declared profiles. Boussinesq uses the existing transform family and its exact persisted matrix groups. Every fixture directory is new; previous outputs are never overwritten. Flux bytes are canonical column-major interleaved complex Float64, little endian: A0 for SQG, Fp/Fm/F0 for Hydrostatic and Boussinesq.

Build against the candidate authoring checkout:

```sh
cmake -S Benchmarks/variable-adoption -B /tmp/wvm-variable-screen \
  -DWVM_SOURCE=/absolute/path/to/checkout \
  -DCMAKE_BUILD_TYPE=Release -DWV_ENABLE_ACCELERATE=ON \
  -DWV_RUNTIME_FFTW_ROOT=/absolute/path/to/native-neon-pthreads
cmake --build /tmp/wvm-variable-screen --parallel 6 --target wv-variable-adoption
python run_screening.py /absolute/path/to/fixtures/manifest.json \
  /tmp/wvm-variable-screen/wv-variable-adoption /new/evidence/folder
```

The runner needs NumPy. It runs four fresh-process blocks, alternating the order of `frozen`, `pruned`, `streamed`, and `pruned-streamed`. Each process uses two warmups and four samples, with twelve prepared outer horizontal workers for pruned selections. `--smoke` uses one block and one sample solely to check fixture loading, output encoding and scientific comparisons; its timings establish no performance result.

Every run preserves stdout, stderr and the complete binary flux payload. The runner checks fixture hashes, exact payload lengths, finiteness and both independent MATLAB and frozen-output errors. Receipts include source hashes/diff, binary/provider identities, raw times, retained capacities and fresh-process peak RSS. The kernel ledger excludes opaque provider internals; process RSS is a distinct measurement. These measurements cover direct complete nonlinear flux only. Prepared/raw event paths are covered by focused C++ tests, not timed here. MATLAB-loaded and complete-model continuation/output qualification remain required before adopting a new default.

The implementation keeps shared-matrix direct interleaved ZGEMM and exact Boussinesq packing. This screen measures horizontal composition and movement reductions. It does not settle split two-DGEMM versus direct-interleaved execution or grouped-matrix scheduling.

## Derivative-resource qualification

`wv-variable-derivative` measures horizontal passive-scalar advection with a manufactured scalar containing low, high and Nyquist frequencies, depth-varying amplitude and prescribed horizontal velocity. Its analytical oracle checks full-grid frequency coverage; no antialias filter or vertical advection participates. It reports exact input preservation, storage stability and scalar error outside the timed region. This workload tests #463's derivative resources and is separate from complete-model/event qualification.

`run_derivative_qualification.py` compares independent production and prior-pruned controls with the bounded-resource candidate for both complete flux and this scalar workload. Build both workers from the same harness directory, setting `WVM_SOURCE` separately to the frozen and candidate source trees. Pass the fixture manifest, both source directories, both build directories and a new evidence directory to the runner; `--smoke` checks the harness without making a timing claim. The [combined goal protocol](../../.github/planning/issue-455-complete-variable-adoption.md) declares the trial counts, selections and payload-retention policy before measurement. Sources, fixtures, harness, binaries and provider identities retain distinct provenance.

## MATLAB source isolation

Run qualification in a fresh MATLAB batch. Before adding dependencies or the authoring checkout, use `restoredefaultpath; clear classes`; then add the intended v4 tools, call `configureCIEnvironment`, and change into that checkout. Assert `which` for `WVModel`, `WVTransform`, `WVNonlinearAdvection`, `WVForcingType` and the relevant transform/geometry classes resolves inside the selected v4 root. Changing directories alone does not unload classes cached from a v5 startup path. Scientific fixtures retain their original generator/source identities; isolated oracle rechecks are separate receipts.
