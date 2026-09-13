# Prepared nonlinear pipeline benchmark

This benchmark loads one Hydrostatic or Boussinesq checkpoint, constructs one
native kernel, and keeps it resident while timing repeated
`nonlinearFluxAndFields` calls. Every call owns a fresh
`beginStateEvaluation`/`endStateEvaluation` scope. The persisted forcing
schedule is reported but is not executed.

The frozen schedule is compact split views, streaming pruned tile 16,
`tiledNonlinear`, `fusedDerivativeAdvection`, and `sharedFieldGradients`, with
12 horizontal workers, 8 pointwise workers, and 1 Hydrostatic or 8 Boussinesq
vertical-group workers.

Configure and build outside a timing window:

```sh
cmake -S Benchmarks/prepared-pipeline -B /tmp/wv-prepared-build \
  -DWVM_SOURCE="$PWD" -DCMAKE_BUILD_TYPE=Release \
  -DWV_RUNTIME_FFTW_ROOT=/absolute/path/to/pinned/fftw
cmake --build /tmp/wv-prepared-build --target wv-prepared-pipeline
```

The worker uses JSON lines on standard input. It emits a `ready` record after
checkpoint loading and kernel preparation, accepts any number of commands such
as the following, and exits on `quit`:

```json
{"command":"run","id":"candidate-0","warmups":2,"samples":8,"payload":"/tmp/candidate-flux.bin"}
```

The campaign driver can keep baseline and candidate executables resident and
alternate measurement blocks without paying checkpoint and kernel preparation
for each block:

```sh
python3 Benchmarks/prepared-pipeline/run_prepared_pipeline.py \
  --input /absolute/path/to/source.nc \
  --worker baseline=/absolute/path/to/baseline/wv-prepared-pipeline \
  --worker candidate=/absolute/path/to/candidate/wv-prepared-pipeline \
  --sequence baseline,candidate,candidate,baseline \
  --warmups 2 --samples 8 --output /absolute/path/to/new-results
```

Warmup and sample times are separate. To keep validation from dominating the
resident screen, output scanning and replay comparison run after the last
warmup and last sample in each block and are reported as `validationSeconds`.
They use `atol=1e-12` and `rtol=1e-10`. Flux payloads are written once per
executable. Physical fields are written once per executable, compared across roles at the same tolerances, and removed after successful comparison with their hashes retained. NumPy, when available, accelerates chunked payload comparisons; a standard-library fallback is supported. The
receipt records input, executable, provider-library, worker-source, CMake,
tracked-source-diff, and source-status hashes, along with source revision,
compiler identity, kernel storage, and actual producer counters. Pass
`--expected-provider-root` to reject a worker linked outside the frozen FFTW
installation.

Sample timing includes beginning the state scope, nonlinear flux/field production, and ending the scope. It excludes loading, preparation, validation, damping, integration and output. These are exploratory RHS timings, not full-model acceptance. Each block requires at least one warmup to establish replay evidence.

Freeze and build the same benchmark sources for every role. The driver rejects stale embedded worker/CMake hashes, differing compiler/provider/execution options, differing producer counts, numerical disagreement, or modified frozen files. It journals completed blocks and explicit failure receipts. Only one worker receives work at a time. Check and record host activity externally; this driver does not assert that the host is idle. Resident roles require enough RAM for both prepared kernels, reported separately from kernel-owned storage.

See [the resident screens and completed paired qualification](RESULTS.md) for the selected assembly/projection candidate, measured benefits, host/order limitations and adoption decision.
