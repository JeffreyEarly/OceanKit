# Variable full-model adoption harness

This author-only harness qualifies complete SQG, Hydrostatic, and Boussinesq `WVModel`
workloads across an independent frozen source build and the service-enabled
candidate source build. It restores a copied MATLAB-authored model-output file,
prepares the stored output graph, injects the selected variable kernel services,
advances fixed RK4, and closes the output. Both graphs include coefficient and
Eulerian output, off-step dense `u`, and a full-grid scalar tracer. SQG marks
that tracer `isXYOnly`, so it uses horizontal advection while preserving the
XYZ storage shape. Hydro uses three-dimensional tracer advection and also
includes two moving particles. The SQG fixture does not include particles
because its transform has no `w` field.

The three candidate policies are:

- `frozen`: scalar vertical matrices, full FFT, one horizontal and pointwise
  worker, established interleaved spectral storage, and nonstreamed nonlinear
  evaluation.
- `interleaved`: Accelerate vertical matrices, pruned tile-16 horizontal work,
  caller-selected horizontal and pointwise workers, streamed nonlinear
  evaluation, and established interleaved spectral storage.
- `compact`: the same optimized topology with compact split fused views.

All policies keep FFTW internal threads and general vertical matrix work at one.
The old source worker has no service-injection API and accepts only `frozen`.
Its report records the model's frozen default boundary plus the identifier from
the same scalar factory. Candidate reports record every backend instance seen by
the injected factory without wrapping matrix calls.

Generate the small immutable correctness fixtures from the WVM authoring root:

```matlab
writeVariableModelAdoptionFixtures("/private/tmp/wvm-variable-model-small",profileSet="small")
```

The generator changes to its own source root before constructing transforms and
runs an 8-by-6-by-9 nonlinear preflight for each family before the declared
fixture. The small grid is exactly 32-by-24-by-33 with four fixed steps. The
large grid is exactly 256-by-256-by-129 with one fixed step. Large generation
has an 8 GiB free-space guard. It refuses instead of substituting a smaller
grid. Boussinesq 256 uses its separately saved exact modal source, so the
fixture does not repeat expensive vertical-mode construction:

```matlab
writeBoussinesqVariableModelAdoptionFixture( ...
    "/private/tmp/wvm-variable-model-bouss256", ...
    "/private/tmp/wvm455-variable-bouss-256-fixtures-v2", ...
    numericalSourceRoot="/path/to/frozen-candidate-a3acbf4d")
```

Run MATLAB after `restoredefaultpath; clear classes`, configure the frozen
candidate source, assert `which("WVModel")` is below that source root, and `cd`
to the source root before calling either generator. The Boussinesq manifest
records the numerical and generator source identities separately, the upstream
modal-source hash, actual output sizes, host physical memory, and free-space
preflight. The 512 grid is the capacity exclusion: measured Boussinesq
scientific arrays alone require 20.29 GiB before model-output copies.

Build the same worker source against distinct WVM checkouts:

```sh
cmake -S Benchmarks/variable-adoption/model-adoption -B /private/tmp/wvm-model-control \
  -DWVM_SOURCE=/path/to/frozen-wvm
cmake --build /private/tmp/wvm-model-control --target wv-variable-model-adoption -j 4

cmake -S Benchmarks/variable-adoption/model-adoption -B /private/tmp/wvm-model-candidate \
  -DWVM_SOURCE=/path/to/candidate-wvm
cmake --build /private/tmp/wvm-model-candidate --target wv-variable-model-adoption -j 4
```

Run the correctness qualification with the calibrated outer worker counts:

```sh
python3 Benchmarks/variable-adoption/model-adoption/run_variable_model_adoption.py \
  --manifest /private/tmp/wvm-variable-model-small/manifest.json \
  --control-worker /private/tmp/wvm-model-control/wv-variable-model-adoption \
  --candidate-worker /private/tmp/wvm-model-candidate/wv-variable-model-adoption \
  --control-source /path/to/frozen-b1579883 \
  --candidate-source /path/to/candidate-a3acbf4d \
  --fftw-base-library /path/to/pinned/lib/libfftw3.3.dylib \
  --fftw-thread-library /path/to/pinned/lib/libfftw3_threads.3.dylib \
  --output /private/tmp/wvm-variable-model-qualification \
  --horizontal-workers 12 --pointwise-workers 8 \
  --host-physical-memory-bytes 137438953472 --mode qualification
```

The driver verifies immutable fixture hashes, exact fixed-step/RHS work, tracer,
particle, dense-output and NetCDF work, actual injected backend construction,
FFTW topology, and every NetCDF group, variable, dimension, attribute and
payload against MATLAB and the independent control. It rotates process order.
The campaign is bound to frozen commits `b1579883` and `a3acbf4d`; source-tree,
binary, and worker-source identities are checked before and after all runs.
Qualification mode makes no timing claim. Final mode accepts the frozen
SQG/Hydro large manifest and the exact Boussinesq-256 manifest. It uses two
excluded warmup blocks plus eight measured blocks per selected profile. Run
profiles as separate campaigns so retained evidence can be archived between
them:

```sh
python3 Benchmarks/variable-adoption/model-adoption/run_variable_model_adoption.py \
  --manifest /private/tmp/wvm-variable-model-large/manifest.json \
  --profile-id variable-stratified-qg-composite-large \
  --control-worker /private/tmp/wvm-model-control/wv-variable-model-adoption \
  --candidate-worker /private/tmp/wvm-model-candidate/wv-variable-model-adoption \
  --control-source /path/to/frozen-b1579883 \
  --candidate-source /path/to/candidate-a3acbf4d \
  --fftw-base-library /path/to/pinned/lib/libfftw3.3.dylib \
  --fftw-thread-library /path/to/pinned/lib/libfftw3_threads.3.dylib \
  --output /private/tmp/wvm-variable-model-sqg-final \
  --horizontal-workers 12 --pointwise-workers 8 \
  --host-physical-memory-bytes 137438953472 --mode final \
  --retention first-block
```

With `first-block`, all four first-block payloads remain available. Every later
payload is hashed, compared with MATLAB and the retained independent control,
and deleted immediately after successful comparisons; failures remain on disk.
After a profile campaign finishes, archive the four retained payloads before
starting the next profile:

```sh
python3 Benchmarks/variable-adoption/model-adoption/archive_variable_model_outputs.py \
  /private/tmp/wvm-variable-model-sqg-final
```

The archiver verifies each decompressed byte count and SHA-256 against the
original worker receipt before deleting the uncompressed payload. It writes a
separate `payload-archive.json` receipt and leaves the completed run receipt
unchanged.

This is a standalone C++ service-injection qualification. The public MATLAB
loaded-variable MEX boundary does not expose these services and is outside this
harness.
