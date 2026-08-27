# Benchmarks

This folder contains authoring-only performance tools and historical profiling scripts. It is not included on the WaveVortexModel runtime package path.

## Reproducible benchmark suites

`runWaveVortexBenchmark` is the canonical performance and memory entry point. It measures a state-advanced `nonlinearFlux()` call while retaining ordinary production caches. State changes use the public `t`, `Ap`, `Am`, and `A0` setters; the runner never clears the variable cache explicitly.

```matlab
addpath("Benchmarks")
results = runWaveVortexBenchmark(suites="core-v1")
```

The registered suites are:

- `smoke-v1`: small cases for every transform family; no score.
- `core-v1`: the canonical constant-stratification nonlinear-advection score.
- `scaling-standard-v1`: standard horizontal and vertical scaling across transform families.
- `scaling-large-v1`: fixed large-memory scaling cases.
- `transform-layout-v1`: an unscored diagnostic of full-complex WV/DFT mapping expressions. It compares the current WV-sorted linear mapping with DFT-sorted, two-dimensional-row, and per-plane alternatives while preserving production behavior.

The runner accepts more than one suite and a subset of case IDs. Suite definitions are versioned in `waveVortexBenchmarkSuites`; changing a case matrix or score definition requires a new suite version. Backends are selected independently from transform families through `waveVortexBenchmarkBackends`.

Every scored case is normalized against the builtin reference registered for its suite in `results/catalog.json`. The runner reads the catalog entry directly; it does not infer a reference from a machine model, MATLAB release, or backend name. Case score 100 matches the reference. Family and suite scores use geometric means, with transform families weighted equally. Same-host builtin-to-candidate speedups are reported separately.

Ordinary artifacts are written beneath the ignored `results/runs` directory. Immutable reference artifacts live beneath `results/reference`. Memory measurements use a fresh MATLAB process and report baseline, persistent, and peak-observed resident memory; these process measurements include MATLAB runtime and allocator behavior.

Reference generation is an explicit authoring operation. Supply `shouldCreateReference=true` together with a temporary or candidate `referenceDirectory`; the runner writes `<referenceDirectory>/<suiteId>` and never edits the catalog automatically. Review and catalog approval are separate steps.

## Raw and published benchmark data

The detailed MATLAB runner artifact is implementation-specific. It retains construction time, cache diagnostics, scoring fields, failures, and other information useful to benchmark authors. Future raw artifacts use schema `1.1.0`, which adds the WaveVortexModel package version and a human-readable processor name without changing the measured operation.

Public results use the language-neutral schema in `schemas/published-benchmark-v1.schema.json`. One published dataset represents one suite, implementation, backend, platform, and run. It contains the timing samples, median runtime, correctness result, and process-memory measurements needed by the benchmark website without requiring consumers to understand MATLAB runner internals. Missing implementation coverage is recorded as `unavailable`; it is never treated as failed or zero performance.

The source-only compiled preview has a stricter paired runner:

```matlab
result = runCompiledPreviewBenchmark
```

Each MATLAB/compiled backend, case, and repeat runs in its own fresh MATLAB process. The raw artifact records active-backend identity, native FFTW identity, exact retained application storage, isolated operation RSS, lifecycle balance, and the availability decision. `publishedWaveVortexBenchmarksFromCompiledPreviewArtifact` converts that one paired artifact into comparable MATLAB and C++ `published-benchmark-v1` datasets. Memory is reported but does not gate preview availability.

`runThreeInterfaceBenchmarkComparison` is the matched complete-workflow companion. It runs one medium `[256 256 129]` nonhydrostatic constant-stratification case in a 150 km by 150 km by 1300 m domain and compares MATLAB builtin, MATLAB with the compiled core, and standalone C++ execution for fixed RK4 and MATLAB-compatible RK3(2), RK5(4), and RK8(7) integration. The deterministic initial condition combines GM energy level 1 with a first-baroclinic red geostrophic spectrum: the geostrophic energy rolls on below mode 4, follows a `k^(-5/3)` range through mode 16, and transitions to `k^(-3)`. The geostrophic component is rescaled to a 0.15 m/s maximum horizontal speed; GM(1) is not rescaled. Default anti-aliasing remains enabled.

The fixed RK4 step is the largest power of two no greater than the frozen-state CFL=0.25 estimate (128 s for the canonical state). The adaptive initial guess is the frozen-state CFL=0.5 estimate (approximately 295.794 s), and no `MaxStep` is supplied, leaving the MATLAB-compatible default in force. Fifty-six RK4 steps give a 7168 s integration, long enough to collect stable method-work statistics without extending the run unnecessarily. A coefficient-only workload records the endpoints. The composite workload stores one restart record and then four records over the first fixed step, including three interior dense-output points, for `u`, two 3-D particles tracking `u`, a 3-D tracer, and a source-linked mooring. One frozen MATLAB-authored model file is shared by all interfaces for each workload. MATLAB builtin uses production MATLAB transforms; both compiled interfaces use the same validated native FFTW provider. Public presentation uses only integration runtime and total process-tree peak RSS sampled during the integration phase. Startup, construction, planning, parsing, and cleanup are excluded from the primary metrics.

Publication requires exact requested/actual integrator and provider identity, matched adaptive step bounds and tolerances, phase-scoped memory evidence, method-appropriate numerical tolerance, unchanged accepted endpoint trajectories when dense output is added, and complete NetCDF output-graph agreement. Integrator-study output payloads use the matched method relative tolerance plus the recorded base component absolute tolerance; structural metadata, schedules, dimensions, names, and record counts must still match exactly. A completed worker matrix that was rejected only by an obsolete stricter numerical policy may be revalidated conservatively from stored category maxima: every failing payload category must fit one of those method tolerances, while any structural difference remains a failure. Compact diagnostics retain accepted/rejected steps, RHS and FSAL work, dense-extension evaluations, secondary RSS measurements, exact standalone persistent/workspace/history/lazy-extension byte ledgers, and every state-sized buffer's producer and last consumer. MATLAB solver storage, allocator/COW behavior, and FFT/provider storage remain explicitly opaque where exact attribution is unavailable.

The runner completes every cross-interface graph comparison and endpoint/dense trajectory gate for one process repeat before releasing that repeat's validated temporary NetCDF outputs. Per-repeat comparison evidence remains in the raw artifact and is combined after all three repeats. This bounds temporary storage to one repeat without changing worker order, timing, RSS sampling, numerical coverage, or the published medians.

Run the canonical medium case with:

```matlab
results = runThreeInterfaceBenchmarkComparison;
```

The detailed raw result and RSS samples are compressed beneath the external sibling archive `../wave-vortex-model-benchmark-artifacts/three-interface/`. They are never committed or copied into the generated website. The source tree contains only the compact normalized record and presentation; its provenance stores fixture hashes, the raw-artifact hash, and the external archive filename, SHA-256, compressed byte count, and location. The author-only standalone kernel worker is built only when `WV_RUNTIME_BUILD_BENCHMARKS=ON`; it is not part of the package or ordinary `wave-vortex-run` interface.

Normalize a MATLAB artifact with explicit platform identity and repository-relative provenance:

```matlab
dataset = publishedWaveVortexBenchmarkFromMatlabArtifact(rawPath,suiteId="scaling-standard-v1",platformId="m5-max",platformName="Apple M5 Max",provenancePath="Benchmarks/results/reference/example/benchmark.json");
```

Legacy `1.0.0` artifacts do not contain a package version or useful processor name, so normalization also requires `implementationVersion` and `processorName`. The original raw artifact is read-only and is never rewritten.

The normalizer returns an ordinary MATLAB structure in canonical field order. Authors can inspect it or encode it with `jsonencode`; Issue #140 owns the reviewed process for writing and registering public datasets.

Published dataset IDs use `<suite>--<implementation>-<backend>--<platform>--<UTC timestamp>`. Approved datasets are listed in `results/catalog.json`; the catalog entry points to the dataset instead of duplicating its machine or implementation metadata. The benchmark website work in Issue #139 owns the comparison rules needed for its plots and tables.

The catalog intentionally excludes transform-layout, storage, retirement, and other engineering gates. Those artifacts continue to support implementation decisions but are not public performance datasets.

The large suites can require substantial time and memory. A partial result is valid, but a family or suite receives no aggregate score unless every required case completes.

The transform-layout suite uses the same artifact entry point but does not participate in benchmark scores or fresh-process memory measurement:

```matlab
results = runWaveVortexBenchmark(suites="transform-layout-v1")
```

It measures extraction, primary insertion, conjugate insertion, combined insertion, and complete horizontal forward/inverse calls. Each strategy owns and reuses a persistent full-complex buffer. Array setup and mapping construction are excluded, while indexing, allocation, conjugation, reshape, transpose, and MATLAB copy-on-write behavior inherent to each expression remain timed. The strict winner has the smallest median; the production `wv-sorted-linear` strategy remains preferred when it is within 3% of that median. MATLAB pointer and copy state is reported as unavailable because no supported API exposes it for these expressions. Whole-process memory comparisons remain separate from this suite.

Issue #70 moved the builtin adapter to a row-oriented Fourier-storage layout. Its integration gate compares the production adapter—not a stand-alone approximation—with the immutable issue #69 medians:

```matlab
results = runWVFourierStorageLayoutIntegrationBenchmark
```

The 3% regression threshold applies only to `[256 256 65]` and `[512 512 129]`, with antialiasing both disabled and enabled. Smaller cases remain descriptive because the row layout was selected as the single production representation even where MATLAB timing noise or fixed overhead can make a legacy expression faster. The integration artifact predates removal of the vertically replicated compatibility properties; the historical suite now constructs equivalent expanded indices during untimed setup through `indicesFromWVGridToDFTGrid`, while its immutable canonical artifact remains unchanged.

The v4.2.1 release audit repeated this gate against the final production source. Its immutable M5 Max/R2026a builtin result is stored under `results/reference/transform-layout-v4.2.1-release-m5-max-r2026a-builtin`.

## Builtin transform storage

`runWaveVortexBuiltinStorageBenchmark` reports exact application-owned transform arrays and repeated externally sampled process RSS without assuming that MATLAB allocation behavior can be inferred from source code:

```matlab
addpath("Benchmarks")
results = runWaveVortexBuiltinStorageBenchmark
```

`runWaveVortexRetirementBenchmark` compares archived `v4.2.1` and candidate source snapshots in three fresh processes per `core-v1` case. Each worker rotates implementation order, records the required 7/3 within-process samples, compares the final numerical outputs, and proves that the builtin adapter executed. The same command includes the generic storage/RSS benchmark in its retirement artifact.

The ledger covers compact Fourier mappings, the reused builtin inverse buffer, dense vertical transform matrices, and known forward/inverse result arrays. MATLAB-internal FFT work storage remains explicitly opaque. Each case runs in three fresh MATLAB processes by default while ordinary production caches stay warm.

`WVTransformConstantStratificationSpeedTest`, `ProfileableSpeedTest`, and `ForcingSpectralMaskPerformanceTest` remain historical investigation scripts. Deterministic correctness checks belong in `UnitTests`; mixed scientific investigations belong in `DeveloperExperiments`.
