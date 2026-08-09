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

Every scored case is normalized against its matching committed builtin reference. Case score 100 matches the reference. Family and suite scores use geometric means, with transform families weighted equally. Same-host builtin-to-candidate speedups are reported separately.

Ordinary artifacts are written beneath the ignored `results/runs` directory. Immutable reference artifacts live beneath `results/reference`. Memory measurements use a fresh MATLAB process and report baseline, persistent, and peak-observed resident memory; these process measurements include MATLAB runtime and allocator behavior.

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

The 3% regression threshold applies only to `[256 256 65]` and `[512 512 129]`, with antialiasing both disabled and enabled. Smaller cases remain descriptive because the row layout was selected as the single production representation even where MATLAB timing noise or fixed overhead can make a legacy expression faster. The integration artifact also confirms that the vertically replicated compatibility indices remain unallocated until a caller explicitly requests them.

The v4.2.1 release audit repeated this gate against the final production source. Its immutable M5 Max/R2026a builtin result is stored under `results/reference/transform-layout-v4.2.1-release-m5-max-r2026a-builtin`.

`WVTransformConstantStratificationSpeedTest`, `ProfileableSpeedTest`, and `ForcingSpectralMaskPerformanceTest` remain historical investigation scripts. Deterministic correctness checks belong in `UnitTests`; mixed scientific investigations belong in `DeveloperExperiments`.
