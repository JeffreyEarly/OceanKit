# Hydrostatic composite workload audit

## Recommendation

Optimize horizontal tracer derivatives using the model’s retained Fourier representation and shared compiled reconstruction infrastructure first. The hydrostatic composite benchmark loses the compiled nonlinear kernel's advantage in tracer derivatives. The output writer is a small cost. A diagnostic substitution of MATLAB `diffX` and `diffY`, with all other compiled operations retained, reduced the profile from **60.815 to 42.607 seconds** and passed complete saved-output comparisons. This is a single exploratory profile pair, not adoption qualification or replacement publication data.

The measured horizontal-derivative difference predicts roughly **26% lower whole-workload time** if the compiled implementation reaches the diagnostic MATLAB throughput. Use about 25% as an experimental hypothesis, not a promised speedup. First evaluate reuse of the prepared retained horizontal transforms, derivative multipliers and worker infrastructure. Physical-grid storage does not imply that tracer derivatives should retain frequencies excluded from the model representation. Preserve the bounded-memory objective of #463; its historical full-frequency requirement is not the intended scientific contract specified by the owner. Real vertical calculus is the second target. A broad tracer, adapter, solver or output rewrite is not justified by this audit.

## Scientific contract correction after the audit

The owner clarified that public model derivatives should not support frequencies outside the model’s retained coefficient set. The initial audit recommendation incorrectly promoted an existing implementation/test requirement from #463 into a scientific requirement. The timing measurements describe that implementation and remain valid; the proposed acceptance criteria are corrected here.

For this benchmark, the initial tracer is in retained horizontal modes and `WVTracer.fluxAtTime` projects its antialiased tendency through `transformFromSpatialDomainWithFourier` and `transformToSpatialDomainWithFourier`. Runge–Kutta combinations therefore remain in that horizontal space up to floating-point error. A grid-valued tracer requires a forward transform, whereas nonlinear reconstruction starts with coefficients, but this difference does not require separate downstream differentiation kernels. Both can use retained horizontal spectra, wavenumber multipliers and the same inverse-transform infrastructure; this does not imply projecting the tracer onto the velocity’s vertical modal basis.

Current public MATLAB `diffX`/`diffY` use the full DFT spectrum, and `addTracer` accepts initial grid values without explicitly projecting them. These are existing behaviors to audit against the intended model semantics, especially for arbitrary initial data and antialias-disabled configurations. Do not preserve excluded-frequency derivatives merely to satisfy legacy compiled tests, and do not silently change public behavior as part of a timing-only experiment. Make projection semantics explicit and verify retained inputs, excluded-mode handling, initialization, tendencies and trajectories before native qualification. The temporary builtin substitution demonstrated the cost of the current path; it did not test this corrected retained-spectrum implementation.

## Source and fixture provenance

Audit base: v4 `39b2a2253976c75bbf746c55a94675155f9bc0ea`, isolated branch `audit/hydrostatic-composite`. The primary authoring checkout is on v5 and was untouched. At audit start PR #519 was still awaiting protected auto-merge. Its constant-stratification changes do not affect the hydrostatic kernel. Hydrostatic runtime, adapter, scientific MATLAB, integrator and worker sources are unchanged from publication source `024add423bf19c682ab5fd793821f398c5e7170f`; the [provenance receipt](../.github/ci-evidence/hydrostatic-composite-audit/provenance-summary.json) records explicit blob comparisons.

The qualified archive `20260913T134854458Z-three-interface-benchmark.json.gz` is present and matches its recorded hash. Both original model fixtures are missing locally: `matched-hydrostatic-exponential-coefficient-endpoint-model.nc` and `matched-hydrostatic-exponential-composite-dense-output-model.nc`. The unchanged scientific recipe was reconstructed with [hydroCompositePrepare](hydroCompositePrepare.m). Both reconstructed files match **all 29 recorded initial-condition values exactly**, but each is one byte larger and has a different hash. They are scientifically reconstructed fixtures, not byte-identical originals. Their exact hashes, original hashes and byte sizes are in the receipt.

Apple M4 Max, MATLAB R2025b Update 4, MATLAB thread budget 16. The audit build reuses the publication's exact FFTW library binaries. Active hydrostatic policy is `compact-native-accelerate`, FFT threads 1, horizontal workers 12, pointwise workers 8. The audit MEX SHA-256 is `8fbd7e590b0b0274483540e7ceb966c8b24b20c9ac396f06764d95b2a3cb8c8a`. The qualified MEX was rebuilt for the isolated checkout; its path-bound identity was validated before execution.

Grid 256 × 256 × 129, 85 vertical modes and 11,439 retained horizontal modes; domain 150000 × 150000 × 1300 m; N²(z) = 2e-5 exp(2z/1300) s⁻². The GM(0.5) plus deterministic red geostrophic state, ode78 final time 7168 s, initial step 269.39985873487035 s, RelTol 1e-3, component-scaled base AbsTol 1e-6 and default maximum step 716.8 s follow publication.

## What “dense output” measures

The composite adds an evolving full-grid tracer, two particles, Eulerian `u` and mooring output. Dense records occur at 0, 32, 64, 96 and 128 s; coefficient records occur at 0 and 7168 s. It is a composite observer workload, not an isolated dense-interpolation or sustained file-throughput benchmark.

Published endpoint medians were MATLAB 36.984 s, MATLAB compiled 13.232 s and standalone 10.448 s. Composite medians were 55.812, 57.582 and 52.523 s respectively. Thus the additional workload costs approximately 18.827, 44.350 and 42.075 s. All **18 historical samples** are retained in the provenance receipt. In particular, compiled composite samples were 57.436675, 57.716320 and 57.581742 s; standalone samples were 52.124302, 52.523112 and 53.330458 s.

Both workloads took 15 accepted steps and no retries. The endpoint has 195 RHS calls; composite has 199 because the continuous extension requires four additional RHS evaluations. Step-count differences cannot explain the slowdown.

## Fresh diagnostic profiles

Each row is one fresh MATLAB process. [hydroCompositeProfile](hydroCompositeProfile.m) uses the benchmark worker's loading, integrator and output contract with `profileCodeHotspots` and before/after native counters. These are exploratory profiles, not an idle-host timing qualification. Normal desktop services remained active; a process snapshot during the compiled composite showed Dropbox using roughly one core. A ten-second native sample was also taken during that profile. No competing experiment was interrupted.

| Workload | MATLAB builtin | MATLAB compiled | Compiled with MATLAB horizontal derivatives |
| --- | ---: | ---: | ---: |
| Coefficient endpoint | 37.208 s | 12.961 s | Not run |
| Composite | 56.836 s | 60.815 s | 42.607 s |

Composite boundaries below are inclusive except the explicitly marked self times. Tracer child rows come from the tracer's own call edges: builtin derivatives also serve nonlinear dynamics, so their global totals would overcount tracer work.

| Boundary | MATLAB builtin | MATLAB compiled | Diagnostic substitution |
| --- | ---: | ---: | ---: |
| Nonlinear coefficient flux | 34.021 s | 14.035 s | 13.070 s |
| Tracer flux, 199 calls | 9.021 s | 30.420 s | 13.461 s |
| Tracer `diffX`, 199 calls | 1.197 s | 9.942 s | 1.115 s |
| Tracer `diffY`, 199 calls | 1.489 s | 8.461 s | 1.353 s |
| Tracer `diffZF`, 199 calls | 2.233 s | 9.366 s | 8.589 s |
| Tracer antialias Fourier pair | 3.070 s | 1.372 s | See receipt |
| Particle flux/sampling | 0.125 s | 2.144 s | 1.820 s |
| RHS array packing, self | 5.860 s | 6.291 s | 5.643 s |
| ode78 arithmetic, self | 4.505 s | 4.715 s | 4.550 s |
| Output delivery | 1.465 s | 0.480 s | 1.655 s |

The compiled nonlinear kernel saves 19.986 s while tracer flux loses 21.399 s. Particle sampling adds another 2.019 s relative to builtin. This explains the near-parity despite much faster coefficient-only dynamics. Solver and RHS packing costs also rise substantially when a full-grid tracer is integrated: compiled endpoint self times are only 1.340 s and 1.378 s respectively.

In the compiled profile, continuous-extension setup costs 1.368 s, including its four RHS calls; interpolation costs 0.215 s. Output delivery includes 0.299 s in NetCDF `putVar`. These costs are too small to explain the 47.853 s profile increment over coefficient-only integration. Do not sum continuous-extension RHS work again with tracer/nonlinear rows. Output timing varies across these individual runs and is not a file-I/O regression qualification.

## Why the derivative path is expensive

[WVRetainedHorizontalOperator](../CompiledKernel/src/WVRetainedHorizontalOperator.cpp) prepares a single full-spectrum horizontal plane when a streaming retained plan exists. `spatialDerivative` then serially copies, transforms, multiplies and copies back each of 129 planes. Each derivative performs 129 forward and 129 inverse **2-D** FFT executions. At 199 RHS calls, `diffX` plus `diffY` account for **102,684 plan executions**, inferred from the executed source path and call counts. The full-grid derivative loop does not use the retained transform's twelve-worker schedule; its FFT plans have one internal thread. The builtin [x](../FastTransforms/@WVFastTransformDoublyPeriodicMatlab/diffX.m) and [y](../FastTransforms/@WVFastTransformDoublyPeriodicMatlab/diffY.m) methods transform only the differentiated axis.

The native sample corroborates this boundary: after merging alternate labels for the same OS thread, the interpreter has 3,646 wall-stack observations. Horizontal derivative leaves comprise 708 FFTW and 527 other arithmetic/copy observations, 33.9% together. Vertical calculus comprises 404 observations: 300 arithmetic/other and 104 BLAS/GEMM. Samples include waits, cover only part of integration, and are not CPU-time fractions or speedup measurements. Separate worker-thread totals must not be added to that denominator. The [sample receipt](../.github/ci-evidence/hydrostatic-composite-audit/compiled-native-sample-summary.json) preserves thread labels and caveats.

Vertical calculus also has a clear secondary cost. Builtin `diffZF` forms a real physical derivative matrix and applies it once. The compiled [stratified calculus](../CompiledKernel/src/WVStratifiedVerticalCalculus.hpp) clears and packs complex buffers, applies projection and reconstruction, computes the zero imaginary lane, normalizes and scatters. It operates in chunks sized to retained horizontal modes, despite receiving real physical columns.

Seven warm isolated samples give these medians: compiled `diffX` **44.663 ms**, builtin `diffX` **3.564 ms**; compiled `diffZF` **44.086 ms**, builtin `diffZF` **11.574 ms**. For the same vertical input already in column order, a prepared real derivative-matrix multiply takes **4.047 ms**, while two real factorized multiplies take **6.050 ms**. These narrower multiplication timings exclude layout changes and public boundaries. Maximum absolute vertical difference against builtin is 2.07e-17 for compiled and 2.19e-17 for real factorization. Every primitive sample is retained in the [receipt](../.github/ci-evidence/hydrostatic-composite-audit/primitives.json).

## Counters and copy limits

Compiled composite integration executes 995 public primitives: 199 each of x derivative, y derivative, vertical derivative, Fourier forward and Fourier inverse. It adds 199 nonlinear producers, 402 total producers, 800 physical reconstructions and 597 cache hits, with **zero duplicate executions**. End totals and deltas are in the [profile receipt](../.github/ci-evidence/hydrostatic-composite-audit/profile-evidence.json).

State-input copying increases by **9.521 GB** during composite integration, versus **0.700 GB** at the coefficient endpoint. Returned output accounting increases by **108.469 GB**, versus **9.101 GB**. Returned bytes are not a measurement of distinct `memcpy` operations or elapsed allocation time. The MEX's 43.813 s self time includes numerical execution, so it must not be labeled adapter overhead. Field-service storage grows by 270.533 MB as the composite's physical fields are materialized; kernel, engine and state-storage capacities remain unchanged. Three field requests per RHS account for 597 plan preparations. The profile supports later inspection of particle field requests and state traffic, but they are secondary to derivative execution.

## Prior investigations and changed assumptions

- [#463](https://github.com/JeffreyEarly/wave-vortex-model/issues/463), adopted through #465/#466, deliberately replaced full-depth derivative scratch with one plane to bound memory while preserving every resolved frequency. Keep the bounded-memory objective; do not undo it by blindly restoring full-volume staging. Its original full-frequency acceptance requirement is historical and is superseded for the proposed work by the owner’s retained-spectrum clarification above. Its original QG 3.385% timing miss remains historical evidence.
- [#338](https://github.com/JeffreyEarly/wave-vortex-model/issues/338) rejected a generic MATLAB tracer horizontal-sharing rewrite because vertical calculus and antialiasing dominated that earlier MATLAB profile. The changed assumption here is measured: the later compiled variable-family derivative path is serial and dominates horizontal tracer cost. This does not justify reviving broad tracer fusion.
- Spectral-kernel-benchmarks commits `230344b`, `bc69d7a`, `c9ef7db` and `49e3e54`, and closed issues #19–25, establish retained 2-D FFT pruning, tile-16 streaming, direct family views and persistent workers. No prior axis-only full-grid `diffX`/`diffY` production experiment was found. Reuse the selected infrastructure; do not repeat rejected inverse zero-fill/preserved-input experiments or the packing crossover.
- WVM [assembly](HYDROSTATIC-ASSEMBLY.md), [matrix/advection](MATRIX-ADVECTION-PIPELINE.md), [tiled advection](TILED-ADVECTION-PIPELINE.md) and [prepared pipeline](prepared-pipeline/RESULTS.md) reports optimize the coefficient graph. They do not establish that the raw tracer derivative path is equally optimized.

## Diagnostic result, verification and next experiment

The two-line diagnostic temporarily routed only public `diffX`/`diffY` to builtin MATLAB. Its patch is archived and **the runtime source was restored**. Tracer time fell 16.959 s, and the measured horizontal boundaries fell 15.935 s. The latter alone predicts 44.880 s from the baseline profile, a 26.2% reduction. The observed 42.607 s is 29.9% lower; roughly 2.3 s of that difference is outside the targeted boundaries and must not be credited to the proposed native optimization. Eliminating both original horizontal derivatives entirely gives an upper bound of 30.3% less total time, or 1.43× speedup, at that measured boundary.

All three complete graph comparisons and all three endpoint/composite coefficient checks passed the existing combined tolerance **1e-13 + 1e-12 × reference scale**. Saved fields, tracer, particles, metadata and output times agree. All composite profiles use identical control/tolerance hashes, 15 accepted steps, zero rejected steps, 199 RHS calls and record counts 2/5/5/5. The largest absolute hybrid difference is 1.819e-12. Raw relative errors can exceed 1e-12 for small coefficients while satisfying the unchanged combined tolerance; see [numerical comparisons](../.github/ci-evidence/hydrostatic-composite-audit/numerical-comparisons.json).

Verification ledger: all five profiles completed; the primitive timing completed; final output/control checks passed; Code Analyzer reports only three inherited benign unused-value notices in the profiling wrapper. Preliminary script-name, missing-module, profiler-export and verification-structure failures are retained separately, including outputs already produced. A late native sample missed the completed process and provides no additional sampling evidence. No timings were repeated or discarded. No runtime change is proposed for adoption, and no new CI or publication campaign was run. Canonical/generated website files are unchanged from the successful checks at the audit base.

For implementation, first test retained-spectrum tracer derivatives using the existing horizontal transform and reconstruction infrastructure, sharing the forward tracer transform where useful. Make the public projection contract explicit. Test both directions, derivative orders 1–4, odd/even grids, retained modes and the intended rejection or projection of excluded modes, real/complex public behavior where supported, aliases/strides, warmed allocation and lifecycle rules, and downstream tracer/control/output agreement. Audit arbitrary tracer initialization and antialias-disabled configurations separately from this retained, antialiased benchmark. Qualify the reviewed native candidate with frozen builds and three fresh processes per coefficient/composite configuration on an idle host, preserving all samples and investigating changes above 3%. The shared kernel also serves standalone and other stratified families; measure those consumers before claiming their benefit. A real vertical operator can follow as a separate increment.

Raw evidence is retained under `OceanKitRepositories/wave-vortex-model-benchmark-artifacts/hydrostatic-composite-audit-20260913`, including fixtures, every output, profiles, native sample, scripts, temporary patch, build/provider identities and failure logs. Existing frozen archives and published website numbers were untouched.
