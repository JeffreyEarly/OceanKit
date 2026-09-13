# Constant-stratification MATLAB vertical calculus

Investigation for [issue #513](https://github.com/JeffreyEarly/wave-vortex-model/issues/513), based on v4.4.0 `9345b79efa441d0ba097f931e1bebf7bb42cc089`. This is follow-up evidence, not a replacement for the frozen #501 publication campaign.

## Provenance and limits

The qualified source `024add423bf19c682ab5fd793821f398c5e7170f` and v4.4.0 have identical kernel, adapter, scientific MATLAB, and benchmark-worker sources. The original receipt archive is present and has SHA-256 `2ef32106c44a607bbe9f0d6aa94f4ba51a06b49cfa1fac1dbe2f4c606ae045a0`. Its original NetCDF fixture payloads, `matched-constant-nonhydrostatic-coefficient-endpoint-model.nc` and `matched-constant-nonhydrostatic-composite-dense-output-model.nc`, are missing from the workspace and `/private/tmp`; the receipt supplies no recoverable absolute payload path. The unchanged recipe reconstructs fixtures with exactly the recorded byte sizes and all recorded initial-condition diagnostics. New fixture hashes are endpoint `12db1772e88177348753e4c6b7fe40bcb1a1a4434f56c67e704f8ecf80cc5ab5` and composite `fc7514fbfa5df87f28c689a7c30ada39085c08f79bc01e1f842de8205e1cd9b1`; these are scientifically reconstructed fixtures, not claimed byte-identical originals.

The original compiled composite samples remain **114.204873, 55.298439, 55.641446 s**. Their pre-integration RSS is 5.518–5.521 GB and peak integration RSS 13.404–13.407 GB, so their timing spread does not coincide with a larger retained working set. Original FFTW plan descriptions and host contention measurements were not archived; no retrospective causal assignment is possible.

## Measured baseline profile

Apple M4 Max, MATLAB R2025b Update 4, native FFTW 3.3.11 NEON/pthreads, FFT budget 16, horizontal/pointwise workers 12/12. The reconstructed composite uses `[256 256 129]`, domain `[150000 150000 1300]` m, N² = 2e-5 s⁻², latitude 45°, the qualified deterministic physical state, ode78 to 7168 s with initial step 418.31529239621887 s, RelTol 1e-3 and component-scaled AbsTol 1e-6. Dense outputs are at 0:32:128 s. MATLAB profile and a separate ten-second native sample are diagnostic evidence, not timing qualification.

The profiled integration took 55.435990 s, 12 accepted steps, zero rejected steps and 160 RHS calls. Output record counts are coefficients 2, dense 5, particles 5, tracers 5.

| Boundary | Calls | Total seconds | Interpretation |
| --- | ---: | ---: | --- |
| Coefficient nonlinear flux | 160 | 19.212 | Includes compiled nonlinear kernel |
| Tracer flux | 160 | 22.336 | Includes the next four primitive rows |
| Vertical `diffZF` | 160 | 17.722 | 16.476 s in callees, 1.246 s wrapper self time including permutation |
| Horizontal `diffX` | 160 | 1.142 | Compiled full-grid derivative |
| Horizontal `diffY` | 160 | 1.211 | Compiled full-grid derivative |
| Tracer antialias forward/inverse | 160 each | 1.243 | Horizontal Fourier pair |
| Particle flux/sampling | 160 | 2.921 | Includes field requests and interpolation |
| RHS array packing and assembly | 160 | 4.714 self | `fluxArray` outside its callees |
| ode78 | 1 | 3.966 self | Solver arithmetic outside callees |
| Output delivery | 5 | 0.508 | Includes 0.244 s NetCDF `putVar` time |

Inclusive rows overlap; do not sum the tracer row with its derivative/antialias rows. The MEX gateway's combined self time is 41.515 s across 1774 calls, but this includes numerical execution and cannot be labeled allocation/copy time.

End-of-run counters: 800 primitives (160 each vertical derivative, x derivative, y derivative, forward antialias and inverse antialias); 160 nonlinear producers; 329 total producers (324 added during integration); 668 physical reconstructions (656 added); 480 cache hits; zero duplicate executions. State input copying totals 7.841 GB (7.701 GB during integration). Output bytes total 87.400 GB (87.265 GB during integration); this is a returned-byte ledger, not a claim that every byte is a separate memcpy. Persistent kernel/engine/field-service storage does not grow during integration.

The native MCR interpreter thread has 6488 sampled wall stacks. Inclusive attribution finds 1654 inside vertical calculus (25.5%), 1074 inside nonlinearFluxImpl (16.6%), 594 in memmove (9.2%), and 52 beneath mxCreate (0.8%). These categories can overlap. Within vertical calculus, 1124 samples are in FFTW execution and 530 in other column-loop work; no `fftw_spawn_loop` frame occurs under vertical calculus in this sample. FFTW synchronization occurs elsewhere. Samples include waits and are not CPU-time measurements; worker-thread totals are not added to the interpreter denominator.

## Prior work and bounded hypothesis

Reviewed spectral-kernel-benchmarks history, closed issues #8, #13 and #19–25, and their published decision records. The existing constant nonlinear kernel already uses retained batched DCT-I/DST-I with qualified endpoint/normalization rules. Prior WVM #94/#160 rejected extrapolating primitive FFT or explicit packing wins to complete workloads. The present target is the later MATLAB primitive column loop, separate from the already optimized nonlinear flux graph; no new provider, generic tracer rewrite, or full-volume staging is justified.

An exact-layout native probe at Nz=129 (DST length 127) shows that a requested budget of 16 can select `rdft-thr-vrank>=1-x2/1` for a single complex column's two real lanes. For 65536 same-kind round trips, exploratory DCT/DST timings are 0.691/0.689 s with that threaded plan, versus 0.069/0.073 s with one requested thread. At 256 columns per batch and budget 16, the corresponding times are 0.0307/0.0309 s. These repeated-buffer diagnostics exclude MATLAB boundaries and most calculus arithmetic. A separate baseline integration then took 120.206415 s. A diagnostic MEX inspected the provider's already-created plans through `FFTW_WISDOM_ONLY`, without measured replanning: the DCT was `rdft-thr-vrank>=1-x2/1`, while the DST was serial `rdft-buffered-127-x2/2-7`. The provider's `threads/rdft-vrank-geq1.c` confirms that `x2` prints the actual worker count. Of 7027 interpreter-thread wall stacks, 4908 (69.8%) were within vertical calculus, 3497 (49.8%) within its `fftw_spawn_loop`, and 2894 (41.2%) within its condition-variable wait. Those wait frames are nested subsets, not additive categories. This directly demonstrates short-plan synchronization in a newly reproduced slow execution. The fast profile had no such vertical-calculus dispatch frames. Plan selection and synchronization explain a major part of the newly observed bimodality; the original archive lacks the plan evidence needed to assign its 114.205 s sample retrospectively.

Eliminating all 17.722 s of `diffZF` would cap the observed whole-workload improvement at 1.47× (32% less time). Eliminating only its 16.476 s callee boundary caps it at 1.42×. The prototype reuses the existing interleaved-complex type-I plans and arena in batches of at most 256 columns, with zero-padded private tails. Scalar-column projection/reconstruction plans and coefficient nonlinear execution remain unchanged. Two plans are added; no full-grid scratch is added. A full-grid derivative makes 512 plan invocations instead of 131072. Execution metrics count actual plan invocations, consistently with other batched kernel operations.

Go/no-go: numerical/control/output agreement at existing tolerances, at least 10% lower composite median in three fresh processes per frozen configuration, and no coefficient endpoint regression over 3%. Compare every sample and investigate any regression; do not use this profile or the primitive probe as adoption qualification.

## Candidate profile and numerical verification

The candidate profile took 49.245032 s. `diffZF` fell from 17.722 to 12.611 s (28.8%); tracer flux fell from 22.337 to 16.923 s. Nonlinear flux remained 19.412 s. Whole-profile time fell 11.2%, consistent with the measured target share. Solver work, tolerance hashes, primitive/producer/copy counters and output counts are unchanged. Explicit kernel/engine storage increases by 104 bytes; the existing arena is reused. FFTW's internal plan allocation remains opaque.

The combined native suite passed 59/59 tests with Apple Clang. Eight focused MATLAB tests passed, including first through fourth F/G derivatives and integrals, real and complex broadband volumes, constant hydrostatic/nonhydrostatic layouts with antialiasing enabled/disabled, matrix-column boundary cases 0/1/255/256/257, and downstream tracer evolution. Native tests also cover the arena-limited 119/120/121-column boundary. All 25 full-output, repeatability and endpoint/composite-control comparisons passed the existing combined absolute/relative tolerance (1e-13 + 1e-12 × reference scale); the largest absolute difference was 2.843e-14. Metadata and output times agree exactly. No tolerance was loosened.

Code Analyzer found only three benign unused-variable warnings in the diagnostic wrapper (a label consumed through `evalc` and explicit owner cleanup). The 52 CI-policy tests, repository/package/generated-artifact checks and `docs:check` passed. Local GCC 14 cannot parse this host's Apple SDK `_Alignof` headers, so Linux GCC CI is a separate pending portability gate; the local GCC attempt is retained and is not reported as passing.

## Candidate qualification

The baseline and candidate were frozen at `88c472639b29fe5a7761ee72ee6201d75f6a7259` and `c3067ad5e294c0bff401a3c2895b75eecc771f6f`. Source aggregates are `3d565dff3d0e6881f47c206b853efc301fd4a0fd82976e3c89068268214a6123` (identical to the original qualified baseline) and `7cfa4cc82f7f1daec6afcc0b870be8210e4ba6aa4656f15998c79ef2adeea92a`. Frozen MEX hashes are `4ec1f366662a9aabcb712c7cdcec191128d6abee2b50f9eba0987e47df3d3a3d` and `be0140d6edad9cf7dd6f77f46b45ad4aa1187185b7fcb77e96e304c05232a720`. Both use the same qualified FFTW provider library hashes. Each integration is a fresh MATLAB process; workloads alternate order within each configuration. Configuration blocks were baseline then candidate, so this is not a randomized crossover campaign.

The twelve initial timing receipts also contain a **post-cleanup provider-inspection failure**: after deleting its compiled owners, the wrapper restores the runner's entry path, exposing a legacy v5 MEX. This did not execute during integration. Creation validates the cached provider, activation verifies its exact path and SHA, and the wrapper asserts the active object's validated native identity before timing. The later inspection reports that the expected cached module is still loaded and fails only the subsequent `which`-path check. Original failure receipts are retained. Supplemental inspections retain the qualified module directory on the path and validate normally. Neither the v5 checkout nor its legacy binary was changed.

Pre-run process snapshots show no competing MATLAB, native benchmarks or builds. Normal desktop services remained active. The initial third candidate pair is retained separately because its snapshots show macOS media analysis consuming a core. It is excluded from the idle-host qualification set for that reason, decided before the supplemental runs, rather than because of its timing. An interval `top` check before the supplement measured 96.97% host idle; both supplemental pre-run snapshots had no process above 40% of one core.

| Configuration/workload | Selected fresh-process samples (s) | Median (s) |
| --- | --- | ---: |
| Baseline coefficient endpoint | 14.896217, 14.739576, 14.650868 | 14.739576 |
| Candidate coefficient endpoint | 14.460604, 14.659385, 14.626509 | 14.626509 |
| Baseline composite | 117.510934, 116.893989, 117.710799 | 117.510934 |
| Candidate composite | 49.466165, 48.584709, 48.829551 | 48.829551 |

Retained candidate repeat 3: endpoint **15.048656 s**, composite **48.805859 s**. Including it leaves the composite conclusion unchanged. The selected composite median falls 58.45% (2.407×) relative to the reproduced slow baseline and 12.24% relative to the historical fast median of 55.641446 s. The coefficient median falls 0.77%; even the excluded third endpoint is only 2.72% slower than its baseline repeat, below the 3% regression threshold. There is no measured coefficient regression requiring an algorithm change. Three selected candidate composite samples span 48.585–49.466 s (1.81% of the median); baseline slow-mode samples span 116.894–117.711 s (0.70%). All fourteen samples have identical per-workload solver controls, tolerance hashes, work counts, output counts and integration counter deltas. All 25 output/control comparisons pass.

The fourteen retained timing processes consumed 18.10 minutes of process lifetime, of which 10.85 minutes were measured integration. Profiles, preparation, tests and comparison runs are separate. An exact active-time split between implementation, review and waiting was not measured.

## Durable evidence and reproduction

The compact [qualification receipt](../.github/ci-evidence/issue-513/qualification.json) retains all sample times, solver/counter records, source and binary hashes, inspection failures and selection reasons. The [output comparisons](../.github/ci-evidence/issue-513/output-comparisons.json), [profile functions](../.github/ci-evidence/issue-513/profile-functions.json), [plan wisdom](../.github/ci-evidence/issue-513/baseline-plan-wisdom.json), and [verification ledger](../.github/ci-evidence/issue-513/verification.md) are committed with this report.

Raw payloads, profiles, sampled stacks, logs, exact diagnostic/runner scripts, source tarballs and binaries are frozen in `wave-vortex-model-benchmark-artifacts/issue-513-vertical-calculus/20260913-v4.4.0/` relative to the OceanKitRepositories workspace. The [archive index](../.github/ci-evidence/issue-513/archive-index.json) records its manifest hash and file count. This is a new local raw archive; it is not hosted in the source repository. The original publication archive remains unchanged.

`issue513Prepare` reconstructs the two named fixtures and refuses to overwrite existing files. `issue513Profile` copies one fixture into a fresh output, validates the compiled provider, and captures integration controls/counters; its optional profiler uses `OceanKit/tools/profiling`. The archive's `runSamples.py`, `runIdleSupplement.py`, `qualifyPairs.m`, `qualifySupplement.m`, `issue513PlanWisdom.cpp`, and sampling extractors preserve the exact executed protocol. Use isolated v4 source roots and the same qualified provider; keep its module directory on the MATLAB path when inspecting identity after wrapper cleanup. The supplied binary archives record the original installation paths and are not relocatable distributions.

## Recommendation

Adopt the bounded vertical-calculus batching change after required portability/integration CI. It addresses measured short-transform scheduling overhead while preserving the existing numerical rules and qualified provider. The conservative expected composite benefit is about 11–12% relative to the already-fast baseline, with much larger improvement when the baseline selects a threaded short DCT. Do not describe that larger ratio as a universal adapter speedup or as a replacement publication result. Three candidate processes cannot guarantee the absence of all future plan/runtime variability.

The residual approximately 49 s composite cost is still far above the original 20.995 s standalone median. Even eliminating the entire measured vertical derivative cannot close that gap. Tracer arrays, MATLAB solver/RHS assembly, repeated primitive boundaries and substantial returned-byte traffic remain the next evidence-led targets if a larger refactor is separately authorized. Output delivery is too small to justify prioritizing it here.
