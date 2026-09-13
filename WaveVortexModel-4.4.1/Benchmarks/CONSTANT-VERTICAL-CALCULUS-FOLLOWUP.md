# Constant vertical calculus: multiplier and real-transform follow-up

This second increment for [issue #513](https://github.com/JeffreyEarly/wave-vortex-model/issues/513) extends [PR #519](https://github.com/JeffreyEarly/wave-vortex-model/pull/519). It compares against the already-batched implementation, not the original slow short-transform mode. The [first investigation](CONSTANT-VERTICAL-CALCULUS.md) preserves the original workload breakdown, prior spectral-kernel-benchmarks decisions, original samples and batching qualification. Neither increment replaces the frozen publication results.

## Hypothesis and measured stages

After batching, the composite profile still spent 12.611 s in 160 public `diffZF` calls, versus 1.085 s in `diffX`. The mathematical transform lengths alone do not explain that gap: DCT-I length 129 and DST-I length 127 both have logical FFT length 256. The vertical implementation also rearranges MATLAB arrays, packs columns, normalizes coefficients, applies derivative factors and dispatches batched transforms. The horizontal primitive follows a different prepared layout.

An instrumented build of the batched baseline measured the following native stages over the same 160-call composite workload. Timers surround each batch stage; their overhead and the profiler make this a diagnostic, not qualification.

| Native stage | Seconds | Share of measured native stages |
| --- | ---: | ---: |
| Pack columns | 0.724 | 6.2% |
| Forward transform | 2.346 | 20.1% |
| Forward normalization | 0.413 | 3.5% |
| Derivative/integral multipliers | 5.179 | 44.4% |
| Inverse normalization | 0.281 | 2.4% |
| Inverse transform | 2.403 | 20.6% |
| Unpack columns | 0.317 | 2.7% |

The multiplier loop recomputed the same mode-dependent factors for every column. Its 5.179 s is comparable to both transforms together, 4.749 s. The transform scratch also carried an imaginary lane filled with zeros even though the native primitive accepts a real field. These are two bounded sources of avoidable work. Removing all of the baseline public derivative would cap the 49.245 s profile's gain at 1.34×; removing only the measured multiplier stage would cap it at approximately 1.12×, before accounting for diagnostic overhead.

The native wall-stack sample of the batched diagnostic includes 641 FFTW synchronization/child-work frames among 1544 vertical-calculus frames. Those include waits for useful parallel work and are not all avoidable overhead. The first increment separately demonstrated pathological two-worker dispatch for individual short columns. This increment retains the batched threading policy and does not repeat that experiment.

## Bounded implementation and exploratory results

The candidate computes the original factors once per primitive call, preserving the original `pow` arithmetic, signs and mode indices. It uses the existing qualified FFTW provider to prepare real DCT-I/DST-I batches, with at most 256 columns and the existing scratch arena. No full-grid allocation is added. The per-call factor vector is O(Nz), with allocation failure handled by the kernel status contract. Endpoint rules, normalization order, integral bottom subtraction and the public complex-input split remain unchanged. Coefficient nonlinear execution, scalar projection/reconstruction and other transform families are unchanged.

The new [primitive diagnostic](issue513VerticalPrimitives.m) measures one untimed warmup and five retained calls on the same real `wvt.u` fixture. Each variant uses the same qualified provider libraries. These are isolated within-process medians, not complete-workload speedups.

| Variant | Public `diffZF` median (s) | Change from batched baseline |
| --- | ---: | ---: |
| Batched baseline | 0.080278 | — |
| Multipliers computed once | 0.052408 | −34.7% |
| Real-only transform batches | 0.068649 | −14.5% |
| Both changes | 0.039470 | −50.8% |

For the combined candidate, `diffX` is 0.007271 s. Input/output MATLAB `permute` medians are separately 0.005894/0.001072 s. They reorder `[Nx Ny Nz]` into vertical-column storage and back; they do not permute the physical coordinates or change the calculus. These separate timings cannot be added to profiled integration costs as an exact decomposition. Public `diffZF` remains about 5.4× slower than `diffX` in this diagnostic, down from about 10.5× in the isolated batched baseline.

The combined whole-composite profile falls from 49.245 to 43.717 s (11.2%). Public `diffZF` falls from 12.611 to 6.431 s (49.0%), including 1.241 s of MATLAB wrapper self time. The full output comparison passes with maximum absolute difference 1.3323e-15, confined to tracer roundoff; coefficients, Eulerian fields, moorings, particles, metadata and times agree exactly. Profiling and native sampling are disabled for qualification.

| Profile boundary | Batched baseline (s) | Combined candidate (s) |
| --- | ---: | ---: |
| Coefficient nonlinear flux, inclusive | 19.435 | 19.221 |
| Tracer flux, including derivatives | 16.923 | 11.084 |
| Vertical `diffZF`, included in tracer flux | 12.611 | 6.431 |
| Particle flux/sampling, inclusive | 2.839 | 2.942 |
| RHS `fluxArray`, self | 4.562 | 4.736 |
| `ode78`, self | 3.755 | 3.955 |
| Output-file delivery, inclusive | 0.362 | 0.392 |
| NetCDF `putVar`, included in output | 0.165 | 0.177 |

Inclusive rows overlap. The isolated profiles show small increases in several unchanged boundaries, so the whole-profile reduction is smaller than the derivative reduction alone. They are not three-process regression tests. MEX self time falls from 36.518 to 30.273 s, but includes the numerical kernels and cannot be treated as allocation/copy time. State-copy and returned-byte counters remain exactly equal: the optimization reduces repeated arithmetic and transform work without reducing MATLAB's full-volume traffic. End-of-run counters remain 800 primitives, 160 nonlinear producers, 329 total producers, 668 physical reconstructions, 480 cache hits and zero duplicate executions. Integration still adds 7.701 GB of state copies and 87.265 GB of returned output bytes. The latter is a returned-byte ledger, not a count of separate memcpy operations. Full counters, self/inclusive boundaries and native sample classifications are in the committed evidence.

## Numerical verification

The combined native suite passes 59/59 tests with a warning-clean Clang build. Fourteen final focused MATLAB methods pass, covering derivatives one through four, F/G integrals, real/complex inputs, supported constant hydrostatic/nonhydrostatic layouts, antialiasing on/off, batch boundaries and downstream tracer evolution. Twelve unchanged successful methods were reused; only the changed/new two methods were rerun. No tolerance was loosened.

A newly attempted even-grid low-mode fourth-derivative stress case exposed a pre-existing difference from MATLAB's matrix derivative: approximately 1.43e-12 for F and 1.89e-12 for G, identically on the frozen batched baseline and candidate. All failed diagnostics are retained. A separate direct low-mode probe finds bitwise equal baseline/candidate outputs. Existing odd-grid test inputs are preserved. The added even-grid regression instead compares both implementations with independently differentiated retained modes 6 and 3 at the existing 1e-12 tolerance; even-grid matrix-integral coverage is also retained. This avoids presenting the pre-existing low-mode discrepancy as either a new regression or a passing test.

The source-selection kernel fingerprint is refreshed, correcting the stale receipt found by the first PR's MATLAB CI. Its focused source-hash contract passes. Code Analyzer reports only two benign unused-variable warnings for explicit cleanup in the diagnostic harness. All 52 CI-policy tests pass. The first increment's successful documentation-generation check is reused because this increment changes no canonical generated documentation.

## Frozen qualification protocol

The baseline is commit `81692f49a98ceb7b5ef75f345fb5567494eb8b86`, runtime aggregate `7cfa4cc82f7f1daec6afcc0b870be8210e4ba6aa4656f15998c79ef2adeea92a`, MEX SHA-256 `be0140d6edad9cf7dd6f77f46b45ad4aa1187185b7fcb77e96e304c05232a720`. The combined candidate is commit `8a1840a3129a3aed62efd8e4ac6277bdf063e873`, aggregate `209a99eeda112d7b4106821c859eba619a8f2745f44489dc4ed6fb2d8af81700`, MEX SHA-256 `414f573ce55b427dbb23d48f32b7f5bf986a9c921a8783d31986ebf83100e189`. Both use the identical qualified FFTW base/thread libraries, FFT budget 16 and horizontal/pointwise worker counts 12/12 on the same Apple M4 Max and MATLAB R2025b Update 4.

The original qualified NetCDF endpoint/composite fixture payloads remain missing. This increment reuses the exact scientifically reconstructed fixtures and hashes documented in the first report; it does not reconstruct them again or claim byte identity with the missing originals. The v5 checkout, other experiments, original archives and published website numbers remain untouched.

The adoption criterion is at least 10% lower composite median against the batched baseline, with no coefficient endpoint regression above 3% and unchanged numerical/control/output gates. Each workload/configuration receives three fresh MATLAB processes, alternating baseline/candidate order across repeats and reversing workload order on repeat two. Every process has a retained interval host check requiring at least 94% idle and no nonexempt process over 40% of a core before launch. All failed checks and post-run snapshots are retained.

An initial baseline endpoint integration completed in 14.690308 s, then its post-run inspector failed after harness cleanup restored the legacy MEX path. The frozen identity was validated before integration. This extra sample and all failure receipts are retained separately. The corrected runner explicitly restores the qualified module directory before post-run inspection; timed code is unchanged. The final qualification set requires successful pre/post identity checks and process completion.

The initial twelve processes completed successfully, with valid pre/post source and MEX identity. A subsequent host audit found post-run activity for baseline composite repeats 2 and 3 (91.36% idle with photoanalysisd at 97.3% of one core; 86.48% idle with backup/cloud/system activity), and a borderline 93.99% idle reading after baseline endpoint repeat 3. Their times are retained below. Before running any supplements, those three were excluded from the stricter idle-host cohort based on the host snapshots, not their timing. Three baseline-only supplements require three consecutive good prechecks and a good postcheck, using the same 94% idle / 40%-of-one-core thresholds. All passed. The selected candidate samples all passed the same final-interval postcheck thresholds. These are boundary snapshots, not continuous monitoring of contention throughout integration.

| Workload | Selected baseline samples (s) | Candidate samples (s) | Median baseline → candidate |
| --- | --- | --- | --- |
| Coefficient endpoint | 14.609002, 14.912715, 14.376724 | 14.378576, 14.731212, 14.282991 | 14.609 → 14.379 |
| Composite | 49.203848, 48.845695, 49.119515 | 41.478405, 41.769560, 41.825933 | 49.120 → 41.770 |

The composite median falls 14.96% (1.176×); the endpoint median falls 1.58%. No corresponding selected sample has a coefficient regression above 3%. The selected composite spans are 0.73% of the baseline median and 0.83% of the candidate median. The endpoint spans are 3.67% and 3.12%, respectively; keeping each observation matters more than reporting extra decimal places in the medians.

The retained original baseline endpoint repeat 3 is **14.694530 s**; composite repeats 2 and 3 are **49.171966 and 49.344479 s**. Including them in the original twelve-process cohort gives approximately 15.11% lower composite time and 2.15% lower endpoint time. Thus the conclusion does not depend on their exclusion. Together with the initial 14.690308 s inspector-failure sample, all sixteen integration observations are retained. No candidate sample was discarded. Initial order was alternating; the final qualified cohort includes baseline-only supplements and is not a randomized crossover campaign.

All 29 full-output, repeatability and endpoint/composite-control comparisons pass at the existing combined tolerance, 1e-13 + 1e-12 × reference scale. Metadata and output times agree exactly. The maximum absolute difference, 1.81899e-12, occurs in particle positions and also appears between candidate repeats. The maximum relative difference, 1.06219e-12, is a coefficient comparison with only 1.23512e-14 absolute difference; it passes the unchanged combined tolerance. This is separate from the direct profiled-pair comparison above. Complete per-workload solver controls, tolerance hashes, work counts, runtime counters and counter deltas are exactly equal across all successful timing runs; the initial inspector-failure sample also has matching controls and counter deltas.

The sixteen timing processes consumed 16.61 minutes of process lifetime, including 8.12 minutes of integration. Idle-host checks/waits, profiling, implementation, validation and receipt refresh are separate. An exact active-time breakdown of those activities was not measured.

## Recommendation and evidence

Adopt the combined bounded change in PR #519 after required CI. It meets the predeclared 10% composite / 3% endpoint criterion, with approximately 15% less composite time against the already-batched implementation on this fixture and host. The isolated derivative is approximately 51% faster. Three candidate processes do not establish universal performance or eliminate all future FFTW/runtime variability.

The remaining `diffZF` profile cost is 6.431 s out of 43.717 s. Even removing it entirely would cap an additional whole-profile speedup at 1.17×. The MATLAB wrapper has about 1.24 s of self time, including permutations; removing that entire boundary would save only 2.8% of this profile and cannot close the standalone gap. Keep the present two changes narrow; a further tracer/solver/adapter redesign would need its own measured objective and authorization.

The [qualification receipt](../.github/ci-evidence/issue-513-real/qualification.json) retains all times, host selection, exact source/module/provider identities, controls and counters. The [profile evidence](../.github/ci-evidence/issue-513-real/profile-evidence.json), [stage timings](../.github/ci-evidence/issue-513-real/diagnostic-stages.json), [output comparisons](../.github/ci-evidence/issue-513-real/output-comparisons.json), [focused MATLAB results](../.github/ci-evidence/issue-513-real/final-focused-tests.csv) and [verification ledger](../.github/ci-evidence/issue-513-real/verification.md) preserve the diagnostic and numerical findings.

Linux release and sanitizer builds pass for the frozen runtime in [CI run 34774291513](https://github.com/JeffreyEarly/wave-vortex-model/actions/runs/34774291513). MATLAB smoke passes on R2025b and R2026a. Later full/sanitized shards identify a stale constant-kernel digest in the committed forward-integration receipts. A fresh native-capable runner/probe from commit 8a1840a3 passes the existing six lifecycle methods with both reference and native providers (twelve cases, 254.85 s summed method duration); the previously failing receipt contract then passes. The existing collector writes twelve new content-addressed fragments and refreshes the two current receipt pointers. Prior fragments and backed-up receipt pointers remain retained. No validator or runtime change was needed. The final evidence/receipt commit will receive its normal CI rerun. An initial lifecycle-helper invocation failed before tests because of its hyphenated filename; the retry used a valid MATLAB script name. That invocation log was overwritten by the retry, so the archive includes a clearly labeled reconstruction from the tool transcript. All timing and numerical failure receipts remain retained.

The new local raw archive is `wave-vortex-model-benchmark-artifacts/issue-513-vertical-calculus/20260913-v4.4.0-second-stage/`, relative to the OceanKitRepositories workspace. The committed [archive index](../.github/ci-evidence/issue-513-real/archive-index.json) and [manifest](../.github/ci-evidence/issue-513-real/archive-manifest.sha256) bind all raw outputs, profiles, failed diagnostics, native samples, source snapshots, frozen modules and exact executed scripts. It contains 326 files totaling 8,934,152,661 bytes; manifest SHA-256 is `ac05175ad0e5ed2a4aa5a1f79a3f5b14e383a53272febe36f197c92419871071`. Both source tarballs were checked against all 319 manifest source files, and the archived modules/provider libraries match the qualified hashes. The first-stage archive and original publication archive remain unchanged. This local archive is not hosted in the source repository.

Reproduction uses the frozen commits/fixtures/providers above, `issue513VerticalPrimitives` for isolated diagnostics, and the archived `runQualificationV2.py` / `runIdleSupplement.py` for the timing protocol. `compareQualification.m` and `compareSupplements.m` retain all 29 comparisons, and `aggregateQualification.py` checks complete source/provider/control/counter invariants before computing qualification. The archive retains the earlier failed runner as well. `refreshForwardIntegrationReceipts.m` records the separate required lifecycle refresh. Archived binaries retain their original installation paths and are not relocatable distributions.
