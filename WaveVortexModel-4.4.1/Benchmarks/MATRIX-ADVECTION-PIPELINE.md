# Matrix scheduling and immediate FFT advection consumption

This optimization follows the shared-gradient pipeline in PR #476. [Issue #477](https://github.com/JeffreyEarly/wave-vortex-model/issues/477) tracks adoption qualification. The owner requested both directions. The implementation, final paired performance campaign and source-linked consumer qualification are complete. Required hosted checks gate merging. The experiment does not change MATLAB science, checkpoint formats, tracer algorithms or a running experiment.

## Changes

Hydrostatic and Boussinesq can consume each completed inverse FFT derivative plane immediately in the existing x/y/z advection accumulation order. The native retained FFT implementation invokes a synchronous consumer before releasing its plane task; fallback implementations deliver normalized physical output. Cached derivative hits use the existing pointwise path. Derivative storage remains available, capture still follows successful production, and callbacks touch only private flux scratch before publication. This removes a separate dispatch and improves locality without adding a physical-grid buffer. Compact runner execution selects it for Hydrostatic and Boussinesq only; direct callers can select `fusedDerivativeAdvection` explicitly.

Boussinesq now projects the wave-F inertial contribution only at its exact retained `(0,0)` column. Other projection outputs are unused by the coefficient assembly. The public complete projection remains available, and `inertialOnlyProjection=false` preserves the oracle. Exact retained identities and matrix membership remain authoritative; inconsistent identity falls back to the complete projection.

Prepared vertical operators can execute independent matrix groups through one persistent kernel-owned worker pool. Scalar and Accelerate backends declare support for concurrent calls. Packing buffers, when required, are separate for each worker; the measured model matrices already use direct views. Workspace and pool reentry are rejected. The compact Boussinesq runner selects up to eight workers, bounded by available performance cores. Direct construction defaults to one. Accelerate and its environment remain vendor controlled; no global thread settings are changed.

The spectral-kernel-benchmarks history had already accepted outer-group split-real Accelerate scheduling with a single-thread vendor cap. This experiment tests its transfer to the current complete model while leaving the vendor environment unchanged. Previously rejected packing bridges, zero-reset schemes and multi-pass vDSP coefficient assembly were not repeated.

## Isolated screens

Frozen baseline: `wave-vortex-model-benchmark-artifacts/shared-gradient-pipeline-20260911/qualified-candidate-wave-vortex-run`. Raw archive: `OceanKitRepositories/wave-vortex-model-benchmark-artifacts/matrix-advection-20260911`.

| Change | Workload | Measured reduction |
| --- | --- | ---: |
| Immediate FFT advection consumption | EddyTide integration | 8.35% |
| Immediate FFT advection consumption | Larger Hydrostatic integration | 3.71% |
| Immediate FFT advection consumption | Boussinesq integration | 5.51% |
| Exact inertial-column projections | Boussinesq integration | 10.34% |
| Eight matrix workers versus one, after inertial optimization | Complete Boussinesq nonlinear RHS | 26.14% |

Integration screens used two alternating baseline/candidate pairs per profile. Matrix workers were screened in order `1,2,4,8,8,4,2,1`, with two warmups and eight RHS samples per process. The median of process medians was 123.24 ms, 106.23 ms, 94.51 ms and 91.03 ms for one, two, four and eight workers. Kernel preparation stayed near 2.35 seconds. Matrix-group producer count was 23,635 per RHS, corresponding to 47,270 real Accelerate GEMM calls, down from the previous 60,762. The worker screen used twelve horizontal and twelve pointwise workers; the production runner uses twelve horizontal and eight pointwise workers on this host. The combined screen tests that actual runner policy.

These reductions are measured at different boundaries and must not be added. All output/flux comparisons passed. Flux differences in the worker sweep were at most approximately 5.1e-16 relative L2; repeated serial controls also showed final-bit variation.

## Combined complete-model screen

Four alternating baseline/candidate pairs per fixture used the frozen runner and identical integration controls. All sixteen scientific output comparisons and state/step decisions passed, with zero duplicate evaluator executions. Timings below measure integration, excluding checkpoint loading, preparation and output. No warmup campaign or production-default qualification is claimed.

| Workload | Integration-time reduction | Added owned peak |
| --- | ---: | ---: |
| EddyTide | 7.80% | 366,240 bytes |
| Constant-stratification control | 1.42% | 32 bytes |
| Larger Hydrostatic | 3.38% | 212,384 bytes |
| Boussinesq | 38.67% | 584,280 bytes |

The constant-stratification control is unchanged scientifically; its small timing shift is treated as ordinary measurement variation. EddyTide owned storage grows 0.068%, with no added physical-grid derivative buffer. Generic vertical operators retain an exact column-to-matrix lookup. Boussinesq selects eight vertical workers in the actual runner and retains its established eight pointwise workers.

Photos background activity was present immediately before the run; the 36 recorded host samples during the campaign show it had settled. Preserve that caveat and repeat the formal idle-host campaign before adoption. The independent screens and the combined result agree on both optimization directions. Compact results and hashes are in [the evidence directory](../.github/ci-evidence/matrix-advection-exploration/combined-summary.json).

## Verification ledger

- Combined native Release suite: 58/58 passed. Subsequent test-only lifetime correction reran only the two affected kernels successfully.
- ASan/UBSan: seven affected kernel, runtime, FFT and policy tests passed after a test-only correction. Initial H/B tests declared producer counters after the kernel owners, causing use-after-scope during instrumented plan destruction. Counters now outlive the plans; failed and corrected logs are retained. No runtime fix was needed.
- GCC 14: operator and both family kernel tests passed. The local runner target cannot compile the installed Apple SDK Mach header assertions with GCC; actual Linux runner CI remains required.
- MATLAB: all seven affected Hydrostatic/Boussinesq scientific parity methods passed, including compact consumers and concurrent Boussinesq matrix groups. Four integration/source checks passed initially; the stale native-engine identity check passed on focused rerun after updating source-selection metadata. No MATLAB source changed.
- Independent review covered consumer lifetime, derivative cache hits, normalization, density correction, partial inertial output, exact group mapping, worker scratch, failure publication, overflow, reentry and the combined conflict resolution.
- Prerequisite PR #476 cache-capacity probe corrections were verified under Release and sanitizers and pushed separately. PR #476 merged as `5a41977b` and issue #475 is closed. The incoming main tree was byte-identical to `39afaa13`, already present in this experiment; merge `a16da09c` records that ancestry without changing the experiment tree. The v4 main checkout is synchronized.

## Delivery boundary

The exploratory compiled runtime was frozen at `e8806710`; the [adoption ledger](../.github/planning/issue-477-adoption.md) records subsequent qualification and a narrow invalid-worker-count error-boundary correction. The local branches and raw archive retain separate experiments, failed checks, provider/binary/input hashes, requests and numerical comparisons. The combined screen is exploratory: Photos background services used about one core immediately before it, and a process-activity journal records the measurement period. It does not replace formal idle-host qualification.

The final qualification below supersedes the exploratory timing screen. Required hosted CI gates merging. Retain the exact independent reconstruction/projection controls for future qualification; reprofile the qualified combination before selecting another larger refactor.

## Final adoption qualification

The frozen candidate is `400ed1f1`, compared with the exact preserved qualified PR #476 executable. Two warmup pairs and eight measured alternating pairs ran per fixture under `reuse`, using the existing immutable manifest, controls and scientific tolerances. All 32 measured scientific comparisons (plus eight warmup comparisons), state and integration decisions passed. Duplicate evaluator executions were zero. Executable, fixture, provider-library and source hashes were unchanged after measurement.

| Workload | Integration reduction | Paired bootstrap 95% interval for reduction | Process-lifetime reduction |
| --- | ---: | ---: | ---: |
| EddyTide, 256 × 256 × 28 | 7.00% | 5.81–8.41% | 6.79% |
| Constant-stratification control | 1.05% | −0.23–2.36% | 0.15% |
| Larger Hydrostatic | 3.65% | 3.00–4.29% | 2.98% |
| Boussinesq, 256 × 256 × 129 | 39.54% | 38.72–40.34% | 1.29% |

The unchanged constant control is consistent with measurement variation. Boussinesq's short continuation spends approximately 34 seconds in loading/preparation per process; the integration reduction is the relevant result for longer runs. None of the integration or process-lifetime ratios requires the protocol's greater-than-3% regression investigation. The manifest requires improvement for both EddyTide and Boussinesq, which both pass; the generic protocol's abbreviated prose names only EddyTide.

The host was idle at preflight and no competing build, MATLAB, profile or experiment workload ran. Ten-second process observations captured intermittent system-service activity, including macOS media analysis. These snapshots cannot exclude brief overlap with individual integration windows; preserve this qualification caveat rather than claiming perfect quiescence. The paired confidence intervals and the earlier independent screens support the measured direction and scale. No unrelated process or experiment was changed.

Added owned peak is unchanged from the exploratory screen: 366,240 bytes for EddyTide, 212,384 bytes for larger Hydrostatic, 584,280 bytes for Boussinesq and 32 bytes for the constant control. RSS is reported separately in the machine-readable summary. Low-memory timing and storage are not qualified in this campaign: `candidatePolicies` contains only `reuse`, and `lowMemoryPassed` is null. The generic protocol's low-memory storage text is therefore inapplicable. A separate low-memory replay covers scientific correctness and exact step decisions only.

Final source-linked qualification at `5947ac0d` passed all six families with both providers and the strict committed-receipts catalog check. Candidate `400ed1f1` adds those receipts and verification evidence without changing compiled inputs. The narrow invalid-worker-count factory correction at `5947ac0d` was separately checked under Release, ASan/UBSan and GCC. Subsequent adoption finalization changes documentation, evidence and source-selection metadata only.

See [the final decision](../.github/ci-evidence/matrix-advection-adoption/decision.json), [performance summary](../.github/ci-evidence/matrix-advection-adoption/performance-summary.json) and [archive hashes](../.github/ci-evidence/matrix-advection-adoption/archive.json). Raw outputs, protocols and profiles are retained under `OceanKitRepositories/wave-vortex-model-benchmark-artifacts/matrix-advection-adoption-20260912`.
