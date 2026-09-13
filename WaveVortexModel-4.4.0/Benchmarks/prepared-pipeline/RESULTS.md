# Prepared projection and assembly pipeline results

The selected candidate is commit `ed6049a2504dad833e93c9221b66ea68db1ac404`, based on the merged PR #494 baseline at `f65ebf684a403742ed27775de0a2c2aac4a1e2fe`. It keeps the existing FFT, vertical matrix, coefficient, phase, and cache representations. Boussinesq skips the unused horizontal-divergence modal write when vertical velocity is present, consumes the vertical-velocity divergence contribution in the final wave-coefficient loop, and runs the remaining projection, coefficient-assembly, and wave-plus-balanced addition loops on the existing synchronous pointwise executor. Hydrostatic runs its projected-field coefficient loop on that executor. Every worker range owns disjoint coefficients, and cache readiness is published only after a joined assembly or after both vertical products and the joined wave-plus-balanced addition complete.

Adopt the reviewed candidate. The final complete paired campaign passes the declared numerical, integration-control, producer, storage and runtime gates. EddyTide improves by 6.8%; larger Hydrostatic changes modestly, and Boussinesq has an order-sensitive improvement. Retain the host and timing qualifications below; these are measured incremental results, not perfectly isolated performance guarantees.

## Profile basis

The fresh Boussinesq baseline profile used the accepted native tiled pipeline on the `256 x 256 x 129` no-observer fixture. One resident kernel ran 120 timed RHS evaluations with a median of 85.66 ms while macOS sampling was active. Kernel preparation took 2.71 s after a 9.09 s checkpoint load. The profile attributed main-thread and summed active samples as follows:

| Category | Main thread | Summed active |
| --- | ---: | ---: |
| BLAS library | 26.66% | 30.20% |
| FFT library | 19.16% | 38.33% |
| Field assembly and derivative arithmetic | 13.66% | 2.74% |
| Tiled adapter and advection arithmetic | 11.71% | 23.95% |
| Wait | 7.99% | — |
| Phase preparation | 6.64% | 1.16% |
| Projection pointwise arithmetic | 6.47% | 1.13% |
| Other runtime and control | 3.93% | 1.16% |
| State validation | 2.83% | 0.49% |
| Matrix adapter and packing | 0.95% | 0.84% |

Main-thread samples help identify serial work and waits; summed active samples describe aggregate CPU work across workers. Neither alone establishes the full dependency critical path. They cannot be added or treated as interchangeable speedup budgets. The inclusive `projectSpectralFields` and `reconstruct` shares were 18.58% and 35.97%, respectively, but both include downstream vertical matrix work. The selected change only targets their pointwise portions.

## Resident RHS screens

Each Boussinesq stage used two already prepared native workers and the same fixture, compiler, Accelerate matrix backend, pinned FFTW libraries, one FFTW thread, twelve horizontal workers, eight pointwise workers, and eight vertical-group workers. The Hydrostatic screen used its EddyTide fixture and one vertical-group worker with the other execution settings unchanged. Eight counterbalanced blocks each contained two warmups and 24 timed fresh-scope RHS calls. Ratios pair neighboring opposite-order blocks; values below one favor the candidate.

| Screen | Cumulative change | Geometric mean ratio | Reduction | Four paired ratios | Maximum flux difference | Maximum field difference | Reported kernel bytes |
| --- | --- | ---: | ---: | --- | ---: | ---: | ---: |
| Boussinesq fusion | Remove the dead `U` write and fuse the stored divergence addition into final projection | 0.98193 | 1.81% | 0.9859, 0.9795, 0.9791, 0.9832 | `2.29e-21` | Not compared across roles | 4,222,276,248 |
| Boussinesq parallel | Fusion plus pointwise projection ranges and wave-plus-balanced addition | 0.93943 | 6.06% | 0.9343, 0.9379, 0.9425, 0.9430 | `2.54e-21` | `1.33e-15` | 4,222,276,248 |
| Boussinesq assembly | Parallel stage plus pointwise coefficient assembly | 0.88192 | 11.81% | 0.8881, 0.8750, 0.8799, 0.8847 | `2.54e-21` | `1.55e-15` | 4,222,276,248 |
| Hydrostatic projection | Pointwise projected-field coefficient loop | 0.92522 | 7.48% | 0.9304, 0.9353, 0.9261, 0.9093 | `1.78e-22` | `2.84e-14` | 281,630,917 |

Every cross-role comparison that was performed passed the existing `atol=1e-12`, `rtol=1e-10` scientific gate. Within-role replay differences were zero. Every strengthened screen retained matching storage and producer counters between roles. A 24-sample Boussinesq block recorded 24 state validations, 24 phase preparations, 24 tiled nonlinear executions, 96 coefficient assemblies, 192 vertical preparations, 576 vertical-operator executions, and 567,240 matrix-group executions for each role. The Hydrostatic block likewise retained 24 validations, 24 phases, 24 tiled executions, 96 assemblies, 168 vertical preparations, and 240 vertical-operator executions.

The first fusion-only screen predates the strengthened physical-field payload comparison and driver provenance checks. It compared complete flux payloads and checked each role's fields for finiteness and repeatability, but did not compare fields across roles. The later Boussinesq stages contain that fusion unchanged and compare both complete flux and complete fields, so the final cumulative candidate has the stronger evidence. The baseline runtime source is identical across the old and strengthened harness binaries; their executable identities differ because the strengthened worker writes one field payload per role.

## Benchmark throughput and limits

Resident workers make this staged search practical without changing the timed operation. A Boussinesq screen completed 192 timed RHS calls and 16 warmups in 43–45 s, including about 20.4 s to load and prepare the two workers and about 5.6 s of untimed output validation. Keeping both workers resident avoided six additional checkpoint-load and kernel-preparation cycles of roughly ten seconds each. The Hydrostatic screen completed the same 208 total calls in 4.13 s, including 2.07 s of timed RHS work and 1.30 s of untimed validation. These are harness-throughput observations, not model speedups, and loading, preparation, validation, forcing schedules, integration, and output remain outside the RHS ratios.

The strengthened driver verifies immutable input, executable, provider, worker, CMake, and driver hashes before and after a campaign. It rejects differing compiler, provider, execution options, backend, schedule, producer counts, or numerical outputs. Complete field payloads are compared in chunks and then removed only after success, with their hashes and comparison results retained; flux payloads remain in the archive. The six driver tests pass with NumPy and with the standard-library comparison path.

Host journals record intermittent background work, including antivirus services in several screens, system intelligence activity in the fusion screen, and Time Machine activity in the parallel screen. All four paired ratios in every stage favored the candidate, but these short screens have no confidence interval and are sensitive to host scheduling. The Hydrostatic process lasted only about four seconds and is especially dependent on sparse host observations. The separately sampled baseline's 85.66 ms median includes profiler overhead and must not be compared directly with the unsampled 72–74 ms Boussinesq block means.

## Verification before full-model qualification

- All 50 PortableRuntime Release tests and all eight compiled-kernel tests passed.
- The prepared executor plus Hydrostatic and Boussinesq kernel tests passed under AddressSanitizer.
- The changed Hydrostatic and Boussinesq translation units compiled warning-clean with GCC 14. The known Apple deployment-target warning remains unrelated to this source.
- The stored-divergence arithmetic oracle passed exact-bit comparison under Apple Clang and GCC 14 for both `hasW` paths, including inertial modes, signed zeros, cancellation, extreme finite operands, and odd tails. GCC required the current Xcode macOS 26.5 SDK because its older default SDK headers are incompatible with this host toolchain.
- Focused Boussinesq kernel contracts passed after the fusion, parallel, and assembly stages, and focused Hydrostatic kernel contracts passed after its projection change.

## Alternatives retained as history

The candidate deliberately avoids a new matrix or coefficient layout. Earlier Hydrostatic screens found direct strided and packed vDSP assembly about 3.7 and 3.0 times the scalar loop on EddyTide, with 32 or 128 extra bytes per coefficient; the larger case was slower still. Accelerate has no public variable-size grouped GEMM batch that removes the cost of packing these heterogeneous vertical products, and direct split-real matrix views plus the existing vertical group executor are already selected. Repeating those packing and vDSP routes would add traffic and storage to solve a pointwise scheduling problem.

A coordinated multi-field `u/v/w` assembly could reduce factor traversals, and a phase-adjusted factor layout could expose more vector work, but either would add producer, cache, demand, or retained-storage obligations. The current loops use existing buffers and workers and already remove most of the measured pointwise gap. Vector phase preparation was also screened previously: its remaining gain was small, changed scalar-libm results by several ulps, and required special reference-time handling. None of those larger changes is needed for this increment. FFT, BLAS, and tiled advection remain the dominant summed active work and should be assessed from a new post-adoption profile rather than inferred from these loop screens.

## Initial interrupted campaign (historical)

The frozen campaign ran from 20:22:45 to 20:34:05 UTC on September 12, 2026 against the exact accepted #494 runner. It completed two warmup pairs plus eight measured pairs for EddyTide and larger Hydrostatic, then two Boussinesq warmup pairs. All 22 completed full-model comparisons passed the existing scientific tolerances and exact integration-decision checks. Evaluator/actual-producer records, retained owned bytes and maximum-live owned bytes were identical in every completed pair. This is reuse-policy evidence; no low-memory storage or timing qualification is claimed.

| Fixture | Completed measured pairs | Provisional candidate/baseline integration ratio | Provisional paired bootstrap 95% interval | Status |
| --- | ---: | ---: | --- | --- |
| EddyTide `256 x 256 x 28`, RK78 | 8 | 0.93854 | 0.92945–0.94639 | 6.15% observed reduction; host-contaminated |
| Composite Hydrostatic `256 x 256 x 129`, RK4 | 8 | 0.98628 | 0.85675–1.13232 | Inconclusive; strong order effect and host activity |
| Boussinesq `256 x 256 x 129`, RK78 | 0 | — | — | Two correctness warmup pairs passed; timing incomplete |

The idle preflight passed, but the host journal subsequently recorded unrelated Photos, indexing, backup and other activity above the preflight's 10% CPU threshold. All measured Hydrostatic pairs had at least one such observation within conservative request-to-receipt windows. These windows include output hashing, and `ps` percentages are sampled averages rather than continuous instantaneous measurements. They establish a host-quality concern; they do not quantify precisely how much each integration was affected. The bootstrap intervals summarize those observations and cannot remove the contamination or order effect.

The coordinator stopped the remaining campaign instead of completing expensive repeats that could not establish quiet-host acceptance. Both Boussinesq warmup comparisons were already complete; the newly started first measured baseline was interrupted and remains unqualified. No experiment or system service was stopped. No samples were discarded based on measured speed. The archive retains completed runs, comparisons, host observations, the explicit operator-stop record and the partial run. Candidate/baseline executable and provider hashes passed postflight verification.

The 14 focused MATLAB parity/integration tests and all six forward/restart configurations under both providers passed. Source-linked receipts were collected from the frozen runner/probe; all eight strict catalog tests and the repository/provenance checker passed. Scientific source and compiled inputs remained unchanged after `ed6049a2`.

## Final paired acceptance, September 12

The second frozen campaign used an Apple M4 Max (16 logical cores, 128 GiB RAM), macOS 26.6.1, Apple Clang 21, FFTW 3.3.11 NEON/pthreads and Accelerate. Execution used one internal FFTW thread, twelve horizontal workers, eight pointwise workers, and one Hydrostatic/eight Boussinesq vertical-group workers. It completed all 30 comparisons: two warmup pairs and eight alternating measured pairs for each of the same three fixtures. It finished at 23:33:36 UTC after approximately 23 minutes including copying, startup, validation and hashing. Every process passed an idle preflight; no local build or MATLAB test ran concurrently, and no experiment or system service was altered. Outputs lived under `/private/tmp` and were removed only after successful comparison, retaining hashes and comparison records.

| Fixture | Baseline integration | Candidate integration | Candidate/baseline | Paired bootstrap 95% | Complete-process ratio |
| --- | ---: | ---: | ---: | --- | ---: |
| EddyTide `256 × 256 × 28`, RK78 | 16.862 s | 15.719 s | 0.93219 | 0.92493–0.93880 | 0.93277 |
| Composite Hydrostatic `256 × 256 × 129`, RK4 | 5.374 s | 5.294 s | 0.98515 | 0.97461–0.99586 | 0.98771 |
| Boussinesq `256 × 256 × 129`, RK78 | 4.419 s | 3.938 s | 0.89103 | 0.81961–0.96301 | 0.99137 |

Times are geometric means of the eight measured processes per role. Ratios below one favor the candidate. EddyTide reduces observed integration time by 6.8%; larger Hydrostatic by 1.5%. Their order-stratified ratios agree closely: 0.93202/0.93235 for EddyTide and 0.98554/0.98476 for Hydrostatic. The full paired runtime gates pass, with no point estimate exceeding the 3% regression-investigation threshold for integration or complete process lifetime.

Boussinesq's aggregate reduction is 10.9%, but it remains sensitive to process order: baseline-first pairs give 0.99187 and candidate-first pairs 0.80044. Each complete counterbalanced block favors the candidate (0.93320, 0.89795, 0.87434, 0.86033); this supports adoption alongside the independently measured 11.8% resident RHS improvement, but does not establish a stable 10.9% production speedup. Its roughly four seconds integrating are a small part of complete-process time, whose observed improvement is only 0.9%. Startup, persistent matrix loading and host scheduling remain material to this short fixture.

This was not a perfectly isolated host. Intermittent desktop activity above 10% of one CPU was sampled in 14/16 estimated EddyTide integration windows, 1/16 Hydrostatic windows and 7/16 Boussinesq windows; maximum individual unrelated observations were 127.5%, 13.8% and 100%, respectively. No pairs were excluded or chosen by speed. Estimated windows combine process-launch timestamps and reported phase durations; two-second `ps` averages cannot prove instantaneous isolation. Bootstrap intervals describe the observed pairs and do not remove host effects. The original heavily order-sensitive/contaminated campaign remains historical evidence above; the final decision uses the completed balanced campaign, its order sensitivity and the resident screens together, not a claim of zero background activity.

All 30 scientific comparisons pass `atol=1e-12`, `rtol=1e-10`; maximum absolute differences are `1.28e-17`, `1.11e-16` and `2.08e-17` for EddyTide, Hydrostatic and Boussinesq. Integration controls and accepted-step decisions match exactly under the existing normalized-error comparison policy. Evaluator and actual-producer records are identical. No owned storage grows: reported total retained/peak differences are at most 40, 288 and 128 bytes, while integrator, forcing, dense-history and occurrence-array ledgers are unchanged. Process peak-RSS geometric ratios are 1.00012, 0.99999 and 0.99999. This qualifies the requested reuse policy; low-memory performance/storage was not rerun.

The source, binaries, providers, fixtures and frozen driver hashes pass postflight checks. [Final machine-readable results](../../.github/ci-evidence/prepared-pipeline/final-summary.json), [host attribution](../../.github/ci-evidence/prepared-pipeline/final-host-attribution.json) and [archive hashes](../../.github/ci-evidence/prepared-pipeline/final-archive-manifest.json) retain every pair. Full evidence is in `wave-vortex-model-benchmark-artifacts/prepared-pipeline-20260912/final-acceptance` relative to the workspace root.

## Adoption decision and next increment

Adopt the frozen runtime after the final combined required CI passes. The complete paired results and strengthened resident screens support the existing-buffer, existing-worker change; no larger layout or cache refactor is needed. The integrated #497 guard additionally checks the actual selected runner/probe source identity during refreshed six-family forward qualification. Compiled runtime inputs remain byte-identical to `ed6049a2`; subsequent edits concern verification, evidence and website benchmark presentation.

After adoption, reprofile the accepted candidate before choosing a larger matrix-to-horizontal layout or coordinated `u/v/w` refactor. The present change uses existing workers and buffers and requires no new coefficient representation. Remaining gains should be justified against the now-smaller serial assembly/projection cost.

Compact [resident-screen results](../../.github/ci-evidence/prepared-pipeline/resident-screens.json), [provisional full-model summary](../../.github/ci-evidence/prepared-pipeline/provisional-summary.json), [host attribution](../../.github/ci-evidence/prepared-pipeline/host-attribution.json), and [archive hashes](../../.github/ci-evidence/prepared-pipeline/archive-manifest.json) accompany the source. Large artifacts live under `wave-vortex-model-benchmark-artifacts/prepared-pipeline-20260912` relative to the workspace root.
