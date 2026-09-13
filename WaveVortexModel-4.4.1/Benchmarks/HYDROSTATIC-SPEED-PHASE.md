# Hydrostatic speed reduction and phase preparation

## Scope and baseline

This increment follows the qualified coefficient-assembly implementation in PR #473. It targets the remaining serial maximum-horizontal-speed reduction used by adaptive damping and the sine/cosine phase preparation. The fused `u/v/w` reconstruction proposal remains a separate follow-up in [the assembly investigation](HYDROSTATIC-ASSEMBLY.md).

The fresh EddyTide profile used the frozen PR #473 executable, a disposable copy of its original fixture, 60 RK78 steps of at most 300 s, and a 12-second steady sampling window. Maximum-speed reduction accounts for 9.23% and phase preparation for 8.44% of 8,201 main-thread samples. These inclusive sampled fractions identify opportunities; they are not precise wall-time measurements or projected speedups. The native execution policy remains twelve horizontal workers, eight pointwise workers, one FFTW internal thread and Accelerate matrices. The running experiment was not modified.

The spectral-kernel-benchmarks variable-family handoff explicitly excludes phase construction and diagnostics from its timed boundary. Its selected FFT/matrix algorithms and workers are retained; the new measurements address these WVM-owned producers.

## Selected implementation

- Reduce supplied `u/v` fields through the existing prepared pointwise executor. Each chunk uses the original `std::hypot` and `std::max` expression, then merges one local maximum under a stack mutex. No field scratch, worker pool or retained result is added. This preserves robust norm scaling, zero initialization and the existing treatment of nonfinite contributions. Validation precedes execution and the output scalar is published only after the synchronous dispatch completes.
- Route both Hydrostatic forcing-cache paths and the Hydrostatic diagnostic maximum through that kernel primitive. Existing evaluation keys, cache ownership and successful-producer metrics remain authoritative. The primitive neither reconstructs fields nor introduces another cache.
- Prepare phases with the same `omega * (t - t0)`, `std::cos` and `std::sin` expressions on disjoint worker ranges. The caller waits before reading phases. Existing exact-state validation, reference-time identity, once-per-evaluation preparation and reentrancy guards remain unchanged. No phase recurrence or cache across evaluations is introduced.

This change is limited to Hydrostatic consumers. Boussinesq and constant/QG execution retain their existing producers until separately measured.

## Alternative phase screen

An isolated eight-worker screen used the real 200-byte factor record and EddyTide's 205,902 coefficients. Three process samples with thirty iterations measured approximately 0.307 ms for worker-parallel scalar sine/cosine, 0.208 ms for Apple `simd::sincos(double2)`, and 0.173 ms for vForce `vvcosisin`, including argument packing and conforming output conversion. The SIMD path adds no arrays; vForce needs 4,941,648 bytes, about 1% of EddyTide owned peak. That fits the broad storage budget but is additional capacity and traffic.

Both vector paths differed from scalar libm by up to three ulps on the tested physical-angle distributions. On this host their cosine at zero was one ulp below one; a production route would need to preserve the exact reference-time identity explicitly. The installed vForce SDK also permits result/denormal variation across platforms. These are bounded measurements, not a general accuracy guarantee.

The extra vector gain beyond scalar workers is only about 0.10–0.14 ms per phase evaluation, suggesting roughly half a percent of whole-EddyTide time under a simple estimate. That does not justify the additional numerical/portability work for this increment. Keep the original scalar math with prepared workers; do not interpret the rejected vDSP coefficient-assembly result as evidence against all vector trigonometry.

## Exploratory results

Three alternating fresh-process pairs, without warmups, compared each candidate to the frozen PR #473 executable:

| Candidate | EddyTide integration ratio | Large Hydrostatic composite ratio |
| --- | ---: | ---: |
| Parallel maximum-speed reduction only | 0.91062 | Not measured in this ablation |
| Parallel reduction and scalar phase preparation | 0.84254 | 0.97813 |

Every full-output comparison and state/step-count check passed at existing tolerances. No field/spectral scratch was added; small byte differences between report roles are metadata strings. These screens select the candidate; the frozen paired campaign below supplies qualification.

## Verification ledger

Independent read-only review of reduction ownership/arithmetic and phase lifecycle/tests passed. The native tests compare one and two workers, cover every field/component/derivative, exercise robust speed values at ordinary, 1e308 and subnormal magnitudes, and retain the existing NaN/infinity semantics. Invalid speed inputs cannot publish a result, and prepared reduction performs no C++ allocation. Phase tests compare every value with the scalar expression at repeated, reversed and extreme finite times, including the exact reference time.

All 57 native contracts passed in Release and ASan/UBSan, as did the GCC 14 Hydrostatic kernel/runtime build and tests. Both Hydrostatic MATLAB methods passed across nine provider/layout combinations, together with five source-selection tests. Source-linked receipt and performance results are appended after completion. No MATLAB scientific source or checkpoint change is intended.

Artifacts: `OceanKitRepositories/wave-vortex-model-benchmark-artifacts/hydrostatic-speed-phase-20260911`, including the baseline profile, phase screen, frozen exploratory executables/sources and paired comparisons. One coordinator owns production code and qualification; one reused Sol agent owns the phase feasibility screen and independent review. Timed runs did not overlap builds or tests on the host.

## Frozen qualification

Runtime implementation commit: `72278e2c`. The source-linked forward-integration refresh passed all six configurations with reference and native providers, followed by all fifteen catalog/matrix checks. PR #473 merged at `be656012`; its tree exactly matches the original `3cec040e` prerequisite. Integrating it into this branch changed no tracked file content. Frozen candidate `8c2b816d` therefore includes the reviewed implementation and refreshed evidence without another runtime correction. Existing source tests remain applicable; the runner was rebuilt to record the frozen commit.

The established qualification driver ran two warmup groups and eight measured groups per fixture, alternating baseline/reuse/low-memory order. Baseline is the frozen PR #473 executable, not the older pre-assembly version. All fixture, provider, thread and integration settings were retained. Both executable and input hashes passed postflight checks.

| Workload | Reuse integration ratio (95% paired bootstrap interval) | Low-memory integration ratio | Reuse complete-process ratio |
| --- | ---: | ---: | ---: |
| EddyTide 256 × 256 × 28 | 0.82700 (0.81949–0.83437) | 0.82748 | 0.83010 |
| Constant nonhydrostatic composite 256 × 256 × 129 | 0.97366 (0.92644–1.00601) | 0.97063 | 0.99541 |
| Hydrostatic composite 256 × 256 × 129 | 0.97679 (0.96681–0.98633) | 0.97246 | 0.98457 |

The unchanged constant control has a wider interval, including no improvement; one slow baseline integration sample contributes to the apparent mean gain. No speedup is attributed to the Hydrostatic change in that control. No workload triggers the greater-than-3% runtime investigation gate.

All numerical comparisons passed at existing tolerances, with identical integration controls and accepted/rejected-step decisions. Every default-policy evaluation had zero duplicate producers. The representative Eddy report has 157 state validations/phase preparations (initial state plus 156 RHS evaluations), 156 horizontal-speed reductions, and zero RHS evictions/recomputations/duplicate executions. Existing cache reuse remains intact.

No field/spectral scratch, persistent worker storage or prepared C++ allocation was added. Maximum-live and retained owned storage were no larger than baseline for reuse. Low-memory maximum-live storage was at most baseline, retaining that policy's existing larger-fixture savings. Whole-process timing includes preparation/output; these fixtures do not qualify an output-heavy production schedule.

[Performance summary](../.github/ci-evidence/hydrostatic-speed-phase/performance.json), [verification ledger](../.github/ci-evidence/hydrostatic-speed-phase/verification.json), and [archive hashes](../.github/ci-evidence/hydrostatic-speed-phase/archive.json) preserve the result. `docs:check` passed with zero generated differences. No MATLAB source, package metadata or released snapshot changed. Required hosted CI is the remaining merge gate.
