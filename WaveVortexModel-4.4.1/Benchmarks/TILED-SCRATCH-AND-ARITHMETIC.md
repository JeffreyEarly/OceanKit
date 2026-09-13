# Tiled scratch reuse and arithmetic specialization (#491)

This increment follows the accepted tiled-advection pipeline in PR #490. It reuses native provider storage and exposes disjoint accumulation spans to the compiler. The measured implementation is frozen at `0a819a495101370c7edd174a9fbd62d5a8715541`; the control is the exact qualified #490 executable (`78025904`, merged through `1ec6d768`). Later receipt/report changes do not alter the measured compiled inputs.

## Storage and arithmetic

Four base fields use the existing `16 W M` complex tile. Their maximum requirement is `4 W D M`, where `D <= 4`. Each target-z tile is overwritten with projected flux after its last derivative read. Additional packed capacity falls from `(4 + 2 T) W D M` to `T W D M` complex values. The existing shared-resource guard and synchronous dispatch protect the entire reuse interval. Worker/lane slices remain disjoint. No family API or broader workspace loan is added.

Ordinary and density-z accumulation use separate private helpers. Local compiler alias qualifiers express the separation already guaranteed by owned flux/derivative scratch and distinct physical-field slices. Clang and GCC generate vector instructions. The code preserves axis order, correction-before-multiplication, ordinary `derivative + 0.0`, and compiler contraction semantics under existing flags. No fast-math, vendor change, FFT-count change, or MATLAB scientific change is involved.

## Independently attributable stage screen

The stage driver calls the actual native retained operator using the representative mode sets, deterministic nonzero-density fields, 12 workers and both three- and four-target workloads. Each of four fresh-process blocks rotates/reverses baseline, memory-only and combined roles. Each process warms twice and measures nine calls. Complete physical fields and split target outputs satisfy the unchanged comparison tolerance, with identical FFT counts.

| Nz | Targets | Memory-only / baseline | Combined / memory-only | Combined / baseline |
| ---: | ---: | ---: | ---: | ---: |
| 28 | 3 | 1.01603 | 0.93526 | 0.95025 |
| 28 | 4 | 0.99717 | 0.96104 | 0.95832 |
| 129 | 3 | 1.00446 | 0.93974 | 0.94393 |
| 129 | 4 | 0.99081 | 0.94741 | 0.93870 |

Ratios below one are faster. The memory change is approximately timing-neutral in this screen; the private loops reduce the complete stage by 3.9–6.5% relative to memory-only. These timings exclude vertical MM, assembly, forcing and integration and cannot support a whole-model speed claim by themselves.

## Complete-model qualification

| Workload | Integration reduction | Paired bootstrap 95% interval | Process-lifetime reduction | Owned peak reduction |
| --- | ---: | ---: | ---: | ---: |
| EddyTide, 256 × 256 × 28 | 2.17% | 1.37–2.91% | 2.17% | 46.12 MB |
| Constant-stratification control | 3.63% | -0.25–8.94% | 0.43% | 0.00 MB |
| Larger Hydrostatic | 0.25% | -6.20–6.49% | 0.52% | 35.65 MB |
| Boussinesq, 256 × 256 × 129 | 5.51% | -15.71–23.09% | 0.57% | 40.74 MB |

All 32 measured pairs and eight warmup pairs passed scientific comparisons and identical integration/state decisions. The other actual producer counts are unchanged, and tiled executions/reused columns are positive for the two MM families. No observed mean integration or process-lifetime regression exceeds 3%; the wide intervals on the larger fixtures do not rule out changes of that size. The unchanged constant control should be interpreted as measurement variation. Owned storage and process RSS are separate measurements; full values are retained in the evidence. Provider scratch savings match the formula exactly: 46,122,048 / 35,648,256 / 40,740,864 bytes for EddyTide / larger Hydrostatic / Boussinesq. Full-model accounting includes incidental output/orchestration differences of 0 / 8 / 128 additional bytes; the constant control differs by 8 bytes in output-sink accounting. These small differences are not scientific-storage savings.

EddyTide shows a modest measured improvement. The unchanged control's interval includes zero. Larger Hydrostatic and Boussinesq timing remains inconclusive; the nominal Boussinesq 5.51% reduction is not a supported speedup claim. Slow Boussinesq samples occur in both versions and concentrate in the second execution of a pair. Alternation balances that order effect, but eight pairs leave a wide interval.

Each campaign starts after an idle-host preflight, with no concurrent builds or MATLAB jobs. Intermittent operating-system services appear in the retained activity journal, so these are ordinary local-host measurements with recorded interference, rather than a guarantee of an entirely silent system. No outlier was removed from the timing statistics.

The preliminary two-pair screen is preserved separately: EddyTide was 0.86% slower and larger Hydrostatic 2.17% faster. Boussinesq had one slow baseline outlier, so its apparent 21% gain was not used as an adoption claim. The final campaign uses the established extended windows, two warmup pairs and eight alternating measured pairs per profile, plus an unchanged constant-stratification control. Host process observations are retained; no experiment or unrelated process is modified.

The initial campaign completed the first three profiles, then ran out of disk space during Boussinesq's first warmup output sync. All completed comparisons were retained. After hashing this task's redundant generated payloads and retaining their comparison/failure reports, 14.1 GB was freed. Only Boussinesq was resumed, with two fresh warmups and eight measured pairs, the unchanged driver/controls, and source/binary/provider/fixture postflight checks. The combined summary uses the original aggregation and bootstrap code. [Recovery evidence](../.github/ci-evidence/tiled-scratch-simd/resumption.json) and the archive record identify the two campaign segments. No successful workload or scientific suite was repeated.

## Verification and boundaries

Fresh Clang and GCC provider/operator plus four family tests passed. Release and ASan/UBSan combined suites each passed 58/58. The new checks cover exact reduced capacity, split/interleaved views, both target counts, depth and worker tails, alternating shared-resource calls, zero warmed allocations, signed zeros and cancellation-sensitive density arithmetic. Independent review found no remaining alias/lifetime/arithmetic blocker. Scientific tolerances, integration controls and evaluator counts remain acceptance gates.

All 12 affected MATLAB/scientific-source checks and all six forward/restart family cases with both providers passed. All eight strict catalog tests passed after collecting the new fragments; every fragment matches the executed frozen runner and probe hashes.

The first new-test builds found two mixed-type auto declarations; they were corrected before qualifying fresh executables. An ungated command had run stale tests after the build failed; those results are explicitly excluded. An initial macOS sanitizer invocation requested unsupported leak detection and aborted before running tests. The supported address/undefined-behavior rerun passed. Failed logs remain in the archive. No source change followed these successful checks.

Representative Boussinesq output, output hashes/comparisons, frozen binaries, disassembly and failed/successful logs remain in `OceanKitRepositories/wave-vortex-model-benchmark-artifacts/tiled-scratch-simd-20260912`. Compact evidence and hashes are linked below. Low-memory timing is outside this increment; existing fallback correctness remains covered. Tracer work, a wider real-workspace loan, RK78 combination loops and Boussinesq assembly/projection batching remain separate candidates.

See [performance results](../.github/ci-evidence/tiled-scratch-simd/performance-summary.json), [producer and memory checks](../.github/ci-evidence/tiled-scratch-simd/producer-and-memory.json), [host observations](../.github/ci-evidence/tiled-scratch-simd/host-observations.json), and [archive hashes](../.github/ci-evidence/tiled-scratch-simd/archive.json). Raw artifacts are retained under `OceanKitRepositories/wave-vortex-model-benchmark-artifacts/tiled-scratch-simd-20260912`. Required hosted checks gate merging.
