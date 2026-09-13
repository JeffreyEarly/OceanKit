# Hydrostatic coefficient assembly optimization

The subsequent [speed/phase increment](HYDROSTATIC-SPEED-PHASE.md) builds on this qualified implementation. Fused `u/v/w` assembly remains a separate recorded follow-up.

## Decision and boundary

Adopt contiguous coefficient ranges on the existing prepared pointwise executor. Each coefficient retains its original field/component selection, complex multiplication/addition order, derivative scaling, and split/interleaved output. Dispatch joins before the vertical transform reuses modal storage. A one-worker executor executes inline. No new full-grid arrays, workers, cache keys, MATLAB behavior, or checkpoint fields are introduced. Small-grid speed tuning is outside this increment.

Baseline: v4 main `e5134612d5a339707eca8243c17432691cb4c378`, with the unchanged qualified runtime from PR #471. Frozen baseline executable SHA-256: `f1a642a4fe53176f257b8a5bd8b8bf9ad3cfa5481fdc4d9c892b19d13cdaa3e9`. The executable predates merge; the compiled-input equivalence is documented in [the preceding profile](POST-CACHE-PROFILE.md).

The investigation addresses coefficient construction, explicitly excluded by section 1.4 of the spectral-kernel-benchmarks [handoff](https://github.com/JeffreyEarly/spectral-kernel-benchmarks/blob/6f2e4cd9e2f092433a584645e573d818ecaca42b/docs/wvm-compiled-core-integration-handoff.md). Closed issues #22 and #25 already establish the prepared worker and compact pipeline baseline. No FFT or matrix-kernel survey was repeated.

## Investigation

A temporary timing probe isolated assembly from vertical multiplication, normalization and horizontal inverse transforms. For the frozen EddyTide continuation it measured 2,028 reconstructions: assembly 1.1150 s, vertical multiplication 0.2803 s, derivative normalization 0.0256 s, horizontal inverse 0.7888 s, against 4.1505 s integration. Assembly accounts for about 27%; post-matrix normalization is less than 1%. The instrumented run matched its paired baseline numerically and had an integration ratio of 1.0012. Instrumentation is retained only in the artifact archive.

Short exploratory campaigns used three alternating fresh-process pairs, unchanged fixtures, no warmups, and complete output comparisons. They screen candidates; final qualification is recorded separately.

| Candidate | EddyTide integration ratio | Large Hydrostatic ratio | Decision |
| --- | ---: | ---: | --- |
| Existing executor, contiguous coefficient ranges | 0.88280 | 0.94019 | Select the simple implementation |
| Executor partitioned by complete horizontal mode | 0.88798 | 0.96103 | Reject; avoiding index division did not improve complete workload |
| Field switch specialized outside worker loop | 0.87442 | 0.95504 | Reject mixed gain; extra template dispatch does not clearly beat the simpler candidate |

All comparisons passed at existing tolerances, with identical state/step counts. The selected candidate retained the same owned storage. The first prototypes had a small-size serial threshold; the final implementation uses the existing executor for every extent, allowing the native small fixtures to exercise worker partitioning. One worker remains inline; performance at tiny sizes was not tuned.

## Accelerate and vDSP screen

Accelerate already supplies the selected vertical matrix multiplication. Its vDSP vector arithmetic was tested separately for the measured coefficient-assembly bottleneck. The scratch benchmark includes both interleaved and production-hot split destinations, the actual 200-byte mode-factor record, and simultaneous wave/geostrophic support. It covers the `u`, all-component, value expression; it is not full-model qualification of a vDSP implementation.

Direct strided vDSP uses six arithmetic passes and 32 additional bytes per coefficient. Packing uses four complex-input conversions, one factor-gather pass and six arithmetic passes, retaining 128 additional bytes per coefficient. Split output avoids a final conversion. The tests include these recurring costs; allocation occurs before timing. Ordinary finite synthetic data agrees within about 4.5e-16 scaled error. Compiler-specialized C++ is exact on this fixture, but its unconditional zero-factor arithmetic needs separate extreme-value/selection review before production use.

At EddyTide’s `18 * 11,439 = 205,902` coefficients, three process samples with thirty iterations measured 0.5587 ms for the original scalar split loop, 0.4440 ms for specialized C++, 2.0716 ms for strided vDSP and 1.6732 ms for packed vDSP. Paired ratios are 0.7943, 3.7074 and 2.9947 respectively. Extra vDSP scratch is 6,588,864 bytes strided or 26,355,456 bytes packed.

At the large coefficient count, five process samples with twelve iterations gave strided vDSP 6.41 times and packed vDSP 3.33 times the original scalar split-output loop. The compiler-specialized loop was 0.81 times the scalar loop. These are single-producer microbenchmarks, not expected whole-model speedups. The extra scratch and slower arithmetic screen out straightforward vDSP adoption on the current layout. The successful production change retains the original conditional arithmetic and adds no vDSP dependency.

## Larger pipeline options

1. **Fuse velocity value assembly into existing modal slots.** When the compiled evaluation plan requests `u/v/w` together, one factor traversal could fill the three existing modal slots, then reconstruct those fields before assembling `eta`. This could reduce four factor traversals to two without new full-volume storage. It needs an explicit batched producer, proof that the projected-flux slots are dead at that point, and demand-sensitive handling so a damping-only `u/v` request does not acquire unnecessary work. This is the preferred assembly-pipeline experiment if the next profile still warrants it.
2. **Reuse phase-adjusted coefficients within an evaluation.** An explicit internal producer could form the two phase-adjusted wave arrays once and feed several field reconstructions. This would remove repeated phase multiplications across fields, but adds `32 * Nj * Nkl` bytes if both arrays are retained. Reassociating `(factor * coefficient) * phase` as `factor * (coefficient * phase)` changes rounding and extreme-value behavior. It needs its own numerical and adaptive-control qualification, evaluator ownership, and low-memory accounting. Measure its remaining critical-path share after this increment before implementing it.
3. **Prepare compact hot factors and fuse coefficient operations.** The current 200-byte records contain much more than a particular field needs. A structure-of-arrays or tiled hot-factor representation could improve vectorization and reduce factor traffic. Evaluate fused compiler SIMD first; the vDSP screen does not prove that a new layout cannot work, but provides no reason to pay for a library-oriented rewrite now. Include preparation, extra retained factors, rare component paths and multi-field consumers in any follow-up screen.
4. **Reuse or batch modal fields across derivatives.** Horizontal derivatives share a modal field before their wavenumber multiplier. Retaining or batching those fields could eliminate assembly and some vertical work, but increases intermediate lifetimes and complicates the deliberately streamed derivative path. Fusing a larger reconstruction/vertical/FFT pipeline is an architectural experiment, not a free loop change. Require a measured complete-model benefit and preserve both memory policies before adopting it.
5. **Prioritize full-grid tracer vertical calculus for the large tracer workload.** The prior profile attributes about 28% of main-thread samples to its packing/scaling/scattering loop. That is a separate, larger measured opportunity on that workload. Start with prepared workers or a fused real-column pipeline, preserving all resolved frequencies and the existing operator. Do not infer an EddyTide improvement from a tracer-specific optimization.

These options are independent of cache correctness. Phase factors are already prepared once per evaluation; the potential reuse above concerns subsequent arithmetic and representation. Boussinesq has an analogous assembly loop, but should receive its own representative profile and qualification before extending this change.

## Verification and evidence

Investigation artifacts are retained under `OceanKitRepositories/wave-vortex-model-benchmark-artifacts/hydrostatic-assembly-20260911`: frozen sources/executables, protocols, paired reports, output comparisons, timing probe and vDSP source/logs. The experiment's own executable and files were untouched. One coordinator owns production edits and timing; one reused Sol agent owns the bounded vDSP screen and independent review.

Verification ledger before source freeze: independent read-only implementation and test review passed; all-field/component/derivative serial versus two-worker equality passed for compact split and established interleaved storage; GCC 14 warning-clean kernel build/test passed (the existing Apple deployment warning remains); 57 native Release contracts passed; all 22 focused MATLAB/provenance methods passed, including Hydrostatic parity across nine providers/layouts. Sanitizer and paired qualification results follow below. Documentation consistency passed with zero generated differences; final repository/provenance and whitespace checks passed. No generated website or production MATLAB source changes are required.

### Frozen qualification

Candidate runtime commit: `3e907c10`. Two warmup groups and eight measured groups per workload ran the frozen baseline, candidate reuse and candidate low-memory executables, reversing order on alternate groups. The existing qualification driver, numerical tolerances, fixtures, integration controls, native FFTW provider and thread policy were unchanged. Variable-family execution uses twelve horizontal outer workers, eight pointwise workers, one FFTW internal worker and Accelerate vertical multiplication.

| Workload | Reuse integration ratio (95% paired bootstrap interval) | Low-memory integration ratio | Reuse complete-process ratio |
| --- | ---: | ---: | ---: |
| EddyTide 256 × 256 × 28 | 0.89242 (0.88429–0.89913) | 0.88487 | 0.89534 |
| Constant nonhydrostatic composite 256 × 256 × 129 | 1.00498 (0.98945–1.01743) | 0.98853 | 1.00064 |
| Hydrostatic composite 256 × 256 × 129 | 0.95529 (0.94238–0.96613) | 0.95371 | 0.96757 |

Every output comparison passed at the existing tolerances; integration controls and accepted/rejected-step decisions were identical. Default-policy duplicate execution counts remained zero. No full-volume or spectral scratch was added. Reported owned-memory ratios differ slightly from one because report/request strings differ between roles; the largest low-memory maximum-live ratio was 1.000000064, well within the 1.03 gate. The larger fixtures retain the existing low-memory savings; those savings are not caused by this loop patch. No workload triggered the greater-than-3% runtime investigation gate. Whole-process results include preparation and output; this is not qualification of a normally output-heavy production run.

The combined sanitizer suite also passed: six script checks passed initially, and all 51 instrumented binaries passed on retry after disabling unsupported Apple LeakSanitizer startup. ASan and UBSan stayed enabled. The initial startup failures are preserved; Linux required CI retains its leak checks. No source correction or repeated performance campaign was needed.

[Performance summary](../.github/ci-evidence/hydrostatic-assembly/performance.json), [verification ledger](../.github/ci-evidence/hydrostatic-assembly/verification.json), and [archive hashes](../.github/ci-evidence/hydrostatic-assembly/archive.json) bind the evidence. The artifact archive includes exact input copies and a separate replay manifest; the originally executed protocol remains unchanged. PR #473 merged at `be656012` after required hosted CI passed.
