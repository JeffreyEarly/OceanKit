# Native tracer derivatives

Issue [#520](https://github.com/JeffreyEarly/wave-vortex-model/issues/520) follows the [hydrostatic composite audit](HYDROSTATIC-COMPOSITE-AUDIT.md). The v4 candidate accelerates public derivatives used by real tracer volumes while preserving their MATLAB definitions. Horizontal derivatives use axis-only real FFT pairs shared by compiled transform families. Variable-stratification Hydrostatic and Boussinesq vertical derivatives and integrals apply prepared balanced F/G matrices directly to the volume.

The frozen final candidate is commit `e9cba05235ce1a2f7988f3d402d1a56537067faf`, compiled source aggregate `91487ae2e6311d95e0d41aaa591247bba420c0a7ca943cee0e664c88873454e1`, and MEX SHA-256 `a16e7d1cc83707d028f1665dbdd5cff9f411e0c99072d9e7334b3cbeb25796fa`. Repeated-process qualification passed; release verification is tracked in the verification ledger.

## Horizontal derivatives

For each physical field, `diffX` performs batched real forward transforms along x, multiplies each frequency by `(ik)^order`, and performs the corresponding inverse transforms. `diffY` applies the same operation along y. The inverse normalization uses the transformed-axis length. Even-grid Nyquist coefficients are zero for odd derivative orders and retained for even orders, matching existing public behavior.

The native implementation reuses the retained provider's persistent worker pool and x-row plans. It prepares two y-axis real plans and bounded scratch once. Each warmed call processes all resolved input modes independently of the retained wave-vortex coefficient map. This preserves arbitrary real inputs and configurations with antialiasing disabled; the existing tracer tendency projection remains unchanged. Reference providers continue through the generic implementation.

Independent axis pairs require four separable transform passes for `diffX` and `diffY` together. A shared complete two-dimensional forward transform followed by two complete inverse transforms would require six. Existing compact retained-transform tiling remains unchanged. This is distinct from the broader MATLAB tracer-fusion approach rejected in #338.

## Hydrostatic and Boussinesq vertical calculus

The direct-volume vertical path handles real full-volume inputs on `WVTransformHydrostatic` and `WVTransformBoussinesq`. Matrix inputs retain the existing MATLAB-facing path; real vertical-first matrices use the prepared matrices internally through a packed path, while complex inputs retain the established complex calculation. Constant-stratification and QG implementations are unchanged by this vertical extension. Boussinesq uses its balanced F/G matrices for tracer calculus; the wave-dependent reconstruction matrices are unchanged.

Setup prepares five elementary real matrices: the first F derivative, first G derivative, second G derivative, F integral and G integral. The public higher derivatives retain their existing compositions:

| Input family | Order 1 | Order 2 | Order 3 | Order 4 |
| --- | --- | --- | --- | --- |
| F | `F1` | `G1*F1` | `G2*F1` | `G1*G2*F1` |
| G | `G1` | `G2` | `G1*G2` | `G2*G2` |

The MEX receives the original `[Nx Ny Nz]` real volume. Treating its horizontal points as rows, the kernel evaluates the equivalent of `A*D.'` directly, without either MATLAB permutation. The packed path used for vertical-first matrices processes at most `Nx*Ny` columns per chunk. Both paths reuse existing real workspace slots, and higher orders reuse the same elementary formulas rather than introducing a new numerical operator. Internal in-place Boussinesq reconstruction retains the established complex helper because the real matrix backend requires disjoint input and output buffers.

## Frozen qualification and recommendation

The three-fresh-process campaign on the Apple M4 Max passed all twelve runs, with profiling disabled and unchanged source, MEX, provider and fixtures. Every saved output graph passed the original absolute 1e-13 plus relative 1e-12 tolerance. All coefficient runs used 195 RHS evaluations; all composite runs used 199, with 15 accepted steps, zero retries and identical output counts.

| Workload | Configuration | All integration samples (s) | Median (s) |
| --- | --- | --- | ---: |
| coefficient-endpoint | baseline | 12.503877, 12.415890, 12.461762 | 12.461762 |
| coefficient-endpoint | candidate | 12.437828, 12.456795, 12.495818 | 12.456795 |
| composite-dense-output | baseline | 53.732419, 53.756837, 53.965720 | 53.756837 |
| composite-dense-output | candidate | 32.101743, 31.884297, 32.023581 | 32.023581 |

The composite median decreases **40.4%**, from 53.756837 to 32.023581 s (1.679× faster). The coefficient endpoint decreases 0.04%, effectively unchanged. No configuration regresses by more than 3%. Composite sample spans are 0.233 s baseline and 0.217 s candidate; coefficient spans are 0.088 and 0.058 s. No timing sample was excluded. The earlier profiled comparison showed a larger 47% total reduction; it remains diagnostic rather than the basis for the release claim. The original constant-stratification 114.205 s outlier belongs to a different, preserved campaign and was not reclassified or discarded.

**Recommendation: adopt** the combined narrow optimization for v4.4.1. It removes the measured tracer bottleneck without changing the scientific algorithm or undertaking a broad tracer/adapter refactor. Boussinesq derivative diagnostics and downstream correctness support the shared implementation; these Hydrostatic integration medians do not establish a Boussinesq whole-workload speedup. Further marshaling or particle optimizations are outside this change.

See [all qualification samples and comparisons](../.github/ci-evidence/issue-520-tracer/qualification-summary.json) and the [freeze record](../.github/ci-evidence/issue-520-tracer/freeze-combined.json). Raw evidence is archived under `wave-vortex-model-benchmark-artifacts/tracer-derivatives-20260913`. Release-only metadata and generated documentation changes do not alter the frozen compiled inputs.

## Development measurements

The reconstructed hydrostatic exponential fixture has a 256 × 256 × 129 grid and 85 retained vertical modes. The original published NetCDF model is unavailable. The reconstruction matches all 29 archived initial-condition diagnostics exactly but has a different file identity. Frozen archives and published benchmark rows remain unchanged.

| Profile boundary | Frozen compiled audit | Horizontal-only candidate | Combined direct-volume candidate |
| --- | ---: | ---: | ---: |
| Composite integration | 60.8147 s | 42.6575 s | 32.3362 s |
| Tracer flux | 30.4200 s | 12.4569 s | 5.1794 s |
| Tracer `diffX` + `diffY` | 18.4030 s | 1.4429 s | 1.4455 s |
| Tracer `diffZF` | 9.3660 s | 8.6276 s | 1.4242 s |

These profiles were collected separately and are developmental diagnostics, not an idle-host paired qualification. Untargeted timing differences cannot be attributed to the derivative changes. The direct run used adaptive `ode78`, accepted 15 steps with no rejected steps, and made 199 right-hand-side evaluations. Its source aggregate was `48aad8e068fe2c950c98db6d7a52a8ed18e5b2ae2664c64f2e274a729e405748` and its MEX SHA-256 was `55c084de8dabb722fd73bba875100642c43136e61bebeb14c11a074cdff7c491`; it predates the final Boussinesq extension and final frozen identity.

The historical Hydrostatic seven-sample warm primitive diagnostic reports these medians:

| Primitive | Compiled | MATLAB |
| --- | ---: | ---: |
| `diffZF` | 7.7548 ms | 10.3846 ms |
| `diffX` | 3.5208 ms | 3.4713 ms |

The full direct-run output graph passed the existing comparison limits across 70 variables and 64 records. Maximum relative error was `7.888164409073408e-13`; maximum absolute error was `1.8189894035458565e-12`. Coefficients, Eulerian fields, moorings, particles and the tracer all passed.

The final Hydrostatic prepared kernel storage increases by 12,246,304 bytes (11.68 MiB) at this grid, including horizontal worker scratch and the five small real vertical matrices. Kernel and engine ledgers both include this allocation; their increases must not be added. No extra full spatial-volume scratch is allocated for vertical calculus. This is an application-owned storage ledger, not measured process peak RSS.

## Boussinesq derivative diagnostic

A seven-sample warm diagnostic of the final source at 256 × 256 × 129 measured compiled `diffZF` median 7.8654 ms versus MATLAB 11.1760 ms, and compiled `diffX` 2.9789 ms versus MATLAB 3.5030 ms. Maximum absolute differences were 5.0191e-18 and 2.7105e-19 respectively, passing the existing tolerance. All samples are in [the diagnostic receipt](../.github/ci-evidence/issue-520-tracer/boussinesq-primitives.json). The input was read from the retained Boussinesq repeat-zero baseline archive; an older repeat-minus-two request points to a removed output, so it was not used. This diagnostic measures derivative calls, not Boussinesq integration performance.

## Verification status

Native tests cover horizontal derivative orders 1–4, nonsquare odd and even grids, full-spectrum and Nyquist inputs, padded layouts, immutable inputs, stable prepared storage and allocation-free warmed execution. Ownership tests cover failure of either y-plan acquisition, destruction and retry. Raw primitive and kernel tests cover Hydrostatic and Boussinesq F/G vertical derivatives through fourth order, both integrals, a chunk-tail extent, input preservation, stable storage and warmed allocation behavior.

The final native suite first passed 61 of 62 tests and exposed an in-place Boussinesq reconstruction alias. The final guard sends aliased internal reconstruction through the established helper; all five affected native tests then passed, with unaffected earlier passes retained. Four final MATLAB checks passed: Hydrostatic and Boussinesq variable-stratification parity, public layout parity, and downstream compiled-state consumers. This covers real-volume dispatch plus complex-volume and matrix fallbacks. Independent review passed after the alias correction.

All six representative forward-integration configurations passed for reference and native providers, including coefficient, tracer and particle evolution and continuation. The source-bound receipts were refreshed and their committed contract passed. Code Analyzer reported zero blocking findings across 187 production MATLAB files.

A hosted horizontal-only CI run passed its GCC Release job before that campaign was canceled to combine the vertical work. It predates both variable-stratification vertical extensions and does not qualify the final combined source. Final combined hosted CI and release publication follow the successful local qualification.

The qualification protocol freezes baseline and candidate builds, runs three fresh processes for each coefficient and composite configuration on an idle host, compares every complete output graph, and retains timing and numerical receipts. Any regression exceeding 3% requires investigation.

Evidence: [verification ledger](../.github/ci-evidence/issue-520-tracer/verification.md), [historical direct primitive samples](../.github/ci-evidence/issue-520-tracer/primitives-combined-direct.json), [historical direct profile](../.github/ci-evidence/issue-520-tracer/direct-composite.json), [historical direct output comparison](../.github/ci-evidence/issue-520-tracer/direct-comparison.json), [earlier horizontal primitive samples](../.github/ci-evidence/issue-520-tracer/exploratory-primitives.json), [earlier horizontal composite profile](../.github/ci-evidence/issue-520-tracer/exploratory-composite.json), [earlier horizontal output comparison](../.github/ci-evidence/issue-520-tracer/axis-comparison.json).
