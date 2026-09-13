# Shared Hydrostatic and Boussinesq field/gradient preparation

This increment follows [speed/phase preparation](HYDROSTATIC-SPEED-PHASE.md). It shares modal and horizontal-spectrum intermediates across field values and gradients in the two matrix-based models. No MATLAB scientific source, general full-grid calculus, tracer algorithm, FFT primitive or matrix backend changes.

## Execution design

An explicit kernel evaluation owns a bounded table keyed by registered state-view generation, normalized field and requested component. Each entry owns modal coefficients and a horizontal spectrum after vertical reconstruction. Aliases normalize before lookup: density uses eta, surface velocities use u/v, and surface displacement uses pi. Pressure remains distinct from pi to retain its existing scaling order. Composite horizontal vorticities consume their derivative operands. Additional state-view registration/removal cannot resurrect cached values when coefficient buffers are reused.

Modal assembly and base-spectrum preparation have independent success flags. A failed downstream inverse leaves successful upstream preparation available for retry. Evaluation teardown clears every validity flag; prepared capacity remains reusable. Standalone operations start fresh and clear on exit. Derived coefficient tendencies bypass the primary-state table. The existing evaluator continues to own physical-field/derivative caching and consumer scheduling; this table supplies its shared numerical intermediates.

Four prognostic entries are preallocated so normal nonlinear RHS evaluations add no allocation. Wider diagnostic/component requests can acquire additional entries on first use, bounded by five views, five components and nine normalized fields. They are retained through the evaluation without eviction. Additional capacity is included in the kernel's reported storage and reused subsequently. This is not a claim that first execution of every arbitrarily broad diagnostic plan is allocation-free.

`WVVariableExecutionOptions::sharedFieldGradients` defaults to true. Setting it false retains independent reconstruction for C++ numerical qualification. Both runtime physical-storage policies use the same optimized algorithm. The user explicitly prioritized this shared pipeline over preserving the previous low-memory growth target; no second low-memory derivative implementation is introduced.

## Operations and arithmetic

For horizontal derivatives, prepare the modal field and vertical transform once, then multiply the horizontal spectrum by ik or il into disposable scratch before inverse transformation. Each required physical field still gets its inverse FFT. Cached spectra never alias mutable projection scratch, and horizontal transforms preserve their input.

Hydrostatic vertical derivatives reuse the base modal coefficients, retaining the complementary F/G basis and the original equivalent-depth and N2 normalizations.

Boussinesq uses the existing finite-resolution v4 vertical operators:

- F derivative: projectF, reconstructG, then multiply by -N2/g.
- G derivative: projectG, divide modal coefficients by h0, then reconstructF.

These real vertical matrices and scalings are independent of horizontal mode. Consequently the discrete vertical operator commutes with horizontal synthesis. Applying that same operator to the retained spectrum before inverse FFT avoids packing every physical-grid column through repeated vertical-calculus batches. It does not substitute an analytical wave-basis derivative. Arbitrary caller grids and tracers still use their original full-grid calculus, preserving frequencies outside the model's retained spectrum.

Moving horizontal multipliers and the Boussinesq discrete vertical operation changes floating-point order. Qualification compares against the independent C++ path and MATLAB at the existing tolerances; equality of integration controls/step decisions is checked separately. Mathematical commutation is exact; floating-point equality is not assumed.

For one nonlinear RHS with all prognostic fields, shared assembly produces four fields in either model. Hydrostatic uses seven reconstruction operators plus three projection operators. Boussinesq uses eight base reconstruction operators, eight discrete-z operators and eight flux projection operators. Actual wrapper counts, separate from logical field requests, verify those totals. The original Boussinesq path additionally ran its discrete-z operators in batches covering all physical horizontal columns.

## Evidence and verification

Artifacts are retained in `OceanKitRepositories/wave-vortex-model-benchmark-artifacts/shared-gradient-pipeline-20260911`. The preceding EddyTide reprofile is in the sibling `hydrostatic-post-speed-phase-20260911` archive. The Boussinesq fixture is a scientific-variable-preserving copy of the retained 256 × 256 × 129 composite source, excluding dense, particle and tracer observer groups. Source hashes and the derivation script are retained. The production experiment was not changed.

The Boussinesq steady-state baseline profile attributes 40.5% of main-thread samples to vertical-calculus self work, with no active tracer. A first window sampled matrix loading/preparation and is excluded; both traces remain archived. This profile motivated applying the same discrete operator before horizontal synthesis. Sampling percentages are attribution estimates, not end-to-end speedup forecasts.

The first Hydrostatic exploratory screen passed all scientific output comparisons and identical state decisions; three integration ratios had geometric mean 0.87117 against the qualified speed/phase executable. Additional owned peak was 33,678,880 bytes. Final qualification supersedes this short screen.

Independent lifecycle/storage review and common-cache tests cover view identity, teardown, failure recovery, bounded growth and hot-path allocation. Family tests compare shared and independent value/x/y/z reconstruction across components, layouts and worker settings, plus scoped and standalone nonlinear producer counts. An inverse failure after successful vertical preparation verifies retry reuse. Vertical backends return void noexcept, so a backend-status failure cannot be injected without redesigning that interface; partial cache readiness is tested directly.

## Frozen qualification

The qualified candidate is `aefc6484`, integrating the merged PR #474 prerequisite. Its compiled implementation is unchanged from the combined verification commit `60516944`; the intervening commit refreshes source-linked JSON receipts. The baseline is the frozen qualified speed/phase executable from `8c2b816d`. Protocol, executable/input/source hashes, individual comparisons and raw measurements are archived. The final campaign ran once, for 16 minutes 24 seconds, with two warmup and eight measured pairs per fixture, alternating execution order, native FFTW and one execution thread. No competing builds or benchmarks ran during measurement.

Ratios below are paired geometric means, candidate/baseline; smaller is faster. Intervals are paired bootstrap 95% intervals. Process lifetime includes startup and output as well as integration.

| Workload | Integration ratio (95% interval) | Integration time change | Process lifetime ratio | Additional owned peak |
| --- | --- | --- | --- | --- |
| EddyTide 256 × 256 × 28 | 0.8775 (0.8718–0.8833) | 12.3% less | 0.8830 | 33,679,008 bytes; 6.69% |
| Constant nonhydrostatic composite | 1.0150 (0.9836–1.0421) | 1.5% more | 1.0029 | 80 bytes; metadata only |
| Larger variable Hydrostatic composite | 0.9521 (0.9419–0.9609) | 4.8% less | 0.9686 | 90,820,768 bytes; 4.09% |
| Boussinesq 256 × 256 × 129, no observers | 0.2885 (0.2869–0.2901) | 71.1% less; 3.47× throughput | 0.9570 | 126,893,408 bytes; 2.45% |

All scientific comparisons passed at the existing 1e-10 relative / 1e-12 absolute tolerances. Integration controls and step decisions were identical, and frozen files were unchanged after the campaign. No integration or complete-lifetime ratio exceeded the 3% investigation threshold. The constant control's interval includes no change; it is not evidence of a speedup. Boussinesq median integration fell from 2.169 to 0.625 seconds. Its roughly 34 seconds of matrix loading/preparation dominate this short continuation's process lifetime, so the integration gain must not be reported as whole-process acceleration.

Both physical-field policies use the shared numerical pipeline. A separate four-fixture low-memory correctness smoke passed against the retained final reuse outputs, including scientific values and step decisions. Low-memory timing and the former 3% storage-growth target are not qualified in this increment; `lowMemoryPassed: null` in the timing summary records that omission explicitly.

## Verification ledger

- Native Release: 58/58 tests, 21.19 seconds; ASan/UBSan: 58/58, 35.70 seconds. Apple LeakSanitizer is disabled; allocation tests and Linux leak checks remain enabled.
- GCC 14: both MM kernel and runtime tests passed. One raw-flux test compared separately rounded sums exactly; its one-ulp discrepancy is now bounded by eight machine epsilons times the summed absolute advection contributions. Borrowed-field equality remains exact and scientific tolerances are unchanged.
- MATLAB parity: 12 focused methods across both MM families, supported providers and layouts.
- Source-linked forward integration: six configurations under reference and native FFTW, plus 15 catalog/compatibility methods. Receipts refreshed once after runtime review.
- Documentation check: 2,026 files, 4,145 routes, zero failures and zero generated differences. No MATLAB scientific source or released package payload changed.
- Independent lifecycle and family review used two Sol agents with the coordinator owning builds, numerical design and qualification. No unresolved correctness findings remain.

The repository retains compact [performance](../.github/ci-evidence/shared-gradient-pipeline/performance.json), [verification](../.github/ci-evidence/shared-gradient-pipeline/verification.json) and [archive hashes](../.github/ci-evidence/shared-gradient-pipeline/archive.json). Large fixtures, profiles, binaries and logs remain in the local benchmark archive. Documentation/evidence-only commits after the frozen candidate do not invalidate these numerical results. Required hosted CI remains a merge gate.
