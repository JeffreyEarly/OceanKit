# Post-cache C++ profile, 2026-09-11

Follow-up: the [Hydrostatic assembly investigation](HYDROSTATIC-ASSEMBLY.md) implements and qualifies the recommendation below, records the vDSP screen, and identifies the next pipeline experiments. This page remains the pre-optimization profile.

This exploratory profile selects the next bounded optimization after PR #471. It follows [the C++ optimization workflow](CPP-OPTIMIZATION.md). No runtime or MATLAB scientific source was changed, and this is not a new speedup qualification.

## Source and workload

The inspected v4 main revision is `0512a2d2a6238d31f91f60daef3427145a6b569d`. The executable is the frozen qualified candidate from #470, SHA-256 `f1a642a4fe53176f257b8a5bd8b8bf9ad3cfa5481fdc4d9c892b19d13cdaa3e9`. Every C++ source/header and CMake input recorded by the final campaign was checked against merged main before execution. The embedded build commit predates merge; the exact compiled-input comparison establishes its relevance to main.

The original frozen fixtures were copied into disposable continuations, with their original SHA-256 verified before and after execution. The integration interval was extended to obtain a steady sampling window, and primary grid output moved to the extended endpoint. Original initial dense/passive schedules were retained. The production experiment was not changed. These runs do not represent normal output-heavy production timing.

- EddyTide: 256 × 256 × 28, adaptive RK78, 300-second maximum step, 120 accepted steps, no rejection, 1,560 RHS calls; nonlinear advection plus adaptive damping.
- Constant nonhydrostatic composite: 256 × 256 × 129, fixed RK4, dt 0.5, end time 40, with tracer and particles.
- Variable Hydrostatic composite: 256 × 256 × 129, fixed RK4, dt 0.5, end time 12, with tracer and particles.

All use default `reuse`. The variable-family policy reports compact split views, Accelerate vertical multiplication, 12 horizontal outer workers and eight pointwise workers, with one FFTW internal worker. Thus request `threads: 1` does not mean the complete variable-family execution uses only one CPU. The constant fixture retains its qualified requested FFT thread count of 16.

## Findings

macOS `sample` requested 1 ms intervals over 12-second windows. The table uses exclusive samples on the main thread, separating a function's own work from its callees. The second EddyTide window ended when the process completed and has fewer samples. Waiting worker stacks are not counted as useful CPU work, and sample fractions are not precise elapsed-time measurements.

| Workload | Main-thread samples | Finding |
| --- | ---: | --- |
| EddyTide | 8,326 | `WVTransformHydrostaticKernel::reconstruct` self work: 28.93%; prepared vertical multiplication including callees: 7.24% |
| EddyTide confirmation | 4,903 | Reconstruction self work: 29.21%; prepared vertical multiplication including callees: 7.71% |
| Constant composite, steady window | 7,875 | Main-thread condition-variable waits: 17.66%; scalar-derivative preparation self work: 7.50%; nonlinear-flux routine self work: 6.78%. Worker stacks prominently execute FFTW. |
| Hydrostatic composite | 8,770 | `verticalCalculus` self work: 27.78%; reconstruction self work: 10.95%; prepared vertical multiplication including callees: 9.12% |

The EddyTide report attributes 94.25% of integration time to the RHS. The first continuation completed in about 43.27 seconds of integration, but sampling perturbs execution and this number is descriptive, not an uninstrumented benchmark. All completed runs report zero default-policy duplicate executions. The profile does not indicate a need for another cache architecture change.

In [Hydrostatic reconstruction](../CompiledKernel/src/WVTransformHydrostaticKernel.cpp), the modal coefficient loop is serial and performs field selection, component selection, complex coefficient/phase arithmetic, derivative scaling, and writes into prepared compact storage. Some derivative normalization also follows the vertical multiplication. The function-level profile makes this a strong candidate; a short targeted measurement should distinguish those loops before changing them. It does not prove every self sample belongs to one source line.

The large tracer case exposes a separate cost in [Hydrostatic vertical calculus](../CompiledKernel/src/WVTransformHydrostaticKernel.cpp): clearing and gathering full physical columns into complex spectral storage, scaling, and scattering around the prepared matrix operators. This is a follow-up for tracer-heavy workloads, not a reason to expand the EddyTide experiment.

## Recommended next increment

Optimize Hydrostatic modal coefficient assembly, beginning with the measured EddyTide reconstruction path. First measure the assembly and post-matrix scaling loops separately. Then test branch specialization and/or distributing independent coefficient ranges over an existing prepared executor. Preserve per-coefficient arithmetic order, phase/component semantics, prepared storage and both evaluation policies. Avoid adding full-volume buffers or revisiting FFT/matrix algorithms.

Use one implementation owner and, when delegation is authorized, one independent reviewer/test worker. A 30–45 minute investigation checkpoint should establish whether the candidate merits qualification. A 5–10% complete-EddyTide integration improvement is a reasonable screening objective, not a forecast. Halving a cost accounting for roughly 29% of main-thread samples suggests about 14–15% potential integration-time reduction under a simple critical-path approximation; actual scheduling and bandwidth effects must be measured. Extend to Boussinesq only after confirming the analogous loop is material there.

Run affected compiler, kernel/runtime, component/derivative, policy and MATLAB parity checks before freezing a candidate. Then use the established EddyTide and larger fixtures for uninstrumented paired timing and memory qualification. Preserve their existing numerical tolerances and investigate regressions above 3%. If assembly is not materially improved, record the result and move to the tracer vertical-calculus path rather than broadening the experiment indefinitely.

The existing spectral-kernel-benchmarks [handoff](https://github.com/JeffreyEarly/spectral-kernel-benchmarks/blob/6f2e4cd9e2f092433a584645e573d818ecaca42b/docs/wvm-compiled-core-integration-handoff.md), section 1.4, explicitly excludes phase and coefficient construction from the variable-stratification benchmark boundary. Its selected matrix representation, pruned tile-16 FFT algorithm and worker guidance remain the baseline. The retained [#455 adoption plan](../.github/planning/issue-455-complete-variable-adoption.md) already excludes repeated primitive packing surveys; this proposed WVM assembly experiment does not repeat them.

## Evidence and limitations

Raw profiles, requests, protocols, reports, disposable outputs and the two profiling scripts are retained in `OceanKitRepositories/wave-vortex-model-benchmark-artifacts/cpp-post-cache-20260911`. The compact `summary.json` SHA-256 is `e3388df97a996acfd4a3c3e0dbd2a97d1b1fff0b3e1cf8fb37f866967b37e3a9`.

The first constant sample overlapped preparation and is retained but excluded from the steady-state findings; its replacement delayed sampling eight seconds. The first window of the second EddyTide run briefly overlapped the preceding continuation and is also excluded. Its selected `sample-idle.txt` was attached after the other process exited; the retained stacks cover integration. The second run's complete-process time is not used for a comparison. Original fixtures remained unchanged, reports completed successfully, and the numerical implementation is the already qualified executable. No new before/after numerical or speedup claim is made.

Verification ledger: compiled input and frozen fixture hashes matched; all disposable continuations completed and reported zero default duplicate executions; sampled call-tree accounting had no negative exclusive counts; local Markdown links, repository/compatibility checks and whitespace review passed. `docs:check` passed with no generated changes. No numerical implementation changed, so native numerical suites and another paired performance campaign were not repeated for this documentation increment.
