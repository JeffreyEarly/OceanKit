# WaveVortex portable runtime

This directory contains the optional MATLAB-independent portable runtime. MATLAB remains WaveVortexModel's primary interface. The qualified `WVModel`, field, observer, NetCDF, and restart composition supports constant stratification in hydrostatic and nonhydrostatic modes, `WVTransformBarotropicQG`, `WVTransformStratifiedQG`, `WVTransformHydrostatic`, and `WVTransformBoussinesq` behind transform-neutral state and service contracts. Barotropic QG retains MATLAB's compact `A0`-only `kl` layout through integration, output, append, and restart; it never allocates compatibility `Ap` or `Am` state. The runtime is a source-only checkpoint-to-checkpoint tool for advanced users and the v1 source-level surface for statically linked C++ observing-system, output-schedule, and forcing implementations.

The runtime supports:

- fixed-step RK4, MATLAB `ode23`-compatible Bogacki--Shampine integration with a third-order accepted solution and second-order embedded estimate, MATLAB `ode45`-compatible Dormand--Prince integration with a fifth-order accepted solution and fourth-order embedded estimate, and MATLAB `ode78`-compatible Verner integration with an eighth-order accepted solution and seventh-order embedded estimate; all three adaptive methods include method-owned continuous output, with RK78 using Verner's seventh-order extension;
- all twelve stable supplied forcing identities on their MATLAB-applicable baseline transforms, including explicit antialiasing, horizontal/vertical damping, vertical diffusivity and constant-stratification narrow-band forcing;
- the five built-in observer records (`WVCoefficients`, `WVEulerianFields`, `WVMooring`, `WVLagrangianParticles`, and `WVTracer`) on the transform-specific layouts documented below;
- MATLAB-compatible checkpoint and time-series NetCDF data for the documented transform, observer, and diagnostic contracts.

Arbitrary MATLAB forcing or observing-system subclasses are not supported. Multi-file, named-group observing-system output is authored in MATLAB and executed by the command-line program without reproducing that scientific configuration in a second format.

`WVModel` is the thin, move-only runtime façade. It retains the frozen extension catalog together with the resolved forcing, observers, numerical system, integrator, output evaluation, driver, and sink; `WVModelState` separately owns canonical coefficients and dynamic observer state. The same façade backs the standalone runner and the production MEX right-hand-side path. It adds no numerical algorithm or state-sized copy.

Library clients may configure output with the stable source-level C++ `WVModelOutputFile` and `WVModelOutputGroup` builders. The builders emit canonical records and are discarded; `WVModelOutputConfiguration::compile` then applies the same record-preserving path used by `WVModel::createFromModelOutputFiles` for a restored MATLAB-authored sibling set. The resulting descriptor is the sole compiled graph shared by the plan, output evaluator, and NetCDF sink. Destination remapping changes only paths keyed by stable file identifier and preserves group membership, restart designation, complete observer configuration, declared schemas, and schedule identity, version, configuration, and typed cursor.

Output schedules are resolved by exact source-linked identity and version. `WVOutputPlan` retains one immutable schedule per group, while `WVOutputDriver` owns one small continuation cursor and one cached occurrence per group. Continuation cursors are distinct from sink-owned destination progress, which records committed record counts, time-last evidence, full typed cursors, and physical and committed ragged-axis offsets. A continuation advances only after its route is committed, so a failed coincident route can be replayed with the exact cursor without repeating successful routes. Delivery records expose only the latest selected event and retry attempts; cumulative history is aggregated in metrics, keeping storage independent of the event-window length. The legacy evenly spaced schedule retains its existing NetCDF scalars and ordinal semantics; new algorithmic providers persist a typed configuration and a bounded cursor envelope. `WVOutputSchedulePayloadSchema` resolves finite-real, integer, and Boolean occurrence fields once, and `WVOutputSchedulePayload` carries their values in fixed 4 KiB text-free storage addressed by numeric slot. State-triggered schedules are reserved but unsupported.

## Source API v1

The [variable evaluation lifecycle](VARIABLE-EVALUATION.md) describes shared RHS/output scopes, the default `reuse` and optional `low-memory` policies, and the scoped field-query API.

`wave-vortex-portable-source-api-v1`, version 1.0, identifies the C++ extension and embedding surface under development: source-API identity constants; explicit catalog construction and capability discovery; `WVObservingSystem`, `WVOutputSchedule`, and `WVForcing`; their construction/planning records and data-only `WVObservationSchema`/`WVObservationBatch` boundary; `WVModelOutputConfiguration`; `WVModel`; and `runWaveVortex()`. Legacy MATLAB-encoding adapters, declarations under `PortableRuntime/src`, `detail` namespaces, and undocumented implementation declarations are not extension APIs.

Consumers select one WaveVortexModel checkout and compile the runtime, application-owned catalog, statically linked extensions, and runner together. Rebuild the complete application whenever that checkout changes. The currently unused C++ API may change without backward source compatibility; it promises no binary ABI across commits, compilers, or build configurations. MATLAB APIs and persisted scientific behavior retain backward compatibility. It provides no dynamic discovery, separately loadable plug-in, or distributed runtime/extension binary.

`WVObserverFactoryRegistration` has exactly five public constructor inputs: type identity, contract version, factory, optional data-only configuration resolver, and optional data-only output-plan resolver. No legacy operation callback or persistence metadata is exposed. `WVExtensionCatalogBuilder` is the only mutable registration owner: register everything before freezing, treat a rejected registration as invalidating the builder, and freeze exactly once. Frozen catalogs are immutable shared owners retained by descriptors, output configurations, models, MEX handles, and runners. The builder and the caller's original catalog handle may then be destroyed, and independent catalogs may bind the same identities differently in one process.

Source API versioning is independent from exact data contracts. `wave-vortex-portable-pair-v1`, each paired type identity/version, the portable observer-graph schema, every schedule identity/version/configuration/payload schema/cursor, every observation-schema identity/version, `wave-vortex-run-request-v1` or `wave-vortex-run-request-v2`, and the compiled-kernel contract are validated separately. Matching one never relaxes another.

The portable reference runtime and source-linked consumers are qualified on Ubuntu with GCC or Clang and on macOS with AppleClang. The optimized native FFTW runner is Apple-silicon-only. Windows/MSVC source-linked builds are unsupported. The [current AlongTrackSimulator consumer](https://github.com/satmapkit/AlongTrackSimulator/commit/511aa6af9c60353b3de4d371dfe0df43159027cf) supplies the external application-owned catalog and reusable-runner proof. Source API qualification recompiles that consumer against the explicitly selected WaveVortexModel checkout. The older ATS #4 commit recorded in source-selection provenance is a historical baseline, not the current consumer source.

## Standard compatibility matrix

The [generated compatibility matrix](COMPATIBILITY.md) assembles the MATLAB-authoritative forcing, variable, integration, observer, sampling, output, and persistence contracts. It links exact support or rejection records to test fixtures and keeps execution evidence separate from catalog completeness. Final standard-parity qualification is recorded under [issue #307](https://github.com/JeffreyEarly/wave-vortex-model/issues/307).

## Forward-integration support and evidence

The [forward-integration slice](INTEGRATION.md) records the six transform configurations and five execution profiles: explicit and CFL-selected RK4, RK23, RK45, and RK78. It links each support row to method/configuration witnesses and identifies the representative complete-model continuation cases. Source references establish test coverage; executed qualification receipts separately record outcomes, provider, tolerances, and source provenance. This slice is an input to the complete compatibility matrix under #306.

Execution is forward in time. Backward integration, execution of arbitrary MATLAB function handles in C++, and experimental adaptive-cell integration are outside this supported surface. MATLAB-authored opaque `N2Function` payloads remain supported scientific provenance and are preserved through restart; the C++ runtime does not execute those handles. Source-linked C++ implementations continue to use the documented catalog and integration contracts.

Spline sampling preserves MATLAB’s transform-specific boundary conventions. On small horizontal grids, a circularly shifted query can remain beyond the final grid knot: the three-dimensional transform sampler returns zero, while Barotropic QG extends the end spline polynomial. The MATLAB small-grid behavior is tracked separately in [#419](https://github.com/JeffreyEarly/wave-vortex-model/issues/419).

## Controlled termination

Source-linked applications can pass a `WVIntegrationControl` to the additive `WVModel::advanceToTime` overload. The model state remains caller-owned; the returned `WVModelAdvanceResult` contains status, termination details and a complete metrics snapshot, including work completed while draining a stop request.

```cpp
WVIntegrationControl control;
control.shouldStop = [&stopFlag](const WVIntegrationProgress &) {
    return stopFlag.load(); // for example, an application-owned atomic<bool>
};
auto result = model.advanceToTime(state, requestedFinalTime, initialStep, control);
if (!result) {
    // result.termination distinguishes integration, output and callback failure.
} else if (result.termination.stopped()) {
    // state and result.termination.finalAcceptedTime describe an accepted state.
}
```

The callback is borrowed for this invocation and runs synchronously at the initial accepted state, after an accepted step, or after a complete output occurrence across all coincident routes. It must not mutate the model, integration state or output configuration, or re-enter the driver. No callback runs inside RHS evaluation, individual element loops, a rejected step attempt or an incomplete output occurrence. A throwing callback is a failure, including nonstandard exceptions. Output failures take precedence over a pending stop; retrying the same plan and sink completes only the failed route, preserves committed cursors and retains the stop request.

A stop requested by a generic sink or control first drains every output occurrence covered by the current accepted step. Interpolated output is never returned as the accepted model state. With an attached model-file sink, termination additionally waits for an authored checkpoint occurrence that contains all required dynamic observer state in **each destination**, possibly in coincident sibling groups. The driver may take additional steps to that boundary and clips the last such step. This keeps coefficients, particles, tracers and schedule cursors immediately readable by MATLAB from any individual file and by C++ from the complete graph. Authored schedules and the existing NetCDF encoding are preserved. Empty observer collections carry MATLAB's existing `AnnotatedClassArray` marker.

Create/replace destinations have no committed records merely because their source cursor points at an earlier checkpoint: an immediate stop may drain to their first actually persisted authored occurrence. An append destination with a complete current checkpoint can stop before taking a step. Existing floating-point time-coincidence tolerances apply to accepted endpoints and scheduled times.

Checkpoint spacing therefore limits graceful-stop latency. A stop request is latched at its original boundary and time; the actual accepted stop time can be later. A bounded lookahead uses temporary schedule cursors without changing execution cursors or retaining an occurrence history. If no common complete checkpoint occurs at or after the accepted state and before the requested final time, the run returns an explicit **output failure**; it never exceeds the requested target or claims a restartable controlled stop. Failure retains the accepted in-memory state and already committed records, with the same retry/failure limitations as other output errors.

The existing overloads remain source compatible. Built-in controlled integration preserves each method's step-selection behavior, including `ode45`/`ode78` final-step handling. Calls without a callback retain their original numerical paths. Control adds no state-sized workspace; the output driver retains bounded scalar termination bookkeeping, and stop-only lookahead storage is reported in `controlledStopWorkspaceMaximumLiveBytes`. Source-linked custom integrators that do not implement the additive controlled overload continue to compile; requesting output-free controlled integration on one returns `unsupportedOperation`. Their existing accepted-step interface already supports control through `WVOutputDriver`.

The standalone executable maps **SIGINT** to this same stop path. The first interrupt requests graceful termination; a second SIGINT uses normal process termination. SIGTERM retains its ordinary abort behavior. The reusable `runWaveVortex` API installs no signal handlers and accepts optional application-owned control. Successful controlled runs return exit code 0 and JSON `status: "stopped"`; the `termination` object records reason, requested boundary/time, actual accepted time and callback count. Callback failures use exit code 7, output failures use 6, and numerical/integration failures use 5, with distinct structured failure stages and reasons. Coefficient-only output writes the final accepted checkpoint; explicitly scheduled checkpoint output contains only the committed scheduled occurrences through the accepted endpoint.


## Build

A portable reference build requires CMake 3.20, a C++17 compiler, and NetCDF C:

```sh
cmake -S PortableRuntime -B build/portable -DCMAKE_BUILD_TYPE=Release
cmake --build build/portable --parallel
ctest --test-dir build/portable --output-on-failure
```

On Apple silicon, the optimized runner can be built with:

```sh
PortableRuntime/buildWaveVortexRun.sh build/native
```

The script verifies the pinned FFTW 3.3.11 archive, builds it in the ignored `.compiled-backend-cache`, and writes the example executable to `build/native/wave-vortex-run`. Omitting the build-directory argument uses `.compiled-backend-cache/runtime-build`. WaveVortexModel distributes no FFTW archive, library, MEX file, or executable. Redistributing a locally linked executable requires compliance with FFTW's GPL license.

For a new CMake cache, `WV_RUNTIME_ENABLE_COMPACT_VARIABLE_POLICY` defaults to `ON` only when configuring on Apple silicon with native FFTW and Accelerate enabled. It defaults to `OFF` for portable, reference-only and other ineligible configurations. An explicit cached `ON` or `OFF` remains authoritative; an ineligible explicit `ON` is rejected. For native SQG, Hydrostatic and Boussinesq execution, the compact policy selects split matrices, pruned horizontal transforms and streamed physical tendencies. Horizontal workers follow available performance cores up to twelve, pointwise workers use up to eight, and FFTW internal threads stay at one. The stages execute in sequence, and the runner does not change process-wide BLAS settings. The `variableKernelPolicy` report records the selection and effective topology.

An explicit `--threads` value greater than one retains its meaning as FFTW internal threads and selects the established full-FFT/interleaved path with Accelerate matrices. Omit that option to use the compact policy when enabled. Reference builds and non-variable transforms retain their existing policies. This runner selection does not expand the public MATLAB compiled preview's supported transform families; C++ applications can separately supply `WVVariableKernelServices` when constructing a persisted model.

## Density diagnostic foundation

`WVNoMotionProfile` provides density, inverse material height and APE for an explicitly supplied strictly monotone density profile. It does not fit a profile to a three-dimensional model state. Include [`WVNoMotionProfile.hpp`](include/WaveVortexRuntime/WVNoMotionProfile.hpp) and link the existing portable runtime library. Successful construction owns the knots and normalized shape-preserving cubic coefficients; failed construction preserves the previous object. Scalar queries allocate no workspace and leave their caller-owned result unchanged on failure. Retained storage scales with the number of vertical knots, not the number of parcel queries.

Density must strictly decrease as height increases. Plateaus have no unique inverse and are rejected. Heights remain inside the profile domain; density inversion permits and clamps only the same eight-ulp endpoint allowance as MATLAB. Nonfinite inputs and finite inputs whose normalization or cubic arithmetic cannot be represented fail explicitly. An uninitialized object rejects queries.

For total parcel density, `inverse` returns material height $s$, giving displacement $\eta_{\mathrm{true}}=z-s$. `availablePotentialEnergy` retains the derivative-integral form

$$\mathrm{APE}(z,s)=\frac{g}{\rho_0}\int_s^z(r-z)\rho'_{\mathrm{nm}}(r)\,dr.$$

It integrates the local cubic over only the crossed intervals. This avoids subtracting nearby densities or large hydrostatic-pressure terms, preserving tiny-displacement energy. With heights in metres, density in kg m-3 and gravity in m s-2, the result is energy per reference mass in m2 s-2.

Run-request v2 accepts an optional runtime selection inside `execution`:

```json
"densityDiagnostics": {
  "contract": "wave-vortex-density-diagnostics-v1",
  "reference": "actual"
}
```

The object and either member may be omitted. Omission selects `actual`: the corrected current-distribution `rho_nm` contract, whose recovery algorithm is fixed to `dampedLeastSquares` by this contract version. Explicit `reference: "initial"` selects the coherent initial-profile approximation, which needs no current-profile recovery. Unknown contract versions, reference names and solver fields are rejected; legacy MATLAB solver choices are never inferred from restart files. V1 request syntax is unchanged and rejects this v2 extension. MATLAB's existing request writer needs no change to obtain the corrected default; an explicit alternative is authored in the v2 JSON request.

The runner evaluates full-natural-grid `rho_nm`, `eta_true`, `ape` and `apv` through the existing field and observer services for the four wave-bearing transform configurations, with either antialias setting. Its `densityDiagnosticContract` report identifies the selected reference, recovery policy and default/request provenance, with `outputEvaluation: "available"`. Actual `rho_nm` requests always invoke current-density recovery, including under an initial-reference selection. `profileRecovery` reports dampedLeastSquares when recovery ran; `selectedProfileRecovery` separately identifies whether the selected displacement reference requires recovery. The `densityEvaluation` report records actual recovery, inversion, APE/APV work, reuse and event storage. The selection is runtime-only and must be supplied again on a later run to preserve an explicit alternative; no MATLAB restart schema or saved flag is changed.

`WVNoMotionProfileRecovery`, declared in [`WVNoMotionProfileRecovery.hpp`](include/WaveVortexRuntime/WVNoMotionProfileRecovery.hpp), provides the separate current-distribution recovery primitive. `recover` borrows a `WVRealVolumeConstView` in MATLAB's x-fastest `(Nx,Ny,Nz)` layout, bottom-to-top integration weights, domain depth and the original reference profile. It returns an exactly horizontally uniform stable profile directly. Otherwise it extracts the current density extrema and normalized raw moments, then fits strictly ordered density nodes using the original normalized shape as the starting guess. Weights are divided by depth without further renormalization. The recovered profile can be passed to `WVNoMotionProfile::create` with the matching height grid.

`moments` exposes the compact distribution summary, and `fitMoments` accepts supplied moments and current extrema for qualification or reuse. A stable summary has a nonempty `stableProfile` and empty `moments`, avoiding an unnecessary moment calculation. Nonuniform extraction retains only a 4,096-sample tile and O(Nz) vectors; it never copies the volume. The fit allocates O(Nz²) workspace once and solves its damped augmented system with Householder QR. It does not form normal equations or add a linear-algebra dependency.

`WVNoMotionRecoveryReport` records the algorithm, exit flag, termination reason, iteration/evaluation counts, accepted/rejected steps, costs, gradient and maximum residual. `qualified` is true only when the solver has a positive exit, normalized maximum moment residual at most 1e-8, and finite strictly decreasing physical output at the supplied current endpoints. A small step alone cannot certify the result. Default budgets are 2,000 iterations and 5,000 evaluations, with 1e-12 gradient, step and relative-cost tolerances. Failed recovery preserves the caller's previous profile while returning diagnostics. `workspaceBytes` reports peak owned vector capacity, including unpublished profile staging and partial work before failure, and excluding borrowed input, caller output and object/allocator metadata; it is not process RSS. Optional workspace outputs on `moments` and `WVNoMotionProfile::create` expose their corresponding construction peaks without changing the numerical operations.

The private `WVDensityEventEvaluation` helper joins recovery and profile calculus inside the existing public field evaluator. Field plans accept full-natural-grid `rho_nm`, `eta_true`, `ape` and `apv`; interpolation and reductions remain unsupported. One output event owns reconstructed total density and the geometry it borrows into the helper. Sequential consumers share completed recovery, profile construction and material-height inversion, even if their individual plans are destroyed between calls. Direct evaluation creates an invocation-local scope; subsequent invocations always rebuild derived density state, including when coefficient addresses and times are unchanged. Ordinary density plans request the same event lifetime without requiring forcing or moving observers. Prepared NetCDF payloads survive failed sibling routes without retaining the numerical helper or repeating density work. An explicit density contract passed to field-plan/model-output construction is preserved by observer rebinding.

`rho_nm` always denotes the actual recovered profile. The explicit initial reference selects the profile used for displacement, APE and APV: mixed requests can recover actual `rho_nm` while evaluating displacement/APE/APV against the initial profile. An initial-reference request without `rho_nm` does not recover a profile; a profile-only request does not invert parcels. Reference selection is immutable within a resolved plan and preserved by rebinding. Mixing reference selections within one active event is rejected. All active output views are validated before work, and their values are published only after every requested numerical stage succeeds. Event counters expose recovery, profile construction, inverse, APE and APV passes and reuse; live storage returns to zero when the event or direct invocation ends. Density-only evaluation avoids retaining the unrelated four-volume primitive cache.

APV uses `zeta_z - zeta_x*diffX(eta_true) - zeta_y*diffY(eta_true) - (zeta_z+f)*diffZG(eta_true)`. It shares the selected displacement, applies the existing transform-specific G-grid derivative, and retains one APV result per event. Constant stratification uses the existing FFT/DST/DCT plans with retained vertical modes and zero first derivatives at horizontal Nyquist modes. Event peaks describe owned vector capacities; the runner combines component peaks as a storage upper bound, not process RSS or a measured simultaneous peak.

Matching normalized moments does not imply bitwise identical recovered nodes: the raw-moment inverse problem is ill-conditioned, and floating-point reductions and QR arithmetic affect the final fit. Qualification therefore compares physical density, material-height and APE differences in addition to the solver's residual and termination. The retained JAMES cases use one-part-per-million density-range and depth bounds for C++/MATLAB profile comparisons; independent small-distribution tests retain their separate analytic tolerances. These are cross-implementation checks, not estimates of error relative to an exact continuous sorted state.

## MATLAB-authored run bundles

A portable run has two parts:

- one or more NetCDF files authored by MATLAB, which remain authoritative for model configuration, state, forcing, observers, output groups, schedules, and restart progress;
- a small JSON request, which selects integration, execution, destination paths, and the report location.

For a new in-place continuation, first close the restart-capable MATLAB model-output file, then author the request with the metadata-only writer:

```matlab
WVModel.writePortableRunRequest("run.json","initial-condition.nc",finalTime=86400);
```

After building the native runner from the repository root as shown above, the resulting request runs without editing by invoking the executable at its generated path:

```sh
.compiled-backend-cache/runtime-build/wave-vortex-run --request run.json
```

Omitted v2 controls select MATLAB `ode78` as `adaptive-rk78`, relative tolerance `1e-3`, absolute-tolerance scale `1e-6`, an initial step from the minimum advective/oscillatory CFL `0.5` candidate after state restoration, a maximum step equal to one tenth of the continuation interval, native FFTW, and an automatically hardware-bounded thread count. These reproduce standard `WVModel` defaults rather than custom MATLAB-session settings that were never persisted. Every explicit option overrides its corresponding default.

Run-request v2 also accepts fixed RK4 with a positive explicit `initialStep`; fixed RK4 with positive `cfl` and an `advective`, `oscillatory`, or `min` constraint; and explicit MATLAB `ode23`, `ode45`, or `ode78` controls. The adaptive methods are serialized as `adaptive-rk23`, `adaptive-rk45`, and `adaptive-rk78`. [`examples/portable-run-request-v2.json`](examples/portable-run-request-v2.json) shows the minimal generated document. List the complete sibling NetCDF set and map every stable output-file identifier when using `create` or `replace`. Paths are resolved relative to the request file.

The runner first parses the schema identity and version and routes to the strict matching decoder. Run-request v1 remains unchanged: existing documents retain fixed RK4 and MATLAB `ode23` behavior, reports, and errors, and receive no v2 defaults. The typed method and step policy are resolved once. Raw inspection then reconstructs the complete NetCDF graph as owning data-only records without invoking an extension factory or constructing a runtime implementation. Semantic preflight uses the explicitly supplied frozen catalog to resolve every paired forcing, observer, and output schedule, compares registered data-only output plans with persisted schemas, validates transform CFL capability and the complete destination policy, and compiles the sole output graph. Only then does it construct the FFT provider and state-sized runtime storage. An unavailable native provider fails at that boundary with source-build instructions, before state-sized allocation or output mutation, and never falls back. The JSON cannot add observers, forcings, groups, or schedules. AlongTrackSimulator's `authorAlongTrackPortableRunBundle` is one external package-specific authoring example.

CFL-selected RK4 and the omitted adaptive initial-step default evaluate transform-owned candidates once after restoring the segment-start state. Advective selection combines effective horizontal resolution and `uvMax` with the vertical `dz/w` restriction for three-dimensional transforms. Oscillatory selection uses the highest active wave frequency; transforms without waves report infinity. The selected step is not recomputed during RK stages or output events, and the existing final partial-step behavior is preserved. Reports separate parsing, preflight, provider construction, startup, and CFL evaluation; distinguish requested values from active method, tolerances, steps, provider, and threads; identify the controller, candidates, selected step, and no-fallback status; and retain the existing integration, RHS-work, and exact-storage diagnostics.

`create` and `replace` require a complete destination map and never mutate source files. Their destination progress begins empty even when scheduling resumes from a noninitial source cursor. `append` may use the source destinations with an empty map or a complete remap to an existing compatible file set. Append reconstructs and validates the complete destination set read-only—including graph, schema, record counts, time-last values, full typed schedule cursors, and ragged committed/physical offsets—before any file is reopened with mutation capability. A partial remap, source alias, incompatible graph, unsupported paired implementation, or continuation/destination mismatch fails before output mutation.

## Legacy run and restart

The optimized command-line program requires an explicit FFT provider. Complete-model continuation restores and appends the selected file's supported observing-system graph and schedules:

```sh
wave-vortex-run saved-model.nc \
    --restart-mode model \
    --output-policy append \
    --delta-t 1 --final-time 100 \
    --fft-provider native-fftw

wave-vortex-run saved-model.nc \
    --restart-mode model \
    --output-policy append \
    --integrator adaptive-rk23 \
    --delta-t 1 --initial-step 1 --maximum-step 10 --final-time 100 \
    --relative-tolerance 1e-3 --absolute-tolerance 1e-6 \
    --fft-provider native-fftw

wave-vortex-run saved-model.nc \
    --restart-mode model \
    --output-policy append \
    --integrator adaptive-rk45 \
    --delta-t 1 --initial-step 1 --maximum-step 10 --final-time 100 \
    --relative-tolerance 1e-3 --absolute-tolerance 1e-6 \
    --fft-provider native-fftw

wave-vortex-run saved-checkpoint.nc continued-checkpoint.nc \
    --restart-mode coefficients \
    --output-policy create \
    --integrator adaptive-rk78 \
    --delta-t 1 --initial-step 1 --maximum-step 10 --final-time 100 \
    --relative-tolerance 1e-3 --absolute-tolerance 1e-6 \
    --fft-provider native-fftw
```

For adaptive integration, `--delta-t` remains the backward-compatible initial-step default. `--initial-step` overrides it explicitly. `--maximum-step` defaults to one tenth of the requested continuation interval, matching MATLAB's default bound; name it explicitly when comparing runs or continuing the same controller policy across segments. Use `adaptive-rk23`, `adaptive-rk45`, or `adaptive-rk78` for the corresponding MATLAB controller semantics. RK78 endpoint-only execution retains the 11 state-equivalent base workspace from its accepted solution and embedded estimate and performs no extension RHS evaluation or extension-only allocation. When an accepted step contains an interior output request, the method retains `f1` and `f6` through `f12`, stores the accepted-step initial state in an otherwise-dead base buffer, lazily allocates `f14` through `f17`, evaluates those four stages once, and reuses them for every interior sample in that step. The lazy buffers are released before the next step, interpolated states are written only to the output driver's separate staging state, and changing output times does not change accepted or rejected steps. Reports separate base and extension RHS counts, retained base stages, current and maximum extension workspace, maximum-live state-equivalent arrays, dense evaluations, cache builds and reuse, and wall time. Run-request v1 remains fixed to RK4 and MATLAB `ode23`; run-request v2 adds exact MATLAB `ode45` and `ode78` selection plus explicit or CFL-selected RK4.

Use `--restart-mode coefficients --output-policy create` with positional input and output paths for an explicit reduced checkpoint-only workflow. `create` refuses existing files; `replace` must be named to authorize atomic replacement. The legacy complete-model form accepts only `append` and validates compatibility before mutation. It remains available for scripts that operate on a single output file; `--request` is the preferred complete multi-file boundary and cannot be mixed with legacy semantic flags.

Plans, caches, integrator history, derived forcing operators, and scratch are rebuilt after restart rather than persisted.

## Architecture

`WVIntegrationStateLayout` freezes the transform identity, spatial dimensions, ordered coefficient-family identifiers, each family's natural spectral rank, and any observer-owned state blocks before state allocation. The constant-stratification adapter declares rank-2 `Ap`, `Am`, and `A0` families and retains its stabilized `WVState` views. Transform-neutral integration uses ordered `WVCoefficientFamilyView` values, so a one-family transform owns and advances only that family; it does not allocate dummy `Ap` or `Am` arrays or a state-sized compatibility copy.

`WVBarotropicQGIntegrationSystem` is the concrete one-family implementation. Its transform-specific decoder validates persisted doubly periodic axes and the equivalent depth, mode index, gravity, planetary radius, rotation rate, latitude, antialias flag, and times before publishing spatial shape `[Nx,Ny]` and one rank-1 family `A0[Nkl]`. Its resolved QG forcing engine executes ordinary nonlinear PV advection, adaptive damping, fixed amplitude (including `WVNarrowBandGeostrophicForcing` records), linear and quadratic bottom friction, and beta-plane PV advection in stage, priority, and original-ordinal order. Spatial operations reuse one RHS-scoped primitive-field reconstruction and the kernel's existing 4H+5R scratch; the forcing façade owns no array-sized workspace. The same resolved system supplies `WVModel`, all four integrators, full-grid and position-sampled QG fields, XY particles, rank-2 tracers, transactional multi-file output, append, and restart. Pseudo-topographic wave generation, vertical particles, rank-3 tracers, moorings, and later closures are rejected during allocation-light preflight.

The #281 Donut control compared the unchanged 256-by-256, `j=1`, antialiased nonlinear QG workload with #280 at `e5061a6c`, using the pinned native FFTW provider in three fresh process pairs with two warmups and seven timed RHS evaluations per process. The paired runtime changes were +1.39%, -0.255%, and +0.020% (median +0.020%). Exact retained system storage changed from 6,176,827 to 6,177,214 bytes (+387 bytes, 0.0063%); 118,480-byte compact state storage, 4,734,976-byte 4H+5R scratch, three plans, and zero persistent full-Hermitian storage were unchanged. This is a routing-overhead control, not a comparison with the deliberately slow direct-DFT correctness provider; new forcing workloads perform additional science and are reported separately.

`WVIntegrationSystem` supplies the selected right-hand side, constraints, error scaling, and optional field-evaluation service. `WVTimeIntegrator` advances accepted state, and `WVDenseOutput` evaluates within an accepted step. RK4, RK23, RK45, RK78, and output interpolation traverse the frozen family descriptions rather than assuming three equal rank-2 arrays. Output orchestration depends only on those contracts, so another integrator, coefficient-family layout, or state block does not require a method branch in the driver.

Checkpoint inspection resolves the persisted transform identity before reading transform configuration or inspecting coefficient variables. The selected persistence adapter then publishes its allocation-light spatial and coefficient-rank description; only a successful complete preflight may allocate and load state-sized arrays. Existing WaveVortexModel 4.x constant-stratification NetCDF, run-request v1, observer, forcing, output, restart, and export encodings are unchanged. A transform-neutral owning checkpoint record is available to transform-specific adapters, while each adapter remains responsible for its own encoding. There is no dynamic transform plug-in or native C++ transform-authoring API.

`WVModel` composes these services but does not replace them. High-level restart, RHS, step, integration, output progress, and metrics operations delegate to the same contracts used before the façade. Model-output construction inspects all sibling files together, restores the latest complete compatible state, and compiles the recovered graph directly into the existing plan and sink. It retains shared ownership of the frozen catalog for as long as resolved implementations depend on it.

`WVExtensionCatalogBuilder` is the only mutable extension-registration boundary. Its strongly typed `addObserverFactory`, `addOutputScheduleFactory`, and `addForcingFactory` operations reject duplicate or incomplete registrations. A rejected registration invalidates the builder, mutation after freezing fails, and a builder can freeze only once. `addBuiltInExtensions()` installs the qualified built-ins through this same source-linked path.

Freezing produces one immutable `std::shared_ptr<const WVExtensionCatalog>` with separate `WVObserverCatalog`, `WVOutputScheduleCatalog`, and `WVForcingCatalog` subcatalogs. The caller passes it explicitly through inspection, semantic resolution, output-graph compilation, descriptor construction, `WVModel`, MEX handles, and the runner. Runtime owners retain the shared catalog, so a live model is independent of the builder and the caller's original pointer. Independent catalogs can coexist in one process and can bind the same source identity differently without interference. There is no process-global mutable registry, implicit sealing event, binary plug-in ABI, or distributed compiled extension.

`WVObservingSystem` is the stable source-linked C++ implementation boundary for a paired MATLAB observer. `WVPortableObserverDescriptor` creates a distinct immutable resolved instance for every observer record, even when records share a stateless factory. Each resolved instance owns its typed construction configuration and one declarative execution plan; there are no `records*` or `owns*` observer discriminator calls in runtime consumers. Runtime routes carry resolved-observer pointers and variable ordinals rather than class names. At an event, the resolved observer's `prepareOccurrence` operation fills an evaluator-owned `WVObserverOccurrenceWorkspace` with coordinates, logical extents, sample-time metadata, values, and ragged relationships. Element loops remain in the specialized particle, tracer, and field-evaluation kernels.

`WVFieldEvaluationService::createEventPlan` resolves field identities, dependencies, interpolation operations, and position-set slots during construction. `prepareEventGeometry` accepts each occurrence's coordinate and extent views; `evaluateEvent` samples one prepared occurrence, while `evaluateEventBatch` unions compatible primitive dependencies across coincident occurrences, reconstructs them once, and retains separate interpolation for each geometry. Every catalog field that permits position sampling is supported, including the derived horizontal fields `ssh`, `ssu`, and `ssv`. Primitive reconstruction, interpolation, derived-field evaluation, and scratch remain central; observers receive only resolved sampled values and cannot request whole fields for private interpolation. Version 1 evaluates one model or dense-output state at the scheduled trigger time. Per-sample times are output metadata and do not request additional state interpolation.

`WVObservationOccurrenceIdentity` is destination-independent. During the prepared-event lifetime it borrows the immutable resolved-observer and logical group/schedule records plus the exact typed cursor, payload schema, and bounded payload, so semantic comparison is allocation-free and value-exact across independently compiled or resumed plans. Scalar cursor, payload, geometry, and field-plan fingerprints remain diagnostics, not exact cache keys. The evaluator mints a collision-free owner/generation/slot token as the authoritative in-flight cache key used by the evaluator and sink. A later-route failure therefore retries without preparing or evaluating again and without repeating a committed route. Only compatible semantic occurrences share data, and storage is bounded by immutable plans plus currently in-flight occurrences.

`WVModelOutputNetCDFSink` remains the sole owner of transactional MATLAB-compatible persistence. Its documented create, replace, and append factories require the compiled output plan and complete source/graph preflight before schema discovery or destination mutation. Observer implementations supply a data-only `WVObservationSchema` plus initial and event-scoped `WVObservationBatch` values and never call NetCDF. Schemas cover fixed and unlimited axes, real, complex, integer, Boolean, and text values, coordinate roles, and nested contiguous-ragged row counts or offsets. Batches make borrowed versus owned buffers explicit and carry construction-resolved variable ordinals in the event path. The numerical kernel has no MATLAB, MEX, NetCDF, or Apple API dependency.

The generic NetCDF adapter stores variable-size observations on flat unlimited axes and commits each axis count with the time record. Inspection reconstructs nonlegacy schemas from those declarations; the five legacy fixed-shape observers retain their existing MATLAB names, metadata, dimensions, and restart encoding. Fixed, event-variable, fixed-by-variable, state-coupled, zero-length, and nested-ragged layouts are generic contract capabilities. State-triggered scheduling, multi-state sampling within one occurrence, three-dimensional Barotropic QG particles or tracers, and real ADCP, glider, profiling-float, satellite, shipboard, or additional mooring implementations remain outside source API v1. Arbitrary MATLAB subclass execution, dynamic plug-ins, binary ABI compatibility, and Windows/MSVC source-linked execution are also explicitly excluded; unsupported graphs fail without silent fallback.

Registrations include the exact paired MATLAB/C++ identity and contract version and must be added before the builder is frozen. Observer registrations expose only the exact five v1 inputs described above. Semantic resolution performs no registration discovery, and integration performs no class-name lookup.

The forcing subcatalog maps exact MATLAB identities and versions to construction-time factories and generic persistence schemas. Each factory returns an immutable source-linked forcing implementation that owns typed configuration and derived operators and is called once per forcing stage or constraint pass. Constant-stratification and Barotropic QG implementations remain behind their resolved system boundaries; the engine, integrators, reader, writer, and checkpoint code contain no forcing-class switch, and no string lookup occurs in coefficient or grid loops.

Library and test clients can reuse command-line behavior through `runWaveVortex(argc,argv,catalog)`. The executable's small `main` creates the built-in catalog and delegates to this entry point; a source-linked client can instead freeze and supply its own catalog without changing runner dispatch.

See the website's portable-runtime user and developer pages for the supported compatibility profile and extension boundaries.

## Scientific variable-stratification records

`WaveVortexRuntime/WVStratifiedModalRecord.hpp` exposes the scientific decoder added in #295. It reads shared rigid-lid F/G records written by `WVTransformStratifiedQG` and `WVTransformHydrostatic`, plus the grouped wave records of `WVTransformBoussinesq` (#302). The scientific decoder is independent of model execution. SQG and Hydrostatic model integration are described below.

The internal scientific contract is `wave-vortex-stratified-modal-record-v1`, separate from the extension source API and `wave-vortex-spectral-operators-v1` execution contract. Root `model_version`, `WVTransform` and `AnnotatedClass` identify the v4 transform schema. Numeric arrays use the existing MATLAB property names and units: `x/y/z`, `Lx/Ly/Lz`, `j`, `kl`, `k/l`, `N2`, `rho_nm0`, `dLnN2`, `P0/Q0`, `h_0`, `z_int`, `PF0inv/QG0inv` and `PF0/QG0`, with `rho0`, gravity and rotation parameters. Default stratified `WVTransform.writeToFile` saves and `WVModel` output files automatically include evaluated `N2`, `rho_nm0` and `k/l`, so the scientific reader accepts them without extra export arguments. These additive writer fields remain optional for MATLAB reconstruction: required-property lists, constructors and readers are unchanged, and older MATLAB files remain valid. Explicit property-only saves (`shouldAddRequiredProperties=false`) retain their existing behavior. To use an older file missing numeric inputs with this C++ reader, load it in MATLAB and save a new copy with the default writer. C++ neither evaluates serialized MATLAB functions nor solves a vertical eigenproblem. Backward compatibility is required for MATLAB; the development-only compatibility exception applies solely to C++.

### Inspect and read

`WVStratifiedModalReader::inspect(path, inspection)` validates every scientific payload and returns geometry, explicit retained mode keys and logical matrix-byte counts. It checks double types, named dimensions and matrix orientation, physical extents, periodic horizontal coordinates, increasing vertical coordinates with bottom/surface endpoints, stable stratification, integer retained vertical and horizontal keys, primary Hermitian membership, duplicate keys, the declared antialias mask, positive preconditioners and equivalent depths, quadrature normalization, rigid-lid G zeros and unit-maximum inverse-matrix preconditioning. NetCDF fill values and nonfinite payloads are rejected. The barotropic G row/column is intentionally zero and `h_0(1)=P0(1)=Q0(1)=1` are the MATLAB sentinel conventions. No generic invertibility requirement is imposed on truncated or zero-barotropic matrices.

Named NetCDF dimensions are authoritative for orientation: the C API sees `PF0inv(j,z)` and `PF0(z,j)`, corresponding to MATLAB column-major `[Nz,Nj]` and `[Nj,Nz]`. No guessed transpose is allowed. Vertical endpoints use `1e-10*Lz` tolerance because MATLAB's WKB quadrature inversion returns endpoints to interpolation accuracy; their actual values are retained. Quadrature integrates a constant to `Lz` within `1e-10` relative error. Fourier keys allow floating-point roundoff when converting stored radians/metre to signed integers; this does not merge nearby keys or matrices.

Inspection retains geometry and membership metadata, scaling as `O(Nx+Ny+Nz+Nj+Nkl)` for shared records, plus `O(Nj*Nkl+Nj*nK2unique)` wave depths/preconditioners for Boussinesq, and a duplicate-key set during validation. Matrix payloads are scanned through a fixed 4096-double buffer. All shared and, when applicable, grouped matrix payloads pass inspection before `read(path, shared_ptr<const WVStratifiedModalRecord>&)` allocates their owned copies. The same open file is used for validation and loading; callers must not modify a file during a read. No model state, FFT plan, derived product or output file is allocated or mutated by inspection. Failed inspection/read leaves the previous result unchanged.

### Immutable scientific ownership and prepared views

The record provides const numerical accessors and `matrixView()` for typed, borrowed preconditioned views with exact source/name/action/family identities. Borrowed views remain valid while the record is alive. `groups()` exposes explicit group identities and retained-column membership. Shared QG/Hydrostatic operators use group zero for all retained columns. Boussinesq records expose the exact persisted `K2unique/iK2unique` membership, including discontiguous columns, for their grouped wave matrices. No integer-radius grouping, approximate deduplication or per-mode matrix expansion occurs.

Each successful read creates a new immutable in-process source identity. It is an ownership generation, not a file-path or numerical hash, and is never persisted. Reuse the same record and prepared operators during time evolution. A new read or changed geometry creates new identities and requires new preparation; changing `t`, `t0` or coefficient arrays does not mutate an existing scientific record. This avoids accidental reuse after a file is overwritten with different modal data.

`prepareVertical()` produces #357 operators for unpreconditioned F/G reconstruction, projection and cross-family products. Given the saved diagonal preconditioners, the scientific matrices are

- `Finv = PF0inv * diag(P0)` and `F = diag(1/P0) * PF0`;
- `Ginv = QG0inv * diag(Q0)` and `G = diag(1/Q0) * QG0`;
- `GToF = F * Ginv` and `FToG = G * Finv`.

Layouts name `F-modal`, `G-modal`, `F-grid` or `G-grid` as appropriate and use the record's `modeSetIdentity()`. Raw preconditioned views instead use `PF0-modal` and `QG0-modal`. The caller selects split/interleaved storage, backend and overwrite/add behavior. Derived numerical matrices are prepared once and copied into the execution owner; the scientific record can then be released if no scientific views are needed. Execution workspaces have their own ownership and lifetimes. `horizontalSpecification()` preserves the exact retained `(k,l)` order and WVM forward-unit normalization for an explicit number of vertical planes.

`persistentBytes()` accounts for scientific objects, numerical capacities, group membership and identity capacity. Add prepared-operator and workspace accounting separately when those owners are retained. `inspection.scientificMatrixBytes` is the logical four-matrix payload size, not another allocation to add. String inline storage, allocator/control-block overhead and NetCDF/BLAS/FFTW internal allocations prevent these estimates from being a process-memory bound.

### Verification

`TestWVStratifiedModalRecord` tests malformed and incomplete records, reversed matrix dimensions, late-payload errors beyond the scan buffer, allocation-failure cleanup, immutable ownership, exact views/groups, changed scientific identities and zero prepared application allocations. `TestStratifiedModalRecord` in MATLAB generates authoritative files across profiles, nonuniform vertical grids, odd/even nonsquare horizontal grids, mask settings and retained basis sizes. It compares decoded values and all six prepared operations, independent horizontal FFT output, modal subspace identities, quadrature, and decoded vertical differentiation/integration formulas against MATLAB; it also checks MATLAB reloads for Stratified QG and Hydrostatic and scientific independence from time changes.

For the MATLAB parity test, optionally set `WV_MODAL_DUMP_EXECUTABLE` to a previously built `WVStratifiedModalDump`; otherwise the test builds its own temporary reader. The task verification ledger is `.github/planning/issue-295-verification.md`.

The modal record now implements the pure numerical `WVStratifiedModalSource` interface in `CompiledKernel`. Pass its shared owner directly to `WVTransformStratifiedQGKernel::create` with an FFT engine to use the #296 standalone kernel. The core's [Stratified QG documentation](../CompiledKernel/README.md#stratified-qg-numerical-kernel) describes fields, projections, derivatives, diagnostics and explicit beta behavior. The SQG runtime adapters added in #297 connect this kernel to `wave-vortex-run`, forcing, observers, output and restart. MATLAB APIs and the additive default-save behavior integrated with #295 are unchanged.

## Hydrostatic model runtime

The modal record supplies `WVTransformHydrostaticKernel` (#299) and full model execution through `WVHydrostaticForcingEngine` and `WVHydrostaticIntegrationSystem` (#300). The runner restores MATLAB-authored `WVTransformHydrostatic` models with all three reference-time amplitude families, the persisted F/G operators, numeric stratification and the opaque `N2Function` payload. Execution does not evaluate MATLAB functions or recompute vertical modes. Output preserves the payload so MATLAB can reopen it.

Hydrostatic forcing factories live in the existing extension catalog. They reuse the stable forcing objects and dispatch coarse numerical operations to the Hydrostatic kernel. Spatial closures include the variable-stratification background-density correction, bottom quadrature weights, and pressure response to prescribed topography. Adaptive damping uses the effective horizontal and retained vertical resolution after explicit antialiasing. Fixed-amplitude restoration applies to Ap, Am and A0.

The shared stratified field service provides standard fields, diagnostics, surface variables and fixed, moving and event-based linear/spline interpolation. Integrated observers include 3D tracers with full or horizontal-only advection, and particles with full or fixed-depth XY advection; moorings, Eulerian fields and tracked variables use the existing output providers. Explicit/CFL RK4 and adaptive RK3(2)/RK4(5)/RK7(8) use the existing integration and output machinery, including dense output, multi-file graphs, create/replace/append, restart and transactional preflight.

`TestPortableHydrostatic` compares isolated and composed forcing tendencies, continuation, fields, observers, output and restart against MATLAB. Set `WV_STABLE_FORCING_DUMP` and `WV_STABLE_FORCING_RUNNER` to existing `WVStableForcingDump` and `wave-vortex-run` executables; build `WVStratifiedQGFieldDump` in the same directory for the shared field probe. `WV_STABLE_FORCING_NATIVE=1` adds native FFTW parity. `TestWVHydrostaticRuntime` tests input/output alias rejection, invalid metadata, prepared allocation counts and setup failure cleanup. The focused Hydrostatic workflow runs reference parity and boundary contracts under release and ASan/UBSan builds.

The runtime owns prepared numerical and observer storage; scientific records can be shared across independent models. Physical-field reuse is limited to a single forcing evaluation, so calls after in-place state changes cannot reuse stale fields. The current implementation may reconstruct fields again for observers and projects spatial closures individually. Large-case timing, memory scaling, lifecycle evidence and performance refinement are the #301 handoff. This increment establishes runtime functionality and parity, without an end-to-end performance claim.

## Stratified QG model integration

`WVModel.writePortableRunRequest` accepts default v4 SQG model-output files. The C++ runtime restores only `A0[Nj,Nkl]`, retaining one immutable scientific owner shared by the kernel, model state and output adapters. Model graph inspection retains the scientific matrices; the standalone `WVStratifiedModalReader::inspect` API retains its bounded matrix scan. Restore revalidates the file before sharing the selected owner. No vertical eigenproblem or MATLAB function execution occurs in C++.

SQG supports nonlinear advection, adaptive damping, fixed/narrow-band amplitudes, explicit antialiasing, beta PV advection, linear/quadratic bottom friction and vertical diffusivity. The forcing catalog supplies exact identity/version, stage, priority and incompatibility evidence. Prepared operations use canonical interleaved coefficients and the existing full FFT gather/embed schedule; matrix schedule tuning remains separate.

The field service provides applicable geostrophic volume/surface fields and energy/speed diagnostics, fixed profiles, and three-dimensional fixed, moving and event samples. Periodic horizontal interpolation and nonuniform vertical linear/not-a-knot spline interpolation match MATLAB, including outside-domain zeros. Moorings, Eulerian fields, coefficients, XY particles at fixed depths and three-dimensional tracer grids integrate with the existing output graph. An SQG tracer must explicitly use `WVTracer(...,isXYOnly=true)`: MATLAB's default `addTracer` configuration attempts to read the unavailable `w` field. The runtime rejects vertical particle/tracer advection at preflight. MATLAB behavior is unchanged.

Explicit and CFL RK4, RK23, RK45 and RK78 use the common integrators and continuous-output machinery. Adaptive tolerances use energy per radial bin and vertical mode. Tracer derivatives use all full-grid horizontal Fourier modes before optional projection to the retained set; their even-grid Nyquist derivative is zero. This preserves tracer behavior when the model's retained set is truncated.

Create, replace, append and segmented restart preserve multi-file identities, independent group schedules, dense fields and integrated observer state. Fresh outputs preserve numeric scientific matrices and the opaque uint8 `N2Function` persistence payload so MATLAB can reopen them through `WVModel.modelFromFile`. The payload is copied without deserialization; it is not used as a numerical cache identity. Existing MATLAB validators and save formats remain unchanged; SQG request validation is additive.

`TestPortableStratifiedQG` covers MATLAB-authored forcing pairs/compositions, samples, integrated observers, output policies, restart and preflight rejection with reference and native FFT providers. C++ kernel tests independently check full-grid tracer derivatives, retained projection, invalid buffers and allocation-free prepared calls. See `.github/planning/issue-297-verification.md` for the verification ledger. The [SQG qualification guide](QUALIFICATION.md) describes the complete-model evidence, reproducible reference/native runs, numerical tolerances, lifecycle/storage checks and the fixed-RK4 endpoint distinction.

## Stable forcing and compatibility catalog

`contracts/portable-forcing-compatibility-v1.json` is the forcing slice of `portable-compatibility-matrix-v1`. Its 144 rows cover the twelve stable supplied MATLAB forcings, constant stratification in hydrostatic/nonhydrostatic mode, Barotropic QG, Stratified QG, Hydrostatic and Boussinesq, with transform antialiasing on/off. There are 124 supported pairs and 20 intentional MATLAB incompatibilities. Thermal damping remains explicitly excluded as under development. Later transform integration issues extend this same slice; #306 assembles the full feature catalog.

`WVAntialiasing` persists `Nj` and prepares its radial/vertical tendency mask once. Its spectral stage does not constrain amplitudes. Double antialiasing is rejected during model preflight. On constant-stratification transforms, the MATLAB-recognized `antialias filter` name also determines the effective resolution for adaptive damping, pseudo-topographic damping avoidance and CFL selection. Barotropic and Stratified QG retain MATLAB's geometry-based effective resolution.

Horizontal and vertical damping persist `nu` and `kappa`; vertical diffusivity persists `kappa_z` and `shouldForceMeanDensityAnomaly`. Their spatial forcing stage and priority are preserved, while the constant-stratification implementation evaluates the equivalent diagonal field operator and wave/vortex projection directly in coefficient space. Unequal viscosity/diffusivity retain wave/geostrophic coupling and reference-time phases. This adds no FFTs, physical workspaces or dense matrices. The mean-density-anomaly source is zero for constant stratification and inert for SQG. SQG vertical diffusivity preserves MATLAB's `-f*kappa_z*diffZG(eta,n=3)`, including the intermediate modal projection for variable N2 and truncated bases. Horizontal and vertical damping remain SQG incompatibilities; Barotropic QG rejects all three closures. Narrow-band forcing preserves its own class identity and persisted selected coefficients through the fixed-amplitude implementation.

Constant-stratification factories receive `WVForcingPreparation`, containing adaptive-damping availability and the resolved horizontal/vertical resolution. Registration-owned preflight and preparation hooks run before factory construction; the prepared facts are immutable during execution. Rebuild source-linked C++ extensions for this signature. No MATLAB constructor, required-property list or property-only save format changed.

Regenerate with `generatePortableForcingCompatibility` on the MATLAB tools path. MATLAB constructors and actual attachment determine applicability, stage and priority; unexpected errors stop generation. The generator checks the documented stable inventory against every file in `Forcing/`. The evidence supplement records implementation and test references. Factory availability, scientific applicability and qualified acceptance remain separate; pending rows cannot pass `validatePortableForcingCompatibility(matrix,requireComplete=true)`.

`TestPortableStableForcing` consumes the constant-stratification and Barotropic QG rows; `TestPortableStratifiedQG` consumes the SQG rows for isolated tendencies and MATLAB → C++ append → MATLAB restoration. It also covers odd/even grids, zero/nonzero parameters, both mean-density flags, fixed-amplitude ordering, converted resolution, segmented versus uninterrupted RK4, adaptive RK78/CFL, invalid append targets and prepared native allocation counts. Native FFTW and the direct reference provider share numerical tolerances; the reference provider intentionally allocates per-line scratch. CTest checks registry identity/version/stage/priority and all eighteen incompatible rows. Catalog tests reject missing, duplicate, contradictory or unproven rows and verify deterministic regeneration. Test references are enforced structurally; passing numerical tests remains necessary for qualification.

Build `WVStableForcingDump` and `wave-vortex-run` for the MATLAB suite. Optional environment variables `WV_STABLE_FORCING_DUMP`, `WV_STABLE_FORCING_RUNNER` and `WV_STABLE_FORCING_NATIVE=1` reuse a prepared build and enable both providers; otherwise the test builds a temporary reference executable. Verification evidence and the later-transform handoff are recorded in `.github/planning/issue-290-294-verification.md`.

### Boussinesq scientific records and standalone kernel

The modal decoder validates `PFpmInv`, `QGpmInv`, `PFpm`, `QGpm`, `QGwg`, `Ppm`, `Qpm`, `h_pm`, `K2unique`, and the one-based `iK2unique` membership from default MATLAB files. Wave matrices retain their declared group order and exact column membership; groups are not reconstructed by rounding a radius. Validation rejects missing/unfilled/nonfinite arrays, wrong named dimensions or orientation, invalid/unreferenced groups, inconsistent group depths, nonpositive preconditioners, rigid-lid violations and the singular `N2 == f^2` projection. Shared geostrophic matrices remain common to all columns. Preparation builds normalized grouped operations without changing the persisted values.

Pass this scientific owner to `WVTransformBoussinesqKernel` for the [standalone numerical surface](../CompiledKernel/README.md#boussinesq-numerical-kernel). Issue #302 adds kernel capability; #303 adds complete-model dispatch as described below, with end-to-end qualification in #304. MATLAB scientific methods and save defaults are unchanged.

### Boussinesq model integration

`wave-vortex-run` continues complete default MATLAB-authored `WVTransformBoussinesq` model files through `WVBoussinesqForcingEngine` and `WVBoussinesqIntegrationSystem`. The shared catalog has 23 supported Boussinesq forcing/configuration pairs and one intentional double-antialias incompatibility. Spatial forcing projects all four tendencies, including vertical momentum. Topographic responses use each persisted wave group's Fw boundary values; adaptive coefficient tolerances use wave h_pm energy and balanced h_0 energy. Reference/native FFTs use the scalar matrix backend in complete-model execution.

The model supports linear passive observers, full 3D and XY particles/tracers, fixed/moving/event interpolation, Eulerian fields, moorings, RK4 with explicit/CFL steps, adaptive RK3(2)/RK4(5)/RK7(8), dense output, multi-file schedules, segmentation and restart. Create/replace output writes the selected restart's exact shared and wave matrices, group membership, preconditioners, depths and opaque N2Function payload. Append preserves each destination's existing scientific source. Bundle compatibility compares the full wave payload; unsupported or inconsistent models fail preflight before destination mutation.

MATLAB's portable run-request helper accepts Boussinesq additively and includes wave arrays when checking source-file compatibility. Previously accepted calls and default model saves retain their behavior. `TestPortableBoussinesq` compares all applicable forcing tendencies/continuations, composed closures, field sampling, observer state, integrators, exact persisted arrays and transactional rejection against MATLAB. CTest `boussinesq-runtime` covers aliases, invalid state, allocation-failure cleanup, prepared allocation counts and source-backed wave structures. The focused Boussinesq workflow records release and sanitizer contract/lifecycle evidence. Longer continuation, MATLAB resumption, exact wave-source persistence and repeated ownership/storage qualification are documented in [QUALIFICATION.md](QUALIFICATION.md#boussinesq-qualification). Complete release scope remains an explicit qualification gate; contract-only reports cannot claim full readiness.
