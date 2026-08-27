# WaveVortex portable runtime

This directory contains the optional MATLAB-independent portable runtime. MATLAB remains WaveVortexModel's primary interface. The qualified `WVModel`, field, observer, NetCDF, and restart composition supports constant stratification and `WVTransformBarotropicQG` behind transform-neutral state and service contracts. Barotropic QG retains MATLAB's compact `A0`-only `kl` layout through integration, output, append, and restart; it never allocates compatibility `Ap` or `Am` state. The runtime is a source-only checkpoint-to-checkpoint tool for advanced users and the stable v1 source-level surface for statically linked C++ observing-system, output-schedule, and forcing implementations.

The runtime supports:

- fixed-step RK4, MATLAB `ode23`-compatible Bogacki--Shampine integration with a third-order accepted solution and second-order embedded estimate, MATLAB `ode45`-compatible Dormand--Prince integration with a fifth-order accepted solution and fourth-order embedded estimate, and MATLAB `ode78`-compatible Verner integration with an eighth-order accepted solution and seventh-order embedded estimate; all three adaptive methods include method-owned continuous output, with RK78 using Verner's seventh-order extension;
- the frozen v1 forcing subset (`WVNonlinearAdvection`, `WVAdaptiveDamping`, `WVFixedAmplitudeForcing`, `WVBottomFrictionLinear`, `WVBottomFrictionQuadratic`, `WVPseudoTopographicWaveGeneration`, and `WVBetaPlanePVAdvection`);
- the five qualified built-in observer records (`WVCoefficients`, `WVEulerianFields`, `WVMooring`, `WVLagrangianParticles`, and `WVTracer`) for constant stratification, and the MATLAB-valid coefficient, Eulerian-field, XY-particle, and rank-2 tracer forms for Barotropic QG;
- MATLAB-compatible checkpoint and time-series NetCDF data for the documented constant-stratification and Barotropic QG subsets.

Arbitrary MATLAB forcing or observing-system subclasses are not supported. Multi-file, named-group observing-system output is authored in MATLAB and executed by the command-line program without reproducing that scientific configuration in a second format.

`WVModel` is the thin, move-only runtime façade. It retains the frozen extension catalog together with the resolved forcing, observers, numerical system, integrator, output evaluation, driver, and sink; `WVModelState` separately owns canonical coefficients and dynamic observer state. The same façade backs the standalone runner and the production MEX right-hand-side path. It adds no numerical algorithm or state-sized copy.

Library clients may configure output with the stable source-level C++ `WVModelOutputFile` and `WVModelOutputGroup` builders. The builders emit canonical records and are discarded; `WVModelOutputConfiguration::compile` then applies the same record-preserving path used by `WVModel::createFromModelOutputFiles` for a restored MATLAB-authored sibling set. The resulting descriptor is the sole compiled graph shared by the plan, output evaluator, and NetCDF sink. Destination remapping changes only paths keyed by stable file identifier and preserves group membership, restart designation, complete observer configuration, declared schemas, and schedule identity, version, configuration, and typed cursor.

Output schedules are resolved by exact source-linked identity and version. `WVOutputPlan` retains one immutable schedule per group, while `WVOutputDriver` owns one small continuation cursor and one cached occurrence per group. Continuation cursors are distinct from sink-owned destination progress, which records committed record counts, time-last evidence, full typed cursors, and physical and committed ragged-axis offsets. A continuation advances only after its route is committed, so a failed coincident route can be replayed with the exact cursor without repeating successful routes. Delivery records expose only the latest selected event and retry attempts; cumulative history is aggregated in metrics, keeping storage independent of the event-window length. The legacy evenly spaced schedule retains its existing NetCDF scalars and ordinal semantics; new algorithmic providers persist a typed configuration and a bounded cursor envelope. `WVOutputSchedulePayloadSchema` resolves finite-real, integer, and Boolean occurrence fields once, and `WVOutputSchedulePayload` carries their values in fixed 4 KiB text-free storage addressed by numeric slot. State-triggered schedules are reserved but unsupported.

## Source API v1

`wave-vortex-portable-source-api-v1`, version 1.0, is a source-compatibility promise for the documented extension and embedding surface: source-API identity constants; explicit catalog construction and capability discovery; `WVObservingSystem`, `WVOutputSchedule`, and `WVForcing`; their construction/planning records and data-only `WVObservationSchema`/`WVObservationBatch` boundary; `WVModelOutputConfiguration`; `WVModel`; and `runWaveVortex()`. Legacy MATLAB-encoding adapters, declarations under `PortableRuntime/src`, `detail` namespaces, and undocumented implementation declarations are not extension APIs.

Consumers select one WaveVortexModel checkout and compile the runtime, application-owned catalog, statically linked extensions, and runner together. Rebuild the complete application whenever that checkout changes. An incompatible change to the documented v1 surface requires a new source-API major version; v1 promises no binary ABI across commits, compilers, or build configurations. It provides no dynamic discovery, separately loadable plug-in, or distributed runtime/extension binary.

`WVObserverFactoryRegistration` has exactly five public constructor inputs: type identity, contract version, factory, optional data-only configuration resolver, and optional data-only output-plan resolver. No legacy operation callback or persistence metadata is exposed. `WVExtensionCatalogBuilder` is the only mutable registration owner: register everything before freezing, treat a rejected registration as invalidating the builder, and freeze exactly once. Frozen catalogs are immutable shared owners retained by descriptors, output configurations, models, MEX handles, and runners. The builder and the caller's original catalog handle may then be destroyed, and independent catalogs may bind the same identities differently in one process.

Source API versioning is independent from exact data contracts. `wave-vortex-portable-pair-v1`, each paired type identity/version, the portable observer-graph schema, every schedule identity/version/configuration/payload schema/cursor, every observation-schema identity/version, `wave-vortex-run-request-v1` or `wave-vortex-run-request-v2`, and the compiled-kernel contract are validated separately. Matching one never relaxes another.

The portable reference runtime and source-linked consumers are qualified on Ubuntu with GCC or Clang and on macOS with AppleClang. The optimized native FFTW runner is Apple-silicon-only. Windows/MSVC source-linked builds are unsupported. The [AlongTrackSimulator ATS #4 integration](https://github.com/satmapkit/AlongTrackSimulator/commit/ba57981f336ad5bbbc0907dcd74fcd4fcd137708) is the external application-owned catalog and reusable-runner proof; final source API qualification recompiles that consumer against the explicitly selected WaveVortexModel checkout.

## Build

A portable reference build requires CMake 3.20, a C++17 compiler, and NetCDF C:

```sh
cmake -S PortableRuntime -B build/portable -DCMAKE_BUILD_TYPE=Release
cmake --build build/portable --parallel
ctest --test-dir build/portable --output-on-failure
```

On Apple silicon, the optimized runner can be built with:

```sh
PortableRuntime/buildWaveVortexRun.sh
```

The script verifies the pinned FFTW 3.3.11 archive, builds it in the ignored `.compiled-backend-cache`, and writes the executable to `.compiled-backend-cache/runtime-build/wave-vortex-run`. WaveVortexModel distributes no FFTW archive, library, MEX file, or executable. Redistributing a locally linked executable requires compliance with FFTW's GPL license.

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
    --fft-provider native-fftw --threads 18

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

`WVIntegrationSystem` supplies the selected right-hand side, constraints, error scaling, and optional field-evaluation service. `WVTimeIntegrator` advances accepted state, and `WVDenseOutput` evaluates within an accepted step. RK4, adaptive RK3(2), and output interpolation traverse the frozen family descriptions rather than assuming three equal rank-2 arrays. Output orchestration depends only on those contracts, so another integrator, coefficient-family layout, or state block does not require a method branch in the driver.

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
