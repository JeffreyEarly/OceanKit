# Portable transform qualification

The [forward-integration slice](INTEGRATION.md) records current method support across all six transform configurations and links representative cross-feature continuation evidence. Its execution receipt separates new qualification from source-bound inherited results. The transform-specific reports below retain their original methods, workloads, revisions, and tolerances; their historical case descriptions are not expanded when a later method is qualified.

## Historical evidence and current inventories

The SQG, Hydrostatic and Boussinesq validators default to `evidenceScope="current"`. This requires the complete current test inventory and the existing numerical, catalog, continuation and lifecycle checks. The result is `scope="current-inventory-only"` with `currentReadiness=false`: source and execution freshness require independent qualification and are not inferred from a reported commit label. A known historical execution is rejected in this mode, even if new test names are appended. New live qualification continues to use this strict default.

Use `validation = validatePortableHydrostaticQualification(report,evidenceScope="historical")` to inspect the retained Hydrostatic report; the SQG and Boussinesq validators use the same option. The result identifies `scope="historical-workload"` and `currentReadiness=false`. It preserves the source commit and binds the original report artifact and its independent test inventory. This mode accepts only the exact retained report or its existing reference-only/contracts-only projection. Each projection comes from the complete successful original artifact, so failed or missing original results cannot be hidden by choosing a subset. Original numerical results, methods, timings, tolerances and source identities remain unchanged.

The committed [historical inventory](qualification/historical-test-inventory-v1.json) derives literal test methods from the clean Git revisions recorded by those reports and includes the corresponding source-file hashes. Its digest is pinned in `portableQualificationTestInventory`, independent of incoming evidence. With full Git history available, `python3 tools/generate_historical_qualification_inventory.py` reproduces it; shallow CI needs only the committed authority and artifacts. Tampered inventories, changed artifact bytes and altered historical results fail validation. These historical checks do not establish current production-source compatibility or replace the current [forward-integration qualification](INTEGRATION.md).

## Boussinesq qualification

`qualifyPortableBoussinesq(outputPath,runner=runner,native=true)` records complete MATLAB-C++-MATLAB continuation after `configureCIEnvironment`. A supplied runner needs matching `WVStableForcingDump`, `WVStratifiedQGFieldDump`, `WVBoussinesqLifecycleProbe`, `WVStratifiedModalDump`, `WVBoussinesqKernelDump` and `WVHydrostaticKernelDump` probes beside it. With no runner it builds the reference configuration. Native qualification runs both reference and native FFTW. Complete model execution uses scalar matrix multiplication; standalone kernel tests separately cover scalar/Accelerate parity.

`contracts/boussinesq-qualification-cases-v1.json` declares six 600-second cases spanning exponential and thermocline stratification, nonuniform vertical grids, odd/even horizontal grids, modal truncation, explicit/transform antialiasing, linear/nonlinear dynamics and full/XY observer advection. Larger horizontal domains exercise the hydrostatic limit; smaller domains exercise vertical motion. Initial vertical kinetic fractions measure the distinction, with a ceiling of `1e-3` in the limit cases and a floor of `0.01` in the nonhydrostatic cases. Constrained initial states prevent initial filtering from satisfying the nonlinear evolution requirement of `1e-4` relative coefficient change.

Each provider runs uninterrupted C++ continuation and segmented C++ append/restart, compared with MATLAB throughout the same interval. A third path runs C++ for the first half and MATLAB for the second half, verifying executable restoration of the model graph. RK4 covers explicit and CFL-selected steps; RK3(2) and RK7(8) cover explicit and default adaptive controls. Output is every 100 seconds with dense fields every 25 seconds. Comparisons include each coefficient family, 3D fields and vorticity, energy/enstrophy, tracer and particle state, moorings, output schemas, schedules and restored observer identities. Wave matrices, group membership, preconditioners, depths and opaque N2Function bytes must remain exactly unchanged across each C++ append, before MATLAB's writable reload.

The lifecycle probe constructs, warms up, advances and destroys six complete models per fixture. Reference/native qualification measures 8×6×9, 12×10×13 and 20×16×17 grids; native FFTW additionally measures 64×48×49 with 24 modes. Reference-only hosted jobs use the two smaller grids to bound direct-DFT cost. Measurements cover sixteen prepared RK4 steps, RHS and forcing-service work counts, setup/runtime and retained capacities. Every cycle must release its scientific source and retain exactly the prepared storage; native execution must make zero prepared C++ application allocations. Output and restart costs are recorded separately by the continuation reports. Capacity accounting excludes allocator metadata and opaque provider/NetCDF internals; timings describe the recorded host and workload.

`wave-vortex-boussinesq-qualification-v1` links executed tests to all 23 supported Boussinesq forcing/configuration pairs plus the intentional double-antialias rejection. `validatePortableBoussinesqQualification` rejects missing, duplicate, stale, contradictory and out-of-tolerance evidence, including altered regime, restart, wave-payload and work measurements. `scope="contracts"` deliberately omits the six longer cases and emits `wave-vortex-boussinesq-contracts-v1`; it cannot claim complete continuation qualification. The focused workflow runs release and ASan/UBSan contracts and lifecycle checks independently of optional Full CI, uploading machine-readable reports. Manual dispatch with `complete=true` additionally runs complete release qualification. Linux instrumentation includes leak detection; local Apple Silicon uses ASan/UBSan without unsupported LeakSanitizer.

This qualification covers the supported Boussinesq transform surface. It does not add runtime eigenmode construction, new third-party forcing/observer types or runner selection of Accelerate. Broader feature-catalog and remaining non-transform readiness work remain separate. MATLAB scientific methods, save defaults, existing API behavior and package dependencies are unchanged.

The recorded [Apple Silicon Boussinesq report](qualification/boussinesq-apple-silicon-v1.json) contains 44 passing tests, 47 catalog/provider rows, twelve longer continuation results and seven lifecycle results. Maximum coefficient-family error is `1.58e-07`; maximum particle-position error is `3.37e-09 m`. Every lifecycle reports zero retained growth, zero prepared allocations and released ownership. The largest native fixture retains `202,313,652` bytes, with a median `8.713 s` per sixteen RK4 steps. Its short horizontal domain produces MATLAB projection-conditioning warnings; this case qualifies finite state and lifecycle/storage behavior, while the declared smaller fixtures establish numerical parity. The report records its tested source commit before the artifact itself is added.

Standard diagnostic catalog/evaluation and forcing-tendency diagnostics were completed in #314, #305 and #315. The completed compatibility matrix and [final source-bound parity decision](../.github/ci-evidence/issue-307-qualification/standard-parity-decision.json) are recorded through #306 and #307. This historical transform qualification supplies its family evidence; the final decision receipt combines it with subsequent changes and current-source qualification.

## Hydrostatic qualification

`WVTransformHydrostatic` supports MATLAB-authored model continuation with Ap/Am/A0 state, variable stratification, all twelve stable forcing identities where MATLAB permits them, full or horizontal-only tracers and particles, moorings, surface and volume fields, dense output, multiple output files, append and restart. MATLAB reopens the resulting model graphs directly. The scientific modal source and opaque `N2Function` payload are preserved. Create/replace writes the selected restart's immutable source to every new destination; append retains each existing file's payload. Independently serialized MATLAB function handles need not have identical opaque bytes. No MATLAB execution or eigensolver is needed during C++ continuation; MATLAB constructors, scientific methods and save defaults are unchanged.

After `tools/configureCIEnvironment`, run:

```matlab
qualifyPortableHydrostatic("hydrostatic-reference.json");
qualifyPortableHydrostatic("hydrostatic-native.json",runner="/path/to/build/wave-vortex-run",native=true);
```

The supplied runner needs matching `WVStableForcingDump`, `WVStratifiedQGFieldDump` (the shared stratified field probe), `WVHydrostaticLifecycleProbe`, `WVStratifiedQGLifecycleProbe`, `WVStratifiedModalDump` and `WVHydrostaticKernelDump` executables beside it. Without a supplied runner the tool builds the reference configuration. Native qualification includes the reference provider and native FFTW, plus scalar/Accelerate matrix parity in the numerical kernel tests. Complete model runs currently use the scalar matrix backend; this evidence does not claim Accelerate selection through the runner.

The `wave-vortex-hydrostatic-qualification-v1` report links executed tests to the existing forcing catalog, including all 23 supported Hydrostatic forcing/antialias pairs and the intentional antialias-on-antialias incompatibility. It records the source revision, working-tree status, MATLAB release, provider, timing, storage and measured continuation errors. `validatePortableHydrostaticQualification` rejects missing, duplicate, stale, contradictory or out-of-tolerance evidence. The focused Hydrostatic workflow runs reference release and ASan/UBSan contract/lifecycle coverage on pull requests and main pushes, uploading reports independently of optional Full CI. A manual dispatch with `complete=true` additionally runs full release continuation qualification for milestone acceptance. `qualifyPortableHydrostatic(...,scope="contracts")` runs the 39 reference tests covering kernel parity, every forcing pair, short all-integrator continuation, dense output, output policies/restart, malformed inputs and lifecycle ownership. Its `wave-vortex-hydrostatic-contracts-v1` report deliberately excludes the six longer trajectories; the validator rejects treating it as complete continuation qualification. Full readiness uses the complete release/native reports together with the instrumented contract report. Long trajectories run in complete scope, avoiding repeated direct-DFT timing under instrumentation.

### Numerical coverage

`contracts/hydrostatic-qualification-cases-v1.json` declares six 600-second trajectories using exponential and thermocline stratification, nonuniform vertical grids, odd/even horizontal grids, modal truncation, transform and explicit antialiasing, linear/nonlinear dynamics, and full/XY observer advection. RK4, RK3(2) and RK7(8) cover explicit, CFL-selected and default adaptive stepping. Each uninterrupted trajectory and two-segment append/restart trajectory is compared against MATLAB. Primary output is every 100 seconds and dense output every 25 seconds, including the initial state.

The coefficient comparison bounds each of Ap, Am and A0 separately. Further comparisons cover fields, energy, enstrophy, tracers, all three particle coordinates, tracked fields, moorings, dense fields, forcing identities, observer ordering, output variable schemas, schedules and ordinals. Explicit/CFL cases use a `2e-7` relative maximum-error ceiling; default adaptive cases use `5e-4`, below MATLAB's default `1e-3` relative tolerance. Particle positions have a `0.02 m` ceiling. Nonlinear cases must change the coefficient state by at least `1e-4` in relative maximum norm. The existing focused Hydrostatic tests additionally exercise non-prefix modal subsets, malformed metadata, exact forcing tendencies and continuations, and transactional multi-file create/replace/append behavior.

The CFL fixture requests the CFL number that selects a five-second step at each segment's initial state. This keeps MATLAB and C++ on a shared step grid while checking Hydrostatic advective/oscillatory step selection. It does not change the existing endpoint policies: C++ shortens a final step, while MATLAB can take a full step and interpolate. Fixtures avoid the existing MATLAB adaptive-damping degeneracy when explicit antialiasing leaves only one positive vertical mode.

### Runtime, memory and lifecycle scope

The shared stratified lifecycle probe constructs, warms up, advances and destroys six complete models per case. Both providers run 8×6×9, 12×10×13 and 20×16×17 grids; native FFTW additionally runs 64×48×49 with 24 modes. Large cases are native-only because the reference implementation uses a direct DFT. Reference-only runs use the 8×6×9 and 12×10×13 lifecycle fixtures and rely on the dedicated SQG workflow for its lifecycle regression. The declared `referenceOnlyLifecycleGrids` bounds the unoptimized direct-DFT workload on hosted runners; all six continuation fixtures remain required in complete scope, and both scopes retain kernel/forcing parity cases. Native-plus-reference runs retain all seven lifecycle/provider cases described above. Each lifecycle records setup time, sixteen measured RK4 steps, 64 RHS evaluations, reconstruction/projection counts and retained capacities. It verifies modal-source ownership release after destruction, unchanged prepared capacities and zero native C++ application allocations. These prepared-step measurements occur between output events; the continuation reports separately include real output/restart work.

Retained bytes count integration-system, integrator, state and output capacities. Scientific bytes are also reported for inspection. Allocator metadata and opaque FFTW/NetCDF internals are outside this accounting; process RSS in the runtime report is diagnostic, not an exact live-array ledger. Timing is descriptive, with no machine-dependent CI speed threshold.

The recorded [Apple Silicon Hydrostatic report](qualification/hydrostatic-apple-silicon-v1.json) contains 41 passing tests, 47 forcing/provider results, twelve continuation results and seven lifecycle results. The maximum coefficient-family discrepancy is `5.19e-9`; the maximum particle-position discrepancy is `2.62e-6 m`. All lifecycle cases report zero retained growth and zero prepared allocations. Representative native measurements are:

| Grid / modes | Retained capacities | Median time for 16 RK4 steps |
| --- | ---: | ---: |
| 8×6×9 / 4 | 372,498 bytes | 0.00271 s |
| 12×10×13 / 6 | 856,354 bytes | 0.0116 s |
| 20×16×17 / 8 | 2,543,002 bytes | 0.0486 s |
| 64×48×49 / 24 | 66,545,802 bytes | 2.91 s |

These samples ran on a shared development host; they establish descriptive baselines rather than isolated speed comparisons. The forcing-service counters record 128 physical-bundle reconstructions and 192 spatial projections per 64 RHS evaluations. They do not include every field reconstruction performed by the separate sampling/output services.

Hydrostatic spatial closures still project separately, and integrated observers may reconstruct fields after the forcing RHS. The qualification records this schedule before optimization. Backend selection and safe per-RHS field reuse remain possible measured follow-ups. The Boussinesq kernel (#302) and complete-model integration (#303) are implemented. Boussinesq runtime qualification is described above.

## Stratified QG qualification

`WVTransformStratifiedQG` supports complete MATLAB-authored portable model continuation with A0-only state, the nine applicable stable forcings, Eulerian and surface fields, three-dimensional fixed/moving/event sampling, fixed-depth XY drifters, rank-three XY tracers, moorings and linear passive observers. Model graphs retain multi-file/group identities, schedules, dense output, create/replace/append and restart. MATLAB reopens the output directly. See [README.md](README.md#stratified-qg-model-integration) for the execution surface and scientific persistence contract.

A three-dimensional SQG tracer must use `WVTracer(...,isXYOnly=true)`. This describes its horizontal advection on `[Nx,Ny,Nz]`; it does not flatten the tracer. MATLAB SQG has no `w`, so vertical tracer advection is an intentional preflight incompatibility. Hydrostatic models have their own qualification above. Existing MATLAB constructors, readers, required properties and save behavior remain unchanged.

## Run and inspect the evidence

From the authoring repository, configure the MATLAB dependencies with `tools/configureCIEnvironment`, then run:

```matlab
addpath("tools");
qualifyPortableStratifiedQG("sqg-reference.json");
```

The tool builds a temporary reference runner and probes, executes the qualification and existing regressions, writes JSON evidence, then validates it. A failed or incomplete run writes `status="incomplete"` and raises an error. Native qualification requires an explicitly built native runner and matching probes beside it:

```matlab
qualifyPortableStratifiedQG("sqg-native.json",runner="/path/to/build/wave-vortex-run",native=true);
```

Build targets are `wave-vortex-run`, `WVStableForcingDump`, `WVStratifiedQGFieldDump`, `WVStratifiedQGLifecycleProbe` and `WVStratifiedModalDump`. The native run executes both reference and native FFTW providers. No fallback is permitted. The authoring tool and tests are not required to execute an exported model bundle.

The JSON uses `wave-vortex-sqg-qualification-v1`. It records source revision and dirty status, MATLAB release/platform, the exact digest of the existing forcing catalog, executed test outcomes, catalog-row/provider results, continuation errors, selected stepping controls, runtime work, reported storage and repeated lifecycle measurements. Intentional forcing incompatibilities are provider-independent preflight results; they are recorded once under the reference run. `validatePortableStratifiedQGQualification` rejects missing, duplicate, stale-catalog, contradictory, failed or out-of-tolerance evidence. The evidence is a result of executing the existing compatibility slice, not another independently maintained support matrix.

The dedicated `Stratified QG qualification` workflow runs reference-provider release and ASan/UBSan qualification on Linux and uploads separate JSON artifacts. It runs independently of the optional Full suite. Native FFTW qualification is performed on the supported local Apple Silicon build; the report identifies the actual provider/library. An artifact demonstrates its recorded configuration, rather than certifying every compiler/platform or an exhaustive cross-product.

The recorded [Apple Silicon reference/native qualification](https://github.com/JeffreyEarly/wave-vortex-model/blob/main/PortableRuntime/qualification/stratified-qg-apple-silicon-v1.json) contains 34 passing tests, 41 forcing/provider results, 12 continuation cases and six lifecycle results. Its source revision is recorded before the artifact itself is added to the repository.

## Numerical coverage and tolerances

The shared forcing catalog contains 144 rows, including SQG's 17 supported and 7 intentionally incompatible rows across both transform-antialias configurations. Existing exact-pair tests compare every applicable RHS and append continuation, with odd/even grids and representative closure orderings. Their MATLAB-authoritative applicability and evidence links remain in the same catalog used by C++ contracts.

`contracts/stratified-qg-qualification-cases-v1.json` declares six longer complete-model cases. They span odd/even grids, one through six retained modes, explicit and transform antialiasing, fixed amplitudes, nonlinear closure compositions, linear passive observers, RK4/RK23/RK78, explicit steps, CFL-selected RK4 and default adaptive initial steps. Each run covers 600 seconds with primary output every 100 seconds and dense output every 25 seconds; its uninterrupted and two-segment trajectories are compared against MATLAB. Nonlinear cases must change A0 by at least `1e-4` in relative maximum norm, preventing nearly stationary trajectories from passing as substantive qualification.

Comparisons include coefficients, fields, energy, enstrophy, tracer state, particle positions and tracked fields, mooring records, dense fields, observer ordering, schedule times and committed ordinals. Explicit/fine-CFL trajectories use a `2e-7` relative maximum-error ceiling, while default adaptive segmented trajectories use `5e-4`, below MATLAB's default `1e-3` relative tolerance. Particle positions have a `0.02 m` ceiling. Reports retain the measured errors, not merely Boolean outcomes. Existing finer exact-pair tests continue to impose their tighter operation-level tolerances.

Fixed RK4 has an explicit endpoint distinction: C++ shortens the final step to land on the requested time; MATLAB's `WVArrayIntegrator` takes a full step and interpolates back. Segment boundaries can therefore change the discrete trajectory. The CFL case checks that reducing the step by ten decreases the endpoint coefficient difference by at least one hundred, in addition to meeting the fine-step parity ceiling. This qualification preserves both established behaviors and does not claim bitwise endpoint equality when the step grids differ.

## Runtime and storage

The lifecycle probe constructs, advances and destroys six complete models at each of three grid sizes. It checks that the scientific modal owner expires, prepared retained capacities remain constant over sixteen measured RK4 steps after warmup, and native prepared steps make zero application C++ allocations. A qualification-discovered fix reuses RK4's existing coefficient/block view vectors instead of allocating two metadata vectors per step; no state-sized buffer or numerical formula changed.

Runtime reports retain RHS counts, field reconstruction/projection counts, timing, integrator storage-ledger agreement and RSS diagnostics. SQG closures and distinct event occurrences currently reconstruct independently. Further reuse or backend tuning should follow measurements of the same workload. Retained-capacity accounting excludes allocator metadata and opaque NetCDF/FFTW internals; RSS is a process measurement rather than proof of exact live array storage. Timing is descriptive, with no machine-dependent CI speed threshold.

The authoring repository retains the [issue-298 verification ledger](https://github.com/JeffreyEarly/wave-vortex-model/blob/main/.github/planning/issue-298-verification.md). Issue #306 consumes this qualification alongside the same expanded forcing slice when assembling the full feature catalog.

## Catalog provenance after additive transform support

The recorded SQG/Hydrostatic measurements retain their original source revisions, test results and catalog SHA-256. Their `catalogPath` now points to an exact archived copy of the catalog used for those runs. Validation checks those bytes against the original digest and requires the transform's current configurations, complete rows, inventory and referenced evidence definitions to match exactly. Adding Boussinesq rows therefore preserves historical measurements without claiming a fresh long qualification run. Changes to a recorded transform's own compatibility semantics still invalidate its evidence. Current qualification tools continue to record the current complete catalog.
