# Variable evaluation lifecycle

The C++ execution policy applies to both constant-stratification configurations, Hydrostatic, Boussinesq, Barotropic QG and Stratified QG. Each RHS evaluation and each output event owns an immutable-state scope. Nonlinear advection, forcing, tracers and particles share the RHS scope; coincident output destinations share the output scope. New stages, rejected-step retries, restored states and dense-output states receive new scopes.

The default `reuse` policy retains shared results for that scope. A dependency key identifies the variable or internal producer together with its component, derivative, density reference and applicable forcing stage or prefix. Grid quantities have no sampling-geometry key, so different destination geometries can share the same reconstruction. Single-use derivatives retain their streamed storage when the prepared consumer schedule proves that the storage can be reused.

`low-memory` uses the same producers and dependency keys. At a consumer boundary it can explicitly evict an unpinned intermediate; a later request can recompute it. This is a memory/runtime trade-off, and the execution report records both eviction and recomputation. Allocation failure returns an error; it does not change the selected policy.

## Selecting the policy

Run-request v2 accepts the optional execution property:

```json
"execution": {
  "variableEvaluationPolicy": "low-memory"
}
```

Omission selects `reuse`. Request v1 and checkpoint formats are unchanged. The policy controls C++ execution storage and is not MATLAB scientific state.

C++ clients can set `WVVariableKernelServices::variableEvaluationPolicy` during model construction, or call `WVModel::setVariableEvaluationPolicy` before execution. Standalone forcing engines and `WVFieldEvaluationService` provide the corresponding setter. Changing policy while an evaluation is active is rejected.

## Scoped field queries

A standalone `WVFieldEvaluationService::evaluate` call starts a fresh scope. To share work across several prepared field plans, keep one `WVFieldEvaluationSession` alive:

```cpp
WVFieldEvaluationSession session;
auto status = fields.beginEvaluationSession(state, session);
if (!status) return status;
status = fields.evaluate(firstPlan, state, firstOutputs, firstOutputCount);
if (!status) return status;
return fields.evaluate(secondPlan, state, secondOutputs, secondOutputCount);
```

The session closes on destruction, including an error return. The service and borrowed state must outlive it. The caller must not modify coefficients, time or additional state until the session closes. Output arrays remain caller-owned; their contents do not depend on the session remaining active.

An owner and a monotonically advancing evaluation generation establish cache identity. Pointer, time and shape checks reject foreign views inside that explicit scope. Equal times or reused buffer addresses never establish unchanged state across scopes. Component views are registered explicitly with the owning kernel; arbitrary same-time arrays cannot join its scope.

Kernel constraints and projections reject writes into active borrowed coefficient arrays. This protects mutations performed through the C++ API; it cannot intercept a caller writing directly into its own buffer. Nested evaluation scopes, dependency cycles and incompatible extents are errors. Failed producers publish no cache entry, and later requests can retry them.

Component coefficient views are explicitly unregistered before low-memory eviction. Scoped cleanup releases pins on both success and failure, including retries within the same session. Prepared complex buffers retain their capacity after logical eviction; subsequent requests recompute their values without repeating allocation. Event teardown closes kernel registrations before clearing the evaluation arena.

The constant-stratification `evaluateRightHandSideWithContext` API returns borrowed advection fields and therefore requires an explicit `beginStateEvaluation`/`endStateEvaluation` scope. The integration system owns this scope through the last tracer and particle consumer. An expired or foreign RHS context is rejected, and its field accessor returns an empty view.

## Metrics and qualification

The runner's `variableEvaluation` report records the requested/effective policy and separate RHS/output execution ledgers. Each ledger reports contexts, producer executions, cache hits, evictions, recomputations, duplicate executions, live bytes and high-water bytes. `outputArena` reports prepared and peak arena capacity. `kernelProducers` independently counts state validation, phase preparation, derived-tendency validation/reconstruction and field/component/derivative reconstruction at the numerical producers. Named reduction counters include horizontal speed, vertical speed and energy loops, including work initiated by public field queries. Forcing metrics also count nonlinear production. Grid-calculus counters identify the physical field (`u`, `v`, `w`, `eta`) and the numerical formulation (`Dxx`, `Dyy`, F-grid `Dzz`, G-grid `Dzz`); constant-stratification Laplacian counters distinguish horizontal and vertical transforms. Reconstruction field indices follow the selected kernel's field enum. These independent counters detect work that bypasses the cache ledger.

Ledger bytes describe live cached values. The existing `storageBytes` and `livenessBytes` reports include allocated adapter buffers, plan metadata and model storage; process peak RSS additionally captures allocator and library storage. These quantities have different scopes and should not be added together. `diagnosticEvaluation.additionalTransientHighWaterBytes` reports allocations beyond retained service storage; the total peak-storage estimate adds this value once, without adding logical diagnostic/event bytes already backed by the arena. Prepared workspace capacity can remain allocated between events while all cached values and identities are invalidated. Moving formerly transient allocations into preparation changes retained capacity without necessarily changing peak storage; the low-memory 3% qualification therefore compares maximum live owned storage and reports retained capacity separately.

[The qualification driver](../Benchmarks/variable-evaluation/run_qualification.py) runs frozen baseline/reuse/low-memory executables against copied model fixtures. It preserves the source files, compares complete NetCDF outputs at the manifest's scientific tolerances and checks integration decisions exactly, and reports paired timings and retained/peak memory. `--smoke` runs one correctness comparison and explicitly cannot qualify performance. The default comparison is bitwise exact. A manifest may select existing scientific tolerances when baseline repeatability evidence shows native backend variability; that evidence and the observed differences must accompany the result. Final qualification uses two warmup and eight measured pairs on an idle host, including the EddyTide `256 × 256 × 28` nonlinear-plus-damping case and representative larger composite fixtures. The active development results and outstanding gates are recorded in [the verification ledger](../.github/planning/issue-470-verification.md).
