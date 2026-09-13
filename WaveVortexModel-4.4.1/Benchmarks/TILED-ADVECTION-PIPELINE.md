# Tiled horizontal reconstruction and advection

This increment integrates the positive benchmark-only screen from `de3d6682` into v4 Hydrostatic and Boussinesq. [Issue #489](https://github.com/JeffreyEarly/wave-vortex-model/issues/489) tracks adoption. The accepted baseline is main `60287e19` (PR #478); the frozen implementation is `78025904a58e0195f28540f7bddf51e50b59525b`. MATLAB scientific behavior and checkpoint formats are unchanged.

## Pipeline and scope

The existing modal assembly and vertical matrix operators prepare u/v/w/eta spectra and each advected target's z-derivative spectrum. One native worker dispatch then reconstructs physical values, consumes x/y/z derivatives in the existing accumulation order, and projects the accumulated flux horizontally. The existing vertical projections and coefficient formulas finish the RHS. Each worker retains inverse-y columns for only one active plane; multiplication by the x wave number reuses those columns. Depth tiles are capped at four and at the worker's actual partition size.

This changes scheduling and memory traffic around the accepted FFTW transforms. It does not repeat the rejected full-depth inverse-y cache (#481), change the MM schedule, or introduce a separate scientific formula. The full physical bundle remains available for adaptive damping and other consumers. Actual FFT execution counts and reused columns are reported alongside evaluator producer metrics. For eligible ordinary nonlinear work, Hydrostatic uses 10 inverse-column batches per depth plane instead of 13, and Boussinesq uses 12 instead of 16; row inverses remain 13 and 16 respectively. The full-stage benchmark also benefits from fewer dispatches and immediate derivative consumption.

The compact native runner enables the path for ordinary Hydrostatic/Boussinesq nonlinear evaluation under `reuse`. Direct C++ construction exposes `tiledNonlinear`, default false. The forcing engine prepares the capability only when its initial workload declares a nonlinear consumer and its initial policy is `reuse`. All four physical evaluator nodes must be empty: `evaluateGroup` publishes them together after successful horizontal execution and vertical projection. Public forcing accumulation follows success. Already or partly available fields, diagnostic/capture calls, borrowed-field APIs, unsupported layouts/providers and low-memory evaluation use the established path.

Capability is fixed at construction. Switching an initially low-memory engine to reuse does not retroactively allocate tiled storage; switching an initially reuse engine to low-memory disables execution but retains its prepared capacity. There is no allocation-failure policy switch. Unsupported retained Nyquist layouts fall back; other preparation failures are reported.

## Correctness and ownership

Base spectra remain immutable. Each target's z spectrum shares its retained flux destination; a tile is gathered before its exact rows are overwritten. Worker depth partitions are disjoint. Hydrostatic requires one additional retained-spectrum slot; Boussinesq reuses existing slots. Native resource and workspace guards cover the complete dispatch. Flux destinations are rejected if they overlap any registered immutable view, including a view other than the one passed to the call.

Focused tests cover both families, split/interleaved layouts, worker partition/tail cases, density correction, unsupported providers, ready and partly ready physical nodes, repeated speed requests, low-memory fallback, same-time changed buffers, failure publication and restored-state recovery. Physical fields agree exactly with the established path. Flux, integration and MATLAB comparisons retain existing tolerances. Producer tests verify the actual reduced column schedule and no warmed allocations; derivative capture remains on its established storage path.

## Qualification

| Workload | Integration reduction | Paired bootstrap 95% interval | Process-lifetime reduction | Added owned peak |
| --- | ---: | ---: | ---: | ---: |
| EddyTide, 256 × 256 × 28 | 12.22% | 11.37–13.23% | 12.14% | 96.46 MB |
| Constant-stratification control | 1.36% | -0.06–3.02% | 0.04% | 0.37 MB |
| Larger Hydrostatic | 5.15% | 3.00–7.34% | 4.79% | 89.98 MB |
| Boussinesq, 256 × 256 × 129 | 10.66% | 10.24–11.05% | 0.71% | 90.71 MB |

All 32 measured pairs and eight warmup pairs passed scientific comparisons and identical integration/state decisions. The other actual producer counts are unchanged, and tiled executions/reused columns are positive for the two MM families. No integration or process-lifetime regression exceeds 3%. The unchanged constant control should be interpreted as measurement variation. Owned storage and process RSS are separate measurements; full values are retained in the evidence.

The formal campaign uses the exact preserved qualified PR #478 runner and a freshly linked candidate at `78025904`, with two warmup pairs and eight measured alternating pairs per fixture. Integration is measured separately from process lifetime, loading, preparation and output. Fixtures, provider libraries, binaries and source identities are frozen. The paired bootstrap 95% reduction intervals are positive for all three affected workloads. Boussinesq process lifetime improves only slightly because loading and preparation dominate its short continuation. The benchmark host is observed without changing any unrelated process or experiment. The preflight was idle; intermittent macOS and security/management service activity occurred during measurement. Ten-second process snapshots are retained, and cannot exclude shorter overlap with individual integration windows. Low-memory performance/storage qualification is outside this increment; focused low-memory correctness remains covered.

Verification ledger:

- Native Release combined suite: 58/58 passed; subsequent constructor workload gate received focused Hydrostatic/Boussinesq runtime reruns.
- Complete ASan/UBSan Debug suite: 58/58 passed, including the final gate and failure recovery tests.
- GCC 14: affected provider/operator and four family kernel/runtime tests passed.
- MATLAB: all 12 selected scientific parity/source-boundary checks passed; compact C++ dump probes exercise the tiled implementation.
- Forward/restart: all six family cases passed with both reference and native FFTW providers (12 retained fragments), including dense output and controlled stop/resume.
- Committed forward-receipt catalog: all eight tests passed, including source and artifact identity validation.
- Independent provider and lifecycle review: no remaining blocker after adding the complete immutable-view output guard.

The runner and stop probe each executed a provenance preflight and reported the exact implementation revision. Their hashes are matched against the binaries recorded in every forward receipt. Existing compact forward reports omit embedded source identity; the separate executed preflight supplies that check without expanding this increment into #488's provenance-tooling work. Later receipt and report commits do not change compiled inputs, so successful numerical and performance evidence remains applicable.

## Follow-up

Reprofile the qualified combined runtime before choosing another speed target. A separate storage-lifetime audit could reuse kernel workspace that remains idle during tiled execution and reduce the added scratch. That is a prospective memory improvement, not part of the measured result.

See [performance results](../.github/ci-evidence/tiled-advection-adoption/performance-summary.json), [producer and memory checks](../.github/ci-evidence/tiled-advection-adoption/producer-and-memory.json), [host observations](../.github/ci-evidence/tiled-advection-adoption/host-observations.json), and [archive hashes](../.github/ci-evidence/tiled-advection-adoption/archive.json). Raw artifacts are retained under `OceanKitRepositories/wave-vortex-model-benchmark-artifacts/tiled-advection-adoption-20260912`. Required hosted checks gate merging.
