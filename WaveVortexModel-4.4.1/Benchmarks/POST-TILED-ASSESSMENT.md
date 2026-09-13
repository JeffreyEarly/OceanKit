# Post-tiled profile and workspace assessment

The frozen implementation in PR #490 (`78025904`, report/evidence head `6550c166`) passed local qualification with integration reductions of 12.22% on EddyTide, 5.15% on larger Hydrostatic, and 10.66% on Boussinesq. Required hosted CI remains the adoption gate. This follow-up profiles that exact executable and assesses the next increment; it changes no runtime or MATLAB source.

## Current costs

Each disposable continuation uses one 20-second macOS `sample` window at 1-ms intervals, after loading/preparation. The original fixtures and provider/binary hashes are verified before and after execution. Continuations are extended and primary output is moved to the endpoint on the copies; other observer schedules remain. All runs completed, reported the expected embedded revision, and recorded zero duplicate RHS/output evaluator executions. These are exploratory CPU-attribution samples, not new speedup measurements. Process snapshots retain background activity. No experiment was changed or competing compute scheduled.

| Summed non-waiting samples | EddyTide Hydrostatic | Boussinesq |
| --- | ---: | ---: |
| FFTW | 38.12% | 36.58% |
| Accelerate BLAS | 2.02% | 29.73% |
| Fused horizontal worker, excluding sampled FFTW callees | 38.32% | 27.29% |
| Field assembly/derivative arithmetic | 8.41% | 2.06% |
| Phase preparation | 3.94% | 1.12% |

The fused-worker category includes gather/scatter, normalization, plane arithmetic and control. It must not be interpreted entirely as movement or removable overhead. FFT/BLAS alone account for about 40% of EddyTide active samples and 66% of Boussinesq. Together with the fused numerical pipeline, they account for about 78% and 94%. Most remaining work is scientific arithmetic and its supporting movement, rather than duplicate evaluation. These categories are approximate; worker activity is not the same denominator as elapsed integration time.

The main thread provides complementary information. Its inclusive tiled-worker shares are 43.11% (EddyTide) and 33.13% (Boussinesq); vertical execution shares are 11.44% and 36.04%. Main-thread assembly/derivative self categories are 7.28% and 11.62%, projection arithmetic 8.44% and 6.14%, phase 2.85% and 6.82%, and state validation 5.16% and 2.75%. EddyTide's RK78 routine itself accounts for 8.79% of main samples; the broader integrator/state category is 13.70% and includes supporting control. Inclusive shares overlap other categories and must not be added.

The larger Hydrostatic composite has active tracers: 74.70% of main samples occur within scalar advection, including 41.05% within vertical calculus. That explains why a common fluid-pipeline improvement has less impact there. Tracer optimization remains outside the owner's requested scope. The complete profile is retained as context, not used to redirect this increment toward tracers.

## A small memory reduction is available

The safest first change stays inside the native provider. Existing `RetainedFFTWResources::tile` is unused throughout `advectionTask`, but already owns `16 W M` complex values. Four base tiles need only `4 W D M`, and the selected depth satisfies `D <= 4`. Reuse that tile for the four bases.

Each gathered target-z tile also becomes dead after that target's axis-2 derivative is consumed. Its projected flux can overwrite the same tile. Other lanes and targets occupy disjoint slices; final scatter still occurs after completing the tile. The shared-resource active guard excludes another transform for the whole synchronous worker dispatch.

Together these reduce additional packed storage from `(4 + 2 T) W D M` to `T W D M` complex values. With 16 bytes per complex value, the saving is `16 W (4 + T) D M` bytes. Here `T=3` for Hydrostatic and `T=4` for Boussinesq; `W=12`; EddyTide has `D=3, M=11439`, and both deeper fixtures have `D=4, M=6631`.

| Fixture | Provider-local saving | Remaining added owned peak, predicted |
| --- | ---: | ---: |
| EddyTide | 46,122,048 bytes | 50,339,192 bytes |
| Larger Hydrostatic | 35,648,256 bytes | 54,335,208 bytes |
| Boussinesq | 40,740,864 bytes | 49,968,504 bytes |

These are lifetime/capacity calculations, not measurements of an implemented candidate. They require no family API or public alias-contract change. Test alternating calls through plans sharing the resource, split/interleaved input, worker tails, density correction, exact fields/flux tolerances, producer counts, failure recovery and warmed allocation/storage accounting.

A second, wider change could lend each kernel's otherwise idle `real_` workspace to the provider's two real planes. The kernels own six physical volumes; the provider needs only `2 W Nx Ny` doubles, and `W <= Nz`. At these grids it would save a further 12,582,912 bytes. It requires a validated per-call borrowed span and should not be retained by preparation. Do not reinterpret real-vector storage as complex objects or overwrite field caches, phase values or spectra still needed by projection/diagnostics. Start with provider-local reuse rather than this wider ownership change.

## Next speed screen

The emitted fused accumulation loop is scalar, including a per-element branch for the invariant eta/z density-correction case. The disassembly at `advectionTask + 6376` lies in that loop. Sampling coalesces leaf offsets, so the loop's exact share cannot be isolated from these reports; the complete worker's 38%/27% active share is only a broad upper bound.

Screen a private ordinary/density-z loop specialization with local compiler alias qualifiers. Flux and derivative are disjoint halves of private provider scratch; velocity and eta are different field slices. Those ownership facts can be made visible to the compiler without introducing a public caller obligation. Start with compiler vectorization under Clang and GCC, not manual intrinsics or another general vDSP sweep.

Preserve x/y/z accumulation order, the addition of correction before velocity multiplication, the existing multiply-subtract contraction, and the ordinary `derivative + 0.0` form. Test signed zeros and nonzero density correction. Do not split the correction into separately accumulated products, reassociate axes, enable fast-math, or change FFT/producer counts. Confirm vectorization and a complete three-/four-target stage improvement before a short whole-model screen. Treat any speedup as unproven until then.

These two provider-local changes can form one bounded follow-up, with separate memory and speed evidence before one combined qualification. If the loop screen loses, retain only a correctness- and memory-qualified workspace reduction without claiming a speed gain. A serial RK78 combination-loop improvement is a subsequent EddyTide-specific candidate. Boussinesq assembly/projection batching is a larger refactor and ranks after this smaller common screen.

## History, evidence and verification

The spectral-kernel-benchmarks history and closed #24 conclusion were rechecked. Active reset and preserved-input alternatives lost, despite reducing nominal movement; generic split/vDSP and page-strided direct access had also failed to displace the accepted path. They are not being reopened. The new proposals exploit storage made idle by #490 and an observed scalar loop in its newly fused producer, rather than repeating those experiments.

One coordinator ran three sequential profiles; one existing reviewer independently audited lifetimes and the proposed loop contract. No rebuild, MATLAB rerun or broad numerical suite was needed because compiled inputs are unchanged. Review covered the shared-resource guard, worker completion, tile last reads, protected live caches, aliasing, correction order and sampling limitations. Repository/provenance, local links and whitespace checks apply to this report/evidence increment.

Compact categories, protocols, scripts and hashes are in [the evidence directory](../.github/ci-evidence/post-tiled-assessment/summary.json). Raw samples, disassembly, requests and outputs remain in `OceanKitRepositories/wave-vortex-model-benchmark-artifacts/tiled-advection-adoption-20260912/post-adoption-profile`. Sampling windows do not establish a hardware-counter utilization measurement or predict a whole-model speedup.
