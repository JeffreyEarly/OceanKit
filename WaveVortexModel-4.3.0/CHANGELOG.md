# Version History

## [Unreleased]

## [4.3.0] - 2026-08-27

### Compiled execution and portable runtime

- Added an explicit, source-only compiled preview for ordinary constant-stratification nonlinear flux. MATLAB remains the default; native support must be built explicitly, and unsupported transforms or forcing are rejected without fallback.
- Added an optional standalone constant-stratification runtime using the same MATLAB-independent numerical core, with fixed RK4, MATLAB `ode23`-compatible adaptive RK3(2), MATLAB-compatible NetCDF restart, and a command-line runner for versioned MATLAB-authored run requests.
- Added a shared allocation-free adaptive Runge--Kutta driver and MATLAB `ode45`-compatible Dormand--Prince integration, preserving the existing `ode23` controller and parity fixtures while reporting exact method-specific workspace, dense-history, error-policy, diagnostic, and stage-buffer liveness records.
- Added endpoint-only MATLAB `ode78`-compatible Verner 8(7) integration as `adaptive-rk78`, with MATLAB controller/work-count fixtures, explicit 11-buffer stage liveness, exact memory diagnostics, restart reconstruction, and no allocation or evaluation of the continuous-extension stages reserved for #284.
- Added MATLAB `ode78`-compatible seventh-order continuous output with four extension RHS stages allocated and evaluated lazily once per accepted step, exact endpoint recovery, cache/work metrics, transactional output retry, interpolated restart, and output-schedule-independent accepted trajectories. Endpoint-only RK78 execution retains its original 11-buffer workspace and performs no extension work.
- Added exact run-request v2 integration contracts for explicit or MATLAB-CFL-selected fixed RK4 and MATLAB `ode23`, `ode45`, and `ode78`, with typed version-first decoding, unchanged v1 behavior, pre-allocation validation, transform-owned CFL candidates, and lifecycle, selection, no-fallback, work, and storage reporting.
- Made a minimal run-request v2 reproduce standard `WVModel` execution defaults: `ode78`, relative tolerance `1e-3`, absolute-tolerance scale `1e-6`, post-restoration CFL `0.5` initial-step selection, a maximum step of one tenth of the continuation interval, and native FFTW with hardware-bounded threads. Reports now separate requested and active integrator/provider settings; explicit options override defaults independently, native unavailability fails actionably before state-sized allocation or output mutation, and v1 plus the positional CLI remain unchanged.
- Generalized the portable integration state boundary to transform-selected spatial and coefficient-family ranks, preserving the constant-stratification `Ap`/`Am`/`A0` source and NetCDF contracts while allowing one-family systems to integrate and checkpoint without dummy state-sized arrays.
- Added a MATLAB-independent Barotropic QG numerical kernel and focused integration system over that transform-neutral boundary, with compact QGPV-normalized `A0` storage, odd/even and nonsquare grids, `j=0/1`, optional transform-level antialiasing, MATLAB-matched fields, derivatives, invariants, nonlinear PV advection, adaptive tolerance vectors, and endpoint-only RK78 integration, plus portable reference and native Apple-silicon FFTW providers. End-to-end `WVModel`, observing-system, output, and restart composition remains separate work.
- Routed every stable Barotropic QG forcing through a construction-resolved QG forcing engine, including adaptive damping, fixed and narrow-band amplitudes, linear and quadratic bottom friction, and beta-plane PV advection, while preserving compact `A0`, accepted-state constraints, generic integrators, shared 4H+5R scratch, and MATLAB numerical parity. Pseudo-topographic wave generation and later closures remain explicitly unsupported for this transform.
- Added MATLAB-compatible portable observing and output support for coefficients, Eulerian fields, moorings, Lagrangian particles, and tracers, including shared field evaluation, continuous and lazy scheduled output, multi-file and multi-group transactional NetCDF persistence, restart, continuation, and append.
- Stabilized the C++17 source-linked extension API around an application-owned catalog, immutable resolved observer and forcing implementations, event-dependent observation geometry, and fixed, variable-length, and ragged output batches. The AlongTrack and linear-bottom-friction extensions prove the workflow; legacy observer adapters and direct `WVModel` internal-service access are no longer public, and no binary-plugin ABI is provided.

### MATLAB model and forcing

- Added persistent, resolution-convertible `WVNarrowBandGeostrophicForcing` and `WVBottomFrictionLinear`, with matching portable implementations. The former fixed-amplitude narrow-band helper remains as a silent deprecated 4.x delegate.
- Verified and documented the mathematical contracts for pseudo-topographic wave generation, beta-plane QGPV advection, and vertical diffusivity, and restored the variable-stratification mean-density-anomaly source disabled by an obsolete class-name check.

### Spectral internals and performance

- Retired the experimental fine-grained WaveVortex FFTW backend and removed the FFTWTransforms dependency; MATLAB's builtin transforms remain the supported default. This is separate from the explicit source-only compiled preview.
- Consolidated Fourier storage around compact `WVFourierStorageLayout` mappings, updated Fourier-at-position reconstruction, removed obsolete expanded-mapping and pseudo-radial/frequency-binning APIs, and added backend-neutral retained-storage and fresh-process RSS measurement.

### Documentation, validation, and benchmarks

- Reorganized transform and forcing documentation, clarified spectral grids, energy, flow components, geometry, density profiles, and projection and reconstruction, improved nonlinear-integration examples, normalized physical-unit metadata, corrected squared-wavenumber and coefficient-tendency dimensions, and adopted ClassDocumentation 1.3.2 through the immutable OceanKit release workflow v0.1.2.
- Reorganized the compiled MATLAB preview and standalone runtime under a concise top-level Compiled execution guide, preserved their public URLs, added an explicit compatibility matrix and minimal MATLAB-to-NetCDF/request-to-runner workflow, and linked detailed architecture to the developer contracts.
- Added a catalog-driven Benchmarks page with deterministic runtime and process-memory charts, accessible tables, platform and MATLAB-toolchain comparisons, downloadable results, and measurements from three physical machines. Empty portable typed-record encoding and decoding were hardened and are covered directly by Linux ASan/UBSan CI.
- Froze the v4.3 matched-interface presentation to issue #312's accepted Donut `[256 256 129]` integration-only runtime and total-peak-process-memory record; final release qualification does not rerun that performance matrix or publish startup as a primary metric.

## [4.2.1] - 2026-08-09

- Raised the minimum MATLAB release to R2025b, aligning the supported runtime, native MPM verification, and release engine.
- Pinned releases to the immutable OceanKit release pilot, added native-MPM clean-install and exported-package gates, and removed authoring tests from the installed runtime path.
- Replaced vertically replicated horizontal Fourier mappings with compact two-dimensional row mappings while preserving the canonical WaveVortex coefficients. On the Apple M5 Max/R2026a release-candidate run, complete builtin forward transforms were 1.05x–2.00x faster and inverse transforms were 1.81x–3.52x faster across the `[256 256 65]` and `[512 512 129]` gate cases with antialiasing on and off; numerical results remained equivalent within the benchmark tolerance. See the `transform-layout-v4.2.1-release-m5-max-r2026a-builtin` artifact for the machine-specific measurements.
- Rewrote the README, homepage, installation instructions, and user guidance around the current `WVTransform` and `WVModel` workflows; completed the NetCDF conventions guide; linked the canonical wavevortexmodel.org site prominently; and adopted the Avenir-first website typography while retaining search.
- Reorganized the generated API reference around user tasks, moved implementation machinery into Developer Topics, replaced release-status jargon with direct descriptions of capabilities and limitations, and restored authored scientific context, coefficient-occupancy tables, and design rationale.
- Corrected the `WVTransform` and `WVModel` API reference, promoted the stored and time-evaluated wave-vortex coefficients, and classified low-level projection and reconstruction arrays as Developer reference.
- Made documentation generation clean, deterministic, transactional, and reviewable with an exact ClassDocumentation 1.3.0 authoring dependency, canonical build/check tasks, generated hierarchy validation, and case-sensitive internal-route checks.
- Quarantined disconnected legacy implementations, corrected low-level adaptive-integrator and geometry spelling, retained a narrow 4.x file fallback for `shouldExludeConjugates`, and repaired the retained barotropic FINUFFT development path.
- Stabilized model-output scheduling and NetCDF restart across multiple files and groups; preserved linear dynamics and shared observing systems; and made initialization, writing, restoration, and handle ownership exception-safe.
- Made operation registration, replacement, and removal atomic and identity-based; corrected multiple-output lookup; and made cache invalidation follow each operation's declared time and coefficient dependencies.
- Made forcing registration identity-based, deterministic, and atomic; preserved vertical-diffusivity, fixed-amplitude, and explicit-antialias configurations through resolution changes and NetCDF restoration; and rejected vertical diffusivity on barotropic transforms.
- Corrected observing-system ownership and composition, linear-time flux evaluation, particle tolerance and tracked-field propagation, and periodic one-based mooring indexing.
- Added compact invariant coverage across every transform family on even and odd grids; corrected parity-aware Nyquist and Hermitian bookkeeping, vector mode/index mappings, coefficient-preserving resolution conversion, and barotropic spatial energy and enstrophy diagnostics.
- Generalized `summarizeDegreesOfFreedom` across every transform family with deterministic grid metadata and primary-component mask counts.
- Implemented deterministic total-flow analytical-solution lookup across primary components and repaired barotropic mode-index conversion used by that lookup.
- Implemented `hasMeanPressureDifference` as an MDA-only boundary diagnostic with a `1e-5` relative pressure tolerance; other flow components and transforms without an MDA component return `false`.
- Completed the vertical-calculus contract: first through fourth derivatives now use ordinary MATLAB order and grid-layout validation, and all three-dimensional transforms provide first antiderivatives with documented G-space projection and bottom-zero conventions.
- Restricted public field and particle interpolation to periodic `linear` and `spline` methods and removed the dead `exact` and off-grid dispatch.
- Documented the v4.2.1 capabilities and limitations, including the five transform families, latitude domain of `5 <= abs(latitude) <= 85`, MATLAB builtin FFTs, fixed and adaptive integration, and periodic `linear` and `spline` interpolation.
- Enforced the supported latitude domain and restored a deterministic, order-independent green test baseline with isolated random state and temporary files.
- Corrected combined random-flow initialization to preserve conjugate inertial `Ap` and `Am` coefficients.

## [4.2.0] - 2026-07-30
- Added `WVPseudoTopographicWaveGeneration` with Darwin-symbol or custom-frequency initialization, adaptive-damping-aware spectral masking, resolution conversion, restart persistence, and deterministic Goff abyssal-hill generation.
- Made transform restoration read-only by default, eliminated hidden duplicate NetCDF handles, and made one-output restoration close its file before returning.
- Changed model restart to complete read-only restoration before acquiring its single writable output handle, and made model output closure explicitly release that handle.
- Fixed horizontal and vertical Laplacian damping to preserve configured viscosity and diffusivity across transform resolution changes, and made vertical damping support constant-stratification transforms by using the exact zero stratification-gradient correction.

## [4.1.1] - 2026-07-26
- Added `waveModeVerticalStructureAtIndex` with optimized wave F and vertical-derivative G endpoint factors for constant-stratification, hydrostatic, and Boussinesq geometries; fixed constant-stratification resolution and explicit-antialias conversions to preserve `N0` and `isHydrostatic`.

## [4.1.0] - 2026-07-24
- registered flow components now automatically expose supported standard and sea-surface variables.
- added flow-component diagnostics for `pi`, `ssh`, `ssu`, and `ssv`, and corrected component pressure to use the selected component rather than total pressure height.
- added periodic two-dimensional interpolation support for surface fields while preserving three-dimensional interpolation behavior.
- added focused tests for composed-component closure and linear two- and three-dimensional interpolation across periodic boundaries.

## [4.0.7] - 2026-05-06
- raised the `InternalModes` dependency floor to `1.3.0` and the `SplineCore` dependency floor to `^2.2.0`.
- added `chebfun` as a direct dependency for the model's direct chebfun helper usage.
- migrated the `eta_true` spline fit to the `SplineCore` 2.x `BSpline` constructor, knot, and basis-matrix APIs.
- modernized WaveVortexModel's `InternalModes` call sites to the lowerCamel mode and quadrature APIs, including the new `modesAtQuadraturePoints` helper.

## [4.0.6] - 2026-05-06
- raised the `ClassAnnotations` and `NetCDF` dependency floors so OceanKit installs resolve the NetCDF function-handle serialization support required by annotated persistence.
- fixed `WVAdaptiveDamping` for `WVTransformConstantStratification` by routing `effectiveJMax` through the constant-stratification geometry superclass.
- fixed constant-stratification `rhoFunction` and `N2Function` handles so they preserve the shape of caller-provided `z` arrays.
- added `zeta_x` and `zeta_y` as known variables for `WVTransformConstantStratification`, matching the hydrostatic and Boussinesq transforms.
- restored `fluxForForcing` on `WVTransformConstantStratification` so diagnostics can compute forcing flux summaries.
- added constant-stratification `F_g` and `G_g` spectrum/cross-spectrum helpers needed by diagnostics.

## [4.0.5] - 2026-05-06
- changed `shouldUseTrueNoMotionProfile` to be a post-initialization setting that is not accepted by constructors or persisted through NetCDF round trips, while preserving transform-copy behavior and invalidating only the `rho_nm` cache when toggled.

## [4.0.4] - 2026-04-22
- added `shouldUseTrueNoMotionProfile` to the constant-stratification, hydrostatic, and Boussinesq transforms so `eta_true` uses `rho_nm0` by default and `rho_nm` when explicitly requested.
- added an advisory `EtaTrueOperation:OptimizationToolboxUnavailable` warning when `shouldUseTrueNoMotionProfile=true` and `rho_nm` is being computed without Optimization Toolbox support.
- added focused tests covering `eta_true` profile selection, warning behavior, transform-copy preservation, and NetCDF round-trip persistence of `shouldUseTrueNoMotionProfile`.
- changed `WVTransform.version` to read from `resources/mpackage.json` instead of a hard-coded version string so runtime metadata stays aligned with the package manifest.

## [4.0.3] - 2026-04-22
- pinned `InternalModes` to `1.0.1` so MPM resolves a `SplineCore`-compatible dependency set for `WaveVortexModel` installs.

## [4.0.2] - 2026-01-16
- new moment based algorithm for computing rho_nm.
- eta_true now computed from high order spline representation of rho_nm; recores SplineCore package.
- added placeParticlesOnIsopycnals function to all WVStratification subclasses.

## [4.0.1] - 2025-12-11
- initial mpm ci release

## [4.0.0](https://github.com/Energy-Pathways-Group/GLOceanKit/releases/tag/v4.0) - 2025-08-05
- Added support for `WVForcing`, a mechanism for adding arbitrary forcing to the model which also also the nonlinear fluxes to be automatically diagnosed.
- Added `WVObservingSystems`, a mechanism for adding user defined observing systems to the model, such as drifters, mooring, along-track altimetry.

## [3.0.0](https://github.com/Energy-Pathways-Group/GLOceanKit/releases/tag/v3.0.0) - 2024-07-25
- Full implemented non-hydrostatic forwards model, passing all unit tests.

## [2.0.0](https://github.com/Energy-Pathways-Group/GLOceanKit/commit/fffd3e6d7822c3844349ec6d65bd262f36f93122) - 2021-11-05
- Fully non-linear forwards model for hydrostatics.
- Never formally declared version 2.0, but on this date we successfully modeled the Cyprus eddy example from the JFM directly in the WaveVortexModel.

## [1.0.0](https://github.com/Energy-Pathways-Group/GLOceanKit/releases/tag/v1.0.0) - 2020-09-18
- All units test passed and declared version 1.0 for the JFM paper.

## [0.9.0](https://github.com/Energy-Pathways-Group/GLOceanKit/commits/master/Matlab/InternalWaveModel/InternalWaveModel.m) - 2018-01-04
- Model expanded to project any state of the fluid onto the complete wave-vortex basis, and time step the state forward linearly. The code was moved to Matlab.

## [0.0.1](https://github.com/Energy-Pathways-Group/GLOceanKit/commits/master/GLOceanKit/GLInternalWaveInitialization.m) - 2014-12-10
- Initially released as a linear internal wave model that time stepped forward non-hydrostatic internal gravity waves in variable stratification using internal wave modes. The code was written in C and Objective-C.
