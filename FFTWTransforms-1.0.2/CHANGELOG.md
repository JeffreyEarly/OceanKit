# Changelog

## Unreleased

- Allocate preserving c2r scratch lazily so construction, r2c, and
  destructive-only c2r retain no spectrum-sized scratch.
- Preserve fixed benchmark history and add authoring/exported-package release gates.

- Extracted the FFTW transform history from GLNumericalModelingKit into a
  dedicated OceanKit authoring repository.
- Added the FFTWTransforms 0.1.0 MPM contract, source-only licensing policy,
  portable CI, documentation source, and shared OceanKit release caller.
- Migrated the validated bundled-FFTW r2c/c2r and DCT-I/DST-I implementations
  into a self-contained production package root.
- Preserved zero-copy allocation, destructive-inverse, alignment, transactional
  build, static eligibility, and structured fallback contracts.
- Quarantined feasibility benchmarks and canonical artifacts under the
  authoring-only `tools/benchmarks` tree.
## [1.0.2] - 2026-08-09
- Allocate preserving c2r scratch lazily so destructive-only consumers retain no spectrum-sized scratch.

## [1.0.1] - 2026-08-08
- Fix local MEX builds from protected MPM installations.

## [1.0.0] - 2026-08-08
- Initial release with zero-copy r2c/c2r, DCT-I/DST-I, and validated MATLAB-bundled FFTW capability detection.
