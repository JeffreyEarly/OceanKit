# FFTWTransforms

FFTWTransforms is an OceanKit MATLAB package for reusable real-to-complex,
complex-to-real, DCT-I, and DST-I transforms backed by the FFTW library shipped
with MATLAB.

The package is source-only. It builds MEX gateways locally against the active
MATLAB installation and does not distribute FFTW libraries, downloaded FFTW
source, or precompiled MEX files.

## Current support

The initial validated backend is MATLAB R2026a on Apple silicon (`maca64`). The
package can be installed on MATLAB R2024b or later, but unsupported releases and
architectures report structured unavailability instead of attempting a build.

```matlab
capabilities = FFTWBackend.capabilities();
if capabilities.build.isRequired && capabilities.build.isPossible
    capabilities = FFTWBackend.build();
end
```

`FFTWBackend.capabilities()` never builds the backend. `FFTWBackend.build()`
uses the active MATLAB installation's bundled FFTW library and returns the same
structured capability record whether the build succeeds or is unavailable.

## Transform classes

`RealToComplexTransform` stores only the nonredundant half spectrum. The final
ordered transform dimension is compressed, so `dims=[2 1]` produces half-x
storage and `dims=[1 2]` produces half-y storage. FFTW inverse results are
unnormalized; multiply them by the transform's `scaleFactor`.

`RealToRealTransform` implements normalized DCT-I and DST-I transforms for real
or interleaved-complex arrays. DCT-I preserves the physical length. DST-I omits
the two physical endpoints from the coefficient array, ignores forward input
endpoints, and restores exact zero inverse endpoints.

Allocating calls return MATLAB-managed arrays. Caller-preallocated outputs must
be reassigned. A preserving complex-to-real inverse copies the half spectrum
once, while the destructive inverse avoids an explicit spectrum copy for a
uniquely owned input; both its spectrum and real output must be reassigned.

See the [transform contracts](Documentation/WebsiteDocumentation/transforms.md)
for examples and ownership details.

## Package status

FFTWTransforms is being prepared for its 1.0 release:

1. establish the authoring repository and MPM contract;
2. migrate and harden the production transform backend;
3. preserve benchmark history and validate an exported package; and
4. publish `FFTWTransforms 1.0.0` through OceanKit's shared release workflow.

Experimental gateways, benchmark harnesses, canonical results, exploratory
scripts, and historical `stash` content are quarantined under
`tools/benchmarks`. They remain available for authoring and history checks but
are not part of the runtime API or exported package.

## History

Development began in
[`GLNumericalModelingKit/Matlab/Spectral/FFTW`](https://github.com/JeffreyEarly/GLNumericalModelingKit/tree/3c5fce2b2df9c892418676ef75ffbc5752216e55/Matlab/Spectral/FFTW).
The subtree's filtered commit history forms the beginning of this repository.
See [Documentation/HISTORY.md](Documentation/HISTORY.md) for the pinned
source revision and commit mapping.

## Licensing

Original FFTWTransforms source is available under the MIT License. The bundled
`fftw3.h` file retains its own permissive header notice and is not relicensed by
the project. See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

FFTW itself is GPL-licensed. If you redistribute a locally linked MEX binary,
you are responsible for determining and satisfying the applicable FFTW license
obligations. This repository distributes neither those binaries nor FFTW source
archives or libraries.
