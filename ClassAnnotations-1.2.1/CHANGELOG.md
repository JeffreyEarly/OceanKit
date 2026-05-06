# Version History

## [1.2.1] - 2026-05-06
- declared the required `NetCDF ^1.0.2` dependency so annotated function-handle persistence resolves the NetCDF serialization support it calls.

## [1.2.0] - 2026-04-10
- added `CAEnableDetailedDescriptionSidecars()` so documentation builders can opt into annotation markdown sidecars without imposing sidecar lookup overhead on ordinary runtime annotation construction
- changed `CAPropertyAnnotation` and `CAMethodAnnotation` so sidecar markdown is loaded only when documentation builders enable it explicitly
- added docs-only `className` and `sizeText` metadata on `CAObjectProperty` so object-property documentation can describe class and shape information when reflected MATLAB validation is absent or incomplete

## [1.1.0] - 2026-04-09
- added a shared annotated post-load hook so subclasses can restore optional persisted state after generic constructor-based reconstruction
- automatically include required dimension annotations when selected numeric properties are written through annotated persistence
- added annotated read/write support for empty object arrays and strengthened round-trip coverage for nested numeric properties

## [1.0.1] - 2025-12-10
- Initial CI release
