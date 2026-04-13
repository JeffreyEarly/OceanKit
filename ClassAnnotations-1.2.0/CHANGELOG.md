# Version History

## [1.2.0] - 2026-04-10
- added `CAEnableDetailedDescriptionSidecars()` so documentation builders can opt into annotation markdown sidecars without imposing sidecar lookup overhead on ordinary runtime annotation construction
- changed `CAPropertyAnnotation` and `CAMethodAnnotation` so sidecar markdown is loaded only when documentation builders enable it explicitly

## [1.1.0] - 2026-04-09
- added a shared annotated post-load hook so subclasses can restore optional persisted state after generic constructor-based reconstruction
- automatically include required dimension annotations when selected numeric properties are written through annotated persistence
- added annotated read/write support for empty object arrays and strengthened round-trip coverage for nested numeric properties

## [1.0.1] - 2025-12-10
- Initial CI release
