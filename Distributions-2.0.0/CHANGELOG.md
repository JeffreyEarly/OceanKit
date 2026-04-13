# Version History

## [2.0.0] - 2026-04-09
- standardized the distribution constructors around required-property name-value arguments, including `AddedDistribution(scalings=..., distributions=...)`
- switched known-type NetCDF reconstruction to the inherited `annotatedClassFromFile(...)` and `annotatedClassFromGroup(...)` path
- kept `Distribution.distributionFromFile(...)` as the polymorphic read entry point while moving optional `rho` restoration onto the shared annotated post-load hook
- preserved optional nonempty `rho` when writing distributions into existing NetCDF groups through `writeToGroup(...)`
- restored the documented no-input defaults for `NormalDistribution()`, `RayleighDistribution()`, and `UniformDistribution()`
- split smoke scripts out of `UnitTests` into a shipped `Examples` folder and replaced the script-based checks with `UnitTests/DistributionsUnitTests.m`
- updated README and website usage examples to the new constructor and persistence APIs

## [1.1.0] - 2026-03-27
- added argument validation across constructors and public API boundaries without intended functional changes
- replaced the external `RunningFilter` dependency in `TwoDimDistanceDistribution` with MATLAB `movmean`
- refactored the package into `@ClassName` folders and tightened the public/internal class surface without intended behavior changes
- added inline API documentation for the distribution classes, methods, and properties
- added the standard OceanKit website scaffolding and generated class reference pages

## [1.0.0] - 2025-11-26
- initial release and MPM package setup
