# Developer experiments

This folder contains manual investigations and scratch work for package authors. It is not an automated-test root and is not included on the WaveVortexModel runtime package path.

`LegacyTests` preserves historical scripts that were previously stored alongside automated tests. These scripts may depend on obsolete APIs, external data, personal paths, interactive plots, or fixed output filenames. Inspect a script and its prerequisites before running it. No compatibility or runnability promise applies to this material.

`LegacyImplementations` preserves superseded or disconnected source that is intentionally excluded from the runtime package path. Its contents are retained for author reference, not as callable package code.

Deterministic automated checks belong in `UnitTests`; repeatable performance measurements belong in `Benchmarks`.
