# Legacy implementations

This folder preserves disconnected and superseded implementations for author reference. It is not on the WaveVortexModel runtime package path, and none of these files carries a compatibility, correctness, or runnability promise.

- `Transforms` contains the superseded constant-stratification transform implementation.
- `OffGrid` contains the removed off-grid and external-wave prototype.
- `Stratification` contains a disconnected free-surface projection experiment.
- `Forcing` contains Internal forcing experiments whose numerical behavior was never established as correct.

Do not add this folder recursively to a production MATLAB path. Revived work should be reviewed and reintroduced through tested runtime code rather than used directly from this archive.
