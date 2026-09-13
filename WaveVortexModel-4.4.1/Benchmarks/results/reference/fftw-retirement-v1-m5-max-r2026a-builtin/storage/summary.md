# Builtin transform storage benchmark

- Status: `complete`
- Source: `d88abddb6409b6b5f8f77e01dfe6fbd60233f749`
- MATLAB: `2026a`
- Architecture: `maca64`

| Case | Known persistent (MiB) | Maximum known live (MiB) | Opaque entries | Persistent RSS (MiB) | Peak RSS (MiB) |
|---|---:|---:|---:|---:|---:|
| constant-nonhydrostatic-256x256x65 | 65.260 | 130.260 | 1 | 1233.406 | 1233.938 |
| constant-hydrostatic-256x256x65 | 65.260 | 130.260 | 1 | 1177.703 | 1178.375 |
| constant-nonhydrostatic-512x512x129 | 517.037 | 1033.037 | 1 | 8456.453 | 8457.000 |
| constant-hydrostatic-512x512x129 | 517.037 | 1033.037 | 1 | 8017.078 | 8017.812 |
