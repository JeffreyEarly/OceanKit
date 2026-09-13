# WaveVortex benchmark

- Status: `complete`
- Run: `20260808T180426Z`
- MATLAB: `2026a`
- Architecture: `maca64`

## Suite scores

| Suite | Backend | Score |
|---|---|---:|
| core-v1 | builtin | 100.000 |

## Family scores

| Suite | Family | Backend | Score |
|---|---|---|---:|
| core-v1 | constant-nonhydrostatic | builtin | 100.000 |
| core-v1 | constant-hydrostatic | builtin | 100.000 |

## Timing and scores

| Suite | Case | Transform | Backend | Median (ms) | Reference score | Same-host speedup | Error |
|---|---|---|---|---:|---:|---:|---:|
| core-v1 | constant-nonhydrostatic-256x256x65 | constant-nonhydrostatic | builtin | 130.103 | 100.000 | 1.000 | 0 |
| core-v1 | constant-hydrostatic-256x256x65 | constant-hydrostatic | builtin | 104.823 | 100.000 | 1.000 | 0 |
| core-v1 | constant-nonhydrostatic-512x512x129 | constant-nonhydrostatic | builtin | 1156.642 | 100.000 | 1.000 | 0 |
| core-v1 | constant-hydrostatic-512x512x129 | constant-hydrostatic | builtin | 997.667 | 100.000 | 1.000 | 0 |

## Construction and cache diagnostics

| Suite | Case | Backend | Construction (s) | First call (ms) | Same-state cache hit (ms) |
|---|---|---|---:|---:|---:|
| core-v1 | constant-nonhydrostatic-256x256x65 | builtin | 1.310 | 207.574 | 94.857 |
| core-v1 | constant-hydrostatic-256x256x65 | builtin | 0.505 | 137.880 | 72.191 |
| core-v1 | constant-nonhydrostatic-512x512x129 | builtin | 1.110 | 1284.977 | 851.836 |
| core-v1 | constant-hydrostatic-512x512x129 | builtin | 1.031 | 1103.964 | 665.475 |

## Memory

| Suite | Case | Backend | Persistent increment (MiB) | Peak increment (MiB) | Provider |
|---|---|---|---:|---:|---|
| core-v1 | constant-nonhydrostatic-256x256x65 | builtin | 1274.234 | 1299.234 | macos-ps-rss |
| core-v1 | constant-hydrostatic-256x256x65 | builtin | 1222.094 | 1247.828 | macos-ps-rss |
| core-v1 | constant-nonhydrostatic-512x512x129 | builtin | 7369.875 | 7554.609 | macos-ps-rss |
| core-v1 | constant-hydrostatic-512x512x129 | builtin | 6865.062 | 7051.859 | macos-ps-rss |
