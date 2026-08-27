# WaveVortex benchmark

- Status: `partial`
- Run: `20260808T190619Z`
- MATLAB: `2026a`
- Architecture: `maca64`

## Suite scores

| Suite | Backend | Score |
|---|---|---:|

## Family scores

| Suite | Family | Backend | Score |
|---|---|---|---:|

## Timing and scores

| Suite | Case | Transform | Backend | Median (ms) | Reference score | Same-host speedup | Error |
|---|---|---|---|---:|---:|---:|---:|
| scaling-large-v1 | constant-nonhydrostatic-256x256x129 | constant-nonhydrostatic | builtin | 255.151 | 100.000 | 1.000 | 0 |
| scaling-large-v1 | constant-nonhydrostatic-512x512x257 | constant-nonhydrostatic | builtin | 2329.690 | 100.000 | 1.000 | 0 |
| scaling-large-v1 | constant-hydrostatic-256x256x129 | constant-hydrostatic | builtin | 218.237 | 100.000 | 1.000 | 0 |
| scaling-large-v1 | constant-hydrostatic-512x512x257 | constant-hydrostatic | builtin | 2152.495 | 100.000 | 1.000 | 0 |
| scaling-large-v1 | hydrostatic-512x512x129 | hydrostatic | builtin | 856.233 | 100.000 | 1.000 | 0 |
| scaling-large-v1 | stratified-qg-512x512x129 | stratified-qg | builtin | 345.178 | 100.000 | 1.000 | 0 |
| scaling-large-v1 | barotropic-qg-2048x2048 | barotropic-qg | builtin | 30.321 | 100.000 | 1.000 | 0 |
| scaling-large-v1 | barotropic-qg-4096x4096 | barotropic-qg | builtin | 112.796 | 100.000 | 1.000 | 0 |

## Construction and cache diagnostics

| Suite | Case | Backend | Construction (s) | First call (ms) | Same-state cache hit (ms) |
|---|---|---|---:|---:|---:|
| scaling-large-v1 | constant-nonhydrostatic-256x256x129 | builtin | 1.297 | 407.072 | 183.336 |
| scaling-large-v1 | constant-nonhydrostatic-512x512x257 | builtin | 1.953 | 2556.639 | 1947.495 |
| scaling-large-v1 | constant-hydrostatic-256x256x129 | builtin | 0.475 | 287.309 | 161.202 |
| scaling-large-v1 | constant-hydrostatic-512x512x257 | builtin | 1.708 | 2180.815 | 1473.488 |
| scaling-large-v1 | hydrostatic-512x512x129 | builtin | 1.638 | 985.977 | 567.364 |
| scaling-large-v1 | stratified-qg-512x512x129 | builtin | 0.620 | 335.157 | 127.158 |
| scaling-large-v1 | barotropic-qg-2048x2048 | builtin | 0.942 | 51.141 | 15.813 |
| scaling-large-v1 | barotropic-qg-4096x4096 | builtin | 3.543 | 176.627 | 57.776 |

## Memory

| Suite | Case | Backend | Persistent increment (MiB) | Peak increment (MiB) | Provider |
|---|---|---|---:|---:|---|
| scaling-large-v1 | constant-nonhydrostatic-256x256x129 | builtin | 2328.812 | 2374.938 | macos-ps-rss |
| scaling-large-v1 | constant-nonhydrostatic-512x512x257 | builtin | 12616.469 | 16148.188 | macos-ps-rss |
| scaling-large-v1 | constant-hydrostatic-256x256x129 | builtin | 2228.344 | 2274.750 | macos-ps-rss |
| scaling-large-v1 | constant-hydrostatic-512x512x257 | builtin | 14166.641 | 14533.219 | macos-ps-rss |
| scaling-large-v1 | hydrostatic-512x512x129 | builtin | 7967.516 | 8178.234 | macos-ps-rss |
| scaling-large-v1 | stratified-qg-512x512x129 | builtin | 4252.141 | 4342.516 | macos-ps-rss |
| scaling-large-v1 | barotropic-qg-2048x2048 | builtin | 988.578 | 1000.062 | macos-ps-rss |
| scaling-large-v1 | barotropic-qg-4096x4096 | builtin | 3486.719 | 3531.609 | macos-ps-rss |
