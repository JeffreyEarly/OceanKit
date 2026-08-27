# WaveVortex benchmark

- Status: `complete`
- Run: `20260808T180708Z`
- MATLAB: `2026a`
- Architecture: `maca64`

## Suite scores

| Suite | Backend | Score |
|---|---|---:|
| scaling-standard-v1 | builtin | 100.000 |

## Family scores

| Suite | Family | Backend | Score |
|---|---|---|---:|
| scaling-standard-v1 | constant-nonhydrostatic | builtin | 100.000 |
| scaling-standard-v1 | constant-hydrostatic | builtin | 100.000 |
| scaling-standard-v1 | hydrostatic | builtin | 100.000 |
| scaling-standard-v1 | boussinesq | builtin | 100.000 |
| scaling-standard-v1 | stratified-qg | builtin | 100.000 |
| scaling-standard-v1 | barotropic-qg | builtin | 100.000 |

## Timing and scores

| Suite | Case | Transform | Backend | Median (ms) | Reference score | Same-host speedup | Error |
|---|---|---|---|---:|---:|---:|---:|
| scaling-standard-v1 | constant-nonhydrostatic-64x64x65 | constant-nonhydrostatic | builtin | 26.006 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | constant-nonhydrostatic-128x128x65 | constant-nonhydrostatic | builtin | 39.273 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | constant-nonhydrostatic-256x256x65 | constant-nonhydrostatic | builtin | 147.021 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | constant-nonhydrostatic-128x128x33 | constant-nonhydrostatic | builtin | 24.790 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | constant-nonhydrostatic-128x128x129 | constant-nonhydrostatic | builtin | 81.329 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | constant-nonhydrostatic-128x128x257 | constant-nonhydrostatic | builtin | 155.925 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | constant-hydrostatic-64x64x65 | constant-hydrostatic | builtin | 14.292 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | constant-hydrostatic-128x128x65 | constant-hydrostatic | builtin | 32.560 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | constant-hydrostatic-256x256x65 | constant-hydrostatic | builtin | 107.475 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | constant-hydrostatic-128x128x33 | constant-hydrostatic | builtin | 20.584 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | constant-hydrostatic-128x128x129 | constant-hydrostatic | builtin | 62.076 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | constant-hydrostatic-128x128x257 | constant-hydrostatic | builtin | 130.456 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | hydrostatic-64x64x65 | hydrostatic | builtin | 17.428 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | hydrostatic-128x128x65 | hydrostatic | builtin | 29.580 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | hydrostatic-256x256x65 | hydrostatic | builtin | 97.785 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | hydrostatic-128x128x33 | hydrostatic | builtin | 16.608 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | hydrostatic-128x128x129 | hydrostatic | builtin | 62.691 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | hydrostatic-128x128x257 | hydrostatic | builtin | 130.280 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | boussinesq-64x64x65 | boussinesq | builtin | 25.506 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | boussinesq-128x128x65 | boussinesq | builtin | 54.721 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | boussinesq-256x256x65 | boussinesq | builtin | 192.131 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | boussinesq-128x128x33 | boussinesq | builtin | 29.429 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | boussinesq-128x128x129 | boussinesq | builtin | 111.942 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | boussinesq-128x128x257 | boussinesq | builtin | 273.177 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | stratified-qg-64x64x65 | stratified-qg | builtin | 5.281 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | stratified-qg-128x128x65 | stratified-qg | builtin | 12.875 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | stratified-qg-256x256x65 | stratified-qg | builtin | 35.782 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | stratified-qg-128x128x33 | stratified-qg | builtin | 6.743 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | stratified-qg-128x128x129 | stratified-qg | builtin | 20.626 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | stratified-qg-128x128x257 | stratified-qg | builtin | 48.354 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | barotropic-qg-128x128 | barotropic-qg | builtin | 3.698 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | barotropic-qg-256x256 | barotropic-qg | builtin | 2.968 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | barotropic-qg-512x512 | barotropic-qg | builtin | 3.776 | 100.000 | 1.000 | 0 |
| scaling-standard-v1 | barotropic-qg-1024x1024 | barotropic-qg | builtin | 8.441 | 100.000 | 1.000 | 0 |

## Construction and cache diagnostics

| Suite | Case | Backend | Construction (s) | First call (ms) | Same-state cache hit (ms) |
|---|---|---|---:|---:|---:|
| scaling-standard-v1 | constant-nonhydrostatic-64x64x65 | builtin | 1.103 | 68.577 | 14.404 |
| scaling-standard-v1 | constant-nonhydrostatic-128x128x65 | builtin | 0.415 | 55.587 | 26.461 |
| scaling-standard-v1 | constant-nonhydrostatic-256x256x65 | builtin | 0.396 | 169.847 | 100.701 |
| scaling-standard-v1 | constant-nonhydrostatic-128x128x33 | builtin | 0.256 | 26.413 | 16.773 |
| scaling-standard-v1 | constant-nonhydrostatic-128x128x129 | builtin | 0.275 | 92.140 | 44.545 |
| scaling-standard-v1 | constant-nonhydrostatic-128x128x257 | builtin | 0.315 | 193.880 | 102.291 |
| scaling-standard-v1 | constant-hydrostatic-64x64x65 | builtin | 0.219 | 25.505 | 9.832 |
| scaling-standard-v1 | constant-hydrostatic-128x128x65 | builtin | 0.241 | 37.686 | 19.870 |
| scaling-standard-v1 | constant-hydrostatic-256x256x65 | builtin | 0.284 | 128.396 | 72.468 |
| scaling-standard-v1 | constant-hydrostatic-128x128x33 | builtin | 0.185 | 20.859 | 13.235 |
| scaling-standard-v1 | constant-hydrostatic-128x128x129 | builtin | 0.229 | 75.473 | 33.390 |
| scaling-standard-v1 | constant-hydrostatic-128x128x257 | builtin | 0.271 | 170.917 | 79.898 |
| scaling-standard-v1 | hydrostatic-64x64x65 | builtin | 0.905 | 66.054 | 7.846 |
| scaling-standard-v1 | hydrostatic-128x128x65 | builtin | 0.475 | 40.392 | 16.695 |
| scaling-standard-v1 | hydrostatic-256x256x65 | builtin | 0.518 | 117.359 | 58.730 |
| scaling-standard-v1 | hydrostatic-128x128x33 | builtin | 0.403 | 17.941 | 10.854 |
| scaling-standard-v1 | hydrostatic-128x128x129 | builtin | 0.435 | 72.699 | 35.156 |
| scaling-standard-v1 | hydrostatic-128x128x257 | builtin | 1.227 | 157.527 | 76.102 |
| scaling-standard-v1 | boussinesq-64x64x65 | builtin | 14.149 | 126.945 | 10.861 |
| scaling-standard-v1 | boussinesq-128x128x65 | builtin | 47.897 | 68.012 | 26.646 |
| scaling-standard-v1 | boussinesq-256x256x65 | builtin | 178.319 | 219.589 | 100.504 |
| scaling-standard-v1 | boussinesq-128x128x33 | builtin | 46.328 | 34.257 | 16.155 |
| scaling-standard-v1 | boussinesq-128x128x129 | builtin | 33.107 | 123.344 | 57.136 |
| scaling-standard-v1 | boussinesq-128x128x257 | builtin | 300.038 | 318.072 | 142.308 |
| scaling-standard-v1 | stratified-qg-64x64x65 | builtin | 0.409 | 31.279 | 2.806 |
| scaling-standard-v1 | stratified-qg-128x128x65 | builtin | 0.305 | 13.951 | 4.519 |
| scaling-standard-v1 | stratified-qg-256x256x65 | builtin | 0.291 | 43.264 | 13.918 |
| scaling-standard-v1 | stratified-qg-128x128x33 | builtin | 0.234 | 7.171 | 3.105 |
| scaling-standard-v1 | stratified-qg-128x128x129 | builtin | 0.209 | 21.589 | 6.619 |
| scaling-standard-v1 | stratified-qg-128x128x257 | builtin | 0.937 | 52.059 | 13.550 |
| scaling-standard-v1 | barotropic-qg-128x128 | builtin | 0.178 | 37.912 | 1.770 |
| scaling-standard-v1 | barotropic-qg-256x256 | builtin | 0.081 | 4.076 | 1.014 |
| scaling-standard-v1 | barotropic-qg-512x512 | builtin | 0.079 | 4.356 | 1.806 |
| scaling-standard-v1 | barotropic-qg-1024x1024 | builtin | 0.198 | 15.293 | 4.420 |

## Memory

| Suite | Case | Backend | Persistent increment (MiB) | Peak increment (MiB) | Provider |
|---|---|---|---:|---:|---|
| scaling-standard-v1 | constant-nonhydrostatic-64x64x65 | builtin | 265.688 | 272.922 | macos-ps-rss |
| scaling-standard-v1 | constant-nonhydrostatic-128x128x65 | builtin | 466.859 | 478.734 | macos-ps-rss |
| scaling-standard-v1 | constant-nonhydrostatic-256x256x65 | builtin | 1272.688 | 1297.828 | macos-ps-rss |
| scaling-standard-v1 | constant-nonhydrostatic-128x128x33 | builtin | 329.797 | 341.891 | macos-ps-rss |
| scaling-standard-v1 | constant-nonhydrostatic-128x128x129 | builtin | 749.906 | 764.016 | macos-ps-rss |
| scaling-standard-v1 | constant-nonhydrostatic-128x128x257 | builtin | 1266.406 | 1291.531 | macos-ps-rss |
| scaling-standard-v1 | constant-hydrostatic-64x64x65 | builtin | 261.344 | 268.422 | macos-ps-rss |
| scaling-standard-v1 | constant-hydrostatic-128x128x65 | builtin | 454.906 | 463.859 | macos-ps-rss |
| scaling-standard-v1 | constant-hydrostatic-256x256x65 | builtin | 1220.266 | 1245.750 | macos-ps-rss |
| scaling-standard-v1 | constant-hydrostatic-128x128x33 | builtin | 323.672 | 334.266 | macos-ps-rss |
| scaling-standard-v1 | constant-hydrostatic-128x128x129 | builtin | 728.031 | 742.312 | macos-ps-rss |
| scaling-standard-v1 | constant-hydrostatic-128x128x257 | builtin | 1214.828 | 1240.266 | macos-ps-rss |
| scaling-standard-v1 | hydrostatic-64x64x65 | builtin | 287.172 | 295.531 | macos-ps-rss |
| scaling-standard-v1 | hydrostatic-128x128x65 | builtin | 478.719 | 487.953 | macos-ps-rss |
| scaling-standard-v1 | hydrostatic-256x256x65 | builtin | 1220.344 | 1250.016 | macos-ps-rss |
| scaling-standard-v1 | hydrostatic-128x128x33 | builtin | 346.219 | 356.797 | macos-ps-rss |
| scaling-standard-v1 | hydrostatic-128x128x129 | builtin | 744.766 | 761.266 | macos-ps-rss |
| scaling-standard-v1 | hydrostatic-128x128x257 | builtin | 1239.438 | 1268.375 | macos-ps-rss |
| scaling-standard-v1 | boussinesq-64x64x65 | builtin | 339.953 | 354.750 | macos-ps-rss |
| scaling-standard-v1 | boussinesq-128x128x65 | builtin | 577.734 | 588.891 | macos-ps-rss |
| scaling-standard-v1 | boussinesq-256x256x65 | builtin | 1512.594 | 1546.156 | macos-ps-rss |
| scaling-standard-v1 | boussinesq-128x128x33 | builtin | 397.203 | 406.484 | macos-ps-rss |
| scaling-standard-v1 | boussinesq-128x128x129 | builtin | 1043.703 | 1062.281 | macos-ps-rss |
| scaling-standard-v1 | boussinesq-128x128x257 | builtin | 2344.031 | 2376.906 | macos-ps-rss |
| scaling-standard-v1 | stratified-qg-64x64x65 | builtin | 222.797 | 226.234 | macos-ps-rss |
| scaling-standard-v1 | stratified-qg-128x128x65 | builtin | 314.609 | 318.016 | macos-ps-rss |
| scaling-standard-v1 | stratified-qg-256x256x65 | builtin | 700.406 | 712.938 | macos-ps-rss |
| scaling-standard-v1 | stratified-qg-128x128x33 | builtin | 259.578 | 262.516 | macos-ps-rss |
| scaling-standard-v1 | stratified-qg-128x128x129 | builtin | 472.125 | 479.469 | macos-ps-rss |
| scaling-standard-v1 | stratified-qg-128x128x257 | builtin | 729.828 | 742.406 | macos-ps-rss |
| scaling-standard-v1 | barotropic-qg-128x128 | builtin | 152.969 | 154.422 | macos-ps-rss |
| scaling-standard-v1 | barotropic-qg-256x256 | builtin | 160.656 | 162.406 | macos-ps-rss |
| scaling-standard-v1 | barotropic-qg-512x512 | builtin | 202.812 | 206.672 | macos-ps-rss |
| scaling-standard-v1 | barotropic-qg-1024x1024 | builtin | 360.141 | 365.125 | macos-ps-rss |
