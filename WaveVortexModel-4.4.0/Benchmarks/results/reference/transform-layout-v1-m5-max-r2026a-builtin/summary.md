# WaveVortex benchmark

- Status: `complete`
- Run: `20260809T014404Z`
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

## Construction and cache diagnostics

| Suite | Case | Backend | Construction (s) | First call (ms) | Same-state cache hit (ms) |
|---|---|---|---:|---:|---:|

## Memory

| Suite | Case | Backend | Persistent increment (MiB) | Peak increment (MiB) | Provider |
|---|---|---|---:|---:|---|

## Transform-layout diagnostic

The strict winner is the smallest median. The current `wv-sorted-linear` path remains preferred whenever it is within 3% of that median.

### Extraction and complete-forward winners

| Case | Antialias | Operation | Current (ms) | Strict fastest | Strict (ms) | Preferred | Current / fastest |
|---|---:|---|---:|---|---:|---|---:|
| full-layout-64x48x17-antialias-0 | 0 | extract | 0.080 | two-dimensional-rows | 0.061 | two-dimensional-rows | 1.306 |
| full-layout-64x48x17-antialias-0 | 0 | forward-complete | 0.253 | two-dimensional-rows | 0.206 | two-dimensional-rows | 1.229 |
| full-layout-64x48x17-antialias-1 | 1 | extract | 0.034 | two-dimensional-rows | 0.031 | two-dimensional-rows | 1.096 |
| full-layout-64x48x17-antialias-1 | 1 | forward-complete | 0.190 | wv-sorted-linear | 0.190 | wv-sorted-linear | 1.000 |
| full-layout-65x63x17-antialias-0 | 0 | extract | 0.098 | two-dimensional-rows | 0.044 | two-dimensional-rows | 2.211 |
| full-layout-65x63x17-antialias-0 | 0 | forward-complete | 0.246 | wv-sorted-linear | 0.246 | wv-sorted-linear | 1.000 |
| full-layout-65x63x17-antialias-1 | 1 | extract | 0.021 | two-dimensional-rows | 0.017 | two-dimensional-rows | 1.221 |
| full-layout-65x63x17-antialias-1 | 1 | forward-complete | 0.220 | wv-sorted-linear | 0.220 | wv-sorted-linear | 1.000 |
| full-layout-256x256x65-antialias-0 | 0 | extract | 5.342 | two-dimensional-rows | 2.762 | two-dimensional-rows | 1.934 |
| full-layout-256x256x65-antialias-0 | 0 | forward-complete | 7.774 | two-dimensional-rows | 5.100 | two-dimensional-rows | 1.524 |
| full-layout-256x256x65-antialias-1 | 1 | extract | 3.118 | two-dimensional-rows | 0.949 | two-dimensional-rows | 3.285 |
| full-layout-256x256x65-antialias-1 | 1 | forward-complete | 3.800 | two-dimensional-rows | 3.157 | two-dimensional-rows | 1.204 |
| full-layout-512x512x129-antialias-0 | 0 | extract | 88.623 | two-dimensional-rows | 26.202 | two-dimensional-rows | 3.382 |
| full-layout-512x512x129-antialias-0 | 0 | forward-complete | 95.419 | two-dimensional-rows | 44.311 | two-dimensional-rows | 2.153 |
| full-layout-512x512x129-antialias-1 | 1 | extract | 29.761 | two-dimensional-rows | 9.076 | two-dimensional-rows | 3.279 |
| full-layout-512x512x129-antialias-1 | 1 | forward-complete | 52.374 | two-dimensional-rows | 27.425 | two-dimensional-rows | 1.910 |

### Insertion and complete-inverse winners

| Case | Antialias | Operation | Current (ms) | Strict fastest | Strict (ms) | Preferred | Current / fastest |
|---|---:|---|---:|---|---:|---|---:|
| full-layout-64x48x17-antialias-0 | 0 | insert-primary | 0.138 | two-dimensional-rows | 0.056 | two-dimensional-rows | 2.459 |
| full-layout-64x48x17-antialias-0 | 0 | insert-conjugate | 0.051 | two-dimensional-rows | 0.032 | two-dimensional-rows | 1.572 |
| full-layout-64x48x17-antialias-0 | 0 | insert-complete | 0.156 | two-dimensional-rows | 0.093 | two-dimensional-rows | 1.678 |
| full-layout-64x48x17-antialias-0 | 0 | inverse-complete | 0.416 | dft-sorted-linear | 0.320 | dft-sorted-linear | 1.299 |
| full-layout-64x48x17-antialias-1 | 1 | insert-primary | 0.081 | two-dimensional-rows | 0.039 | two-dimensional-rows | 2.060 |
| full-layout-64x48x17-antialias-1 | 1 | insert-conjugate | 0.011 | dft-sorted-linear | 0.006 | dft-sorted-linear | 1.862 |
| full-layout-64x48x17-antialias-1 | 1 | insert-complete | 0.081 | two-dimensional-rows | 0.049 | two-dimensional-rows | 1.645 |
| full-layout-64x48x17-antialias-1 | 1 | inverse-complete | 0.296 | dft-sorted-linear | 0.279 | dft-sorted-linear | 1.062 |
| full-layout-65x63x17-antialias-0 | 0 | insert-primary | 0.048 | wv-sorted-linear | 0.048 | wv-sorted-linear | 1.000 |
| full-layout-65x63x17-antialias-0 | 0 | insert-conjugate | 0.010 | dft-sorted-linear | 0.007 | dft-sorted-linear | 1.358 |
| full-layout-65x63x17-antialias-0 | 0 | insert-complete | 0.054 | two-dimensional-rows | 0.045 | two-dimensional-rows | 1.214 |
| full-layout-65x63x17-antialias-0 | 0 | inverse-complete | 0.267 | wv-sorted-linear | 0.267 | wv-sorted-linear | 1.000 |
| full-layout-65x63x17-antialias-1 | 1 | insert-primary | 0.021 | wv-sorted-linear | 0.021 | wv-sorted-linear | 1.000 |
| full-layout-65x63x17-antialias-1 | 1 | insert-conjugate | 0.008 | dft-sorted-linear | 0.006 | dft-sorted-linear | 1.341 |
| full-layout-65x63x17-antialias-1 | 1 | insert-complete | 0.024 | two-dimensional-rows | 0.022 | two-dimensional-rows | 1.077 |
| full-layout-65x63x17-antialias-1 | 1 | inverse-complete | 0.250 | wv-sorted-linear | 0.250 | wv-sorted-linear | 1.000 |
| full-layout-256x256x65-antialias-0 | 0 | insert-primary | 13.716 | two-dimensional-rows | 1.965 | two-dimensional-rows | 6.978 |
| full-layout-256x256x65-antialias-0 | 0 | insert-conjugate | 0.125 | dft-sorted-linear | 0.052 | dft-sorted-linear | 2.396 |
| full-layout-256x256x65-antialias-0 | 0 | insert-complete | 13.962 | two-dimensional-rows | 1.941 | two-dimensional-rows | 7.194 |
| full-layout-256x256x65-antialias-0 | 0 | inverse-complete | 16.835 | two-dimensional-rows | 4.172 | two-dimensional-rows | 4.036 |
| full-layout-256x256x65-antialias-1 | 1 | insert-primary | 4.274 | two-dimensional-rows | 0.753 | two-dimensional-rows | 5.674 |
| full-layout-256x256x65-antialias-1 | 1 | insert-conjugate | 0.077 | dft-sorted-linear | 0.034 | dft-sorted-linear | 2.280 |
| full-layout-256x256x65-antialias-1 | 1 | insert-complete | 4.482 | two-dimensional-rows | 0.765 | two-dimensional-rows | 5.856 |
| full-layout-256x256x65-antialias-1 | 1 | inverse-complete | 6.643 | two-dimensional-rows | 2.977 | two-dimensional-rows | 2.232 |
| full-layout-512x512x129-antialias-0 | 0 | insert-primary | 76.126 | two-dimensional-rows | 15.314 | two-dimensional-rows | 4.971 |
| full-layout-512x512x129-antialias-0 | 0 | insert-conjugate | 0.316 | two-dimensional-rows | 0.129 | two-dimensional-rows | 2.453 |
| full-layout-512x512x129-antialias-0 | 0 | insert-complete | 76.521 | two-dimensional-rows | 15.538 | two-dimensional-rows | 4.925 |
| full-layout-512x512x129-antialias-0 | 0 | inverse-complete | 95.180 | two-dimensional-rows | 34.174 | two-dimensional-rows | 2.785 |
| full-layout-512x512x129-antialias-1 | 1 | insert-primary | 26.782 | two-dimensional-rows | 5.733 | two-dimensional-rows | 4.672 |
| full-layout-512x512x129-antialias-1 | 1 | insert-conjugate | 0.227 | two-dimensional-rows | 0.176 | two-dimensional-rows | 1.290 |
| full-layout-512x512x129-antialias-1 | 1 | insert-complete | 27.066 | two-dimensional-rows | 5.840 | two-dimensional-rows | 4.635 |
| full-layout-512x512x129-antialias-1 | 1 | inverse-complete | 46.071 | two-dimensional-rows | 24.330 | two-dimensional-rows | 1.894 |

### Mapping-array and working-array storage

| Case | Strategy | Mapping arrays (MiB) | Persistent full buffer (MiB) | Real input (MiB) | Full spectrum (MiB) | WV source (MiB) | Timed WV result (MiB) | Timed real result (MiB) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| full-layout-64x48x17-antialias-0 | wv-sorted-linear | 0.200 | 0.797 | 0.398 | 0.797 | 0.384 | 0.384 | 0.398 |
| full-layout-64x48x17-antialias-0 | dft-sorted-linear | 0.392 | 0.797 | 0.398 | 0.797 | 0.384 | 0.384 | 0.398 |
| full-layout-64x48x17-antialias-0 | two-dimensional-rows | 0.012 | 0.797 | 0.398 | 0.797 | 0.384 | 0.384 | 0.398 |
| full-layout-64x48x17-antialias-0 | per-plane-linear | 0.012 | 0.797 | 0.398 | 0.797 | 0.384 | 0.384 | 0.398 |
| full-layout-64x48x17-antialias-1 | wv-sorted-linear | 0.098 | 0.797 | 0.398 | 0.797 | 0.186 | 0.186 | 0.398 |
| full-layout-64x48x17-antialias-1 | dft-sorted-linear | 0.191 | 0.797 | 0.398 | 0.797 | 0.186 | 0.186 | 0.398 |
| full-layout-64x48x17-antialias-1 | two-dimensional-rows | 0.006 | 0.797 | 0.398 | 0.797 | 0.186 | 0.186 | 0.398 |
| full-layout-64x48x17-antialias-1 | per-plane-linear | 0.006 | 0.797 | 0.398 | 0.797 | 0.186 | 0.186 | 0.398 |
| full-layout-65x63x17-antialias-0 | wv-sorted-linear | 0.274 | 1.062 | 0.531 | 1.062 | 0.531 | 0.531 | 0.531 |
| full-layout-65x63x17-antialias-0 | dft-sorted-linear | 0.540 | 1.062 | 0.531 | 1.062 | 0.531 | 0.531 | 0.531 |
| full-layout-65x63x17-antialias-0 | two-dimensional-rows | 0.016 | 1.062 | 0.531 | 1.062 | 0.531 | 0.531 | 0.531 |
| full-layout-65x63x17-antialias-0 | per-plane-linear | 0.016 | 1.062 | 0.531 | 1.062 | 0.531 | 0.531 | 0.531 |
| full-layout-65x63x17-antialias-1 | wv-sorted-linear | 0.098 | 1.062 | 0.531 | 1.062 | 0.186 | 0.186 | 0.531 |
| full-layout-65x63x17-antialias-1 | dft-sorted-linear | 0.191 | 1.062 | 0.531 | 1.062 | 0.186 | 0.186 | 0.531 |
| full-layout-65x63x17-antialias-1 | two-dimensional-rows | 0.006 | 1.062 | 0.531 | 1.062 | 0.186 | 0.186 | 0.531 |
| full-layout-65x63x17-antialias-1 | per-plane-linear | 0.006 | 1.062 | 0.531 | 1.062 | 0.186 | 0.186 | 0.531 |
| full-layout-256x256x65-antialias-0 | wv-sorted-linear | 16.250 | 65.000 | 32.500 | 65.000 | 32.247 | 32.247 | 32.500 |
| full-layout-256x256x65-antialias-0 | dft-sorted-linear | 32.373 | 65.000 | 32.500 | 65.000 | 32.247 | 32.247 | 32.500 |
| full-layout-256x256x65-antialias-0 | two-dimensional-rows | 0.250 | 65.000 | 32.500 | 65.000 | 32.247 | 32.247 | 32.500 |
| full-layout-256x256x65-antialias-0 | per-plane-linear | 0.250 | 65.000 | 32.500 | 65.000 | 32.247 | 32.247 | 32.500 |
| full-layout-256x256x65-antialias-1 | wv-sorted-linear | 5.757 | 65.000 | 32.500 | 65.000 | 11.345 | 11.345 | 32.500 |
| full-layout-256x256x65-antialias-1 | dft-sorted-linear | 11.430 | 65.000 | 32.500 | 65.000 | 11.345 | 11.345 | 32.500 |
| full-layout-256x256x65-antialias-1 | two-dimensional-rows | 0.089 | 65.000 | 32.500 | 65.000 | 11.345 | 11.345 | 32.500 |
| full-layout-256x256x65-antialias-1 | per-plane-linear | 0.089 | 65.000 | 32.500 | 65.000 | 11.345 | 11.345 | 32.500 |
| full-layout-512x512x129-antialias-0 | wv-sorted-linear | 128.999 | 516.000 | 258.000 | 516.000 | 256.994 | 256.994 | 258.000 |
| full-layout-512x512x129-antialias-0 | dft-sorted-linear | 257.496 | 516.000 | 258.000 | 516.000 | 256.994 | 256.994 | 258.000 |
| full-layout-512x512x129-antialias-0 | two-dimensional-rows | 1.000 | 516.000 | 258.000 | 516.000 | 256.994 | 256.994 | 258.000 |
| full-layout-512x512x129-antialias-0 | per-plane-linear | 1.001 | 516.000 | 258.000 | 516.000 | 256.994 | 256.994 | 258.000 |
| full-layout-512x512x129-antialias-1 | wv-sorted-linear | 45.376 | 516.000 | 258.000 | 516.000 | 90.083 | 90.083 | 258.000 |
| full-layout-512x512x129-antialias-1 | dft-sorted-linear | 90.418 | 516.000 | 258.000 | 516.000 | 90.083 | 90.083 | 258.000 |
| full-layout-512x512x129-antialias-1 | two-dimensional-rows | 0.352 | 516.000 | 258.000 | 516.000 | 90.083 | 90.083 | 258.000 |
| full-layout-512x512x129-antialias-1 | per-plane-linear | 0.353 | 516.000 | 258.000 | 516.000 | 90.083 | 90.083 | 258.000 |

### Correctness and observable copy semantics

| Case | Strategy | Maximum error | Source arrays unchanged | Buffer reused | Timed clearing | Copy status |
|---|---|---:|---|---|---|---|
| full-layout-64x48x17-antialias-0 | wv-sorted-linear | 0 | yes | yes | no | unavailable |
| full-layout-64x48x17-antialias-0 | dft-sorted-linear | 0 | yes | yes | no | unavailable |
| full-layout-64x48x17-antialias-0 | two-dimensional-rows | 0 | yes | yes | no | unavailable |
| full-layout-64x48x17-antialias-0 | per-plane-linear | 0 | yes | yes | no | unavailable |
| full-layout-64x48x17-antialias-1 | wv-sorted-linear | 0 | yes | yes | no | unavailable |
| full-layout-64x48x17-antialias-1 | dft-sorted-linear | 0 | yes | yes | no | unavailable |
| full-layout-64x48x17-antialias-1 | two-dimensional-rows | 0 | yes | yes | no | unavailable |
| full-layout-64x48x17-antialias-1 | per-plane-linear | 0 | yes | yes | no | unavailable |
| full-layout-65x63x17-antialias-0 | wv-sorted-linear | 0 | yes | yes | no | unavailable |
| full-layout-65x63x17-antialias-0 | dft-sorted-linear | 0 | yes | yes | no | unavailable |
| full-layout-65x63x17-antialias-0 | two-dimensional-rows | 0 | yes | yes | no | unavailable |
| full-layout-65x63x17-antialias-0 | per-plane-linear | 0 | yes | yes | no | unavailable |
| full-layout-65x63x17-antialias-1 | wv-sorted-linear | 0 | yes | yes | no | unavailable |
| full-layout-65x63x17-antialias-1 | dft-sorted-linear | 0 | yes | yes | no | unavailable |
| full-layout-65x63x17-antialias-1 | two-dimensional-rows | 0 | yes | yes | no | unavailable |
| full-layout-65x63x17-antialias-1 | per-plane-linear | 0 | yes | yes | no | unavailable |
| full-layout-256x256x65-antialias-0 | wv-sorted-linear | 0 | yes | yes | no | unavailable |
| full-layout-256x256x65-antialias-0 | dft-sorted-linear | 0 | yes | yes | no | unavailable |
| full-layout-256x256x65-antialias-0 | two-dimensional-rows | 0 | yes | yes | no | unavailable |
| full-layout-256x256x65-antialias-0 | per-plane-linear | 0 | yes | yes | no | unavailable |
| full-layout-256x256x65-antialias-1 | wv-sorted-linear | 0 | yes | yes | no | unavailable |
| full-layout-256x256x65-antialias-1 | dft-sorted-linear | 0 | yes | yes | no | unavailable |
| full-layout-256x256x65-antialias-1 | two-dimensional-rows | 0 | yes | yes | no | unavailable |
| full-layout-256x256x65-antialias-1 | per-plane-linear | 0 | yes | yes | no | unavailable |
| full-layout-512x512x129-antialias-0 | wv-sorted-linear | 0 | yes | yes | no | unavailable |
| full-layout-512x512x129-antialias-0 | dft-sorted-linear | 0 | yes | yes | no | unavailable |
| full-layout-512x512x129-antialias-0 | two-dimensional-rows | 0 | yes | yes | no | unavailable |
| full-layout-512x512x129-antialias-0 | per-plane-linear | 0 | yes | yes | no | unavailable |
| full-layout-512x512x129-antialias-1 | wv-sorted-linear | 0 | yes | yes | no | unavailable |
| full-layout-512x512x129-antialias-1 | dft-sorted-linear | 0 | yes | yes | no | unavailable |
| full-layout-512x512x129-antialias-1 | two-dimensional-rows | 0 | yes | yes | no | unavailable |
| full-layout-512x512x129-antialias-1 | per-plane-linear | 0 | yes | yes | no | unavailable |
