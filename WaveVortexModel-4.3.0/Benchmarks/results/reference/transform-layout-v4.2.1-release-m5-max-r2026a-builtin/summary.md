# Fourier storage layout integration

- Status: `complete`
- Gate passed: `true`
- Reference: WaveVortex issue #69 (`wv-sorted-linear`)
- Gate: no operation above `1.03x` the frozen median on the 256 and 512 workloads

## Timing and frozen-baseline comparison

| Case | Antialias | Gate | Operation | Production (ms) | Issue #69 (ms) | Relative | Error | Pass |
|---|---:|---|---|---:|---:|---:|---:|---|
| full-layout-64x48x17-antialias-0 | 0 | no | extract | 0.140 | 0.080 | 1.751 | 0 | yes |
| full-layout-64x48x17-antialias-0 | 0 | no | insert-primary | 0.087 | 0.138 | 0.626 | 0 | yes |
| full-layout-64x48x17-antialias-0 | 0 | no | insert-conjugate | 0.056 | 0.051 | 1.105 | 0 | yes |
| full-layout-64x48x17-antialias-0 | 0 | no | insert-complete | 0.160 | 0.156 | 1.031 | 0 | yes |
| full-layout-64x48x17-antialias-0 | 0 | no | forward-complete | 0.358 | 0.253 | 1.417 | 0 | yes |
| full-layout-64x48x17-antialias-0 | 0 | no | inverse-complete | 0.431 | 0.416 | 1.035 | 0 | yes |
| full-layout-64x48x17-antialias-1 | 1 | no | extract | 0.141 | 0.034 | 4.168 | 0 | yes |
| full-layout-64x48x17-antialias-1 | 1 | no | insert-primary | 0.061 | 0.081 | 0.751 | 0 | yes |
| full-layout-64x48x17-antialias-1 | 1 | no | insert-conjugate | 0.053 | 0.011 | 4.733 | 0 | yes |
| full-layout-64x48x17-antialias-1 | 1 | no | insert-complete | 0.092 | 0.081 | 1.138 | 0 | yes |
| full-layout-64x48x17-antialias-1 | 1 | no | forward-complete | 0.388 | 0.190 | 2.043 | 0 | yes |
| full-layout-64x48x17-antialias-1 | 1 | no | inverse-complete | 0.452 | 0.296 | 1.526 | 0 | yes |
| full-layout-65x63x17-antialias-0 | 0 | no | extract | 0.195 | 0.098 | 1.981 | 0 | yes |
| full-layout-65x63x17-antialias-0 | 0 | no | insert-primary | 0.124 | 0.048 | 2.576 | 0 | yes |
| full-layout-65x63x17-antialias-0 | 0 | no | insert-conjugate | 0.051 | 0.010 | 5.255 | 0 | yes |
| full-layout-65x63x17-antialias-0 | 0 | no | insert-complete | 0.119 | 0.054 | 2.192 | 0 | yes |
| full-layout-65x63x17-antialias-0 | 0 | no | forward-complete | 0.397 | 0.246 | 1.610 | 0 | yes |
| full-layout-65x63x17-antialias-0 | 0 | no | inverse-complete | 0.534 | 0.267 | 2.001 | 0 | yes |
| full-layout-65x63x17-antialias-1 | 1 | no | extract | 0.045 | 0.021 | 2.143 | 0 | yes |
| full-layout-65x63x17-antialias-1 | 1 | no | insert-primary | 0.056 | 0.021 | 2.711 | 0 | yes |
| full-layout-65x63x17-antialias-1 | 1 | no | insert-conjugate | 0.039 | 0.008 | 5.108 | 0 | yes |
| full-layout-65x63x17-antialias-1 | 1 | no | insert-complete | 0.083 | 0.024 | 3.474 | 0 | yes |
| full-layout-65x63x17-antialias-1 | 1 | no | forward-complete | 0.246 | 0.220 | 1.116 | 0 | yes |
| full-layout-65x63x17-antialias-1 | 1 | no | inverse-complete | 0.391 | 0.250 | 1.565 | 0 | yes |
| full-layout-256x256x65-antialias-0 | 0 | yes | extract | 2.776 | 5.342 | 0.520 | 0 | yes |
| full-layout-256x256x65-antialias-0 | 0 | yes | insert-primary | 2.482 | 13.716 | 0.181 | 0 | yes |
| full-layout-256x256x65-antialias-0 | 0 | yes | insert-conjugate | 0.117 | 0.125 | 0.941 | 0 | yes |
| full-layout-256x256x65-antialias-0 | 0 | yes | insert-complete | 2.510 | 13.962 | 0.180 | 0 | yes |
| full-layout-256x256x65-antialias-0 | 0 | yes | forward-complete | 4.692 | 7.774 | 0.603 | 0 | yes |
| full-layout-256x256x65-antialias-0 | 0 | yes | inverse-complete | 4.788 | 16.835 | 0.284 | 0 | yes |
| full-layout-256x256x65-antialias-1 | 1 | yes | extract | 1.024 | 3.118 | 0.328 | 0 | yes |
| full-layout-256x256x65-antialias-1 | 1 | yes | insert-primary | 0.965 | 4.274 | 0.226 | 0 | yes |
| full-layout-256x256x65-antialias-1 | 1 | yes | insert-conjugate | 0.052 | 0.077 | 0.670 | 0 | yes |
| full-layout-256x256x65-antialias-1 | 1 | yes | insert-complete | 0.807 | 4.482 | 0.180 | 0 | yes |
| full-layout-256x256x65-antialias-1 | 1 | yes | forward-complete | 3.629 | 3.800 | 0.955 | 0 | yes |
| full-layout-256x256x65-antialias-1 | 1 | yes | inverse-complete | 3.042 | 6.643 | 0.458 | 0 | yes |
| full-layout-512x512x129-antialias-0 | 0 | yes | extract | 24.898 | 88.623 | 0.281 | 0 | yes |
| full-layout-512x512x129-antialias-0 | 0 | yes | insert-primary | 20.497 | 76.126 | 0.269 | 0 | yes |
| full-layout-512x512x129-antialias-0 | 0 | yes | insert-conjugate | 0.190 | 0.316 | 0.599 | 0 | yes |
| full-layout-512x512x129-antialias-0 | 0 | yes | insert-complete | 21.597 | 76.521 | 0.282 | 0 | yes |
| full-layout-512x512x129-antialias-0 | 0 | yes | forward-complete | 47.753 | 95.419 | 0.500 | 0 | yes |
| full-layout-512x512x129-antialias-0 | 0 | yes | inverse-complete | 38.420 | 95.180 | 0.404 | 0 | yes |
| full-layout-512x512x129-antialias-1 | 1 | yes | extract | 8.736 | 29.761 | 0.294 | 0 | yes |
| full-layout-512x512x129-antialias-1 | 1 | yes | insert-primary | 7.916 | 26.782 | 0.296 | 0 | yes |
| full-layout-512x512x129-antialias-1 | 1 | yes | insert-conjugate | 0.191 | 0.227 | 0.841 | 0 | yes |
| full-layout-512x512x129-antialias-1 | 1 | yes | insert-complete | 8.048 | 27.066 | 0.297 | 0 | yes |
| full-layout-512x512x129-antialias-1 | 1 | yes | forward-complete | 31.164 | 52.374 | 0.595 | 0 | yes |
| full-layout-512x512x129-antialias-1 | 1 | yes | inverse-complete | 25.471 | 46.071 | 0.553 | 0 | yes |

## Storage contract

| Case | Strategy | Mapping (MiB) | Persistent full buffer (MiB) | Legacy maps materialized | Legacy map bytes |
|---|---|---:|---:|---|---:|
| full-layout-64x48x17-antialias-0 | two-dimensional-rows | 0.023 | 0.797 | no | 0 |
| full-layout-64x48x17-antialias-1 | two-dimensional-rows | 0.011 | 0.797 | no | 0 |
| full-layout-65x63x17-antialias-0 | two-dimensional-rows | 0.032 | 1.062 | no | 0 |
| full-layout-65x63x17-antialias-1 | two-dimensional-rows | 0.011 | 1.062 | no | 0 |
| full-layout-256x256x65-antialias-0 | two-dimensional-rows | 0.499 | 65.000 | no | 0 |
| full-layout-256x256x65-antialias-1 | two-dimensional-rows | 0.176 | 65.000 | no | 0 |
| full-layout-512x512x129-antialias-0 | two-dimensional-rows | 1.998 | 516.000 | no | 0 |
| full-layout-512x512x129-antialias-1 | two-dimensional-rows | 0.702 | 516.000 | no | 0 |
