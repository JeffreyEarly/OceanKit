# Fourier spectrum layout integration

- Status: `complete`
- Gate passed: `true`
- Reference: WaveVortex issue #69 (`wv-sorted-linear`)
- Gate: no operation above `1.03x` the frozen median on the 256 and 512 workloads

## Timing and frozen-baseline comparison

| Case | Antialias | Gate | Operation | Production (ms) | Issue #69 (ms) | Relative | Error | Pass |
|---|---:|---|---|---:|---:|---:|---:|---|
| full-layout-64x48x17-antialias-0 | 0 | no | extract | 0.147 | 0.080 | 1.830 | 0 | yes |
| full-layout-64x48x17-antialias-0 | 0 | no | insert-primary | 0.095 | 0.138 | 0.688 | 0 | yes |
| full-layout-64x48x17-antialias-0 | 0 | no | insert-conjugate | 0.056 | 0.051 | 1.108 | 0 | yes |
| full-layout-64x48x17-antialias-0 | 0 | no | insert-complete | 0.140 | 0.156 | 0.903 | 0 | yes |
| full-layout-64x48x17-antialias-0 | 0 | no | forward-complete | 0.350 | 0.253 | 1.384 | 0 | yes |
| full-layout-64x48x17-antialias-0 | 0 | no | inverse-complete | 0.503 | 0.416 | 1.209 | 0 | yes |
| full-layout-64x48x17-antialias-1 | 1 | no | extract | 0.145 | 0.034 | 4.279 | 0 | yes |
| full-layout-64x48x17-antialias-1 | 1 | no | insert-primary | 0.060 | 0.081 | 0.736 | 0 | yes |
| full-layout-64x48x17-antialias-1 | 1 | no | insert-conjugate | 0.060 | 0.011 | 5.344 | 0 | yes |
| full-layout-64x48x17-antialias-1 | 1 | no | insert-complete | 0.100 | 0.081 | 1.237 | 0 | yes |
| full-layout-64x48x17-antialias-1 | 1 | no | forward-complete | 0.336 | 0.190 | 1.772 | 0 | yes |
| full-layout-64x48x17-antialias-1 | 1 | no | inverse-complete | 0.433 | 0.296 | 1.462 | 0 | yes |
| full-layout-65x63x17-antialias-0 | 0 | no | extract | 0.174 | 0.098 | 1.766 | 0 | yes |
| full-layout-65x63x17-antialias-0 | 0 | no | insert-primary | 0.138 | 0.048 | 2.870 | 0 | yes |
| full-layout-65x63x17-antialias-0 | 0 | no | insert-conjugate | 0.050 | 0.010 | 5.068 | 0 | yes |
| full-layout-65x63x17-antialias-0 | 0 | no | insert-complete | 0.128 | 0.054 | 2.364 | 0 | yes |
| full-layout-65x63x17-antialias-0 | 0 | no | forward-complete | 0.392 | 0.246 | 1.592 | 0 | yes |
| full-layout-65x63x17-antialias-0 | 0 | no | inverse-complete | 0.517 | 0.267 | 1.938 | 0 | yes |
| full-layout-65x63x17-antialias-1 | 1 | no | extract | 0.045 | 0.021 | 2.145 | 0 | yes |
| full-layout-65x63x17-antialias-1 | 1 | no | insert-primary | 0.074 | 0.021 | 3.612 | 0 | yes |
| full-layout-65x63x17-antialias-1 | 1 | no | insert-conjugate | 0.039 | 0.008 | 5.016 | 0 | yes |
| full-layout-65x63x17-antialias-1 | 1 | no | insert-complete | 0.082 | 0.024 | 3.432 | 0 | yes |
| full-layout-65x63x17-antialias-1 | 1 | no | forward-complete | 0.259 | 0.220 | 1.176 | 0 | yes |
| full-layout-65x63x17-antialias-1 | 1 | no | inverse-complete | 0.447 | 0.250 | 1.790 | 0 | yes |
| full-layout-256x256x65-antialias-0 | 0 | yes | extract | 2.777 | 5.342 | 0.520 | 0 | yes |
| full-layout-256x256x65-antialias-0 | 0 | yes | insert-primary | 2.471 | 13.716 | 0.180 | 0 | yes |
| full-layout-256x256x65-antialias-0 | 0 | yes | insert-conjugate | 0.124 | 0.125 | 0.995 | 0 | yes |
| full-layout-256x256x65-antialias-0 | 0 | yes | insert-complete | 2.542 | 13.962 | 0.182 | 0 | yes |
| full-layout-256x256x65-antialias-0 | 0 | yes | forward-complete | 4.594 | 7.774 | 0.591 | 0 | yes |
| full-layout-256x256x65-antialias-0 | 0 | yes | inverse-complete | 4.701 | 16.835 | 0.279 | 0 | yes |
| full-layout-256x256x65-antialias-1 | 1 | yes | extract | 1.019 | 3.118 | 0.327 | 0 | yes |
| full-layout-256x256x65-antialias-1 | 1 | yes | insert-primary | 0.956 | 4.274 | 0.224 | 0 | yes |
| full-layout-256x256x65-antialias-1 | 1 | yes | insert-conjugate | 0.043 | 0.077 | 0.555 | 0 | yes |
| full-layout-256x256x65-antialias-1 | 1 | yes | insert-complete | 0.804 | 4.482 | 0.179 | 0 | yes |
| full-layout-256x256x65-antialias-1 | 1 | yes | forward-complete | 3.572 | 3.800 | 0.940 | 0 | yes |
| full-layout-256x256x65-antialias-1 | 1 | yes | inverse-complete | 3.098 | 6.643 | 0.466 | 0 | yes |
| full-layout-512x512x129-antialias-0 | 0 | yes | extract | 24.883 | 88.623 | 0.281 | 0 | yes |
| full-layout-512x512x129-antialias-0 | 0 | yes | insert-primary | 20.468 | 76.126 | 0.269 | 0 | yes |
| full-layout-512x512x129-antialias-0 | 0 | yes | insert-conjugate | 0.225 | 0.316 | 0.710 | 0 | yes |
| full-layout-512x512x129-antialias-0 | 0 | yes | insert-complete | 21.498 | 76.521 | 0.281 | 0 | yes |
| full-layout-512x512x129-antialias-0 | 0 | yes | forward-complete | 47.662 | 95.419 | 0.499 | 0 | yes |
| full-layout-512x512x129-antialias-0 | 0 | yes | inverse-complete | 38.733 | 95.180 | 0.407 | 0 | yes |
| full-layout-512x512x129-antialias-1 | 1 | yes | extract | 8.865 | 29.761 | 0.298 | 0 | yes |
| full-layout-512x512x129-antialias-1 | 1 | yes | insert-primary | 8.013 | 26.782 | 0.299 | 0 | yes |
| full-layout-512x512x129-antialias-1 | 1 | yes | insert-conjugate | 0.187 | 0.227 | 0.822 | 0 | yes |
| full-layout-512x512x129-antialias-1 | 1 | yes | insert-complete | 8.076 | 27.066 | 0.298 | 0 | yes |
| full-layout-512x512x129-antialias-1 | 1 | yes | forward-complete | 31.209 | 52.374 | 0.596 | 0 | yes |
| full-layout-512x512x129-antialias-1 | 1 | yes | inverse-complete | 25.485 | 46.071 | 0.553 | 0 | yes |

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
