# Compiled backend preview

- Status: `complete`
- Decision: `PREVIEW-AVAILABLE`
- Provider: native FFTW 3.3.11 NEON/pthreads

| Case | MATLAB (ms) | Compiled (ms) | Speedup | Error | Exact retained ratio | Operation RSS ratio |
|---|---:|---:|---:|---:|---:|---:|
| constant-hydrostatic-256x256x65 | 96.716 | 57.723 | 1.676x | 8.855e-15 | 1.922 | 0.160 |
| constant-nonhydrostatic-256x256x65 | 113.332 | 66.681 | 1.700x | 8.034e-15 | 1.922 | 0.122 |
| constant-hydrostatic-512x512x129 | 692.773 | 392.815 | 1.764x | 2.298e-14 | 1.916 | 0.014 |
| constant-nonhydrostatic-512x512x129 | 868.107 | 526.656 | 1.648x | 1.568e-14 | 1.916 | 0.074 |

Memory is descriptive and does not gate preview availability. Exact retained bytes count reachable application-owned arrays; allocator, MATLAB FFT, and FFTW plan-owned memory remain opaque and are represented by isolated operation RSS.
