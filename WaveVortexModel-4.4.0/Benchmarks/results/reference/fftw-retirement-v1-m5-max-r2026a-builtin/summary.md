# FFTW integration retirement gate

- Decision: `RETIRE`
- Status: `complete`
- Baseline: `9652b116b3ffd4ee3372cc5cdeea9700cd6cbc32` (`v4.2.1`)
- Candidate: `d88abddb6409b6b5f8f77e01dfe6fbd60233f749`
- Gate: error at most `1e-12`; no candidate median over `1.03x` the same-host baseline

| Case | v4.2.1 (ms) | Candidate (ms) | Candidate / baseline | Error | Builtin executed | Pass |
|---|---:|---:|---:|---:|---|---|
| constant-nonhydrostatic-256x256x65 | 120.055 | 116.490 | 0.970 | 0 | yes | yes |
| constant-hydrostatic-256x256x65 | 100.681 | 102.184 | 1.015 | 0 | yes | yes |
| constant-nonhydrostatic-512x512x129 | 935.649 | 901.795 | 0.964 | 0 | yes | yes |
| constant-hydrostatic-512x512x129 | 703.678 | 724.506 | 1.030 | 0 | yes | yes |

## Candidate builtin storage and RSS

| Case | Known persistent (MiB) | Maximum known live (MiB) | Persistent RSS increment (MiB) | Peak RSS increment (MiB) |
|---|---:|---:|---:|---:|
| constant-nonhydrostatic-256x256x65 | 65.260 | 130.260 | 1233.406 | 1233.938 |
| constant-hydrostatic-256x256x65 | 65.260 | 130.260 | 1177.703 | 1178.375 |
| constant-nonhydrostatic-512x512x129 | 517.037 | 1033.037 | 8456.453 | 8457.000 |
| constant-hydrostatic-512x512x129 | 517.037 | 1033.037 | 8017.078 | 8017.812 |
