# RK78 combination results (#493)

Adopt the two-pass implementation. On the 256x256x28 EddyTide workload, retained balanced measurements show **2.7% less integration time** (paired bootstrap 95% improvement interval: 1.6–3.9%). The complete eight-pair result is a 2.3% improvement. The larger Boussinesq point estimate did not trigger the 3% regression-investigation gate. Its interval permits about 2.2% integration regression, so this experiment establishes neither a Boussinesq speedup nor statistical absence of regression.

| Workload | Baseline integration | Candidate integration | Candidate/baseline | Paired bootstrap 95% | Owned peak ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| EddyTide 256x256x28, RK78 | 17.299 s | 16.839 s | 0.9734 | 0.9612–0.9843 | 1.000 |
| Boussinesq 256x256x129, RK78 | 4.422 s | 4.307 s | 0.9738 | 0.9276–1.0223 | 1.000 |

The campaign used the default `reuse` policy. Low-memory qualification was not run; the owned-peak ratios apply only to reuse. Times are geometric means. The standalone complex combinations are 17–44% faster at the representative coefficient-family sizes (205,902 and 563,635 values), under Clang and GCC, with bitwise legacy agreement. Whole-model gains are smaller because FFTs, matrix multiplication and the RHS dominate.

## Protocol and limitations

The planned campaign executed two warmup pairs and eight alternating measured pairs per workload with frozen binaries. Activity logs found unrelated CPU activity during integration in pair 7 of each workload. Exclude the complete final counterbalanced block (pairs 6 and 7), retaining six measured pairs with three runs of each order. A targeted replacement attempted afterward could not pass the idle preflight because antivirus scanning persisted. The exclusions use host activity rather than speed; the complete eight-pair results, failed preflights and all numerical comparisons remain in the archive. This is a declared reduction from the planned valid timing sample count. Boussinesq integration windows are about four seconds and original activity sampling was every ten seconds, so its host attribution and speed estimate are less precise.

The first idle preflight also rejected Photos/Spotlight work before measurement; no system service or experiment was changed. The successful campaign passed its idle preflight. Background services appeared later, chiefly during Boussinesq startup; whole-process timings are retained as descriptive evidence. Boussinesq loads about 1.5 GB of persisted matrix data and spends roughly 34 seconds in startup versus four integrating. Budget future campaigns from those costs and consider a prepared-workload benchmark before another broad search.

Scientific comparisons retain the existing relative 1e-10 and absolute 1e-12 tolerances. Actual step decisions, RHS/producer counts and owned memory match. Maximum scientific differences are 1.10e-17 (EddyTide) and 2.08e-17 (Boussinesq). The identical frozen baseline itself varies by about 1.03e-17 across fresh runs. Exact helper arithmetic remains a separate Clang/GCC gate; exact full-model normalized-error values are not a valid baseline contract.

## Implementation and verification

RK78 stages 3–13 and the accepted candidate retain the first stored weighted product, then combine remaining additions and final affine application in one traversal. Seventy-seven traversals become twenty-four for these twelve combinations. Stage storage, numerical term order, constraints, thirteen RHS evaluations per step, controller, error norms and dense extension remain unchanged. The one-pass prototype was rejected because GCC could reverse first-product contraction. No MATLAB scientific code or checkpoint format changed.

- Combined release suite: 58/58.
- Fresh ASan/UBSan affected integration/output suite: 3/3.
- Clang and GCC exact helper checks plus GCC unified integration: pass. Local GCC's existing Apple Mach-header limitation prevents building the CLI; Linux CI owns that complete compiler gate.
- Affected MATLAB/source checks: 10/10; six forward/restart configurations pass with reference and native providers.
- Strict receipt catalog: 8/8; compatibility assembly: 7/7; Code Analyzer: no findings in the three metadata-only MATLAB edits.
- New helper fingerprint is included in receipt production, validation and synthetic fixtures. The compatibility assembly was regenerated from the final test source.

Frozen runtime: `72fc2bdd2ff02bbc431ed36fcdf6dab2e5015c0c`; receipt-source revision: `66ce4599`. Runner and probe stamps/hashes were independently verified. Production source remained unchanged through qualification; subsequent changes are fingerprints, receipts and reporting. [Machine-readable results](results.json) retain both the complete and filtered timing statistics, executable identities and artifact hashes. Raw evidence is in `OceanKitRepositories/wave-vortex-model-benchmark-artifacts/rk78-combinations-20260912`.

The next production target remains the Boussinesq spectral assembly/projection pipeline. RK78 combination fusion yields a modest, measured improvement without a larger solver redesign.
