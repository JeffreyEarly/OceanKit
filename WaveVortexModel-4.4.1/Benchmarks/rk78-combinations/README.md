# RK78 ordered combination experiment

Issue [#493](https://github.com/JeffreyEarly/wave-vortex-model/issues/493) reduces repeated state traversals in stages 3–13 and the accepted RK78 candidate. The private helper first stores the initial weighted derivative, then fuses the remaining ordered additions and the final affine combination. For these twelve combinations, 77 full-state traversals become 24. This retains the first-product rounding boundary required by GCC and does not change stage storage, RHS calls, controller logic, or dense extension.

`CombinationScreen.cpp` compares the actual helper with separate noinline legacy passes. It checks exact real/complex output at N=2,4,8,9, then times alternating roles with two excluded warmups and nine measured samples. Run separate fresh processes for representative coefficient counts; output is JSONL. Use the normal release compiler options without fast-math. For example, from the repository root:

```sh
c++ -std=c++17 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror -ICompiledKernel/include -IPortableRuntime/src Benchmarks/rk78-combinations/CombinationScreen.cpp -o /tmp/wvm-rk78-combination
/tmp/wvm-rk78-combination 617706
```

The benchmark requires GCC or Clang for its noinline attribute and opaque timing barrier. Local GCC additionally needs the qualified Xcode SDK. The production helper is portable C++17.

Full-model acceptance uses the existing numerical tolerances and exact integration decisions. Native baseline repeats are not bitwise deterministic; do not substitute an impossible exact error-estimate gate for isolated exact arithmetic checks. Evidence and final qualification are recorded in [the verification ledger](../../.github/planning/issue-493-verification.md).

See [the completed results and timing exclusions](RESULTS.md) for the adoption decision.
