# Analytical Validation Summary

This note summarizes the V2 analytical-solution validation tests in
`UnitTestsV2`. The tests now exercise the canonical EVP stack:
`IMEigenvalueProblem`, `IMBoundaryCondition`, `IMInternalModes`,
`IMAnalyticalSolution`, `IMAnalyticalInternalModesBasis`, and
`IMAnalyticalGeostrophicZeroAPVModesBasis`.

## Current Coverage

| Test suite | Analytical solution family | Main purpose |
| --- | --- | --- |
| `IMConstantStratificationValidationTests.m` | `IMConstantStratificationSolution` | Constant-stratification spectra, free-surface branch construction, hydrostatic `F` zero mode, exact zero-APV construction, and spectral coordinate smoke checks. |
| `IMExponentialStratificationValidationTests.m` | `IMExponentialStratificationSolution` | Exponential Bessel-branch construction, free-surface branch construction, hydrostatic `F` null branch, exact zero-APV construction, and unsupported-boundary rejection. |
| `IMAnalyticalGeostrophicZeroAPVModesTests.m` | both concrete families | Independent constant and exponential formulas, canonical endpoint responses, continuous quadratic forms, rotations, and numerical convergence. |

## Solver Validation

The numerical solver surface is `IMSolverSpectral` with
`coordinateKind="z"`, `"wkb"`, or `"density"`, plus
`IMSolverFiniteDifference`. The canonical tests compare assembly rows,
endpoint rows, mode selection diagnostics, and basis-set evaluation for these
solver paths.

## Unsupported Cases

| Area | Case | Expected result |
| --- | --- | --- |
| Constant analytical solution | unsupported boundary condition | `IMConstantStratificationSolution:UnsupportedBoundary` |
| Exponential analytical solution | nonzero surface location or unsupported boundary condition | `IMExponentialStratificationSolution:UnsupportedBoundary` or `IMExponentialStratificationSolution:UnsupportedDomain` |
| Canonical mode selection | assessed positive metric and nonnegative form with raw negative discrete values | negative candidates are ignored during selection |
