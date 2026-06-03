# Profiling Agent Guide

## Profiling tools

- The generic profiling hotspot utilities live in `tools/profiling`. Add that
  folder to the MATLAB path before using them.
- Prefer `profileCodeHotspots` for new hotspot-finding work or for replacing
  one-off profiling scripts that need total time, self time, call counts, and
  line hotspot summaries.
- Set `projectRoots` to the repository or package root you are actually trying
  to optimize.
- Start with `priorityTargets`, `topProjectBySelfTime`,
  `topProjectByNumCalls`, and `topActionableLines`.
- Treat raw `topBy*` tables and `topLines` as global context only. They may
  include the profiling wrapper itself, external dependencies, toolbox code,
  builtins, and project call-site lines whose time is really spent in callees.
- Use `topCallsiteLines` to separate wrapper or fan-out lines from true
  compute-line hotspots.
- Use `compareProfileHotspots` only for same-workload before/after comparisons.
  Do not use it to compare different benchmark shapes and then interpret the
  result as a performance regression.
- For comparisons, start with elapsed wall time and function-level rows in
  `functionDiffs`, `regressions`, and `improvements`.
- Treat `lineDiffs` and line-level comparison targets as secondary context.
  Helper extraction, code motion, or line-number churn can make line-level
  "regressions" misleading even when the workload got faster.
- When profiling a package authoring repository, add the package and dependency
  folders from `resources/mpackage.json` before profiling so MATLAB resolves
  package code in the same layout the package expects.
