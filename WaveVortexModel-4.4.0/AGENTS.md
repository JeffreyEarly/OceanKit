# WaveVortexModel agent guidance

Follow the shared OceanKit workspace instructions when they are available. Do not change released package snapshots or unrelated worktrees.

Before C++ profiling, optimization, benchmark design, or optimization qualification, read [the C++ optimization workflow](Benchmarks/CPP-OPTIMIZATION.md). Apply its small-team, early-verification, and source-freeze rules to the task; do not repeat broad qualification by habit.

MATLAB scientific behavior must remain backward compatible unless the user explicitly authorizes a particular change. The v4 C++ implementation has no equivalent backward-compatibility requirement, but numerical algorithms, tolerances, and checkpoint contracts must remain unchanged unless the task authorizes changing them.
