# Profiling Tools

Use these helpers when you want MATLAB profiler output turned into a short list
of hotspots you can actually work on.

## Quick start

Profile a function handle directly:

```matlab
addpath("path/to/OceanKit/tools/profiling")

analysis = profileCodeHotspots(@() myWorkload(), projectRoots=pwd);
```

Start with:

- `analysis.priorityTargets`
- `analysis.topProjectBySelfTime`
- `analysis.topProjectByNumCalls`
- `analysis.topActionableLines`

Use the raw tables such as `topBySelfTime` and `topLines` as context only.
They can include wrapper code, dependencies, toolbox functions, builtins, and
call-site lines whose time is really spent in callees.

## Analyze existing profiler output

If you already ran the profiler yourself, analyze the saved output directly:

```matlab
profile on
myWorkload()
profile off

info = profile("info");
analysis = analyzeProfileInfo(info, projectRoots=pwd);
```

## Compare two runs

Use `compareProfileHotspots` for same-workload before/after comparisons:

```matlab
before = profileCodeHotspots(@() oldWorkload(), projectRoots=pwd, shouldPrintReport=false);
after = profileCodeHotspots(@() newWorkload(), projectRoots=pwd, shouldPrintReport=false);

comparison = compareProfileHotspots(before, after, projectRoots=pwd);
```

`comparison.functionDiffs` and `comparison.lineDiffs` are the detailed joined
tables. `comparison.regressions` and `comparison.improvements` are the short
headline views.

## Compute lines vs call sites

- `topActionableLines` highlights lines that look like real compute kernels.
- `topCallsiteLines` highlights lines that mostly hand time off to callees.

If a constructor call, helper call, or wrapper line dominates `topLines`, check
`topCallsiteLines` before deciding that the line itself needs rewriting.

## Package authoring setup

In an authoring repository, add the package and dependency folders from
`resources/mpackage.json` before profiling. A typical pattern is:

```matlab
root = "/path/to/OceanKitRepositories";
addpath(fullfile(root, "OceanKit", "tools", "profiling"))

pkgRoot = fullfile(root, "my-package");
metadata = jsondecode(fileread(fullfile(pkgRoot, "resources", "mpackage.json")));
for iFolder = 1:numel(metadata.folders)
    addpath(fullfile(pkgRoot, metadata.folders(iFolder).path));
end
addpath(pkgRoot)
```

Do the same for dependencies in dependency order when the workload uses them.
