# AGENTS.md

This file gives shared instructions for AI-assisted edits in repositories under
`/Users/jearly/Documents/OceanKitRepositories`.

Repository-local `AGENTS.md` files may add more specific instructions. When a
more specific `AGENTS.md` exists in a descendant repository, follow the more
local file for repo-specific rules and this file for shared OceanKit rules.

## Required guides

Before making changes, reread the relevant OceanKit guides and follow them
literally rather than from memory:

- `/Users/jearly/Documents/OceanKitRepositories/OceanKit/Documentation/WebsiteDocumentation/developers-guide/package-design.md`
- `/Users/jearly/Documents/OceanKitRepositories/OceanKit/Documentation/WebsiteDocumentation/developers-guide/matlab-style-guide.md`
- `/Users/jearly/Documents/OceanKitRepositories/OceanKit/Documentation/WebsiteDocumentation/developers-guide/caching-style-guide.md` when designing or editing nontrivial cached derived state
- `/Users/jearly/Documents/OceanKitRepositories/OceanKit/Documentation/WebsiteDocumentation/developers-guide/documentation-style-guide.md`
- `/Users/jearly/Documents/OceanKitRepositories/OceanKit/Documentation/WebsiteDocumentation/developers-guide/tutorial-style-guide.md` when creating or editing tutorials
- `/Users/jearly/Documents/OceanKitRepositories/OceanKit/Documentation/WebsiteDocumentation/developers-guide/annotated-persistence-style-guide.md` when editing classes that persist state through `CAAnnotatedClass` and NetCDF-backed files

When a repository instruction conflicts with a generic MATLAB habit, follow the
OceanKit guides.

## Documentation compliance

When editing inline MATLAB API documentation, use only the structured
`class-docs` tokens supported by the documentation style guide. Use the exact
spellings below, including the leading `- `, and put each token on its own
comment line:

- `- Topic: ...`
- `- Declaration: ...`
- `- Parameter <name>: ...`
- `- Returns <name>: ...`
- `- nav_order: ...`
- `- Developer: true`

Do not replace these with prose headings or unsupported pseudo-tokens such as:

- `Topic:`
- `Declaration:`
- `Parameters:`
- `Returns:`
- `Notes:`
- `See also:`
- `Examples:`

For classes, methods, and properties intended for generated API docs, keep the
source comment order consistent with the guide:

1. one-line summary
2. overview or discussion prose
3. optional example
4. structured token lines

Mathematical APIs should include the governing equations or defining relations
when they materially clarify behavior. Use `$$...$$` for rendered equations and
rendered mathematical aliases.

When updating existing documentation:

- preserve accurate existing documentation
- only change documentation that is incorrect, outdated, or directly affected by
  the code change

Before finishing any documentation edit, do a compliance pass on all touched
files and verify that no edited API doc uses unsupported headings where the
guide requires structured tokens.

## MATLAB style rules

Apply the OceanKit MATLAB style guide during edits. In particular:

- classes use `UpperCamelCase`
- methods and properties use `lowerCamelCase`
- preserve standard scientific notation when it carries real mathematical
  meaning, such as `sigma_n`, `sigma_s`, `Mxx`, `Myy`, `Mxy`, `Lx`, `Ly`, and
  `N2`
- do not rename scientifically standard quantities just to satisfy a generic
  naming rule
- use `@ClassName` folders for class-based APIs
- keep implementation-only methods in separate files inside the same class
  folder when that materially improves readability
- validate inputs at public API boundaries
- use structured, actionable error identifiers and messages

## MATLAB formatting rules

- Keep MATLAB function calls on one line whenever the full call fits within the
  repository line-length limit and each argument is a simple literal, variable,
  or short name-value pair.
- Do not split function calls across multiple lines purely for visual
  alignment.
- Only wrap a function call when the one-line form would exceed the line-length
  limit or when one or more arguments are long expressions, nested calls,
  anonymous functions, or inline comments.
- Do not introduce manual `...` continuations for simple calls that fit on one
  line.

## Scope discipline

- Preserve existing public behavior unless the requested task explicitly changes
  it.
- Do not refactor stable code purely for style.
- Keep example-script modernization minimal unless broader restructuring is
  explicitly requested.
- Do not update website content unless the task explicitly includes website
  work.

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

## Missing local assets

- If a requested task depends on local data files, example assets, generated
  outputs, scripts, or other workspace content that is not present, say so
  explicitly as soon as you discover it.
- Do not let missing local assets remain implicit in a failing command, a vague
  verification note, or a later summary.
- When practical, continue with the parts of the task that can still be done,
  but clearly separate completed work from work blocked by the missing asset.
- In the final response, explicitly name any missing local asset that prevented
  full verification or completion, and give the exact path that was expected
  when that is known.

## Pre-finish checklist

Before finishing, verify the touched files against this checklist:

- Documentation tokens follow the exact `class-docs` syntax.
- Accurate existing documentation was preserved unless it needed correction.
- Simple MATLAB calls that fit on one line are on one line.
- No unnecessary manual `...` continuations were introduced.
- Scientific notation with real meaning was preserved.
- Any class refactor still follows the `@ClassName` folder pattern.
