# AGENTS.md

This file gives shared instructions for AI-assisted edits in repositories under
the `OceanKitRepositories/` workspace root.

Repository-local AGENTS.md files may add more specific instructions. When a
more specific AGENTS.md exists in a descendant repository, follow the more
local file for repo-specific rules and this file for shared OceanKit rules.

## Path conventions

In this file, `OceanKitRepositories/` is the workspace root. Paths are relative
to `OceanKitRepositories/` unless explicitly marked otherwise. Do not rely on a
machine-specific parent path.

## Instruction precedence

System, developer, and tool instructions remain higher priority than repository
guidance. Within repository guidance, follow instructions in this order:

1. The most local applicable `AGENTS.md`.
2. This shared OceanKit `AGENTS.md`.
3. The routed focused agent guides.
4. General coding conventions.

The user's task request controls the work goal unless it conflicts with
higher-priority instructions or safety constraints.

When a repository instruction conflicts with a generic MATLAB habit, follow the
OceanKit guides.

## Guide routing

Before making changes, identify the kind of work requested and read the matching
guides. Read only the guides relevant to the current task. If more than one
trigger applies, read all matching guides.

### Core OceanKit guides

- Read `OceanKit/Documentation/WebsiteDocumentation/developers-guide/package-design.md`
  when changing package structure, public APIs, dependencies, package boundaries,
  examples that define package behavior, or exported package behavior.

- Read `OceanKit/Documentation/WebsiteDocumentation/developers-guide/matlab-style-guide.md`
  when editing MATLAB source, examples, tests, or MATLAB-facing APIs.

- Read `OceanKit/Documentation/WebsiteDocumentation/developers-guide/caching-style-guide.md`
  when designing or editing nontrivial cached derived state.

- Read `OceanKit/Documentation/WebsiteDocumentation/developers-guide/documentation-style-guide.md`
  when editing developer-facing, user-facing, website, package, or API documentation.

- Read `OceanKit/Documentation/WebsiteDocumentation/developers-guide/tutorial-style-guide.md`
  when creating or editing tutorials.

- Read `OceanKit/Documentation/WebsiteDocumentation/developers-guide/annotated-persistence-style-guide.md`
  when editing classes that persist state through CAAnnotatedClass and NetCDF-backed files.

### Focused agent guides

Read the focused guide before doing work in these areas:

- `OceanKit/tools/agent-guides/matlab-docs.md` when editing inline MATLAB API documentation.
- `OceanKit/tools/agent-guides/matlab-style.md` when editing MATLAB source, examples, or tests.
- `OceanKit/tools/agent-guides/package-release.md` when changing package dependencies, package metadata, changelogs, release compatibility, or exported package behavior.
- `OceanKit/tools/agent-guides/profiling.md` when profiling or changing profiling utilities.
- `OceanKit/tools/agent-guides/latex.md` when editing LaTeX, BibTeX, manuscripts, papers, or mathematical notes.

If a referenced guide is missing from the current checkout, say so explicitly and
continue with the parts of the task that can be completed safely.

## Scope discipline

- Preserve existing public behavior unless the requested task explicitly changes it.
- Do not refactor stable code purely for style.
- Keep example-script modernization minimal unless broader restructuring is explicitly requested.
- Do not update website content unless the task explicitly includes website work.

## Package snapshot discipline

Do not edit versioned package snapshots inside `OceanKit/`, such as
`OceanKit/WaveVortexModel-4.0.6`, `OceanKit/NetCDF-1.0.2`, or
`OceanKit/InternalModes-1.0.1`. Those folders are released MPM package
payloads, not authoring repositories.

When a requested change affects an OceanKit package, make the code change in the
primary package authoring repository, such as wave-vortex-model, netcdf, or
internal-modes. Only a tagged release/export process should update the
versioned package snapshots under `OceanKit/`.

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

- Applicable focused agent guides were read and followed.
- Accurate existing documentation was preserved unless it needed correction.
- Scientific notation with real meaning was preserved.
- Relevant dependency, package metadata, release, and snapshot rules were followed.
- Missing local assets, skipped verification, or residual risks are reported explicitly in the final response.
