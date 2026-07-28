---
layout: default
title: Release workflow refactor milestones
parent: Developers guide
nav_order: 8
---

# Release workflow refactor milestones

This roadmap migrates OceanKit package releases from the [current workflow to the future workflow](release-workflow) without replacing the repository model or requiring a single disruptive cutover. Each milestone must be independently releasable, testable, and reversible.

## Objectives and boundaries

The migration must:

- preserve one thin, manually dispatched workflow in each package repository
- centralize shared release behavior in OceanKit
- make changelog updates reliably appear in published version-history pages
- make documentation-only, version-only, distribution-only, and full-release modes behave as declared
- validate the exported package before publication
- make published snapshots immutable by default
- provide an explicit recovery path for partial failures

The migration does not introduce a central cross-repository dispatcher, GitHub App authentication, automatic semantic-version selection, automatic releases on every merge, or a new documentation framework.

## Milestone status

| Milestone | Status | Prerequisites | Completion evidence |
| --- | --- | --- | --- |
| 1. Release contract and preflight | Not started | Current workflow documented | Input and mode tests pass without starting MATLAB on invalid requests |
| 2. Changelog and documentation integrity | Not started | Milestone 1 metadata contract | Pilot version history contains the exact new changelog entry |
| 3. Central release engine | Not started | Milestones 1–2 helpers | InternalModes and a second pilot resolve the same central entry point |
| 4. Publication safety | Not started | Stable central engine | Existing snapshots are protected and partial-failure recovery is demonstrated |
| 5. Verification and ecosystem rollout | Not started | Milestones 1–4 | All callers use a tested immutable workflow tag |

Update this table as milestones begin and complete. Record links to the validating workflow runs or pull requests in the completion-evidence column.

## Milestone 1: Release contract and preflight

### Deliverables

- Change reusable-workflow documentation and package callers to use native boolean inputs for documentation and distribution.
- Use the `inputs` context rather than the string-converting `github.event.inputs` context.
- Transport release notes through an environment variable or temporary file and read them with `getenv` or `fileread` in MATLAB.
- Validate `bump` against `none`, `patch`, `minor`, and `major`.
- Require nonempty notes when `bump` is not `none`.
- Verify the manifest, requested documentation hook, release branch, proposed tag, and proposed snapshot before MATLAB setup.
- Make the MATLAB release function write metadata for every mode, not only distribution.
- Condition export-copy and OceanKit-push steps on the distribution input and produced metadata.
- Emit a GitHub Actions step summary containing the selected mode, current version, proposed version, and planned writes.

### Acceptance criteria

- Multiline notes containing quotes, backticks, percent signs, and blank lines reach MATLAB unchanged.
- Invalid input, a missing manifest, missing requested documentation builder, or an existing tag fails during preflight.
- Documentation-only mode reaches the authoring commit step without requiring an export or OceanKit metadata folder.
- Distribution-disabled mode never modifies the OceanKit checkout.
- Existing full-release behavior remains available behind the same manual entry point.

### Rollback

Keep existing callers on `@main` during development. Reverting the central workflow commit restores the previous contract until the first versioned workflow tag is created.

## Milestone 2: Changelog and documentation integrity

### Deliverables

- Update the changelog on every real version bump; remove the current dependency on nonempty notes after Milestone 1 makes notes mandatory.
- Add one central helper that renders `CHANGELOG.md` into `docs/version-history.md` with standard Just the Docs front matter.
- Run the helper after the package-specific documentation builder so later copy or generation operations cannot overwrite it.
- Validate that the generated page contains the new `## [<version>]` heading and the submitted release-note text.
- Make a requested but unavailable documentation builder a hard failure.
- Remove duplicated changelog rendering from package builders only after the central helper is active for that package.
- Audit GitHub Pages source settings and change InternalModes from `InternalModesEVP/docs` to `main/docs`.

### Pilot

Use `spline-core` as the documentation-heavy pilot because it combines generated tutorials, class documentation, cached tutorial assets, changelog rendering, and Pages publishing. Use `netcdf` as the simpler documentation pilot because its builder primarily copies source pages and generates class reference pages.

### Acceptance criteria

- A pilot patch release adds one changelog entry and the identical entry appears in committed `docs/version-history.md`.
- The deployed Pages site displays that entry after the release commit.
- A documentation-only run republishes a manually corrected changelog without creating a tag or OceanKit snapshot.
- Deleting or renaming a requested builder fails with an actionable preflight error.
- Rebuilding documentation twice from unchanged inputs produces no tracked diff.

### Rollback

Retain the package-local changelog rendering until each pilot passes. If central rendering fails, disable it for the pilot and restore the package builder's previous rendering block without changing release version semantics.

## Milestone 3: Central release engine

### Deliverables

- Rename the shared MATLAB entry point to a unique name such as `oceankit_release_package`.
- Keep version calculation, MPM manifest mutation, changelog update, version-history synchronization, export, metadata writing, and validation in that central entry point or its OceanKit-owned helpers.
- Reserve package `tools` functions for package-specific documentation and smoke-test hooks.
- Remove or rename `internal-modes/tools/ci_release.m` and its duplicate changelog helper after confirming no local caller uses them.
- Ensure manifest changes occur only through `matlab.mpm.Package`.
- Store a canonical thin caller template in OceanKit and add a read-only drift checker that allows only documented per-package defaults.

### Acceptance criteria

- `which oceankit_release_package` resolves to the OceanKit checkout for every pilot.
- No package-local function can shadow the central entry point.
- InternalModes version mutation preserves a valid MPM manifest written by the public API.
- The caller drift checker reports all participating repositories and explains every difference.
- The central engine produces equivalent metadata for documentation-only, version-only, distribution-only, and full-release dry runs.

### Rollback

Keep the old central `ci_release` available as a compatibility wrapper for one migration cycle. The wrapper delegates to the new entry point and is removed only after every caller has migrated.

## Milestone 4: Publication safety

### Deliverables

- Fail by default when OceanKit already contains `<Name>-<Version>`.
- Add an explicit `replaceExistingSnapshot` recovery input, defaulting to `false`.
- Require a replacement reason in the workflow summary and use `rsync --delete` when replacement is authorized.
- Record the post-bump authoring commit SHA in the OceanKit commit.
- Publish the authoring commit first, the OceanKit snapshot second, and the tag and GitHub release last.
- On a non-fast-forward OceanKit push, fetch, rebase the snapshot commit, verify that the target folder still does not conflict, and retry a bounded number of times.
- Report completed writes and the appropriate recovery mode when a later publication step fails.

### Acceptance criteria

- An attempted ordinary replacement of an existing snapshot fails before either repository is modified.
- An authorized replacement removes a file that existed only in the incomplete snapshot.
- Two simulated releases of different packages can resolve an OceanKit push race without losing either snapshot.
- A simulated OceanKit push failure produces an authoring commit but no tag or GitHub release.
- Retrying that partial release with `bump=none` and the recorded version completes the snapshot without producing another version bump.

### Rollback

Publication-order changes remain behind the untagged workflow during testing. If race retry or recovery reporting is unreliable, retain immutable snapshot checks and return to manual OceanKit conflict resolution before creating the first stable workflow tag.

## Milestone 5: Verification and ecosystem rollout

### Deliverables

- Create a temporary MPM repository containing the export and install the exported package on a clean MATLAB path with manifest-compatible OceanKit dependencies.
- Define a fast default smoke check that loads the package and verifies its manifest; allow a package-specific focused smoke-test hook.
- Do not run entire package test suites inside the release workflow unless a package explicitly opts in.
- Exercise the completed workflow with `spline-core` and a package without generated website documentation, such as `geographic-projection`.
- Tag the tested reusable workflow as an immutable release such as `mpm-release-v1.0.0`.
- Update package callers to the immutable tag in dependency-aware batches.
- Verify the authoring commit, changelog, generated version history, Pages deployment, tag, GitHub release, exported manifest, and OceanKit commit after each batch.
- Retain the previous reusable-workflow tag and document caller rollback as a one-line reference change.

### Rollout batches

1. Pilot `spline-core` and `geographic-projection`.
2. Migrate foundational providers: `netcdf`, `class-annotations`, `class-docs`, `distributions`, and `chebfun`.
3. Migrate dependent numerical packages: `internal-modes`, `advection-diffusion-models`, and `wave-vortex-model`.
4. Migrate cross-organization consumers: `wave-vortex-model-diagnostics` and `AlongTrackSimulator`.

Provider packages move before consumers so exported smoke tests resolve against already migrated dependency snapshots.

### Acceptance criteria

- The exported package installs from a temporary repository rather than from sibling authoring checkouts.
- Pilot release and documentation-only runs complete with no manual repository repair.
- Every participating caller references the same immutable workflow tag unless an exception is recorded.
- All requested Pages sites publish from `main/docs` and show their latest released changelog entry.
- The canonical caller drift check is clean after the final batch.
- The previous workflow tag can be restored without modifying the shared implementation.

### Rollback

Roll back one package by changing its caller to the previous reusable-workflow tag. Do not move or rewrite an existing workflow release tag.

## Completion definition

The refactor is complete when all callers use the tested immutable workflow tag, all supported execution modes behave independently, changelog entries are validated in generated and deployed version-history pages, exported snapshots pass a clean-path smoke test, existing snapshots are protected by default, and release recovery no longer requires inferring completed writes from raw logs.

