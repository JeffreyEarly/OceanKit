---
layout: default
title: Release workflow
parent: Developers guide
nav_order: 7
---

# Release workflow

OceanKit packages are authored in independent Git repositories and distributed as versioned snapshots from the central OceanKit repository. This page documents the release workflow in use as of July 2026 and the target design for its incremental modernization. The [release workflow refactor milestones](release-workflow-refactor-milestones) turn the target design into a staged migration plan.

## Current workflow

### Repository roles

Each package authoring repository owns its MATLAB source, tests, package manifest, changelog, documentation source, generated website, and a small GitHub Actions entry point. Fourteen authoring repositories currently call the shared workflow, including `spline-core`, `netcdf`, `internal-modes`, and `wave-vortex-model`; ten of those repositories provide the conventional `tools/build_website_documentation.m` hook.

The OceanKit repository has two release responsibilities:

- `.github/workflows/reusable-mpm-release.yml` implements the shared GitHub Actions job.
- `tools/ci_release.m` and `tools/update_changelog.m` implement the shared MATLAB release operations.
- Versioned directories such as `SplineCore-2.2.0` and `WaveVortexModel-4.1.1` are released MPM package snapshots.

Versioned package directories in OceanKit are distribution artifacts. Package development belongs in the corresponding authoring repository.

### Trigger and inputs

Every participating authoring repository contains `.github/workflows/release-mpm.yml`. This workflow is manually dispatched from that repository and calls `JeffreyEarly/OceanKit/.github/workflows/reusable-mpm-release.yml@main`.

The legacy caller exposes four inputs:

- `bump`: `none`, `patch`, `minor`, or `major`
- `notes`: release notes used for the changelog and GitHub release
- `shouldBuildWebsiteDocumentation`: string-valued `true` or `false`
- `shouldPackageForDistribution`: string-valued `true` or `false`

The WaveVortexModel pilot adds four optional reusable-workflow inputs without changing those callers:

- `documentationPackageSpecifier`: an exact authoring-tool package such as `ClassDocumentation@1.3.0`
- `shouldCheckWebsiteDocumentation`: whether to run a documentation check before mutation
- `documentationCheckTask`: the build-tool task, defaulting to `docs:check`
- `shouldPromoteUnreleased`: whether to promote the authoritative `Unreleased` changelog section

Supplying a documentation package or enabling either pilot boolean activates the strict pilot path. Callers that omit them retain the legacy release behavior.

The package repository supplies its own `GITHUB_TOKEN` for commits, tags, and releases in that repository. It also supplies `MPM_REPO_PAT`, which checks out and pushes to the OceanKit repository.

### Release sequence

The reusable workflow performs these operations in one job:

1. Check out the calling repository into `my-repo` with full history.
2. Check out the OceanKit commit that defines the reusable job into `OceanKitWorkflow`, using `job.workflow_repository` and `job.workflow_sha`. Release tools and dependency snapshots come only from this immutable checkout.
3. When distribution is requested, check out current OceanKit `main` separately into writable `OceanKitPublish` using `MPM_REPO_PAT`.
4. Validate the manifest, requested mode, proposed version, tag, snapshot, documentation hooks, and initial repository state.
5. Set up MATLAB R2025b Update 2 with caching enabled.
6. Register `OceanKitWorkflow` as the MPM repository.
7. Create a `matlab.mpm.Package` for `my-repo` and install it in authoring mode.
8. For an opted-in caller, install and report the exact documentation package and run the requested documentation check before mutation.
9. Call `ci_release`, using runner-temporary files for free-form notes, metadata, release content, and pilot exports.
10. Validate the resulting authoring diff and write canonical JSON metadata for every mode.
11. Commit and push the validated authoring changes.
12. When an export was produced, copy, commit, and push the new snapshot to `OceanKitPublish`.
13. For a real version bump, create the `v<version>` tag and GitHub release last.

The OceanKit commit message records the package folder and the workflow's original calling-repository SHA. The workflow logs and summary record the exact reusable-workflow ref and SHA.

### Version and changelog handling

The central `ci_release` function reads the package through `matlab.mpm.Package`. For `major`, `minor`, and `patch` releases it constructs a new `matlab.mpm.Version` and assigns it to the package, allowing the public MPM API to rewrite `resources/mpackage.json`.

Legacy callers retain the current dispatch-note behavior: when a real bump has nonempty release notes, `update_changelog` inserts a dated `## [<version>]` entry before the existing release entries, while an empty note leaves the changelog unchanged.

An opted-in caller may instead promote exactly one nonempty `## [Unreleased]` section. The workflow creates a fresh empty section, moves the complete reviewed body under the dated release version, and passes that same Markdown to the GitHub release through a file rather than executable MATLAB text.

Free-form dispatch notes are written to runner-temporary storage and read by MATLAB. Multiline text, quotes, backticks, percent signs, and blank lines therefore remain data rather than generated source. A WaveVortexModel 4.1.1 release attempt exposed the previous interpolation defect before this transport was corrected.

### Documentation and online version history

When documentation is requested, `ci_release` calls the package repository's `build_website_documentation` function if MATLAB can find it. The function is responsible for copying hand-authored pages from `Documentation/WebsiteDocumentation`, generating class or tutorial pages as needed, and writing `docs/version-history.md` from `CHANGELOG.md`.

This arrangement makes online version history fragile:

- Changelog-to-version-history logic is duplicated across package repositories.
- A missing builder silently skips the requested documentation work.
- The central workflow does not verify that `docs/version-history.md` contains the new version.
- Repositories with older copy-based builders can retain stale generated files.
- A documentation build can succeed even when GitHub Pages publishes a different branch. InternalModes currently publishes `/docs` from `InternalModesEVP`, while the release workflow updates `main/docs`.

Most package sites use GitHub Pages legacy branch publishing from `main/docs`. After the release commit reaches `main`, GitHub Pages rebuilds the site independently of the release workflow.

WaveVortexModel provides a deterministic `buildtool docs:check` command and requires the authoring-only generator `ClassDocumentation@1.3.0`. The pilot installs that exact version, reports its resolved path, and runs the check before release mutation. Other packages may have a documentation builder without an equivalent check command, so the exact generator and verification contract remains opt-in.

### Package export

When distribution is requested, `ci_release` copies the authoring repository into `dist/<Name>-<Version>` and removes authoring-only directories:

- `.git`
- `.github`
- `docs`
- `tools`
- `Documentation`
- `OceanKit`
- the nested `dist` directory

The export retains runtime source, manifest-declared package folders, tests, examples, `README.md`, and `CHANGELOG.md` unless the package layout excludes them by other means.

The workflow now writes canonical metadata in runner-temporary storage for every mode. A non-distribution run records that no export exists, skips the writable OceanKit checkout and every snapshot-publication step, and remains a complete supported mode. The legacy text metadata remains available to direct `ci_release` callers but no longer controls the shared workflow.

### Repository-specific variation

The release caller is nearly identical across the participating repositories, but defaults differ. Some packages default to documentation or distribution being disabled.

Not every repository provides a documentation builder. Requesting documentation in those repositories currently produces no error and no generated documentation changes.

`internal-modes/tools/ci_release.m` has the same name as the central entry point. Because the package `tools` directory is added last and therefore takes path precedence, InternalModes can execute its older local implementation instead of OceanKit's current implementation. That implementation writes the manifest JSON directly rather than using the public MPM package API.

### Current failure and recovery characteristics

The workflow has successfully produced authoring commits, tags, GitHub releases, documentation updates, and OceanKit snapshots, but several boundaries remain weak:

- Existing callers still use the mutable `@main` version until they migrate to a tested workflow tag.
- The two legacy boolean inputs remain string-valued until the caller-contract migration.
- Invalid input combinations are discovered after the relatively expensive MATLAB setup.
- Documentation generation is not required to prove that the new changelog entry reached the published version-history page.
- No clean-path test installs and exercises the exported OceanKit snapshot before publication.
- Non-opted-in `bump=none` distributions can still modify an already published version folder.
- The final `rsync` does not delete files that disappeared from a replacement export.
- Authoring commits, tags, releases, and OceanKit commits are separate writes and cannot be transactional.
- Concurrent releases from different repositories can race when pushing directly to OceanKit.

Recovery is manual. The operator must determine which writes completed, repair or resume the missing steps, and avoid applying a second version bump when retrying a partially completed release.

## WaveVortexModel release-safety pilot

[OceanKit issue #3](https://github.com/JeffreyEarly/OceanKit/issues/3) provides a limited, backward-compatible pilot needed before WaveVortexModel 4.2.1. Existing callers retain their current contract unless they opt into the pilot. WaveVortexModel requests `ClassDocumentation@1.3.0`, runs `buildtool docs:check` before mutation, promotes its nonempty `Unreleased` section, and consumes the tested workflow through the immutable `mpm-release-v0.1.0` tag. The tagged reusable YAML checks out its colocated OceanKit tools at `job.workflow_sha`, so later changes to OceanKit `main` cannot alter a pinned release.

The pilot also writes metadata in non-distribution modes, conditions OceanKit publication on actual export output, rejects an existing snapshot, and creates the tag and GitHub release only after the authoring commit and OceanKit snapshot succeed. It intentionally does not migrate all callers to native booleans, rename the central release engine, generalize documentation handling across the ecosystem, or implement snapshot-replacement and push-race recovery.

Only this pilot blocks [WaveVortexModel issue #19](https://github.com/JeffreyEarly/wave-vortex-model/issues/19) and the 4.2.1 release. The remaining v1 roadmap is important shared-infrastructure work but does not block that maintenance release.

## Future workflow

### Design principles

The future workflow keeps the current repository model while hardening its boundaries:

- Package repositories retain a small local dispatch workflow so releases use the package's Actions history, permissions, secrets, and `GITHUB_TOKEN`.
- OceanKit owns all shared release behavior through a versioned reusable workflow and uniquely named MATLAB entry point.
- `CHANGELOG.md` is the source of truth for release history, and generated documentation is validated against it.
- Every requested mode either completes independently or fails before publication.
- Versioned OceanKit snapshots are immutable by default.
- Expensive setup and irreversible publication occur only after cheap validation succeeds.
- A release is tested in the same exported form that users install.

### Stable caller contract

After the pilot, each package migrates to one canonical thin caller with native inputs:

- `bump`: choice with `none`, `patch`, `minor`, and `major`
- `buildDocumentation`: boolean
- `publishPackage`: boolean
- `documentationPackageSpecifier`: optional exact authoring-tool package specifier
- `documentationCheckTask`: optional canonical package check, normally `docs:check`
- `replaceExistingSnapshot`: boolean, default `false`, exposed only while an explicit replacement path remains necessary

The caller references an immutable OceanKit workflow release tag such as `@mpm-release-v1.0.0`, not `@main`. OceanKit stores a canonical caller template and a drift check reports undocumented differences across package repositories.

`CHANGELOG.md` is the release-note source of truth. A real version bump requires a nonempty `## [Unreleased]` section, promotes that content to a dated version heading, creates a fresh empty `Unreleased` section, and uses the promoted section as the GitHub release body. Legacy dispatch notes remain a transitional pilot input only and are transported through environment variables or files rather than inserted into MATLAB source.

### Preflight

A cheap preflight runs before MATLAB setup and validates:

- the bump value and boolean input types
- the presence and readability of `resources/mpackage.json`
- a nonempty `Unreleased` section for a real bump
- the package documentation builder when documentation is requested
- the absence of the proposed tag
- the absence of an existing OceanKit snapshot unless replacement was explicitly authorized
- the expected release branch and a clean checkout

Failures identify the exact invalid condition and perform no repository writes.

### Central release engine

OceanKit exposes a uniquely named MATLAB function, such as `oceankit_release_package`, so a package-local helper cannot shadow it. Package repositories provide hooks only for package-specific operations, initially `build_website_documentation` and an optional focused smoke-test function.

The central function:

1. Reads and updates the manifest only through `matlab.mpm.Package`.
2. Promotes the authoritative `Unreleased` changelog section for every real version bump.
3. Runs the requested package-specific documentation build.
4. Regenerates `docs/version-history.md` centrally from the current changelog after the package-specific builder finishes.
5. Verifies that the generated page contains the release version and current changelog entry.
6. Exports the package only when distribution is requested.
7. Always writes release metadata, including which optional outputs were produced.
8. Validates the export before returning control to GitHub Actions.

The central version-history helper removes changelog rendering from individual documentation builders. This establishes one format and one validation point across the package ecosystem.

### Independent execution modes

The reusable workflow conditions later steps on the metadata and native boolean inputs:

| Mode | Version/changelog commit | Documentation commit | OceanKit snapshot | Tag and release |
| --- | --- | --- | --- | --- |
| Documentation only | No | Yes | No | No |
| Version only | Yes | Optional | No | Yes |
| Distribution without bump | No | Optional | New snapshot only | No |
| Full release | Yes | Optional | Yes | Yes |

Distribution without a bump is allowed only when the target snapshot does not exist. Replacing an existing snapshot requires the explicit replacement input and should be reserved for recovery from a provably incomplete publication.

### Documentation integrity

Documentation generation follows a fixed order: update changelog, run the package builder, regenerate version history centrally, then validate. The release fails if requested documentation lacks a builder or if the generated version history is missing the current release.

All participating package sites publish `main/docs`. Pages configuration is checked during migration and recorded alongside the caller-workflow audit. Generated `docs/` remains committed output; `Documentation/WebsiteDocumentation` remains the source for hand-authored pages.

### Export verification and snapshot safety

Before any release publication, a clean MATLAB path registers a temporary MPM repository containing the new export, installs the exported package with manifest-compatible dependencies, and runs a fast package smoke test when one is defined. This check complements authoring-repository tests; it verifies the artifact and dependency graph users receive.

Publishing fails when `<Name>-<Version>` already exists in OceanKit. If an authorized recovery replaces an incomplete snapshot, synchronization uses delete semantics so removed files cannot remain from the earlier export.

### Publication order and recovery

After every generated artifact passes validation, publication proceeds in this order:

1. Commit and push the version, changelog, and documentation changes to the authoring repository.
2. Copy, commit, and push the new snapshot to OceanKit when distribution was requested.
3. Pull and rebase before retrying an OceanKit push that lost a concurrent race.
4. Create the version tag and GitHub release last.

The OceanKit commit records the post-bump authoring commit SHA. Although writes to two repositories still cannot be atomic, this order avoids announcing a GitHub release before its requested OceanKit snapshot exists.

Release metadata and step summaries identify the completed writes and provide an exact recovery command or rerun mode. The previous reusable-workflow tag remains available so callers can be rolled back without reverting the shared implementation.

### Migration

The future workflow is introduced incrementally rather than replacing the working system in one step. The WaveVortexModel compatibility pilot establishes the first safe vertical slice, after which a documentation-heavy package and a simple package serve as ecosystem pilots before callers move in dependency-aware batches. The [release workflow refactor milestones](release-workflow-refactor-milestones) define the implementation order, GitHub tracking issues, acceptance criteria, and rollback boundaries.
