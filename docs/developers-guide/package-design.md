---
layout: default
title: Package design
parent: Developers guide
nav_order: 1
---

# Package design

This document defines the standard authoring-repository pattern for OceanKit packages. The goal is not to make every repository identical; the goal is to make them structurally predictable enough that package authors, release tooling, and AI-assisted edits can move between repositories without relearning the basics each time.

## Required

The following pieces should exist in every OceanKit package authoring repository.

- MATLAB source code at the repository root. Use `@ClassName` folders for class-based APIs, plain `.m` files for standalone functions, and domain folders such as `Operations`, `Integrators`, or `Forcing` when the package needs them.
- `resources/mpackage.json`. This is the MPM contract for the package. Keep the package name, version, id, display name, summary, description, provider, dependencies, `releaseCompatibility`, and `schemaVersion` accurate. Keep the `folders` entries aligned with the real repository layout.
- `README.md`. Provide a concise package summary, a quick start, and any citations or scientific references that a user needs immediately.
- `.github/workflows/release-mpm.yml`. OceanKit packages use the shared reusable release workflow from the `OceanKit` repository. Each authoring repository should expose the same release entry point.

## Recommended

The following pieces should exist in most OceanKit packages even when the package is relatively small.

- `CHANGELOG.md`. Keep release notes close to the source repository so CI can publish both version bumps and user-facing change history.
- `UnitTests/`. Put automated tests in a dedicated folder and keep them runnable without manual setup beyond normal package dependencies.
- `Documentation/README.md`. Explain the source/build relationship for documentation in the same way the other package repositories do.
- `Documentation/WebsiteDocumentation/`. Treat this tree as the canonical source for hand-authored website pages.
- `tools/build_website_documentation.m`. The repository should be able to regenerate its `docs/` output locally and in CI.
- `docs/`. Commit generated GitHub Pages output so the repository can publish the documentation site directly from versioned markdown.

## Optional

These pieces are common and useful, but they are package-dependent rather than universal.

- `Examples/` for user-facing scripts and tutorials
- `figures/` for README or documentation assets
- `Extras/` for analytical notes, experimental tools, or secondary assets
- domain-specific folders such as `FastTransforms`, `ObservingSystems`, or `FlowComponents`
- additional scripts in `tools/` for local authoring, release preparation, or package creation

Optional folders should still have a clear reason to exist. Avoid adding a new top-level folder when the contents belong naturally in an existing domain folder or in the package root.

## OceanKit vs Authoring Repositories

OceanKit itself is the core MPM repository. It stores released package snapshots such as `SplineCore-2.0.0` or `WaveVortexModel-4.0.2` for distribution.

Authoring happens in the package repositories themselves, for example:

- `spline-core`
- `wave-vortex-model`
- `netcdf`

That split is intentional:

- the authoring repository owns source code, tests, docs source, and release automation
- the `OceanKit` repository receives exported release snapshots for MPM installation

Do not treat the package snapshot inside `OceanKit` as the canonical place to edit a package.

## Release and Export Behavior

The current `OceanKit/tools/ci_release.m` export path copies a package root into `dist/<Name>-<Version>` and explicitly excludes repository-level authoring assets such as:

- `.git`
- `.github`
- `docs`
- `tools`
- `Documentation`
- `OceanKit`

That means documentation source, generated GitHub Pages output, release workflow files, and authoring scripts are part of the source repository, but they are not part of the packaged MPM payload.

Design repositories with that split in mind:

- package code and packaged assets belong in the package root
- documentation and release tooling belong at repo level for authors
- generated `docs/` should exist for publishing, but not because runtime code depends on it

## A Good Default Layout

For a typical OceanKit package, a good starting point is:

```text
package-root/
  @ClassName/
  DomainFolder/
  UnitTests/
  Documentation/
  docs/
  resources/
  tools/
  .github/workflows/
  README.md
  CHANGELOG.md
```

Not every repository needs every folder, but new packages should start close to this pattern unless there is a strong reason not to.
