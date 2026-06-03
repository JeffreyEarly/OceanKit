# Package And Release Agent Guide

## Dependency and release discipline

Treat `resources/mpackage.json` as source code. When an OceanKit package directly
calls a class, function, method, or static helper from another OceanKit package,
the caller must declare that package as a direct dependency. Do not rely on
transitive dependencies to make runtime symbols available.

When a change uses a symbol or behavior that was added after the currently
declared dependency floor, raise `compatibleVersions` to the first released
package version that contains that symbol or behavior. Verify that version from
tags, manifests, or OceanKit package snapshots rather than from memory.

Before releasing a package, test at least one clean MATLAB path using OceanKit
package snapshots that satisfy the manifest, preferably the lowest allowed
versions. Authoring-repo tests against sibling `main` checkouts are useful, but
they do not prove the released install graph works.

For multi-package dependency fixes, release and export provider packages first,
then bump dependent package manifests, then export dependents. After exporting,
run the consumer's focused tests from the exported OceanKit snapshots.

Remove manifest dependencies that the package does not actually use. False
dependencies can create install-order problems and obscure the real ownership of
runtime requirements.

When dependency metadata changes, add a changelog note explaining the dependency
floor or dependency removal.

## Package metadata edits

- Never directly edit `mpackage.json` files. When package metadata must change,
  work from the package root through the public `matlab.mpm.Package` API and
  let MATLAB rewrite `resources/mpackage.json`. Inspect the available
  `matlab.mpm.Package` methods and properties to determine the supported
  operation instead of patching JSON by hand.
