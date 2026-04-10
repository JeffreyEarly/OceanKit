---
layout: default
title: Caching style guide
parent: Developers guide
nav_order: 5
---

# Caching style guide

Caching is useful in OceanKit when a class repeatedly needs expensive derived
state from a small amount of canonical scientific state. Typical examples are
basis matrices, factorizations, piecewise-polynomial representations, compiled
constraint systems, and quadrature operators.

The goal is not to cache everything. The goal is to make expensive repeated
work fast without making object state ambiguous or invalidation fragile.

## Core principles

- Cache only derived state, never the canonical scientific state that defines
  the object.
- Keep cached state clearly separate from public scientific state.
- Prefer explicit cache ownership over generic bags of fields.
- Invalidate caches from the mutation points that actually change their
  dependencies.
- Make stale-cache bugs hard to write and easy to test.

The existing `BSpline` implementation is the model to follow. It stores
piecewise-polynomial data in dedicated implementation properties and refreshes
or clears them from narrow change hooks such as
`splineCoefficientsDidChange()` and `tKnotDidChange()`.

## What should be cached

Cache a value when all of the following are true:

- it is materially more expensive to recompute than to store
- it is derived from other object state
- it is reused across multiple calls or iterations
- its dependency set is small and well understood

Good cache candidates include:

- basis matrices evaluated on a fixed grid
- Cholesky, QR, LU, or normal-equation factorizations
- compiled constraint matrices
- quadrature grids and quadrature operators
- piecewise-polynomial conversion tables

Poor cache candidates include:

- cheap scalar formulas
- values used only once in a call path
- public state that callers think of as part of the model itself
- derived values whose dependencies are broad, implicit, or hard to invalidate

## Recommended class pattern

For handle classes, prefer this pattern:

1. Store canonical model state in ordinary properties.
2. Store cached implementation state in private or protected properties.
3. Compute cached values lazily in a getter or helper method.
4. Clear or refresh the affected caches in explicit invalidation methods when
   dependencies change.

Use dedicated properties with specific names such as `basisMatrix_`,
`roughnessMatrix_`, or `ppCoefficients_` rather than a generic
`variableCache` struct.

```matlab
classdef ExampleSpline < handle
    properties (Access = private)
        knotPoints_
        coefficients_
        basisMatrix_ = []
        roughnessMatrix_ = []
    end

    methods
        function value = basisMatrix(self)
            if isempty(self.basisMatrix_)
                self.basisMatrix_ = ExampleSpline.matrixForDataPoints( ...
                    self.samplePoints_, knotPoints=self.knotPoints_);
            end
            value = self.basisMatrix_;
        end

        function knotPointsDidChange(self)
            self.basisMatrix_ = [];
            self.roughnessMatrix_ = [];
        end
    end
end
```

That structure makes the dependency graph visible in source code. When a new
cached value is added, its invalidation point is obvious and local.

## Prefer explicit invalidation over field lists

Avoid cache systems that rely on string lists such as
`kDependentVariables = {'X','Vq','L'}` and repeated calls to `rmfield(...)`.

That pattern looks flexible, but it is fragile in practice:

- every new cached quantity must be added to one or more dependency lists
- forgetting one list silently creates stale results
- field names become part of the invalidation mechanism
- refactors become risky because renaming a cached field is also a dependency
  migration

OceanKit code should prefer explicit invalidation methods such as:

- `knotPointsDidChange()`
- `coefficientsDidChange()`
- `distributionDidChange()`
- `quadratureGridDidChange()`

Each method should clear exactly the cached properties that depend on that
piece of state.

## Public API and visibility

Cache properties are usually implementation details. They should normally be:

- `Access = private`, or
- `GetAccess = public, SetAccess = private` only when developer inspection is
  genuinely useful

Do not make caches part of the ordinary public scientific API just because a
class uses them internally. External callers should not be able to mutate a
cache directly.

When a cached value is useful for debugging or developer-oriented
documentation, prefer documenting it as a developer topic rather than exposing
the entire cache container as a user-facing property.

## Dependent properties

Dependent properties are fine for user-facing derived quantities, but they
should usually be backed by dedicated hidden storage if the computation is
expensive.

For example:

- user-facing property: `basisMatrix`
- hidden storage: `basisMatrix_`

Do not rely on a dependent property alone as the cache. Use the dependent
getter to populate and return hidden storage.

## Keyed caches for parameter sweeps

Some workloads genuinely benefit from keyed subcaches. A common example is
reusing a matrix factorization across many searches over one or two parameters
while the grid stays fixed.

When a keyed cache is justified:

- keep it private
- use it only for a clearly expensive quantity
- make the key include every dependency that changes the result
- clear the keyed cache whenever a non-keyed dependency changes

For a Matérn covariance factor, the key is not just `maternT` and
`maternAlpha` if the evaluation grid can also change. The effective key is the
full dependency tuple:

```text
(grid signature, spline order or quadrature rule if relevant, maternT, maternAlpha)
```

If computing and maintaining a full key is awkward, prefer the simpler and
safer design: keep the keyed cache per instance and clear it whenever the grid
definition changes.

For new code, prefer modern MATLAB containers such as `dictionary` when the
package already targets a release that supports them cleanly. Use
`containers.Map` only when legacy compatibility or required behavior makes it
the better choice.

## Batching updates

If multiple parameter changes would otherwise trigger repeated recomputation,
it is reasonable to provide a narrow batching mechanism.

Good examples:

- `beginUpdate()` and `endUpdate()`
- `willSetMultipleParameters()` and `didSetMultipleParameters()`

Use batching only to suppress redundant recomputation. Do not let batching
become the main source of cache correctness. After the batch ends, the class
should have exactly the same valid state it would have had if the same changes
were applied one at a time.

## Interaction with persistence

For classes that persist scientific state, caches should normally be treated as
rebuildable implementation state rather than serialized state.

In most cases:

- persist canonical state
- do not persist caches
- rebuild caches lazily after construction or load

Only persist a cache when it is part of the actual scientific contract and the
format version explicitly supports it.

## Testing expectations

Every nontrivial cache should have unit tests that cover:

- first access computes the cached value
- repeated access reuses the cached value
- each dependency change invalidates the affected cache
- unrelated state changes do not clear unrelated caches
- batched updates recompute at most once at the end of the batch

If a cache is keyed, add at least one test that proves changing a non-keyed
dependency cannot accidentally reuse an old keyed entry.

## Recommended default for OceanKit

Unless there is a strong reason otherwise, new OceanKit classes should use
this default:

- canonical state in ordinary properties
- cached implementation state in dedicated private properties
- lazy computation in getters or helper methods
- explicit invalidation hooks tied to the actual mutation boundaries
- small keyed caches only for clearly expensive repeated parameter sweeps

This keeps the object model readable, preserves scientific correctness, and
fits the pattern already used successfully in `BSpline`.
