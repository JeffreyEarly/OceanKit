---
layout: default
title: MATLAB style guide
parent: Developers guide
nav_order: 2
---

# MATLAB style guide

This guide focuses on MATLAB class and function design within OceanKit packages. Repository layout, package structure, and documentation-site mechanics are covered separately in the package design guide and the documentation style guide.

When a repository already has a local convention, preserve local consistency unless the change is part of an intentional cleanup.

## Naming

- Class names use `UpperCamelCase`.
- Methods and properties use `lowerCamelCase`.
- Boolean names should read like assertions or switches:
- `isHydrostatic`
- `hasClosure`
- `shouldAntialias`
- Use short repository prefixes for public classes when the class would otherwise live in the global MATLAB namespace and the repository already uses that convention. Examples include `WVTransform` and `CAPropertyAnnotation`.
- Standalone mathematical primitives may omit a repository prefix when the repository context already provides enough namespace, as in `TensorSpline` or `StudentTDistribution`.
- User-defined subclasses outside the package should not reuse a package prefix such as `WV` unless they are intentionally extending that namespace.
- Preserve domain-standard symbols when they carry real scientific meaning. Examples include `N2`, `Lxyz`, `Nxyz`, `xi`, `tKnot`, and `qgpv`.
- Prefer singular names for a single object and plural names for collections.
- Prefer verb phrases for methods and noun phrases for properties.

## Constructors and function signatures

- A class should have one primary constructor and a clear initialization path.
- Additional construction paths should usually be expressed as explicit static factories, for example `waveVortexTransformFromFile` or `geometryFromFile`.
- Use positional arguments for the core identity of the operation or object, and name-value arguments for modifiers, options, and rare inputs.
- Avoid positional boolean flags.
- Prefer modern MATLAB signatures built around:
- `arguments`
- `arguments (Repeating)`
- `arguments (Output)`
- Treat those forms as the default for new public APIs unless they make the code significantly harder to read or create unnecessary boilerplate.
- When a constructor or function forwards options to another API, use `namedargs2cell(options)` rather than rebuilding the name-value list manually.

## Callable form

Choose the callable form that matches the semantics of the operation.

- Use a standalone function for stateless mathematics, transforms, validators, index or shape helpers, and low-level operations that do not depend on object identity.
- Use an instance method for behavior that depends on object state, mutates object state, exposes a core capability of the object, or uses cached or derived state owned by the object.
- Use a static factory for alternate construction paths, reconstruction from serialized state, or constructors that need explicit source naming, for example `waveVortexTransformFromFile`.
- Prefer explicit factory names over generic names such as `fromFile` when multiple object types in the repository may support similar reconstruction paths.

## Class organization and properties

- Put the class definition, property blocks, constructor, and short accessors in `@ClassName/ClassName.m`.
- Put large or specialized methods in separate files inside the same `@ClassName` folder.
- Group properties by access pattern:
- public mutable state
- public read-only state
- dependent properties
- protected or private implementation details
- Default to read-only properties unless external mutation is part of the public API.
- Use `Dependent` for computed values that should not be stored independently.
- Keep cached or implementation-only state clearly separate from public scientific state.
- Use `Hidden` and `protected` access sparingly and only when they communicate a real API boundary.
- Markdown files should soft-wrap. Code, inline comments, and everything else should soft-wrap unless readability clearly benefits from wrapping.

## Method design

- Keep public method names descriptive and action-oriented.
- Use `arguments (Input)` and `arguments (Output)` when they make the interface or return contract clearer.
- For handle classes, mutate in place intentionally. Return `self` only when it materially improves API clarity or chaining.
- Keep method bodies focused. If a method is doing construction, validation, computation, persistence, and plotting, split it.
- Silent reshaping or transposing should be avoided unless it is part of the documented API contract.

### Shape contracts

- Public APIs should document expected input and output shapes.
- Elementwise or pointwise transforms should preserve input shape unless the API explicitly promises normalization.
- One-dimensional grids, coordinate vectors, and flattened observation arrays should return column vectors by default.
- Tensor and grid outputs should follow `ndgrid` ordering consistently.
- If an API intentionally normalizes shape, it must do so consistently and say so explicitly.

## Numerical and scientific computing conventions

- Prefer column vectors for one-dimensional grids and samples unless there is a strong reason otherwise.
- Use `ndgrid`-style dimensional ordering consistently for tensor products and gridded data.
- Keep grid-size and domain-size variables explicit, for example `[Lx Ly Lz]` and `[Nx Ny Nz]`.
- State units in comments or documentation for physical quantities.
- Make tolerances explicit in numerical comparisons.
- Avoid hidden conversions that change the interpretation of physical quantities.
- Do not rename scientifically standard quantities just to satisfy a generic naming rule.

## Errors and warnings

- Prefer structured error identifiers in new code, for example `ConstrainedSpline:InvalidGrid` or `WVTransform:UnknownVariable`.
- Error messages should explain what failed and, when possible, what the caller should do next.
- Use warnings only for recoverable situations.
- Avoid `disp` for ordinary control flow.
- Prefer explicit argument validation and clear `error(...)` calls over silent coercion.
- Reserve `assert` for tests or internal invariants, not for user-facing API validation.

## Persistence and file I/O

- For objects that support serialization, use a clear pair of APIs:
- an instance method such as `writeToFile`
- a static constructor or factory such as `geometryFromFile`, `waveVortexTransformFromFile`, or `annotatedClassFromFile`
- Keep `path` as the first file-related input and put optional behavior in name-value arguments.
- Prefer explicit options such as `shouldOverwriteExisting`, `shouldReadOnly`, and `iTime`.
- Round-trip file I/O should be covered by tests.
- When the serialized representation is NetCDF-backed, keep attributes, dimensions, and variable names stable unless there is a versioned migration.

## MATLAB modernization

This section describes the preferred direction for new OceanKit APIs and opportunistic refactors. It is not a mandate to rewrite stable code wholesale.

- Prefer `string` scalars and string arrays for new public text APIs.
- Prefer `mustBeText` at compatibility boundaries when an API should accept both `char` and `string`.
- Prefer string arrays over cell arrays of char for public collections of names, labels, and identifiers.
- Prefer `dictionary` or `configureDictionary(...)` for new keyed containers in subsystems already using modern MATLAB containers.
- Use `containers.Map` only when needed for legacy compatibility or when a required behavior is not yet matched by the preferred container.
- Prefer name-value syntax such as `foo(bar=...)` over legacy comma-pair calling when writing new code.
- Prefer `join`, `replace`, `startsWith`, `endsWith`, `contains`, and related string operations over manual char-array manipulation for new text-handling code.
- Keep modernization incremental. Do not rewrite stable legacy code just to change container or text types unless the surrounding API is already being refactored.
