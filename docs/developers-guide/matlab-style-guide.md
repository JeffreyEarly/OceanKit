---
layout: default
title: MATLAB style guide
parent: Developers guide
nav_order: 2
---

# MATLAB style guide

This guide covers MATLAB class and function design within OceanKit packages. Repository structure and documentation mechanics belong in the package design guide and documentation style guide instead.

Prefer local consistency unless the change is part of an intentional cleanup.

## Core principles

- Prefer consistency with the local codebase over abstract style purity.
- Make APIs explicit, predictable, and readable rather than clever.
- Preserve scientific meaning in naming, shapes, units, and defaults.
- Use modern MATLAB features when they improve clarity and safety.
- Do not refactor stable code purely for style.

## Naming

- Classes use `UpperCamelCase`.
- Methods and properties use `lowerCamelCase`.
- Boolean names read like assertions or switches: `isHydrostatic`, `hasClosure`, `shouldAntialias`.
- Methods should read like verbs. Properties should read like nouns.
- Prefer singular names for single objects and plural names for collections.
- Preserve standard scientific notation when it carries real meaning, such as `N2`, `Lxyz`, `Nxyz`, `xi`, `tKnot`, and `qgpv`.
- Use short repository prefixes for public classes when the class would otherwise live in the global MATLAB namespace and the repository already uses that convention, as in `WVTransform` or `CAPropertyAnnotation`.
- Standalone mathematical primitives may omit a repository prefix when the repository context already provides enough namespace, as in `TensorSpline` or `StudentTDistribution`.
- Subclasses outside the package should not reuse a package prefix such as `WV` unless they intentionally extend that namespace.

## API design

- Give each class one primary constructor and a clear initialization path.
- Use positional arguments for the core identity of the operation or object, and name-value arguments for modifiers, options, and rare inputs.
- Avoid positional boolean flags.
- Prefer explicit alternate construction paths such as `geometryFromFile` or `waveVortexTransformFromFile`.
- Prefer modern MATLAB signatures built around `arguments`, `arguments (Repeating)`, and `arguments (Output)` for new public APIs unless they materially hurt readability or add unnecessary boilerplate.
- When forwarding options to another API, use `namedargs2cell(options)` instead of rebuilding name-value pairs manually.

## Callable form

- Use standalone functions for stateless mathematics, transforms, validators, index helpers, shape helpers, and low-level operations that do not depend on object identity.
- Use instance methods for behavior that depends on object state, mutates object state, or exposes a core capability of the object.
- Use static factories for alternate construction paths, reconstruction from serialized state, or source-specific constructors.
- Prefer explicit factory names over generic names such as `fromFile` when multiple object types may support similar reconstruction paths.

## Class organization and properties

- Put the class definition, property blocks, constructor, and short accessors in `@ClassName/ClassName.m`.
- Put large or specialized methods in separate files inside the same `@ClassName` folder.
- Group properties by access pattern: public mutable state, public read-only state, dependent properties, then protected or private implementation details.
- Default to read-only properties unless external mutation is part of the public API.
- Use `Dependent` for computed values that should not be stored independently.
- Keep cached or implementation-only state clearly separate from public scientific state.
- Use `Hidden` and `protected` only when they communicate a real API boundary.

## Method design and shape contracts

- Keep public method names descriptive and action-oriented.
- Keep methods focused. If one method is doing construction, validation, computation, persistence, and plotting, split it.
- For handle classes, mutate in place intentionally. Return `self` only when it materially improves API clarity.
- Public APIs should document expected input and output shapes.
- Elementwise or pointwise transforms should preserve input shape unless the API explicitly promises normalization.
- One-dimensional grids, coordinate vectors, and flattened observation arrays should return column vectors by default.
- Tensor and grid outputs should follow `ndgrid` ordering consistently.
- Avoid silent reshaping or transposing unless it is part of the documented contract.
- If an API intentionally normalizes shape, do it consistently and say so explicitly.

## Numerical and scientific computing conventions

- Keep units explicit in code comments or API documentation for physical quantities.
- Make tolerances explicit in numerical comparisons.
- Keep grid-size and domain-size variables readable, for example `[Lx Ly Lz]` and `[Nx Ny Nz]`.
- Avoid hidden conversions that change the interpretation of physical quantities.
- Do not rename scientifically standard quantities just to satisfy a generic naming rule.

## Errors, validation, and robustness

- Validate inputs at public API boundaries.
- Prefer clear, structured error identifiers using a package or class prefix, for example `ConstrainedSpline:InvalidGrid` or `WVTransform:UnknownVariable`.
- Error messages should explain what failed and, when possible, what the caller should do next.
- Use warnings only for recoverable situations.
- Avoid silent coercion and avoid `disp` for ordinary control flow.
- Reserve `assert` for tests or internal invariants, not user-facing API validation.

## Persistence and file I/O

- Prefer explicit persistence pairs such as an instance method `writeToFile` plus a source-specific factory such as `geometryFromFile`, `waveVortexTransformFromFile`, or `annotatedClassFromFile`.
- Avoid hidden I/O side effects in constructors.
- Keep `path` as the first file-related input and put optional behavior in name-value arguments.
- Prefer explicit file options such as `shouldOverwriteExisting`, `shouldReadOnly`, and `iTime`.
- When a format is NetCDF-backed, keep dimensions, variable names, and attributes stable unless there is a versioned migration.

## MATLAB modernization

- Prefer `string` scalars and string arrays for new public text APIs.
- Use `mustBeText` at compatibility boundaries when an API should accept both `char` and `string`.
- Prefer string arrays over cell arrays of char for public collections of names, labels, and identifiers.
- Prefer `dictionary` or `configureDictionary(...)` for new keyed containers in subsystems already using modern MATLAB containers.
- Use `containers.Map` only when legacy compatibility or required behavior makes it necessary.
- Prefer name-value syntax such as `foo(bar=...)` over legacy comma-pair calling in new code.
- Prefer `join`, `replace`, `startsWith`, `endsWith`, `contains`, and related string operations over manual char-array manipulation in new text-handling code.
- Modernize incrementally. Do not rewrite stable legacy code just to change container types or text types unless the surrounding API is already being refactored.

## Summary

Write MATLAB code that is clear, explicit, scientifically faithful, and progressively modern without introducing churn for its own sake.
