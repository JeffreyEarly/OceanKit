---
layout: default
title: Documentation style guide
parent: Developers guide
nav_order: 3
mathjax: true
---

# Documentation style guide

OceanKit documentation has two layers:

1. website pages written directly in markdown under `Documentation/WebsiteDocumentation`
2. inline class, method, and property comments that are parsed by the `class-docs` package into generated API documentation

This guide covers both. The key rule is simple: only document tokens that the current `class-docs` parser actually understands.

## Website markdown pages

Website pages are the source of truth for hand-authored GitHub Pages content. Keep these pages in `Documentation/WebsiteDocumentation` and treat `docs/` as generated output.

### Standard front matter

Use the standard Just the Docs front matter fields used across OceanKit sites:

- `layout` Usually `default`.
- `title` The page title shown in navigation and page heading.
- `parent` The immediate navigation parent for child pages.
- `grand_parent` Use only when a page needs to render two navigation levels down.
- `nav_order` Numeric ordering within the current navigation level.
- `has_children` Use on section index pages that only organize child pages.
- `permalink` Use when a section index should have a stable site path such as `/developers-guide`.
- `mathjax` Set to `true` for pages that contain equations or rendered symbols such as `$$\xi$$` or `$$K$$`.

### Page structure

- Use a section index page when a topic has multiple child pages.
- Keep page titles short and stable.
- Keep the source pages under `Documentation/WebsiteDocumentation` clean and hand-maintained.
- Soft-wrap markdown source files. Do not hard-wrap ordinary prose paragraphs.
- Regenerate `docs/` from source rather than editing `docs/` first.
- Use ordinary markdown for examples, notes, warnings, and cross-links.

## Inline API documentation for class-docs

The `class-docs` package extracts MATLAB metadata plus a small set of structured tokens from the detailed help text. The parser is case-insensitive, but OceanKit should standardize on the capitalized token spellings shown below.

### Supported structured tokens

These are the exact metadata tokens currently supported:

```text
- Topic: ...
- Declaration: ...
- Parameter <name>: ...
- Returns <name>: ...
- nav_order: ...
- Developer: ...
```

Do not promise or rely on other structured tokens. In particular, `class-docs` does not currently expose special structured handling for tokens such as:

- `- Example:`
- `- Examples:`
- `- See also:`
- `- Notes:`

Those can still appear as ordinary prose headings or markdown, but they are not parsed as metadata.

### What each token does

- `- Topic: ...` On class comments, declares topic groups and their order in the generated class index. On method and property comments, assigns that item to one topic path.
- `- Declaration: ...` Supplies the displayed declaration block.
- `- Parameter <name>: ...` Adds a parameter entry to the generated Parameters section.
- `- Returns <name>: ...` Adds a return-value entry to the generated Returns section.
- `- nav_order: ...` Optionally controls ordering within a topic when multiple items appear together. Smaller values sort earlier.
- `- Developer: ...` Marks an item as internal so it renders under Developer Topics.

### Generated API scope and property coverage

`class-docs` generates a class index page plus one page for each public, non-hidden method and property selected by the documentation build.

That means public properties are first-class API documentation targets, not secondary metadata. Document every public, non-hidden property that should appear in the generated API reference, including dependent properties and read-only computed properties.

Use the same comment structure for documented properties as for methods:

1. one-line summary
2. overview or discussion prose
3. optional example
4. metadata token lines

If a dependent property also has a getter or setter with its own help text, that accessor documentation does not replace the property comment. Document the property where it is declared.

Private, protected, hidden, or implementation-only storage slots may still have source comments for local clarity, but they are not the minimum generated-doc baseline. If a public item is useful as internal reference but should live under Developer Topics, keep it documented and mark it with `- Developer: true`.

### Topic authoring and hierarchy

Every documented method or property should normally have exactly one `- Topic:` line. If an item has no topic, the current `class-docs` output places it under `Other`. Use that fallback only intentionally.

Topic metadata appears in two places, and the distinction matters:

- Class-level `- Topic:` lines declare the topic groups and their display order for the class index.
- Method-level and property-level `- Topic:` lines assign each documented item to one topic or topic path.

Put the ordered class-level topic list near the end of the class overview block, after the prose and optional example, and before `- Declaration: ...`.

Nested topics use the em dash character `—`, not a hyphen-minus `-`.

For example:

```text
- Topic: Create a spline
- Topic: Create a spline — Prepare knot sequences
- Topic: Create a spline — Prepare knot sequences — Boundary handling
```

For most classes, stop at one or two topic levels. Use a third level only when it materially improves navigation. Reuse the same topic vocabulary across related methods and properties so the generated index stays coherent.

### Topic naming and ordering

Prefer imperative, task-oriented topic names, effectively using active voice. Use names such as `Create a spline`, `Inspect spline properties`, `Evaluate the spline`, `Transform the spline`, `Prepare knot sequences`, and `Compile constraints`.

- Prefer user-task language over internal implementation labels.
- Keep sibling topic names parallel in grammar and scope.
- Reuse topic names exactly rather than creating near-duplicates.
- Prefer `Create a spline` over `Spline creation` or `Creating splines`.
- Avoid passive or vague names unless the topic is inherently descriptive.

When defining the class-level topic list, a good default order is:

1. construction or creation topics
2. inspection or reference topics
3. core use or evaluation topics
4. transform, analysis, or preparation topics
5. developer or internal topics last

If a class has developer-only groups, still declare those groups in the class-level topic list so their order is stable. Pair the corresponding methods or properties with `- Developer: true`.

### Ordering within a topic

Use `- nav_order:` when deterministic ordering within one topic materially matters. Do not require it by default.

If `- nav_order:` is omitted, item order falls back to current parser behavior and default ordering. That fallback is acceptable when no specific sequence matters, but it should not be treated as an authored semantic guarantee.

### Recommended source-comment order

For classes, methods, and properties that are meant to generate good API pages, use this order inside the help text:

1. one-line summary
2. overview or discussion prose
3. optional example
4. metadata token lines

That gives MATLAB metadata a clean summary line and keeps the detailed discussion readable in source form. Once the metadata token block begins, keep the rest of the comment block limited to token lines only. Do not return to ordinary prose after the tokens.

### Summary vs discussion

The short summary shown on generated pages comes from MATLAB metadata, not from a separate token. In practice:

- the first help-comment line becomes the short description
- the remaining detailed help text becomes the Overview or Discussion section after structured tokens are stripped out

Do not invent a `- Summary:` token. It is not supported.

### Developer-only items

Use `- Developer: true` for internal implementation details that should be documented but not emphasized as primary public API.

The parser currently treats values such as `true`, `1`, and `yes` as truthy, but OceanKit should standardize on:

```text
- Developer: true
```

Use this sparingly. A page becomes more useful when developer topics are clearly internal rather than a second copy of the public API.

### Mathematical exposition

For mathematically meaningful APIs, mathematical exposition should be treated as part of the documentation, not as an optional flourish. This applies to both public items and internal items marked with `- Developer: true`.

- When a class, method, property, or parameter represents a mathematical object, transform, operator, normalization, constraint, discretization, or state relation, the documentation should explain the governing math whenever it materially clarifies behavior.
- Include the defining equation, operator, or relation when practical.
- Define symbols near their first use and connect them back to code identifiers, dimensions, and units when relevant.
- When a code identifier is itself a mathematical symbol, introduce both forms early in the summary or discussion. For example, document `xi` as the coefficient vector `$$\xi$$`, and `K` as the spline order `$$K$$`.
- Use `$$...$$` delimiters consistently for both displayed equations and rendered symbol forms.
- Because generated `class-docs` page titles and headings use the raw MATLAB identifier, put the rendered mathematical alias in the short summary and/or discussion text rather than relying on the page title.
- Use math to explain the API, not as decoration. If there is no meaningful mathematical interpretation, ordinary prose is sufficient.

## Example patterns

The following patterns are compatible with the current parser and match the style used in the spline-core reference classes.

### Class header with ordered topics

```matlab
classdef ExampleSpline < handle
    % Create, inspect, and evaluate a tensor-product spline model.
    %
    % `ExampleSpline` stores one knot vector per dimension together with a
    % coefficient array for fast evaluation on rectilinear grids.
    %
    % ```matlab
    % spline = ExampleSpline(grid, values, S=3);
    % valuesFit = spline(Xq, Yq);
    % ```
    %
    % - Topic: Create a spline
    % - Topic: Inspect spline properties
    % - Topic: Evaluate the spline
    % - Topic: Transform the spline
    % - Topic: Prepare fit inputs
    % - Topic: Solve fit systems
    % - Declaration: classdef ExampleSpline < handle
end
```

### Documented property

```matlab
properties (SetAccess = private)
    % Grid vectors used to define the interpolation lattice.
    %
    % In one dimension `gridVectors` contains one column vector. In higher
    % dimensions it stores one grid vector per tensor axis.
    %
    % ```matlab
    % spline.gridVectors
    % ```
    %
    % - Topic: Inspect interpolation grids
    gridVectors
end
```

### Documented constructor or method

```matlab
methods
    function self = ExampleSpline(grid, values, options)
        % Create a spline from rectilinear grid data.
        %
        % Use this constructor when the grid vectors are already known and
        % the sampled values should define the fitted spline.
        %
        % - Topic: Create a spline
        % - Declaration: self = ExampleSpline(grid,values,options)
        % - Parameter grid: numeric vector in 1-D or cell array of grid vectors
        % - Parameter values: sampled values on the grid
        % - Parameter options.S: spline degree scalar or vector
        % - Returns self: ExampleSpline instance
    end
end
```

### Developer-only item with optional nav_order

```matlab
properties (SetAccess = private)
    % Cached normal-equation matrix used by the fit.
    %
    % This cache stores the matrix assembled from the weighted design matrix.
    %
    % - Topic: Solve fit systems
    % - nav_order: 20
    % - Developer: true
    systemMatrix
end
```

## Practical rules

- Put each structured token on its own comment line.
- Keep all ordinary prose and examples above the token block.
- Keep token spelling consistent even though the parser is case-insensitive.
- Give every documented public method and property a `- Topic:` line unless `Other` is intentional.
- Reuse class-level topic names exactly and keep the class-level topic list in the intended display order.
- Use `- nav_order:` only when you need to pin within-topic order.
- Use ordinary markdown for examples and explanatory prose.
- Use `mathjax: true` on hand-authored pages that include equations or rendered symbols.
- Prefer `$$...$$` when you need rendered equations or rendered mathematical aliases in OceanKit docs.
- Keep topic names and declarations accurate when signatures change.
- Document only tokens that the current parser supports.
