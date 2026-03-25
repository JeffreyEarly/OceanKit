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

- `- Topic: ...` Assigns the class, method, or property to a topic in the generated site.
- `- Declaration: ...` Supplies the displayed declaration block.
- `- Parameter <name>: ...` Adds a parameter entry to the generated Parameters section.
- `- Returns <name>: ...` Adds a return-value entry to the generated Returns section.
- `- nav_order: ...` Controls ordering within a topic when multiple items appear together.
- `- Developer: ...` Marks an item as internal so it renders under Developer Topics.

### Topic hierarchy

Nested topics use the em dash character `—`, not a hyphen-minus `-`.

For example:

```text
- Topic: Create a spline
- Topic: Create a spline — Prepare knot sequences
- Topic: Create a spline — Prepare knot sequences — Boundary handling
```

Use human-readable topic names. Keep them stable across related methods so generated indexes stay coherent.

### Recommended source-comment order

For classes, methods, and properties that are meant to generate good API pages, use this order inside the help text:

1. one-line summary
2. overview or discussion prose
3. optional example
4. metadata token lines

That gives MATLAB metadata a clean summary line and keeps the detailed discussion readable in source form.

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

## Example pattern

The following pattern is compatible with the current parser:

```matlab
function xi = exampleMethod(self, inputValue, options)
% Compute the coefficient vector `xi`, written mathematically as $$\xi$$.
%
% This method solves the linear system
%
% $$
% \mathbf{A}\xi = b,
% $$
%
% where $$\mathbf{A}$$ is the assembled operator, $$\xi$$ is the unknown
% coefficient vector returned by the method, and $$b$$ is the right-hand
% side assembled from `inputValue`.
%
% ```matlab
% xi = obj.exampleMethod(3, mode="fast");
% ```
%
% - Topic: Example topic — Example subtopic
% - Declaration: xi = exampleMethod(self,inputValue,options)
% - Parameter inputValue: scalar numeric input
% - Parameter options.mode: optional operating mode
% - Returns xi: coefficient vector $$\xi$$ that solves the assembled system
```

## Practical rules

- Put each structured token on its own comment line.
- Keep token spelling consistent even though the parser is case-insensitive.
- Use ordinary markdown for examples and explanatory prose.
- Use `mathjax: true` on hand-authored pages that include equations or rendered symbols.
- Prefer `$$...$$` when you need rendered equations or rendered mathematical aliases in OceanKit docs.
- Keep topic names and declarations accurate when signatures change.
- Document only tokens that the current parser supports.
