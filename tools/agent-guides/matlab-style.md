# MATLAB Style Agent Guide

## MATLAB style rules

Apply the OceanKit MATLAB style guide during edits. In particular:

- classes use `UpperCamelCase`
- methods and properties use `lowerCamelCase`
- preserve standard scientific notation when it carries real mathematical
  meaning, such as `sigma_n`, `sigma_s`, `Mxx`, `Myy`, `Mxy`, `Lx`, `Ly`, and
  `N2`
- do not rename scientifically standard quantities just to satisfy a generic
  naming rule
- use `@ClassName` folders for class-based APIs
- keep implementation-only methods in separate files inside the same class
  folder when that materially improves readability
- validate inputs at public API boundaries
- use structured, actionable error identifiers and messages

## MATLAB formatting rules

- Keep MATLAB function calls on one line whenever the full call fits within the
  repository line-length limit and each argument is a simple literal, variable,
  or short name-value pair.
- Do not split function calls across multiple lines purely for visual
  alignment.
- Only wrap a function call when the one-line form would exceed the line-length
  limit or when one or more arguments are long expressions, nested calls,
  anonymous functions, or inline comments.
- Do not introduce manual `...` continuations for simple calls that fit on one
  line.
