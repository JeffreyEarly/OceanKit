# MATLAB Documentation Agent Guide

## Documentation compliance

When editing inline MATLAB API documentation, use only the structured
`class-docs` tokens supported by the documentation style guide. Use the exact
spellings below, including the leading `- `, and put each token on its own
comment line:

- `- Topic: ...`
- `- Declaration: ...`
- `- Parameter <name>: ...`
- `- Returns <name>: ...`
- `- nav_order: ...`
- `- Developer: true`

Do not replace these with prose headings or unsupported pseudo-tokens such as:

- `Topic:`
- `Declaration:`
- `Parameters:`
- `Returns:`
- `Notes:`
- `See also:`
- `Examples:`

For classes, methods, and properties intended for generated API docs, keep the
source comment order consistent with the guide:

1. one-line summary
2. overview or discussion prose
3. optional example
4. structured token lines

Mathematical APIs should include the governing equations or defining relations
when they materially clarify behavior. Use `$$...$$` for rendered equations and
rendered mathematical aliases.

When updating existing documentation:

- preserve accurate existing documentation
- only change documentation that is incorrect, outdated, or directly affected by
  the code change

Before finishing any documentation edit, do a compliance pass on all touched
files and verify that no edited API doc uses unsupported headings where the
guide requires structured tokens.
