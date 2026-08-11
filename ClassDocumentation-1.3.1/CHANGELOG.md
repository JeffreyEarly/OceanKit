# Version History

## [1.3.1] - 2026-08-11
- rendered machine-readable property units as consistent upright MathJax while preserving raw annotation metadata, legacy unit spellings, readable dimensionless and directional prose, and literal fallback for unknown values

## [1.3.0] - 2026-08-08
- added explicit `tutorialOutputCapture(...)` support so tutorials can emit captured console output as fenced markdown blocks with preserved ordering, normalized text cleanup, and regression coverage
- preserved class landing pages when an API member is named `index`
- added configurable method-page grandparents while retaining the existing `Classes` default
- made repeated scalar metadata deterministic by using the first occurrence while preserving repeatable parameter and return metadata

## [1.2.0] - 2026-04-10
- added scoped annotation-sidecar loading support for documentation builders so generated class reference pages keep markdown sidecars without incurring annotation sidecar lookup overhead in ordinary runtime code
- normalized reflected friend-access metadata so friend-only methods and properties stay out of published public API docs
- preserved reflected property summaries, topics, declarations, access, and visibility when merging `CAAnnotatedClass` property annotations, while still enriching reflected property pages with annotation metadata such as units and dimensions
- added property type rendering for documented properties so generated pages show reflected validation class and size metadata and fall back to `CAObjectProperty` `className` and `sizeText` metadata when reflection is absent
- added source-stamped tutorial build reuse so unchanged tutorial pages, figures, movies, and movie posters can be copied from a previous docs build instead of being regenerated, and added `rebuildWebsiteDocumentationFromSource()` to rebuild docs trees while preserving selected generated assets

## [1.1.1] - 2026-03-22
- Preserved existing tutorial image files when regenerated figures are pixel-identical, reducing unnecessary PNG churn in downstream repos

## [1.1.0] - 2026-03-21
- Added `TutorialDocumentation` and `TutorialBuildRuntime` for generating tutorial pages and figures from runnable MATLAB scripts

## [1.0.1] - 2026-01-01
- De-denting docs to fix markdown parsing
