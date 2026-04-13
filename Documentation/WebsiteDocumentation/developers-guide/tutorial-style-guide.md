---
layout: default
title: Tutorial style guide
parent: Developers guide
nav_order: 6
---

# Tutorial style guide

OceanKit tutorials should teach by combining narrative, math, code, and figures. This guide stays intentionally lightweight and points to the existing parser and build machinery where the tutorial format becomes mechanical.

## Core style

- Start from the scientific question or workflow, not from package internals.
- Alternate explanation with runnable code instead of dumping long uninterrupted scripts.
- Keep tutorial code on the scientific path. Do not add tutorial-local defensive error checking, file-existence guards, or input validation around required assets or standard package APIs.
- Use math when it clarifies the scientific model or notation, not as decoration.
- Use backticks only for literal MATLAB or API identifiers. Render mathematical symbols such as $$\kappa$$, $$\Delta t$$, or $$\mathbf{u}$$ with `$$...$$`, not Markdown code quotes.
- Use figures to confirm and explain the main point of a section.
- Let each section teach one main step or idea.
- Prefer self-contained synthetic or lightweight examples before asset-heavy case studies.

## Keep blocks short

- Usually keep a code block to about `10-20` lines or one main step before returning to prose or a figure.
- Usually keep displayed math to one main relation or a short cluster of closely related equations before returning to explanation.
- If a code block or math block starts to feel long, break it up with narrative, interpretation, or a figure.
- Combine adjacent setup and solve steps when they are one easy mental action, such as preparing trajectories and immediately fitting the model from them.
- Avoid long stretches that are only prose, only code, or only equations.

## Standard OceanKit pattern

For new OceanKit tutorials, the standard pattern is a runnable MATLAB script in `Examples/Tutorials` that is rendered into a website page during the documentation build.

- Put tutorial source scripts in `Examples/Tutorials`.
- Begin the script with a `%% Tutorial Metadata` section.
- In that metadata section, use `% Title:`, `% Slug:`, `% Description:`, and optionally `% NavOrder:`.
- Use `%%` headings to define tutorial sections.
- Use single-`%` comment lines for narrative prose.
- Use section headings to mark teaching beats, not tiny implementation mechanics. When a short preparation step exists only to feed the immediately following solve step, keep them in one section.
- Guard optional `tutorialFigureCapture(...)`, `tutorialMovieCapture(...)`, and `tutorialOutputCapture(...)` calls inline at the call site. Do not define top-of-script no-op fallback handles.
- The optional `tutorialFigureCapture(...)`, `tutorialMovieCapture(...)`, and `tutorialOutputCapture(...)` guards are the standard exception to the no-defensive-checking rule above.
- Call `tutorialFigureCapture(...)` when a figure should appear in the rendered page.
- Call `tutorialMovieCapture(...)` when a generated movie should appear in the rendered page.
- Call `tutorialOutputCapture(...)` when command-window output should appear in the rendered page as a fenced text block. Ordinary `disp(...)` and `fprintf(...)` output is not captured automatically.
- Add each new tutorial source file to the tutorial source list in the repository `tools/build_website_documentation.m` script.

Treat the script as the source of truth. The generated markdown page in `docs/` is build output.

## Tutorial assets

Every script-built tutorial owns one page and one asset root:

- page: `docs/tutorials/<slug>.md`
- asset root: `docs/tutorials/<slug>/`

Within that root:

- generated PNG figures, generated poster images, and generated movies live directly in `docs/tutorials/<slug>/`
- curated static tutorial-owned files live only in `docs/tutorials/<slug>/preserved/`

Tutorial builds are reused at the whole-tutorial level when the tutorial source file is unchanged:

- each tutorial writes a source-based build stamp beside `docs/tutorials/<slug>.md`
- when that stamp matches a previous build, the docs build reuses the prior page and asset root without executing the tutorial script
- use `build_website_documentation(rebuildTutorials=true)` when package or dependency changes should refresh tutorial outputs
- whole-tutorial invalidation is source-based in this pass, not dependency-declaration based

For generated movies, use the normal tutorial pipeline rather than a separate manual docs path:

- expensive generated movies should be reused by stamp instead of rebuilt on every docs run
- movie-level stamps still apply inside a forced or source-invalidated rebuild
- keep tracked tutorial movies at or below `25 MB`
- `preserved/` is for curated website assets, not the normal generated tutorial outputs

## Nuts and bolts

The current mechanics live in the implementation, not in this guide:

- [`class-docs/TutorialDocumentation.m`](https://github.com/JeffreyEarly/class-docs/blob/main/TutorialDocumentation.m)
- [`class-docs/TutorialBuildRuntime.m`](https://github.com/JeffreyEarly/class-docs/blob/main/TutorialBuildRuntime.m)
- [`spline-core/Examples/Tutorials/BSplineFoundations.m`](https://github.com/JeffreyEarly/spline-core/blob/main/Examples/Tutorials/BSplineFoundations.m)
- [`spline-core/tools/build_website_documentation.m`](https://github.com/JeffreyEarly/spline-core/blob/main/tools/build_website_documentation.m)

Older live-script or export-driven tutorial formats still exist in some repositories, but they are legacy context rather than the default pattern for new OceanKit tutorials.

## Small advice

- Make figure captions interpret the result rather than restating the axes.
- Keep notation, terminology, and API names consistent with the package reference.
- Use enough narrative that a reader understands why the next code block exists.
- End sections with a visible result when possible: a figure, a recovered parameter, a comparison, or a short interpretation.
- Stay within the MathJax commands already supported by the site. Prefer simple supported notation such as $$\mathbf{z}$$ over unsupported commands such as `\boldsymbol`.
