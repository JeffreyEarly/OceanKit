Documentation
=============

+ The `WebsiteDocumentation` folder contains the markdown files that act as the source of truth for the OceanKit GitHub Pages site.
+ The generated `/docs` folder is built from that source and committed for publishing.
+ The `../tools/build_website_documentation.m` script handles the copy/build step for this repository.

The site is published at [jeffreyearly.github.io/OceanKit](https://jeffreyearly.github.io/OceanKit/) by GitHub Pages from `main:/docs`. Keep the project-site `url` and `baseurl` in the canonical `_config.yml` and use Liquid's `relative_url` filter for site assets and links.

Editable branding lives in `Branding/`; its [README](Branding/README.md) describes the SVG masters and raster exports. Export the branding before building the documentation:

```matlab
addpath("tools")
build_website_documentation(rootDir=pwd)
```

OceanKit has no `buildfile.m` or `docs:check` task. Verify that `docs/` matches `WebsiteDocumentation/` with `diff -rq Documentation/WebsiteDocumentation docs`, run `git diff --check`, and inspect the GitHub Pages build and rendered site after publication. Package authoring repositories may provide their own `docs:check` tasks.
