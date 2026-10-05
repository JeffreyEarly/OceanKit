# OceanKit branding

The project mark is a circular navy O containing one broad teal ocean wave. Generous space above the wave keeps the symbol open, and a clear gap between the wave and ring preserves the shape in monochrome. Use this identity across the OceanKit ecosystem; package projects retain their own identities.

The source and export layout follows [WaveVortexModel's branding folder](https://github.com/JeffreyEarly/wave-vortex-model/tree/main/Documentation/Branding). The original concept was selected in the OceanKit branding conversation, then drawn as native SVG paths with flat colors.

## Source artwork

- `ocean-o-concept.png` preserves the selected imagegen concept at its original 1254 × 1254 resolution. Use the SVG-derived exports for production.
- `ocean-mark.svg` contains the colored mark on a transparent background.
- `ocean-mark-monochrome.svg` uses `currentColor`; set the root SVG's `color` or use it inline to choose a single color. The gap around the wave remains transparent.
- `icon-rounded.svg` places the mark on a pale rounded square with transparent exterior corners.
- `icon-square.svg` has an opaque square background for system-masked shortcuts.
- `icon-small.svg` simplifies the crest and increases the clear space around the wave for tiny displays.
- `logo-horizontal-light.svg` and `logo-horizontal-dark.svg` combine the rounded icon with outlined Avenir Next Medium lettering. The suffix describes the intended background: navy lettering on light surfaces, pale lettering on dark surfaces.
- `social-card.svg` is the editable 1200 × 630 sharing layout.
- `generation-prompt.md` records the imagegen prompt and the vector-production decisions.
- `preview.html` displays the sources and exported sizes together.

The SVGs contain native paths with no embedded raster images or runtime font dependency. Font files are not distributed with the artwork.

## Published exports

The canonical website assets live in `Documentation/WebsiteDocumentation/assets/branding/`. The documentation build copies them into the committed `docs/assets/branding/` directory. Edit sources and rebuild instead of changing the generated `docs` tree.

| File | Dimensions | Use |
| --- | --- | --- |
| `icon-512.png`, `icon-1024.png` | 512 × 512, 1024 × 1024 | Website illustrations and repository avatars; transparent exterior corners |
| `icon-rounded.svg`, `icon-square.svg` | Scalable | Full icon in rounded or opaque square form |
| `icon-small.svg` | Scalable | Simplified small-display artwork |
| `favicon-16.png`, `favicon-32.png`, `favicon-48.png`, `favicon-96.png` | Named pixel sizes | Browser and search icons |
| `favicon.ico` | 16, 32, and 48 px images | Browser compatibility |
| `apple-touch-icon.png` | 180 × 180 | Opaque square shortcut icon |
| `logo-horizontal-light.svg`, `logo-horizontal-dark.svg` | 325 × 100 | Website header and scalable wordmarks |
| `logo-horizontal-light.png`, `logo-horizontal-dark.png` | 650 × 200 | Wordmark raster fallback at twice the SVG dimensions |
| `social-card.png`, `social-card.svg` | 1200 × 630 | Link-sharing preview |
| `ocean-mark.svg`, `ocean-mark.png` | Scalable, 1600 × 1600 | Standalone colored mark |
| `ocean-mark-monochrome.svg`, `ocean-mark-monochrome.png` | Scalable, 1600 × 1600 | Single-color mark; PNG is black |

Use the small icon for favicons and the square icon for Apple-touch shortcuts. Use the social card as the repository's GitHub social preview and the website's sharing image.

## Re-exporting

`export-assets.mjs` copies the vector masters, outlines the wordmark and sharing-card lettering, and renders the PNG and ICO exports. It requires Node.js, `sharp` 0.35.4, `@napi-rs/canvas` 0.1.100, and Avenir Next Medium. These are graphics-authoring tools, not MATLAB runtime dependencies. Pass `--modules-dir` if the Node packages are outside the usual module search path.

```sh
node Documentation/Branding/export-assets.mjs --font '/System/Library/Fonts/Avenir Next.ttc' --modules-dir /path/to/node_modules
```

For a TTC font collection the exporter selects `AvenirNext-Medium` by PostScript name; `--font-face` can select another named face deliberately. For a single-face TTF/OTF, supply the medium font directly. Font glyphs are converted to paths in the exported SVGs. Rasterizing the committed SVGs needs no installed font.

After exporting, run the documentation build and check described in `Documentation/README.md`. The concept image and editable masters stay outside the website source; production assets are copied into the website source by the exporter.

The flat palette uses navy `#052f60`, teal `#119ba0`, and pale `#edf7fa`. The social card uses pale lettering on navy with supporting text in `#b0e6e6`. Preserve the contrast and the clear space around the wave when preparing other layouts. Keep the standalone mark on a light surface or use the monochrome mark in a suitable light color on dark surfaces.
