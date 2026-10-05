// Export the native SVG artwork and outlined wordmarks for the website.
// Usage: node export-assets.mjs --font /path/to/font.ttf [--modules-dir /path/to/node_modules]
// A TTC collection selects AvenirNext-Medium, or the PostScript name passed with --font-face.
import fs from 'node:fs/promises';
import path from 'node:path';
import { createRequire } from 'node:module';
import { fileURLToPath } from 'node:url';

const source = path.dirname(fileURLToPath(import.meta.url));
const destination = path.resolve(source, '../WebsiteDocumentation/assets/branding');
const args = process.argv.slice(2);
function option(name) {
  const index = args.indexOf(name);
  return index === -1 ? undefined : args[index + 1];
}
const font = option('--font');
if (!font) throw new Error('Supply --font with the Avenir Next font used for the wordmarks.');
const modules = option('--modules-dir');
const require = createRequire(modules ? path.join(path.resolve(modules), '_oceankit_export.cjs') : import.meta.url);
const sharp = require('sharp');
const { createCanvas, GlobalFonts, SvgExportFlag } = require('@napi-rs/canvas');
// Canvas loads only the first face of a TTC. Repackage the named face in memory
// so a system font collection cannot silently substitute its first, bold face.
function selectFontFace(buffer, postScriptName) {
  if (buffer.toString('ascii', 0, 4) !== 'ttcf') return buffer;
  for (let face = 0; face < buffer.readUInt32BE(8); face++) {
    const start = buffer.readUInt32BE(12 + 4 * face);
    const count = buffer.readUInt16BE(start + 4);
    const tables = Array.from({ length: count }, (_, index) => {
      const record = start + 12 + 16 * index;
      return { record, tag: buffer.toString('ascii', record, record + 4),
        offset: buffer.readUInt32BE(record + 8), length: buffer.readUInt32BE(record + 12) };
    });
    const names = tables.find(table => table.tag === 'name');
    const strings = names.offset + buffer.readUInt16BE(names.offset + 4);
    let matches = false;
    for (let index = 0; index < buffer.readUInt16BE(names.offset + 2); index++) {
      const record = names.offset + 6 + 12 * index;
      if (buffer.readUInt16BE(record + 6) !== 6) continue;
      const length = buffer.readUInt16BE(record + 8);
      const offset = strings + buffer.readUInt16BE(record + 10);
      const bytes = Buffer.from(buffer.subarray(offset, offset + length));
      const platform = buffer.readUInt16BE(record);
      const name = platform === 0 || platform === 3 ? bytes.swap16().toString('utf16le') : bytes.toString('ascii');
      if (name === postScriptName) matches = true;
    }
    if (!matches) continue;
    let cursor = 12 + 16 * count;
    const output = Buffer.alloc(cursor + tables.reduce((sum, table) => sum + Math.ceil(table.length / 4) * 4, 0));
    buffer.copy(output, 0, start, start + 12);
    let headOffset;
    for (const [index, table] of tables.entries()) {
      const record = 12 + 16 * index;
      buffer.copy(output, record, table.record, table.record + 16);
      output.writeUInt32BE(cursor, record + 8);
      buffer.copy(output, cursor, table.offset, table.offset + table.length);
      if (table.tag === 'head') { headOffset = cursor; output.writeUInt32BE(0, cursor + 8); }
      cursor += Math.ceil(table.length / 4) * 4;
    }
    let checksum = 0;
    for (let index = 0; index < output.length; index += 4) checksum = (checksum + output.readUInt32BE(index)) >>> 0;
    output.writeUInt32BE((0xb1b0afba - checksum) >>> 0, headOffset + 8);
    return output;
  }
  throw new Error(`Font collection does not contain ${postScriptName}.`);
}
const fontBuffer = selectFontFace(await fs.readFile(path.resolve(font)), option('--font-face') ?? 'AvenirNext-Medium');
if (!GlobalFonts.register(fontBuffer, 'OceanKitBrand')) {
  throw new Error('The supplied font could not be loaded.');
}
await fs.mkdir(destination, { recursive: true });

function svgParts(svg) {
  const match = svg.match(/<svg\b([^>]*)>([\s\S]*?)<\/svg>\s*$/);
  if (!match) throw new Error('Expected an SVG document.');
  const viewBox = match[1].match(/viewBox="([^"]+)"/)?.[1];
  const content = match[2].replace(/\s*<(?:title|desc)\b[^>]*>[\s\S]*?<\/(?:title|desc)>/g, '');
  return { viewBox, content };
}

function svgDocument(width, height, title, content) {
  return `<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" role="img" aria-labelledby="title">\n  <title id="title">${title}</title>\n${content}\n</svg>\n`;
}

// Convert glyphs to paths so the exported wordmarks need no installed fonts.
function lettering(text, fontSize, fill) {
  const measure = createCanvas(1, 1).getContext('2d');
  measure.font = `500 ${fontSize}px OceanKitBrand`;
  const metrics = measure.measureText(text);
  const width = Math.ceil(metrics.width + 2);
  const canvas = createCanvas(width, Math.ceil(fontSize * 1.8), SvgExportFlag.ConvertTextToPaths);
  const context = canvas.getContext('2d');
  context.font = measure.font;
  context.fillStyle = fill;
  context.fillText(text, 0, fontSize);
  const { content } = svgParts(canvas.getContent().toString());
  const inkCenter = fontSize + (metrics.actualBoundingBoxDescent - metrics.actualBoundingBoxAscent) / 2;
  return { width, content, inkCenter };
}

const rounded = await fs.readFile(path.join(source, 'icon-rounded.svg'), 'utf8');
const square = await fs.readFile(path.join(source, 'icon-square.svg'), 'utf8');
const small = await fs.readFile(path.join(source, 'icon-small.svg'), 'utf8');
const icon = svgParts(rounded);
function positionedIcon(x, y, size) {
  return `  <svg x="${x}" y="${y}" width="${size}" height="${size}" viewBox="${icon.viewBox}">${icon.content}\n  </svg>`;
}
async function png(svg, filename, width, height = width) {
  await sharp(Buffer.from(svg), { density: 144 }).resize(width, height).png().toFile(path.join(destination, filename));
}

for (const filename of ['icon-rounded.svg', 'icon-square.svg', 'icon-small.svg', 'ocean-mark.svg', 'ocean-mark-monochrome.svg']) {
  await fs.copyFile(path.join(source, filename), path.join(destination, filename));
}
await png(rounded, 'icon-512.png', 512);
await png(rounded, 'icon-1024.png', 1024);
await png(square, 'apple-touch-icon.png', 180);
await png(small, 'favicon-96.png', 96);
for (const size of [16, 32, 48]) await png(small, `favicon-${size}.png`, size);

// ICO directory entries contain PNG images at each browser display size.
const icoImages = await Promise.all([16, 32, 48].map(async size => ({
  size, data: await fs.readFile(path.join(destination, `favicon-${size}.png`)),
})));
const icoHeader = Buffer.alloc(6 + 16 * icoImages.length);
icoHeader.writeUInt16LE(1, 2);
icoHeader.writeUInt16LE(icoImages.length, 4);
let offset = icoHeader.length;
for (const [index, image] of icoImages.entries()) {
  const entry = 6 + 16 * index;
  icoHeader[entry] = image.size;
  icoHeader[entry + 1] = image.size;
  icoHeader.writeUInt16LE(1, entry + 4);
  icoHeader.writeUInt16LE(32, entry + 6);
  icoHeader.writeUInt32LE(image.data.length, entry + 8);
  icoHeader.writeUInt32LE(offset, entry + 12);
  offset += image.data.length;
}
await fs.writeFile(path.join(destination, 'favicon.ico'), Buffer.concat([icoHeader, ...icoImages.map(image => image.data)]));

for (const [variant, color] of [['light', '#052f60'], ['dark', '#edf7fa']]) {
  const text = lettering('OceanKit', 50, color);
  const width = text.width + 108;
  const logo = svgDocument(width, 100, 'OceanKit', [
    positionedIcon(0, 10, 80),
    `  <g transform="translate(100 ${50 - text.inkCenter})">${text.content}\n  </g>`,
  ].join('\n'));
  const filename = `logo-horizontal-${variant}`;
  await fs.writeFile(path.join(source, `${filename}.svg`), logo);
  await fs.writeFile(path.join(destination, `${filename}.svg`), logo);
  await png(logo, `${filename}.png`, width * 2, 200);
  console.log(`${filename}: ${width} × 100 SVG; ${width * 2} × 200 PNG`);
}

const heading = lettering('OceanKit', 60, '#edf7fa');
const description = lettering('MATLAB tools for ocean science', 30, '#b0e6e6');
const address = lettering('jeffreyearly.github.io/OceanKit', 23, '#b0e6e6');
const social = svgDocument(1200, 630, 'OceanKit — MATLAB tools for ocean science', [
  '  <rect width="1200" height="630" fill="#052f60"/>',
  positionedIcon(64, 123, 384),
  `  <g transform="translate(484 ${285 - heading.inkCenter})">${heading.content}\n  </g>`,
  `  <g transform="translate(488 ${355 - description.inkCenter})">${description.content}\n  </g>`,
  `  <g transform="translate(488 ${402 - address.inkCenter})">${address.content}\n  </g>`,
].join('\n'));
await fs.writeFile(path.join(source, 'social-card.svg'), social);
await fs.writeFile(path.join(destination, 'social-card.svg'), social);
await png(social, 'social-card.png', 1200, 630);

for (const filename of ['ocean-mark', 'ocean-mark-monochrome']) {
  const svg = await fs.readFile(path.join(source, `${filename}.svg`), 'utf8');
  const bounds = svgParts(svg).viewBox.split(/\s+/).map(Number);
  await png(svg, `${filename}.png`, 1600, Math.round(1600 * bounds[3] / bounds[2]));
}
console.log(`Branding assets exported to ${destination}`);
