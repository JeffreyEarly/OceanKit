---
layout: default
title: Home
nav_order: 1
description: "Discover and install MATLAB tools for ocean models, numerical methods, and scientific data."
permalink: /
---

<div class="ocean-hero">
  <img class="ocean-hero-mark" src="{{ '/assets/branding/ocean-mark.svg' | relative_url }}" width="144" height="144" alt="OceanKit: a navy O with a teal wave">
  <div>
    <p class="ocean-eyebrow">OceanKit</p>
    <h1 id="matlab-tools-for-ocean-science">MATLAB tools for ocean science.</h1>
    <p class="ocean-lede">Explore waves and vortices, build numerical models, and work with scientific data. A collection of reusable packages, ready to install with MATLAB Package Manager.</p>
  </div>
</div>

<div class="ocean-actions">
  <a class="btn btn-primary" href="{{ '/installation.html' | relative_url }}">Get started</a>
  <a class="btn" href="{{ '/packages' | relative_url }}">Browse packages</a>
</div>

## A toolkit for your research

<div class="ocean-topics">
  <section>
    <h3 id="ocean-dynamics">Ocean dynamics</h3>
    <p>Model rotating, stratified flows, compute internal modes, and estimate advection and diffusion.</p>
    <a href="{{ '/packages' | relative_url }}#ocean-models">Explore ocean models &rarr;</a>
  </section>
  <section>
    <h3 id="numerical-methods">Numerical methods</h3>
    <p>Interpolate with splines, compute Fourier transforms, and model probability distributions.</p>
    <a href="{{ '/packages' | relative_url }}#numerical-tools">Explore numerical tools &rarr;</a>
  </section>
  <section>
    <h3 id="data-and-geography">Data and geography</h3>
    <p>Read and write NetCDF files, project coordinates, and simulate satellite ground tracks.</p>
    <a href="{{ '/packages' | relative_url }}#data-and-geography">Explore data tools &rarr;</a>
  </section>
</div>

## Install the tools you need

Clone the repository once:

```bash
git clone https://github.com/JeffreyEarly/OceanKit.git
```

Register it in MATLAB, then choose a package:

```matlab
mpmAddRepository("OceanKit","path/to/OceanKit")
mpmsearch(Repository="OceanKit")
mpminstall("WaveVortexModel")
```

MATLAB Package Manager is available in R2024b and newer. Individual packages may require a newer release or a supported platform; see [installation requirements]({{ '/installation.html' | relative_url }}#requirements).

## Build with OceanKit

OceanKit distributes released packages and collects shared guidance for their authors. To contribute to a package, work in its own source repository. The [developers guide]({{ '/developers-guide' | relative_url }}) covers package design, MATLAB conventions, documentation, and releases.
