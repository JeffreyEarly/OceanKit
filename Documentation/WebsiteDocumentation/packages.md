---
layout: default
title: Packages
nav_order: 2
permalink: /packages
description: "Browse OceanKit packages for ocean models, numerical methods, scientific data, and MATLAB package development."
---

# Packages

Install any of these packages after [registering OceanKit]({{ '/installation.html' | relative_url }}). Use the MPM name shown below, for example:

```matlab
mpminstall("InternalModes")
```

Run `mpmsearch(Repository="OceanKit")` in MATLAB to see available releases. Follow each package's source link for its requirements, examples, and contribution workflow.

## Ocean models

| MPM name | What it does | Learn more |
| --- | --- | --- |
| `WaveVortexModel` | Model rotating, stratified flows using wave and vortex modes. | [Documentation](https://wavevortexmodel.org) · [Source](https://github.com/JeffreyEarly/wave-vortex-model) |
| `WaveVortexModelDiagnostics` | Analyze model output, including energy fluxes from forcing and triads. | [Source](https://github.com/Energy-Pathways-Group/wave-vortex-model-diagnostics) |
| `InternalModes` | Compute vertical modes for arbitrary stratification. | [Source and usage](https://github.com/JeffreyEarly/internal-modes) |
| `AdvectionDiffusionModels` | Integrate and estimate advection-diffusion models. | [Source](https://github.com/JeffreyEarly/advection-diffusion-models) |

## Numerical tools

| MPM name | What it does | Learn more |
| --- | --- | --- |
| `SplineCore` | Interpolate and smooth data with B-spline classes. | [Source](https://github.com/JeffreyEarly/spline-core) |
| `Distributions` | Model probability distributions and stochastic processes. | [Source and usage](https://github.com/JeffreyEarly/distributions) |
| `FFTWTransforms` | Reuse FFTW real-to-complex and real-to-real transforms. | [Source and platform requirements](https://github.com/JeffreyEarly/fftw-transforms) |
| `chebfun` | Compute numerically with functions. | [Documentation](https://www.chebfun.org/docs/guide/) · [Source](https://github.com/chebfun/chebfun) |

## Data and geography

| MPM name | What it does | Learn more |
| --- | --- | --- |
| `NetCDF` | Read and write NetCDF files through a MATLAB class interface. | [Source](https://github.com/JeffreyEarly/netcdf) |
| `GeographicProjection` | Project geographic coordinates onto transverse Mercator coordinates. | [Source](https://github.com/JeffreyEarly/geographic-projection) |
| `AlongTrackSimulator` | Sample numerical data with satellite altimetry ground-track patterns. | [Documentation](https://satmapkit.github.io/AlongTrackSimulator/) · [Source](https://github.com/satmapkit/AlongTrackSimulator) |

## Package infrastructure

| MPM name | What it does | Learn more |
| --- | --- | --- |
| `ClassAnnotations` | Annotate MATLAB properties and methods for persistence and documentation. | [Source](https://github.com/JeffreyEarly/class-annotations) |
| `ClassDocumentation` | Generate class references from MATLAB metadata and structured help. | [Source](https://github.com/JeffreyEarly/class-docs) |

For package authoring and releases, start with the [developers guide]({{ '/developers-guide' | relative_url }}).
