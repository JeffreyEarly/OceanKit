---
layout: default
title: Installation
nav_order: 2
---

# Installation

Clone the OceanKit repository and register it with MATLAB Package Manager.

```bash
git clone https://github.com/JeffreyEarly/OceanKit.git
```

```matlab
mpmAddRepository("OceanKit","path/to/OceanKit")
```

Once the repository is registered, search or install packages from it:

```matlab
mpmsearch(Repository="OceanKit")
mpminstall("WaveVortexModel")
```

## Common MPM Commands

Use these commands when working with the OceanKit repository:

```matlab
mpmsearch(Repository="OceanKit")
mpmlist
mpmuninstall("WaveVortexModel")
```

To uninstall every installed package:

```matlab
pkgs = mpmlist;
mpmuninstall([pkgs.Name])
```

## Authoring Workflow

OceanKit is the distribution repository. If you intend to edit and commit changes to a package, clone the package's own Git repository directly and install that package in authoring mode:

```bash
git clone https://github.com/JeffreyEarly/wave-vortex-model.git
```

```matlab
mpminstall("/full/path/to/wave-vortex-model", Authoring=true)
```

That keeps your working copy in the real authoring repository while still allowing MPM to resolve dependencies through OceanKit.
