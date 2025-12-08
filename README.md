# OceanKit-mpm-repo
A MATLAB package repository for a wide variety of oceanography related tools

Clone this repository,
```
git clone https://github.com/JeffreyEarly/OceanKit-mpm-repo.git
```

within Matlab, add this folder as an MPM repository,
```matlab
mpmAddRepository("OceanKit","path/to/folder/OceanKit-mpm-repo.git")
```

and then install whatever you need,
```matlab
mpminstall("distributions")
```

# Packages in this repository

- [class-docs](https://github.com/JeffreyEarly/class-annotations) Builds better website documentation for Matlab classes
- [class-annotations](https://github.com/JeffreyEarly/class-annotations) Annotates Matlab class properties and methods to enable writing to NetCDF files and better class documentation
- [wave-vortex-model-diagnostics](https://github.com/Energy-Pathways-Group/wave-vortex-model-diagnostics) Diagnostic tools for the WaveVortexModel
- [distributions](https://github.com/JeffreyEarly/distributions) Class for distributions and stochastic modeling
- [core-spline](https://github.com/JeffreyEarly/spline-core) Core b-spline classes.
- [geographic-projection](https://github.com/JeffreyEarly/geographic-projection) Tools for projecting geographic coordinates (latitude, longitude) onto transverse Mercator (x,y).
