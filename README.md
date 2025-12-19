# OceanKit
A MATLAB package repository for a wide variety of oceanography related tools.

Clone this repository, e.g., execute
```
git clone https://github.com/JeffreyEarly/OceanKit.git
```
from the command-line. Within Matlab, add this folder as an MPM repository,
```matlab
mpmAddRepository("OceanKit","path/to/folder/OceanKit")
```
and then install any package you need, e.g.,
```matlab
mpminstall("WaveVortexModel")
```
will install the [WaveVortexModel](https://wavevortexmodel.org) and its dependencies.

Use
```matlab
mpmsearch(Repository="OceanKit")
```
to list all the available packages and
```matlab
mpmlist
```
to list all the packages you've already installed. To uninstall a package
```matlab
mpmuninstall("WaveVortexModel")
```
or to uninstall *all* packages,
```matlab
pkgs = mpmlist;
mpmuninstall([pkgs.Name])
```

### Advanced
Note that if you intend to make changes to package to be committed back to github, you will need to check out the repo of interest directly. For example, you'd call
```
git clone https://github.com/JeffreyEarly/wave-vortex-model.git
```
to clone the repo, and then directly install package
```matlab
mpminstall("local/path/to/WaveVortexModel", Authoring=true); 
```
with authoring enabled.

# Packages in this repository

- [AlongTrackSimulator](https://github.com/satmapkit/AlongTrackSimulator) A ground track simulator for satellite altimetry missions
- [chebfun](https://github.com/chebfun/chebfun) Chebfun: numerical computing with functions.
- [class-docs](https://github.com/JeffreyEarly/class-annotations) Builds better website documentation for Matlab classes
- [class-annotations](https://github.com/JeffreyEarly/class-annotations) Annotates Matlab class properties and methods to enable writing to NetCDF files and better class documentation
- [distributions](https://github.com/JeffreyEarly/distributions) Class for distributions and stochastic modeling
- [geographic-projection](https://github.com/JeffreyEarly/geographic-projection) Tools for projecting geographic coordinates (latitude, longitude) onto transverse Mercator (x,y).
- [internal-modes](https://github.com/JeffreyEarly/internal-modes) Quickly and accurately compute the vertical modes for arbitrary stratification.
- [netcdf](https://github.com/JeffreyEarly/netcdf) A better NetCDF interface for Matlab
- [spline-core](https://github.com/JeffreyEarly/spline-core) Core b-spline classes.
- [wave-vortex-model](https://github.com/JeffreyEarly/wave-vortex-model) Solve the stratified nonhydrostatic equations of motion and decompose into wave-vortex components
- [wave-vortex-model-diagnostics](https://github.com/Energy-Pathways-Group/wave-vortex-model-diagnostics) Diagnostic tools for the WaveVortexModel
