# OceanKit-mpm-repo
A MATLAB package repository for a wide variety of oceanography related tools.

Clone this repository, e.g., execute
```
git clone https://github.com/JeffreyEarly/OceanKit-mpm-repo.git
```
from the command-line. Within Matlab, add this folder as an MPM repository,
```matlab
mpmAddRepository("OceanKit","path/to/folder/OceanKit-mpm-repo.git")
```
and then install any package you need,
```matlab
mpminstall("distributions")
```
from the list below. The dependent packages will automatically be installed.

Note that if you intend to make changes to package to be committed back to github, you will need to check out the repo of interest directly. For example, you'd call
```
git clone https://github.com/JeffreyEarly/distributions
```
to clone the repo, and then directly install package
```matlab
mpminstall("local/path/to/distributions", Authoring=true); 
```
with authoring enabled.

# Packages in this repository

- [class-docs](https://github.com/JeffreyEarly/class-annotations) Builds better website documentation for Matlab classes
- [class-annotations](https://github.com/JeffreyEarly/class-annotations) Annotates Matlab class properties and methods to enable writing to NetCDF files and better class documentation
- [core-spline](https://github.com/JeffreyEarly/spline-core) Core b-spline classes.
- [distributions](https://github.com/JeffreyEarly/distributions) Class for distributions and stochastic modeling
- [geographic-projection](https://github.com/JeffreyEarly/geographic-projection) Tools for projecting geographic coordinates (latitude, longitude) onto transverse Mercator (x,y).
- [wave-vortex-model-diagnostics](https://github.com/Energy-Pathways-Group/wave-vortex-model-diagnostics) Diagnostic tools for the WaveVortexModel
