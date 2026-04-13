Fluids paper case study
==============

This folder contains the worked example scripts used to reproduce the linear-velocity-field estimation workflow from Oscroft, Sykulski, and Early (2021).

The scripts are intended to be run from MATLAB by name. Each script resolves its own inputs and outputs relative to this folder, so the current working directory does not matter.

### Folder layout

- `../../../ExampleData/LatMix2011/` contains the local drifter `.mat` inputs shared by the estimation examples.
- `BootstrapData/` contains generated bootstrap-fit outputs written by `GenerateBootstrapFits.m`.
- `Movies/` contains generated `.mp4` outputs written by the movie scripts.
- `LoadFigureDefaults.m` contains the local figure sizing and typography defaults shared by the figure scripts.

### Suggested workflow

1. Run `GenerateBootstrapFits.m` to estimate the bootstrap ensembles and write the `.mat` outputs in `BootstrapData/`.
2. Run `CompareBootstrapFits.m` to compare the candidate models and identify the preferred fit.
3. Run the figure scripts such as `MakeFigureBestSplineFitSite1.m`, `MakeFigureBestSplineFitSite2.m`, or `PlotSite1And2SpectraForManuscript.m`.
4. Run the movie scripts if you want the animated decomposition outputs in `Movies/`.
