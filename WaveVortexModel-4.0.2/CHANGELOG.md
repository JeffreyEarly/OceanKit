# Version History

## [4.0.2] - 2026-01-16
- new moment based algorithm for rho_nm

## [4.0.1] - 2025-12-11
- initial mpm ci release

## [4.0.0](https://github.com/Energy-Pathways-Group/GLOceanKit/releases/tag/v4.0) - 2025-08-05
- Added support for `WVForcing', a mechanism for adding arbitrary forcing to the model which also also the nonlinear fluxes to be automatically diagnosed.
- Added `WVObservingSystems', a mechanism for adding user defined observing systems to the model, such as drifters, mooring, along-track altimetry.

## [3.0.0](https://github.com/Energy-Pathways-Group/GLOceanKit/releases/tag/v3.0.0) - 2024-07-25
- Full implemented non-hydrostatic forwards model, passing all unit tests.

## [2.0.0](https://github.com/Energy-Pathways-Group/GLOceanKit/commit/fffd3e6d7822c3844349ec6d65bd262f36f93122) - 2021-11-05
- Fully non-linear forwards model for hydrostatics.
- Never formally declared version 2.0, but on this date we successfully modeled the Cyprus eddy example from the JFM directly in the WaveVortexModel.

## [1.0.0](https://github.com/Energy-Pathways-Group/GLOceanKit/releases/tag/v1.0.0) - 2020-09-18
- All units test passed and declared version 1.0 for the JFM paper.

## [0.9.0](https://github.com/Energy-Pathways-Group/GLOceanKit/commits/master/Matlab/InternalWaveModel/InternalWaveModel.m) - 2018-01-04
- Model expanded to project any state of the fluid onto the complete wave-vortex basis, and time step the state forward linearly. The code was moved to Matlab.

## [0.0.1](https://github.com/Energy-Pathways-Group/GLOceanKit/commits/master/GLOceanKit/GLInternalWaveInitialization.m) - 2014-12-10
- Initially released as a linear internal wave model that time stepped forward non-hydrostatic internal gravity waves in variable stratification using internal wave modes. The code was written in C and Objective-C.

