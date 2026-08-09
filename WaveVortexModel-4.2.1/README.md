# WaveVortexModel

WaveVortexModel represents rotating, stratified Boussinesq flow on an energetically orthogonal basis of internal waves, inertial oscillations, geostrophic motions, and mean-density anomalies. It can decompose a fluid state into wave–vortex coefficients, reconstruct physical fields, diagnose the resulting flow, and integrate linear or nonlinear dynamics.

**Complete documentation:** [wavevortexmodel.org](https://wavevortexmodel.org)

## Installation

WaveVortexModel requires MATLAB R2025b or newer. The recommended installation uses [OceanKit](https://github.com/JeffreyEarly/OceanKit) as an MPM repository:

```matlab
mpmAddRepository("OceanKit","local/path/to/OceanKit")
mpminstall("WaveVortexModel")
```

See the [installation guide](https://wavevortexmodel.org/installation) for complete runtime and authoring instructions.

## Quick start

Create a small constant-stratification transform, initialize one internal wave, inspect its velocity field, and advance the state with analytical linear dynamics:

```matlab
wvt = WVTransformConstantStratification( ...
    [40e3 40e3 1000], [16 16 9], N0=5.2e-3, latitude=45);

[omega,k,l] = wvt.initWithWaveModes( ...
    kMode=1, lMode=0, j=1, phi=0, u=0.05, sign=1);
[u,v,w] = wvt.variableWithName('u','v','w');

model = WVModel(wvt,shouldUseLinearDynamics=true);
model.integrateToTime(600,shouldShowIntegrationDiagnostics=false);
```

The transform stores the decomposed state in the wave–vortex coefficients `Ap`, `Am`, and `A0`. Physical variables such as `u`, `v`, `w`, density, pressure, energy, and potential vorticity are reconstructed from those coefficients when requested.

## Documentation

- [User guide](https://wavevortexmodel.org/users-guide/)
- [Capabilities and limitations](https://wavevortexmodel.org/users-guide/supported-features.html)
- [Transform API](https://wavevortexmodel.org/classes/transforms/)
- [WVModel API](https://wavevortexmodel.org/classes/wvmodel/)
- [Mathematical introduction](https://wavevortexmodel.org/mathematical-introduction/)

## Citation

If WaveVortexModel contributes to your work, cite the software and the scientific formulation relevant to your application:

- Early, J. J., Fabre-Lima, L., Remy, B. J., Wortham, C., & Sundermeyer, M. A. (2025). [WaveVortexModel](https://doi.org/10.5281/zenodo.4037401). Zenodo.
- Early, J. J., Lelong, M.-P., & Sundermeyer, M. A. (2021). [A generalized wave-vortex decomposition for rotating Boussinesq flows with arbitrary stratification](https://doi.org/10.1017/jfm.2020.995). *Journal of Fluid Mechanics, 912*, A32.
- Early, J. J., Hernández-Dueñas, G., Smith, L. M., & Lelong, M.-P. (2024). [Available potential vorticity and the wave-vortex decomposition for arbitrary stratification](https://doi.org/10.48550/arXiv.2403.20269). arXiv.

Additional acknowledgements and BibTeX downloads are available at [wavevortexmodel.org/acknowledgements](https://wavevortexmodel.org/acknowledgements).
