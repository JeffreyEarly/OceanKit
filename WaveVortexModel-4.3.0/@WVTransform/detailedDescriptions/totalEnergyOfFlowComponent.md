Compute the energy carried by one flow component.

The calculation applies the component's `maskAp`, `maskAm`, and `maskA0` to the transform coefficients and sums the corresponding energy factors.

```matlab
waveEnergy = wvt.totalEnergyOfFlowComponent(wvt.waveComponent);
geostrophicEnergy = wvt.totalEnergyOfFlowComponent(wvt.geostrophicComponent);
```

- Declaration: energy = totalEnergyOfFlowComponent(flowComponent)
- Parameter flowComponent: component whose active coefficient masks select the energy
- Returns energy: horizontally averaged, depth-integrated energy per unit reference density
