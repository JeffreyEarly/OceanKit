LatMix 2011 sample data
==============

This folder contains the shared LatMix 2011 drifter inputs used by the
estimation examples in this repository.

### Files

- `smoothedGriddedRho1Drifters.mat`
- `smoothedGriddedRho2Drifters.mat`

Each file stores synchronous drifter trajectories with the variables
`t`, `x`, `y`, `lat0`, `lon0`, `n_eff_x`, and `n_eff_y`.

The existing Fluids paper workflows trim the final drifter column before
fitting because it is only a partial time series at both sites:

```matlab
x = x(:, 1:(end-1));
y = y(:, 1:(end-1));
```
