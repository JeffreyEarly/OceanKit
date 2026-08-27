Project quasigeostrophic potential vorticity onto `A0` coefficients.

The input is a gridded QGPV field. The output uses the transform's compact spectral layout and may be assigned to `A0`.

```matlab
A0 = wvt.transformQGPVToWaveVortex(qgpv);
```

- Declaration: A0 = transformQGPVToWaveVortex(qgpv)
- Parameter qgpv: QGPV field with the transform's spatial shape
- Returns A0: zero-frequency coefficients with the transform's spectral shape
