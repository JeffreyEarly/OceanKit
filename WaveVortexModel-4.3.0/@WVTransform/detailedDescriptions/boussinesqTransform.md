Create the corresponding nonhydrostatic Boussinesq transform.

The new transform preserves the domain, stratification, planetary parameters, time, coefficients, forcing configuration, and no-motion-profile choice while changing the dynamical approximation.

```matlab
wvtBoussinesq = wvtHydrostatic.boussinesqTransform;
```

- Declaration: wvt = boussinesqTransform()
- Returns wvt: corresponding `WVTransformBoussinesq` instance
