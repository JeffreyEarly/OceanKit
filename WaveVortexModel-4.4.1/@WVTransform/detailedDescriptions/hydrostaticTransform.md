Create the corresponding hydrostatic wave-vortex transform.

For a stratified QG transform, this method constructs a `WVTransformHydrostatic` with the same domain, stratification, planetary parameters, time, coefficients, and forcing configuration.

```matlab
wvtHydrostatic = wvtQG.hydrostaticTransform;
```

- Declaration: wvt = hydrostaticTransform()
- Returns wvt: corresponding `WVTransformHydrostatic` instance
