Create an explicit-antialiasing transform with the same grid.

This method converts an implicitly dealiased transform into a full-grid transform with `WVAntialiasing` attached. It preserves time, coefficients, and compatible forcing.

```matlab
wvtExplicit = wvt.waveVortexTransformWithExplicitAntialiasing;
```

- Declaration: wvtExplicit = waveVortexTransformWithExplicitAntialiasing()
- Returns wvtExplicit: transform using an explicit antialias filter
