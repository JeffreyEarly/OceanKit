Compute a modal autospectrum using the F-basis transform.

The method transforms a gridded F-space field and applies the horizontal Hermitian and vertical normalization factors so summing the result recovers the corresponding quadratic integral.

```matlab
spectrum = wvt.spectrumWithFgTransform(u);
```

- Declaration: spectrum = spectrumWithFgTransform(field)
- Parameter field: F-space field with shape `[Nx Ny Nz]`
- Returns spectrum: nonnegative autospectrum with shape `[Nj Nkl]`
