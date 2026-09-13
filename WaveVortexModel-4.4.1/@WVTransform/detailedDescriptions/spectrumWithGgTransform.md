Compute a modal autospectrum using the G-basis transform.

The method transforms a gridded G-space field and applies the horizontal Hermitian and vertical normalization factors.

```matlab
spectrum = wvt.spectrumWithGgTransform(eta);
```

- Declaration: spectrum = spectrumWithGgTransform(field)
- Parameter field: G-space field with shape `[Nx Ny Nz]`
- Returns spectrum: nonnegative autospectrum with shape `[Nj Nkl]`
