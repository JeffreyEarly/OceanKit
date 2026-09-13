Return spectral-coordinate arrays in wave-vortex layout.

`[K,L,J] = wvt.kljGrid` returns arrays with `wvt.spectralMatrixSize`. `K` and `L` contain angular wavenumbers in radians per meter, and `J` contains dimensionless vertical-mode indices.

```matlab
[K,L,J] = wvt.kljGrid;
```

- Declaration: [K,L,J] = kljGrid()
- Returns K: x-direction angular wavenumber array
- Returns L: y-direction angular wavenumber array
- Returns J: vertical-mode index array
