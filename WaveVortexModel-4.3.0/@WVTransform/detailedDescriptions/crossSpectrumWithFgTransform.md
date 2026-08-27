Compute a real modal cross-spectrum using the F-basis transform.

The result is the normalized real part of the product of the first transformed field and the complex conjugate of the second.

- Declaration: spectrum = crossSpectrumWithFgTransform(firstField,secondField)
- Parameter firstField: first F-space field with shape `[Nx Ny Nz]`
- Parameter secondField: second F-space field with shape `[Nx Ny Nz]`
- Returns spectrum: real cross-spectrum with shape `[Nj Nkl]`
