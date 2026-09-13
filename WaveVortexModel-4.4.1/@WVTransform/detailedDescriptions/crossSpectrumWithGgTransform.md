Compute a real modal cross-spectrum using the G-basis transform.

The result is the normalized real part of the product of the first transformed field and the complex conjugate of the second.

- Declaration: spectrum = crossSpectrumWithGgTransform(firstField,secondField)
- Parameter firstField: first G-space field with shape `[Nx Ny Nz]`
- Parameter secondField: second G-space field with shape `[Nx Ny Nz]`
- Returns spectrum: real cross-spectrum with shape `[Nj Nkl]`
