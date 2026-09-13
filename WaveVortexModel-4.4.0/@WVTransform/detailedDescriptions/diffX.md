Differentiate a gridded field in the periodic x direction.

The derivative is evaluated spectrally. The input and output retain the same spatial layout, and derivative order `n` defaults to `1`.

```matlab
dudx = wvt.diffX(u);
```

- Declaration: derivative = diffX(field,options)
- Parameter field: gridded field with the transform's spatial shape
- Parameter options.n: positive derivative order; default `1`
- Returns derivative: x derivative with the same shape as `field`
