Differentiate a gridded field in the periodic y direction.

The derivative is evaluated spectrally. The input and output retain the same spatial layout, and derivative order `n` defaults to `1`.

```matlab
dudy = wvt.diffY(u);
```

- Declaration: derivative = diffY(field,options)
- Parameter field: gridded field with the transform's spatial shape
- Parameter options.n: positive derivative order; default `1`
- Returns derivative: y derivative with the same shape as `field`
