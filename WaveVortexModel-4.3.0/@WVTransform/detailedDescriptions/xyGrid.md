Return the two-dimensional spatial coordinate arrays.

For a barotropic transform, `[X,Y] = wvt.xyGrid` returns arrays of shape `[Nx Ny]` formed with `ndgrid(wvt.x,wvt.y)`.

```matlab
[X,Y] = wvt.xyGrid;
```

- Declaration: [X,Y] = xyGrid()
- Returns X: x-coordinate array in meters with shape `[Nx Ny]`
- Returns Y: y-coordinate array in meters with shape `[Nx Ny]`
