Return the three-dimensional spatial coordinate arrays.

`[X,Y,Z] = wvt.xyzGrid` returns arrays of shape `[Nx Ny Nz]` formed with `ndgrid(wvt.x,wvt.y,wvt.z)`.

```matlab
[X,Y,Z] = wvt.xyzGrid;
```

- Declaration: [X,Y,Z] = xyzGrid()
- Returns X: x-coordinate array in meters with shape `[Nx Ny Nz]`
- Returns Y: y-coordinate array in meters with shape `[Nx Ny Nz]`
- Returns Z: vertical-coordinate array in meters with shape `[Nx Ny Nz]`
