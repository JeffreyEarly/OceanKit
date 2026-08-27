Periodic y-coordinate axis in meters.

`y` is an `Ny`-by-1 column vector spanning $$0 \le y < L_y$$. The values `Ly` and `Ny` are set during initialization, and

The y coordinate is periodic, which means that
```matlab
dy = Ly/Ny;
y = dy*(0:Ny-1)';
```

Because the endpoint is omitted, `Ly` equals `y(end)-y(1)+dy`, not `y(end)-y(1)`. This is the standard grid for a periodic Fourier transform.

- Topic: Domain attributes — Grid — Spatial
- nav_order: 2
