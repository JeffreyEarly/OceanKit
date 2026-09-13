Periodic x-coordinate axis in meters.

`x` is an `Nx`-by-1 column vector spanning $$0 \le x < L_x$$. The values `Lx` and `Nx` are set during initialization, and

The x coordinate is periodic, which means that
```matlab
dx = Lx/Nx;
x = dx*(0:Nx-1)';
```

Because the endpoint is omitted, `Lx` equals `x(end)-x(1)+dx`, not `x(end)-x(1)`. This is the standard grid for a periodic Fourier transform.

- Topic: Domain attributes — Grid — Spatial
- nav_order: 1
