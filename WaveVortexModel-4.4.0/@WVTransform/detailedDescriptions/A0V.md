Projects $$v$$ onto $$A_0$$.

These projection coefficients map the $$v$$ state variable onto $$A_0$$. In the historical notation of [Early et al. (2021)](https://doi.org/10.1017/jfm.2020.995), they are the row 3, column 2 entries of $$S^{-1}$$ for the primary internal-gravity-wave and geostrophic solutions in equation C5.

For $$k^2+l^2>0, j>0$$ (from either equation B14 or C5) this is written as,

$$
\textrm{A0V} \equiv - i \frac{k h f_0}{\omega^2}
$$

in the manuscript. In code this is computed with,

```matlab
fOmega = f./omega;
A0V = -sqrt(-1)*self.h.*(fOmega./omega) .* K;
```

With a rigid lid the solution at $$k>0, l>0, j=0$$ is from equation B11,

$$
\textrm{A0V} \equiv -i \frac{f k}{g K^2}
$$

which in code is,

```matlab
A0V(:,:,1) = -sqrt(-1)*(f/g_)*K(:,:,1)./K2(:,:,1);
```

The $$k=l=0, j>=0$$ solution is a mean density anomaly,

```matlab
A0V(1,1,:) = 0;
```

- Topic: Developer — Projection coefficients
- Developer: true
