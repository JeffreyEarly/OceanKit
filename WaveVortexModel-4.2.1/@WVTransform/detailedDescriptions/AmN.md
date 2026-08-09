These projection coefficients map the density-displacement state variable onto $$A_-$$. In the historical notation of [Early et al. (2021)](https://doi.org/10.1017/jfm.2020.995), they are the row 2, column 3 entries of $$S^{-1}$$ for the primary internal-gravity-wave and geostrophic solutions in equation C5.

For $$k^2+l^2>0, j>0$$ this is written as,

$$
\textrm{ApN} \equiv - \frac{g K}{2 \omega}
$$

in the manuscript. In code this is computed with,

```matlab
Kh = sqrt(K.*K + L.*L);
ApN = -g*Kh./(2*omega);
```

There are no $$k^2+l^2>0, j=0$$ wave solutions for a rigid lid,

```matlab
ApN(:,:,1) = 0;
```

The inertial solutions at $$k^2+l^2=0$$ do not contribute to $$\eta$$, so that component remains zero.

- Topic: Developer — Projection coefficients
- Developer: true
