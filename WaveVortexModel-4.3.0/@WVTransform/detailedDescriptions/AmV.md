Projects $$v$$ onto $$A_-$$.

These projection coefficients map the $$v$$ state variable onto $$A_-$$. In the historical notation of [Early et al. (2021)](https://doi.org/10.1017/jfm.2020.995), they are the row 2, column 2 entries of $$S^{-1}$$ for the primary internal-gravity-wave and geostrophic solutions in equation C5.

For $$k^2+l^2>0, j>0$$ this is written as,

$$
\textrm{AmV} \equiv \frac{l \omega + i k f_0}{2 \omega K}
$$

in the manuscript. In code this is computed with,

```matlab
alpha = atan2(L,K);
fOmega = f./omega;
AmV = (1/2)*(sin(alpha)+sqrt(-1)*fOmega.*cos(alpha));
```

There are no $$k^2+l^2>0, j=0$$ wave solutions for a rigid lid,

```matlab
AmV(:,:,1) = 0;
```

The inertial solutions occupy the $$k^2+l^2=0$$ portion of the matrix,

```matlab
AmV(1,1,:) = sqrt(-1)/2;
```

- Topic: Developer — Projection coefficients
- Developer: true
