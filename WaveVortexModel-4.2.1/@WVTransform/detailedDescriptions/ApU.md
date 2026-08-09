These projection coefficients map the $$u$$ state variable onto $$A_+$$. In the historical notation of [Early et al. (2021)](https://doi.org/10.1017/jfm.2020.995), they are the row 1, column 1 entries of $$S^{-1}$$ for the primary internal-gravity-wave and geostrophic solutions in equation C5.

For $$k^2+l^2>0, j>0$$ this is written as,

$$
\textrm{ApU} \equiv \frac{k \omega + i l f_0}{2 \omega K}
$$

in the manuscript. In code this is computed with,

```matlab
alpha = atan2(L,K);
fOmega = f./omega;
ApU = (1/2)*(cos(alpha)+sqrt(-1)*fOmega.*sin(alpha));
```

There are no $$k^2+l^2>0, j=0$$ wave solutions for a rigid lid,

```matlab
ApU(:,:,1) = 0;
```

The inertial solutions occupy the $$k^2+l^2=0$$ portion of the matrix,

```matlab
ApU(1,1,:) = 1/2;
```

- Topic: Developer — Projection coefficients
- Developer: true
