Reconstructs density displacement from $$A_-$$.

These reconstruction coefficients map $$A_-$$ onto the density-displacement state variable. In the historical notation of [Early et al. (2021)](https://doi.org/10.1017/jfm.2020.995), they are the row 3, column 2 entries of $$S$$ for the primary internal-gravity-wave and geostrophic solutions in equation C4.

For $$k^2+l^2>0, j>0$$ this is written as,

$$
\textrm{NAm} \equiv \frac{k h}{\omega}
$$

in the manuscript. In code this is computed with,

```matlab
NAm = Kh.*self.h./omega;
```

There are no $$k^2+l^2>0, j=0$$ wave solutions for a rigid lid,

```matlab
NAm(:,:,1) = 0;
```

The inertial solutions occupy the $$k^2+l^2=0$$ portion of the matrix,

```matlab
NAm(1,1,:) = 0;
```

- Topic: Developer — Reconstruction coefficients
- Developer: true
