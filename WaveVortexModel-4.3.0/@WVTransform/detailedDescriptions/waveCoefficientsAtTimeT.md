Return positive- and negative-frequency coefficients at the current time.

The method winds the stored coefficients from reference time `t0` to `t` using `Omega` and leaves `Ap` and `Am` unchanged.

```matlab
[Apt,Amt] = wvt.waveCoefficientsAtTimeT;
```

- Declaration: [Apt,Amt] = waveCoefficientsAtTimeT()
- Returns Apt: positive-frequency coefficients at `t`, with spectral shape
- Returns Amt: negative-frequency coefficients at `t`, with spectral shape
