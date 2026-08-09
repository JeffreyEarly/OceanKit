`Ap` stores the positive-frequency coefficients $$A_+^{k\ell j}$$ for internal gravity waves and the positive-frequency member of the paired inertial representation. The coefficients have units of velocity and use the transform's spectral layout.

These coefficients multiply the positive-frequency wave solutions described by [Early, Lelong, and Sundermeyer (2021)](https://doi.org/10.1017/jfm.2020.995) and the current [available-potential-vorticity formulation](https://doi.org/10.48550/arXiv.2403.20269). The internal-gravity-wave and inertial-oscillation solutions appear as equations (3.18) and (3.15), respectively, in the 2021 paper.

For a three-dimensional wave-bearing transform, the flow constituents occupy `Ap` schematically as follows:

| Vertical mode | $$K_h=0$$ | $$K_h>0$$ |
|:---:|:---:|:---:|
| $$j=0$$ | inertial oscillation | — |
| $$j>0$$ | inertial oscillation | internal gravity wave, $$+\omega$$ |

The inertial entries in `Ap` are the primary members of conjugate pairs whose partners are stored in `Am`. At nonzero horizontal wavenumber, `Ap` contains the positive-frequency member of each internal-gravity-wave pair. The transform's `inertialComponent.maskAp` and `waveComponent.maskAp` properties are the executable definitions of these regions; the table suppresses geometry-specific conjugate storage, excluded Nyquist modes, and antialiasing details.

The stored phase is referenced to `t0`. Linear evolution does not overwrite `Ap`; use `Apt` for the coefficients evaluated at the current `t`. The wave and inertial primary-flow-component masks identify the active locations. Coefficients outside those masks must remain zero.

Together `Ap` and `Am` obey the transform's Hermitian and inertial conjugacy relationships so the reconstructed physical fields are real. Quasigeostrophic transforms have no active `Ap` content.

- Topic: Wave-vortex coefficients
