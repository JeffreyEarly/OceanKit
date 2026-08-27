Negative-frequency wave–vortex coefficient array.

`Am` stores the negative-frequency coefficients $$A_-^{k\ell j}$$ for internal gravity waves and inertial oscillations. The coefficients have units of velocity and use the transform's spectral layout.

These coefficients multiply the negative-frequency wave solutions described by [Early, Lelong, and Sundermeyer (2021)](https://doi.org/10.1017/jfm.2020.995) and the current [available-potential-vorticity formulation](https://doi.org/10.48550/arXiv.2403.20269). The internal-gravity-wave and inertial-oscillation solutions appear as equations (3.18) and (3.15), respectively, in the 2021 paper.

For a three-dimensional wave-bearing transform, the flow constituents occupy `Am` schematically as follows:

| Vertical mode | $$K_h=0$$ | $$K_h>0$$ |
|:---:|:---:|:---:|
| $$j=0$$ | inertial oscillation | — |
| $$j>0$$ | inertial oscillation | internal gravity wave, $$-\omega$$ |

The inertial entries in `Am` are the conjugate partners of the primary inertial coefficients in `Ap`. At nonzero horizontal wavenumber, `Am` contains the negative-frequency member of each internal-gravity-wave pair. The transform's `inertialComponent.maskAm` and `waveComponent.maskAm` properties are the executable definitions of these regions; the table suppresses geometry-specific conjugate storage, excluded Nyquist modes, and antialiasing details.

The stored phase is referenced to `t0`. Linear evolution does not overwrite `Am`; use `Amt` for the coefficients evaluated at the current `t`. The wave and inertial primary-flow-component masks identify the active locations. Coefficients outside those masks must remain zero.

Together `Am` and `Ap` obey the transform's Hermitian and inertial conjugacy relationships so the reconstructed physical fields are real. In particular, the inertial coefficients satisfy `Am = conj(Ap)` on their masks. Quasigeostrophic transforms have no active `Am` content.

- Topic: Wave-vortex coefficients
