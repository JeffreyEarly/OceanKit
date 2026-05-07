Applies quadratic bottom friction to the flow, i.e., $$\frac{du}{dt} = -\frac{Cd}{dz}\lvert \mathbf{u} \rvert\mathbf{u}(x,y,-D)$$. Cd is unitless, and dz is (approximately) the size of the grid spacing at the bottom boundary.

To compare with linear bottom friction where $$\frac{du}{dt} = -r \cdot u(x,y,-D)$$, note that $$- \frac{L_z}{dz} r = -\frac{C_d}{dz} \lvert \mathbf{u} \rvert$$ and you will find a characteristic velocity $$\lvert \mathbf{u} \rvert$$ of about 11.5 cm/s for Cd=0.002 and r=1/(200 days). If $$C_d=0.001$$, then the damping time scale has to double to 1/(400 days) for and equivalent characteristic velocity.

For barotropic QG we want the scaling to work out similarly, but now have $$-r = -\frac{C_d}{D}\lvert \mathbf{u} \rvert$$ where $$D$$ works out to be $$L_z$$. Thus, we will scale the barotropic QG quadratic bottom drag by 4000 m, to match typical oceanic scales.

Using the notation that

$$
|\mathbf{u}(x,y,-D)| = \sqrt{u^2(x,y,-D) + v^2(x,y,-D)}
$$

is the magnitude of the total velocity at the bottom boundary, both nonhydrostatic and hydrostatic transforms linear bottom drag have the form

$$
\begin{align}
\mathcal{S}_u &= -\frac{C_d}{dz} |\mathbf{u}(x,y,-D)| u(x,y,-D) \\
\mathcal{S}_v &= -\frac{C_d}{dz} |\mathbf{u}(x,y,-D)| v(x,y,-D)  \\
\mathcal{S}_w &= 0 \\
\mathcal{S}_\eta &= 0
\end{align}
$$

and for quasigeostrophic transforms,

$$
\begin{align}
\mathcal{S}_\textrm{qgpv} &= -\frac{C_dd}{dz} \left( \frac{\partial}{\partial x} \left( |\mathbf{u}(x,y,-D)| v(x,y,-D) \right) - \frac{\partial}{\partial y} \left( |\mathbf{u}(x,y,-D)| u(x,y,-D) \right) \right)
\end{align}
$$

### Usage

Assuming there is a WVTransform instance wvt, to add this forcing,

```matlab
wvt.addForcing(WVBottomFrictionQuadratic(Cd=0.001));
```
