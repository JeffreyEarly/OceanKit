Applies linear bottom friction to the flow, i.e., $$\frac{du}{dt} = -r \cdot u(x,y,-D)$$. The parameter $$r$$ has units of $$s^{-1}$$ and thus can be set as an inverse time scale.

The linear bottom friction is scaled such that we actually apply, $$\frac{du}{dt} = -\frac{L_z}{dz} r \cdot u(x,y,-D)$$ and the volume integrated effect of friction remains the same regardless of resolution. $$L_z$$ is the total domain depth and $$dz$$ is the spacing at the bottom grid point.

To compare with quadratic bottom friction where $$\frac{du}{dt} = -\frac{C_d}{dz} \lvert \mathbf{u} \rvert $$, note that $$- \frac{L_z}{dz} r = -\frac{C_d}{dz} \lvert \mathbf{u} \rvert$$ and you will find a characteristic velocity $$\lvert\mathbf{u}\rvert$$ of about 10 cm/s for $$C_d=0.002$$.

For both nonhydrostatic and hydrostatic transforms linear bottom drag

$$
\begin{align}
\mathcal{S}_u &= -\frac{L_z}{dz} r \cdot u(x,y,-D) \\
\mathcal{S}_v &= -\frac{L_z}{dz} r \cdot v(x,y,-D)  \\
\mathcal{S}_w &= 0 \\
\mathcal{S}_\eta &= 0
\end{align}
$$

and for quasigeostrophic transforms,

$$
\begin{align}
\mathcal{S}_\textrm{qgpv} &= -\frac{L_z}{dz} r \cdot \zeta(x,y,-D)
\end{align}
$$

where $$\zeta = \partial_x v - \partial_y u$$.

### Usage

Assuming there is a WVTransform instance wvt, to add this forcing,

```matlab
wvt.addForcing(WVBottomFrictionLinear(r=1/(200*86400)));
```
