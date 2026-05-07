Fixed amplitude forcing at the natural frequency of each mode

The fixed-amplitude forcing maintains the amplitude of the wave or geostrophic features, while allowing energy and enstrophy to flux from the feature.

As a simple example, one can set an internal wave mode with amplitude 1 cm/s, and that mode will continue to oscillate and maintain its amplitude. The wave will participate in all the nonlinear dynamics, but its amplitude will be maintained/restored at each time step.


There are several different ways to write this style of forcing mathematically. The equations of motion, written in the spectral domain, take the following form

$$
\frac{\partial}{\partial t} A^{klj} = \sum_i F_i^{klj}
$$

where $$F_i$$ are the different forces applied. The transform computes the spatial forcing (which includes nonlinear advection), the spectral forcing, followed by the spectral amplitude forcing. The `WVFixedAmplitudeForcing` is a spectral amplitude forcing and is thus comptued last. This forcing thus simply adds back the flux the from the spatial and spectral forcing, so that $$\frac{\partial}{\partial t} A^{klj} =0$$ for the modes in question.

In practice, of course, we simply restore the amplitudes to their desired value at the last step, e.g.,

```matlab
A0(self.A0_indices) = self.A0bar
```

### Notes

- This approach is commonly used in forced-dissipative turbulence to maintain some fixed forcing.
- Every mode that is used in `WVFixedAmplitudeForcing` essentially removes a degree-of-freedom from the model, as that mode is no longer free to fully evolve. Thus when you pass the forcing wave-vortex coefficients, e.g. `A0`, it does not fix the amplitude of coefficients that are small to avoid removing degrees-of-freedom.
- One must also be careful not to forcing in the damping region. If you have some sort of small scale damping enabled, you probably do not want to be forcing at those smallest scales.

### Usage

To setup a geostrophic mean flow,

```matlab
% initialize a transform
wvt = WVTransformHydrostatic([Lx, Ly, Lz], [Nx, Ny, Nz], N2=@(z) N0*N0*exp(2*z/L_gm),latitude=33);

% set a geostrophic mode, with no flow at the bottom boundary
wvt.setGeostrophicModes(k=0,l=5,j=1,phi=0,u=u0);
wvt.setGeostrophicModes(k=0,l=5,j=0,phi=0,u=max(max(wvt.u(:,:,1))));

% pass the vortex coefficients to the forcing
force = WVFixedAmplitudeForcing(wvt,name="geostrophic-mean-flow");
force.setGeostrophicForcingCoefficients(wvt.A0);
wvt.addForcing(force);

```

In practice you can initialize the flow in any way you want with any arbitrary structure, and then pass those coefficients to the forcing. The `WVFixedAmplitudeForcing` looks for coefficients that are small and ignores those.
