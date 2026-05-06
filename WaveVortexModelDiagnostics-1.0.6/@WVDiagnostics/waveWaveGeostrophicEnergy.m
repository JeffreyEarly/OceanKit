function E0 = waveWaveGeostrophicEnergy(wvt,mask)
% Note that.
%
% Note that
% energy = waveWaveGeostrophicEnergy(wvt,1,1);
% should produce the same answer as
% flow = wvt.flowComponentWithName("wave") + wvt.flowComponentWithName("inertial");
% [Fp_w,Fm_w,F0_w] = wvt.nonlinearFluxForFlowComponents(flow,flow);
% [Ep_w,Em_w,E0_w] = wvt.energyFluxFromNonlinearFlux(Fp_w,Fm_w,F0_w);
% sum(E0_w(:))
%
% - Topic: Diagnostics — General — Misc — Fluxes in space, [sparseJWavenumberAxis sparseKRadialAxis]
% - Declaration: E0 = waveWaveGeostrophicEnergy(wvt,mask)
% - Parameter wvt: WVDiagnostics object
% - Parameter mask: input argument `mask`
% - Returns E0: output value `E0`
arguments
    wvt
    mask
end


    Apt = wvt.Apt;
    Amt = wvt.Amt;

    Upm = wvt.UAp.*Apt + wvt.UAm.*Amt;
    Vpm = wvt.VAp.*Apt + wvt.VAm.*Amt;
    Wpm = wvt.WAp.*Apt + wvt.WAm.*Amt;
    Npm = wvt.NAp.*Apt + wvt.NAm.*Amt;

    [U,Ux,Uy,Uz] = wvt.transformToSpatialDomainWithFAllDerivatives(        Apm=mask.*Upm);
    [V,Vx,Vy,Vz] = wvt.transformToSpatialDomainWithFAllDerivatives(        Apm=mask.*Vpm);
    W =                          wvt.transformToSpatialDomainWithG(        Apm=mask.*Wpm);
    [ETA,ETAx,ETAy,ETAz] = wvt.transformToSpatialDomainWithGAllDerivatives(Apm=mask.*Npm);

    uNL = -U.*Ux - V.*Uy - W.*Uz;
    vNL = -U.*Vx - V.*Vy - W.*Vz;
    nNL = -U.*ETAx - V.*ETAy - W.*(ETAz + ETA.*shiftdim(wvt.dLnN2,-2));

    % 2. Transform the flux into the spectral domain. Stolen from,
    % [Fp,Fm,F0] = wvt.transformUVEtaToWaveVortex(uNL,vNL,nNL);
    u_hat = wvt.transformFromSpatialDomainWithFourier(uNL);
    v_hat = wvt.transformFromSpatialDomainWithFourier(vNL);
    n_hat = wvt.transformFromSpatialDomainWithFourier(nNL);

    iK = sqrt(-1)*shiftdim(wvt.k,-1);
    iL = sqrt(-1)*shiftdim(wvt.l,-1);
    n_bar = wvt.transformFromSpatialDomainWithGg(n_hat);
    zeta_bar = wvt.transformFromSpatialDomainWithFg(iK .* v_hat - iL .* u_hat);
    F0 = wvt.A0Z.*zeta_bar + wvt.A0N.*n_bar;

    % 3. Compute its energy, stolen from wvt.energyFluxFromNonlinearFlux
    E0 = 2*wvt.A0_TE_factor.*real( F0 .* conj(wvt.A0) );
end