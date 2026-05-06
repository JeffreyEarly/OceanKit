function Epm = geostrophicGeostrophicWaveEnergy(wvt,mask)
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
% - Declaration: Epm = geostrophicGeostrophicWaveEnergy(wvt,mask)
% - Parameter wvt: WVDiagnostics object
% - Parameter mask: input argument `mask`
% - Returns Epm: output value `Epm`
arguments
    wvt
    mask
end


    A0t = wvt.A0t;

    U0 = wvt.UA0.*A0t;
    V0 = wvt.VA0.*A0t;
    N0 = wvt.NA0.*A0t;

    U = wvt.transformToSpatialDomainWithF(A0=mask.*U0);
    V = wvt.transformToSpatialDomainWithF(A0=mask.*V0);
    ETA = wvt.transformToSpatialDomainWithG(A0=mask.*N0);

    uNL = -U.*wvt.diffX(U) - V.*wvt.diffY(U);
    vNL = -U.*wvt.diffX(V) - V.*wvt.diffY(V);
    nNL = -U.*wvt.diffX(ETA) - V.*wvt.diffY(ETA);

    if isa(wvt,'WVTransformConstantStratification')
        [Fp,Fm,~] = wvt.transformUVEtaToWaveVortex(uNL,vNL,nNL);
    elseif isa(wvt,'WVTransformHydrostatic')
        [Fp,Fm,~] = wvt.transformUVEtaToWaveVortex(uNL,vNL,nNL);
    elseif isa(wvt,'WVTransformBoussinesq')
        [Fp,Fm,~] = wvt.transformUVWEtaToWaveVortex(uNL,vNL,zeros(wvt.spatialMatrixSize),nNL);
    else
        error('WVTransform not recognized.')
    end

    % 3. Compute its energy, stolen from wvt.energyFluxFromNonlinearFlux
    Ep = 2*wvt.Apm_TE_factor.*real( Fp .* conj(wvt.Ap) );
    Em = 2*wvt.Apm_TE_factor.*real( Fm .* conj(wvt.Am) );
    Epm = Ep+Em;
end