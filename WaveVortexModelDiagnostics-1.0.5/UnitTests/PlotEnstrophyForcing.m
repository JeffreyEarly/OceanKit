% The interesting result here is the damping is doing something wildly
% different for the qgpv flux, while ALL other fluxes line up as you'd
% expect. All I can imagine is that this particular operator is severly
% misdiagnosing the amount of enstrophy at small scales.
wvd = WVDiagnostics("run9_icR_iner07_tide014_lat32_geo0065_N0052_hydrostatic_res256.nc");

[Z_quadratic, t] = wvd.quadraticEnstrophyOverTime();
[Z_apv, ~] = wvd.exactEnstrophyOverTime();

%%
summed_exact = zeros(size(Z_apv));
summed_quad = zeros(size(Z_apv));
enstrophy_fluxes_exact = wvd.exactEnstrophyFluxesOverTime();
enstrophy_fluxes = wvd.quadraticEnstrophyFluxesOverTime();
enstrophy_fluxes_exact(1) = [];
enstrophy_fluxes(1) = [];
for i=1:length(enstrophy_fluxes_exact)
    summed_exact = summed_exact + enstrophy_fluxes_exact(i).Z0;
    summed_quad = summed_quad + enstrophy_fluxes(i).Z0;
end

% ddt_Z_apv = diff(Z_apv)./diff(t);
% ddt_Z_apv = diff(Z_apv)./diff(t);
%%
figure
plot(t/wvd.tscale,(Z_apv-Z_apv(1))/wvd.zscale,LineWidth=2), hold on
plot(t/wvd.tscale,(Z_quadratic-Z_quadratic(1))/wvd.zscale,LineWidth=2)
plot(t/wvd.tscale,cumtrapz(t,summed_exact)/wvd.zscale,LineWidth=2)
plot(t/wvd.tscale,cumtrapz(t,summed_quad)/wvd.zscale,LineWidth=2)
legend("Z_{APV}", "Z_{QGPV}", "\int F_{APV}", "\int F_{QGPV}")

%%
wvd.plotExactEnstrophyFluxOverTime;
wvd.plotEnstrophyFluxOverTime;

%%
figure
for iForce = 1:length(enstrophy_fluxes_exact)
    plot(t/wvd.tscale,enstrophy_fluxes_exact(iForce).Z0/wvd.z_flux_scale,LineWidth=2), hold on
end
legend(enstrophy_fluxes_exact.fancyName)
set(gca,'ColorOrderIndex',1)
for iForce = 1:length(enstrophy_fluxes)
    plot(t/wvd.tscale,enstrophy_fluxes(iForce).Z0/wvd.z_flux_scale,LineWidth=1), hold on
end


xlabel("time (" + wvd.tscale_units + ")")
ylabel("flux (" + wvd.z_flux_scale_units + ")")
xlim([min(t) max(t)]/wvd.tscale);

%%
wvt = wvd.wvt;
F = wvt.fluxForForcing();
enstrophyScale = wvt.f*wvt.f/(86400*365);
forcingNames = wvt.forcingNames;
eta_true = wvt.eta_true;
tl = tiledlayout("flow");
for iForce=1:length(forcingNames)
    if false %isa(wvt,"WVTransformHydrostatic")
        [Fu,Fv,Feta] = wvt.spatialFluxForForcingWithName(forcingNames(iForce));
        DF_x = - wvt.diffZF(Fv); % w_y - v_z
        DF_y = wvt.diffZF(Fu);  % u_z - w_x
        DF_z = wvt.diffX(Fv) - wvt.diffY(Fu);  % v_x - u_y
    else

        [Fu,Fv,Fw,Feta] = wvt.spatialFluxForForcingWithName(forcingNames(iForce));
        DF_x = wvt.diffY(Fw) - wvt.diffZF(Fv); % w_y - v_z
        DF_y = wvt.diffZF(Fu) - wvt.diffX(Fw);  % u_z - w_x
        DF_z = wvt.diffX(Fv) - wvt.diffY(Fu);  % v_x - u_y
    end
    Z = 2*wvt.A0_TZ_factor.*real( F{forcingNames(iForce)}.F0 .* conj(wvt.A0) );
    total = sum(Z(:))/enstrophyScale;
    fprintf("Total quadratic enstrophy for " + forcingNames(iForce) + " forcing: " + total + "\n");
    
    TZ_A0_j_kR = wvt.transformToRadialWavenumber(Z);
    
    % F_pv = wvt.diffX(Fv) - wvt.diffY(Fu) - wvt.f*wvt.diffZG(Feta);
    % Z2 = wvt.qgpv .* F_pv;
    % fprintf("Total approximate enstrophy for " + forcingNames(iForce) + " forcing: " + int_vol(Z2)/enstrophyScale + "\n");
    % F0_pv = wvt.transformFromSpatialDomainWithFg(wvt.transformFromSpatialDomainWithFourier(F_pv));


    % zeta_x = wvt.diffY(wvt.w) - wvt.diffZF(wvt.v); % w_y - v_z
    % zeta_y = wvt.diffZF(wvt.u) - wvt.diffX(wvt.w);  % u_z - w_x
    % zeta_z = wvt.diffX(wvt.v) - wvt.diffY(wvt.u);  % v_x - u_y

    G_eta = (wvt.N2Function(wvt.Z)./wvt.N2Function(wvt.Z - eta_true)).*Feta;
    % G_eta = Feta;
    Z_NL = - wvt.zeta_x .* wvt.diffX(G_eta) - wvt.zeta_y .* wvt.diffY(G_eta)- wvt.zeta_z .* wvt.diffZG(G_eta);
    Z_NL = Z_NL - wvt.diffX(eta_true) .* DF_x - wvt.diffY(eta_true) .* DF_y - wvt.diffZG(eta_true) .* DF_z;
    F_pv = wvt.diffX(Fv) - wvt.diffY(Fu) - wvt.f*wvt.diffZG(G_eta) + Z_NL;
    F_pv_bar = wvt.transformFromSpatialDomainWithFg(wvt.transformFromSpatialDomainWithFourier(F_pv));
    apv_bar = wvt.transformFromSpatialDomainWithFg(wvt.transformFromSpatialDomainWithFourier(wvt.apv));
    prefactor = wvt.h_0/2; prefactor(1) = wvt.Lz/2;
    Z2 = wvt.apv .* F_pv;
    Z2_bar = 2*2*prefactor.*real(conj(apv_bar) .* F_pv_bar);
    TZ_APV_j_kR = wvt.transformToRadialWavenumber(Z2_bar);
    total2 = sum(TZ_APV_j_kR(:))/enstrophyScale;
    fprintf("Total nonlinear enstrophy for " + forcingNames(iForce) + " forcing: " + total2 + " ("+ int_vol(Z2)/enstrophyScale + ")\n");
end