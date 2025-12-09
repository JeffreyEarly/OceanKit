basedir = "/Users/Shared/CimRuns_June2025/output/";
basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";

runNumber=9; runName = "hydrostatic: geostrophic + waves";
wvd = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","256") + ".nc");

%%
wvt = wvd.wvt;

[kp,bins_0,bins_pm] = wvd.sparsePseudoRadialAxis;
[omegaAxis,bins_omega] = wvd.sparseOmegaAxis;

valid = ~isnan(bins_0);
S_0 = sparse(find(valid), bins_0(valid), 1, numel(wvt.Ap), numel(kp), nnz(valid));

valid = ~isnan(bins_pm);
S_pm = sparse(find(valid), bins_pm(valid), 1, numel(wvt.Ap), numel(kp), nnz(valid));

mask_0 = false(wvd.wvt.Nj,wvt.Nkl,length(kp));
mask_pm = false(wvd.wvt.Nj,wvt.Nkl,length(kp));
for iK = 1:1:length(kp)
    mask_0(:,:,iK) = (bins_0 <= iK);
    mask_pm(:,:,iK) = (bins_pm <= iK);
end

%%
E0 = WVDiagnostics.waveWaveGeostrophicEnergy(wvt,1);
F_wwg_kp = reshape(E0(:).' * S_0,[],1);

Epm = WVDiagnostics.geostrophicGeostrophicWaveEnergy(wvt,1);
F_ggw_kp = reshape(Epm(:).' * S_pm,[],1);

%%
% geo_var.setValueAlongDimensionAtIndex(geo_jk,'t',outputIndex);

pi_w_wwg_kp = zeros(length(kp),1);
pi_w_wwg_omega = zeros(length(omegaAxis),1);
pi_g_ggw_kp = zeros(length(kp),1);

% fprintf("%d: iTime=%d. k=",timeIndex,self.iTime);
tic
for i=1:length(kp)
    E0 = WVDiagnostics.waveWaveGeostrophicEnergy(wvt,mask_pm(:,:,i));
    pi_w_wwg_kp(i) = sum(E0(:));
    Epm = WVDiagnostics.geostrophicGeostrophicWaveEnergy(wvt,mask_0(:,:,i));
    pi_g_ggw_kp(i) = sum(Epm(:));
end

mask_omega = false(wvd.wvt.Nj,wvt.Nkl,length(omegaAxis));
for iK = 1:1:length(omegaAxis)
    mask_omega(:,:,iK) = (bins_omega <= iK);
end

for i=1:length(omegaAxis)
    E0 = WVDiagnostics.waveWaveGeostrophicEnergy(wvt,mask_omega(:,:,i));
    pi_w_wwg_omega(i) = sum(E0(:));
end
toc


%%
flow = wvt.flowComponentWithName("wave") + wvt.flowComponentWithName("inertial");
[Fp_w,Fm_w,F0_w] = wvt.nonlinearFluxForFlowComponents(flow,flow);
[Ep_w,Em_w,E0_w] = wvt.energyFluxFromNonlinearFlux(Fp_w,Fm_w,F0_w);
sum(E0_w(:))

flow = wvt.flowComponentWithName("geostrophic") + wvt.flowComponentWithName("mda");
[Fp_w,Fm_w,F0_w] = wvt.nonlinearFluxForFlowComponents(flow,flow);
[Ep_w,Em_w,E0_w] = wvt.energyFluxFromNonlinearFlux(Fp_w,Fm_w,F0_w);
sum(Ep_w(:)+Em_w(:))

%%

% /Applications/MATLAB_R2025b.app/bin/matlab -nojvm -nodisplay -nosplash
% path = "/Users/Shared/CimRuns_June2025/output/run1_icR_iner0_tide0_lat32_geo0065_N0052_hydrostatic_res512.nc";
% path = "/Users/Shared/CimRuns_June2025/output/run9_icR_iner07_tide014_lat32_geo0065_N0052_hydrostatic_res512.nc";
% path = "/Users/Shared/CimRuns_June2025/output/run18_icR_iner07_tide014_lat32_geo0065_N0052_boussinesq_res512.nc";
% wvd = WVDiagnostics(path);
% wvd.create1DMirrorFluxes();

[kp,omegaAxis] = wvd.diagfile.readVariables("kp","omegaAxis");
[pi_w_wwg_kp_t, F_wwg_kp_t, pi_g_ggw_kp_t, F_ggw_kp_t, pi_w_wwg_omega_t] = wvd.diagfile.readVariables("pi_w_wwg_kp", "F_wwg_kp", "pi_g_ggw_kp", "F_ggw_kp", "pi_w_wwg_omega");

timeIndices = 1:51;
pi_w_wwg_kp = mean(pi_w_wwg_kp_t(:,timeIndices),2);
F_wwg_kp = mean(F_wwg_kp_t(:,timeIndices),2);
pi_g_ggw_kp = mean(pi_g_ggw_kp_t(:,timeIndices),2);
F_ggw_kp = mean(F_ggw_kp_t(:,timeIndices),2);
pi_w_wwg_omega = mean(pi_w_wwg_omega_t(:,timeIndices),2);

%%

figure
tiledlayout(2,2);

nexttile
plot(kp,cumsum(F_wwg_kp)/wvd.flux_scale), hold on
plot(kp,pi_w_wwg_kp/wvd.flux_scale)
xscale('log')
title('wwg'), legend('geostrophic','wave'), xlabel('k')

nexttile
plot(kp,cumsum(F_ggw_kp)/wvd.flux_scale), hold on
plot(kp,pi_g_ggw_kp/wvd.flux_scale)
xscale('log')
title('ggw'), legend('wave','geostrophic'), xlabel('k')

nexttile
plot(omegaAxis,pi_w_wwg_omega/wvd.flux_scale), hold on
xscale('log')
title('wwg'), xlabel('omega')