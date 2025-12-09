basedir = "/Users/Shared/CimRuns_June2025/output/";
basedir = "/Users/jearly/Dropbox/CimRuns_June2025/output/";

runNumber=18; runName = "hydrostatic: geostrophic + waves";
wvd = WVDiagnostics(basedir + replace(getRunParameters(runNumber),"256","512") + ".nc");

%%
wvt = wvd.wvt;

k = wvd.sparseKRadialAxis;
j = wvd.sparseJWavenumberAxis;

k_edges  = [-Inf; 0.5*(k(1:end-1) + k(2:end)); +Inf];
j_edges  = [-Inf; 0.5*(j(1:end-1) + j(2:end)); +Inf];

% [K,J] = ndgrid(k,j);
% figure
% scatter(K(:),J(:),'filled')
% xscale('log'), yscale('log')


[J_ll,K_ll] = ndgrid(j_edges(1:end-1),k_edges(1:end-1));
[J_ur,K_ur] = ndgrid(j_edges(2:end), k_edges(2:end));

jWavenumber = 1./sqrt(wvt.Lr2);
jWavenumber(1) = 0;
J_0 = repmat(jWavenumber,[1 wvt.Nkl]);
Kh = wvt.Kh;
J_pm = 1./sqrt(wvt.g*wvt.h_pm/wvt.f/wvt.f);
J_pm(1,:) = 0; % barotropic mode is a mean?

bins_0 = zeros(wvt.spectralMatrixSize);
bins_pm = zeros(wvt.spectralMatrixSize);
for i=1:length(K_ll(:))
    bins_0(Kh >= K_ll(i) & Kh < K_ur(i) & J_0 >= J_ll(i) & J_0 < J_ur(i)) = i;
    bins_pm(Kh >= K_ll(i) & Kh < K_ur(i) & J_pm >= J_ll(i) & J_pm < J_ur(i)) = i;
end

valid = ~isnan(bins_0);
S_0 = sparse(find(valid), bins_0(valid), 1, numel(wvt.Ap), numel(K), nnz(valid));

valid = ~isnan(bins_pm);
S_pm = sparse(find(valid), bins_pm(valid), 1, numel(wvt.Ap), numel(K), nnz(valid));

mask_0 = false(wvd.wvt.Nj,wvt.Nkl,numel(K));
mask_pm = false(wvd.wvt.Nj,wvt.Nkl,numel(K));
for iK = 1:1:numel(K)
    mask_0(:,:,iK) = (bins_0 <= iK);
    mask_pm(:,:,iK) = (bins_pm <= iK);
end

%%
E0 = WVDiagnostics.waveWaveGeostrophicEnergy(wvt,1);
F_wwg_kp = reshape(reshape(E0(:).' * S_0,[],1),size(K));

Epm = WVDiagnostics.geostrophicGeostrophicWaveEnergy(wvt,1);
F_ggw_kp = reshape(reshape(Epm(:).' * S_pm,[],1),size(K));


%%
% geo_var.setValueAlongDimensionAtIndex(geo_jk,'t',outputIndex);

pi_w_wwg_js_ks = zeros(size(K));
pi_g_ggw_js_ks = zeros(size(K));

% fprintf("%d: iTime=%d. k=",timeIndex,self.iTime);
tic
for i=1:numel(K)
    E0 = WVDiagnostics.waveWaveGeostrophicEnergy(wvt,mask_pm(:,:,i));
    pi_w_wwg_js_ks(i) = sum(E0(:));
    Epm = WVDiagnostics.geostrophicGeostrophicWaveEnergy(wvt,mask_0(:,:,i));
    pi_g_ggw_js_ks(i) = sum(Epm(:));
end

toc

%%
[M_wwg, F_wwg, ks, js] = wvd.quadraticEnergyMirrorTriadFluxes2D(timeIndices=51:251);

figure
pcolor(ks,js,mean(F_wwg,3)), shading flat, xscale('log'), yscale('log')
clim(0.2*[-1 1]*max(max(abs(mean(F_wwg,3)))))