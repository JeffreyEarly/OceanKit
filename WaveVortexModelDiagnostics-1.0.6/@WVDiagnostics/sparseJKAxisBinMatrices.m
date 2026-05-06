function [S_0, S_pm, mask_0, mask_pm] = sparseJKAxisBinMatrices(self)
% Bin matrices for sparseKRadialAxis and sparseJWavenumberAxis
%
% Used by create2DMirrorFluxes to create sparse binning matrices for the (j,k) axes.
%
% - Topic: Transformations — Axes
% - Declaration: [S_0, S_pm, mask_0, mask_pm] = sparseJKAxisBinMatrices(self)
% - Parameter self: WVDiagnostics object
% - Returns S_0: output value `S_0`
% - Returns S_pm: output value `S_pm`
% - Returns mask_0: output value `mask_0`
% - Returns mask_pm: output value `mask_pm`
arguments
    self
end

wvt = self.wvt;

ks = self.sparseKRadialAxis;
js = self.sparseJWavenumberAxis;
N = length(js)*length(ks);

k_edges  = [-Inf; 0.5*(ks(1:end-1) + ks(2:end)); +Inf];
j_edges  = [-Inf; 0.5*(js(1:end-1) + js(2:end)); +Inf];

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
mask_0 = false(wvt.Nj,wvt.Nkl,N);
mask_pm = false(wvt.Nj,wvt.Nkl,N);
for i=1:length(K_ll(:))
    bins_0(Kh >= K_ll(i) & Kh < K_ur(i) & J_0 >= J_ll(i) & J_0 < J_ur(i)) = i;
    bins_pm(Kh >= K_ll(i) & Kh < K_ur(i) & J_pm >= J_ll(i) & J_pm < J_ur(i)) = i;
    mask_0(:,:,i) = Kh < K_ur(i) & J_0 < J_ur(i);
    mask_pm(:,:,i) = Kh < K_ur(i) & J_0 < J_ur(i);
end

valid = ~isnan(bins_0);
S_0 = sparse(find(valid), bins_0(valid), 1, numel(wvt.Ap), N, nnz(valid));

valid = ~isnan(bins_pm);
S_pm = sparse(find(valid), bins_pm(valid), 1, numel(wvt.Ap), N, nnz(valid));

end