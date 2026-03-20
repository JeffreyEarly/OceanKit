function [omegaAxis,bins_omega] = sparseOmegaAxis(self)
% Sparse Omega Axis.
%
% This function is used by create1DMirrorFluxes to create an efficient omega axis.
%
% The bin assignments, `bins_omega`, can be used to create sparse matrices for efficient binning operations.
%
% ```matlab
% valid = ~isnan(bins_0);
% S_0 = sparse(find(valid), bins_0(valid), 1, numel(wvt.Ap), numel(kp), nnz(valid));
% F_wwg_kp_val = reshape(E0(:).' * S_0,[],1);
% ```
%
% - Topic: Transformations — Axes
% - Declaration: [omegaAxis,bins_omega] = sparseOmegaAxis(self)
% - Parameter self: WVDiagnostics object
% - Returns omegaAxis: output value `omegaAxis`
% - Returns bins_omega: output value `bins_omega`
arguments
    self
end

wvt = self.wvt;

omega = self.omegaAxis;
dlog = (log10(2)-log10(1))/4;
omegaAxis = omega(1)*(10.^(log10(1):dlog:log10(omega(end)/omega(1))));

mid    = 0.5*(omegaAxis(1:end-1) + omegaAxis(2:end));
edges  = [-Inf, mid, +Inf];

bins_omega = discretize(wvt.Omega, edges);
end