function [kePeAxis,bins_kepe] = sparseKePeAxis(self)
% Sparse kenetic-potential energy ratio Axis.
%
% This function is used by create1DMirrorFluxes to create an efficient ke-pe axis.
%
% The bin assignments, `bins_kepe`, can be used to create sparse matrices for efficient binning operations.
%
% ```matlab
% valid = ~isnan(bins_0);
% S_0 = sparse(find(valid), bins_0(valid), 1, numel(wvt.Ap), numel(kp), nnz(valid));
% F_wwg_kp_val = reshape(E0(:).' * S_0,[],1);
% ```
%
% - Topic: Transformations — Axes
% - Declaration: [kePeAxis,bins_kepe] = sparseKePeAxis(self)
% - Parameter self: WVDiagnostics object
% - Returns kePeAxis: output value `kePeAxis`
% - Returns bins_kepe: output value `bins_kepe`
arguments
    self
end

wvt = self.wvt;

kePeAxis = reshape(self.kePeAxis,1,[]);
kePeFraction = wvt.A0_KE_factor./(wvt.A0_KE_factor+wvt.A0_PE_factor);

mid    = 0.5*(kePeAxis(1:end-1) + kePeAxis(2:end));
edges  = [-Inf, mid, +Inf];

bins_kepe = discretize(kePeFraction, edges);
end