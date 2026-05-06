function [kp,bins_0,bins_pm] = sparsePseudoRadialAxis(self)
% Create Sparse Pseudo-Radial Wavenumber Axis and Bin Assignments.
%
% This function is used by create1DMirrorFluxes to create an efficient pseudo-radial wavenumber axis.
%
% The bin assignments, e.g. `bins_0` and `bins_pm`, can be used to create sparse matrices for efficient binning operations.
%
% ```matlab
% valid = ~isnan(bins_0);
% S_0 = sparse(find(valid), bins_0(valid), 1, numel(wvt.Ap), numel(kp), nnz(valid));
% F_wwg_kp_val = reshape(E0(:).' * S_0,[],1);
% ```
%
% - Topic: Transformations — Axes
% - Declaration: [kp,bins_0,bins_pm] = sparsePseudoRadialAxis(self)
% - Parameter self: WVDiagnostics object
% - Returns kp: output value `kp`
% - Returns bins_0: output value `bins_0`
% - Returns bins_pm: output value `bins_pm`
arguments
    self
end

wvt = self.wvt;

% Create kPseudoRadial *without* explicit de-aliasing
j2 = 1./wvt.Lr2;
j2(1) = 0; % barotropic mode is a mean?
[kj,kr] = ndgrid(sqrt(j2),wvt.kRadial);
Kh = sqrt(kj.^2 + kr.^2);
allKs = unique(reshape(abs(Kh),[],1),'sorted');
deltaK = max(diff(allKs));
kAxis_ = 0:deltaK:(max(allKs)+deltaK/2);
kPseudoRadial  = reshape(kAxis_,1,[]);    % 1×nK

% this is an expensive calculation so we need to be thoughtful about which
% bins we compute. The idea we use here is that we want to resolve three
% points, log-spaced, in between the damping bin and the maximum bin. So we
% compute that interval, then work backwards.
max_bin = length(kPseudoRadial);
mid    = 0.5*(kPseudoRadial(1:end-1) + kPseudoRadial(2:end));
edges  = [-Inf, mid, +Inf];
if wvt.hasForcingWithName("adaptive damping")
    svv = wvt.forcingWithName("adaptive damping");
    damp_bin = discretize(sqrt(j2(round(svv.j_damp))+svv.k_damp^2),edges);
else
    damp_bin = round(2*max_bin/3);
end

dlog = (log10(max_bin)-log10(damp_bin))/4;
search_bins = cat(2,flip(round(10.^(log10(damp_bin):-dlog:log10(1)))), round(10.^((log10(damp_bin)+dlog):dlog:log10(max_bin))) );
search_bins(1:(find(diff(search_bins) == 0,1,'last')+1)) = [];
search_bins = cat(2,1:(search_bins(1)-1),search_bins);

%% We will use this 'sparse' kPseudoRadial axes from here on out
kp = kPseudoRadial(search_bins);

% Assign each Kh to a radial bin---geostrophic part
Kh = sqrt(j2 + wvt.Kh.^2);
mid    = 0.5*(kp(1:end-1) + kp(2:end));
edges  = [-Inf, mid, +Inf];
bins_0 = discretize(Kh, edges); % 1..nK or NaN for out-of-range

% Assign each Kh to a radial bin---wave part part (which is different for
% non-hydrostatics)
j2 = 1./(wvt.g*wvt.h_pm/wvt.f/wvt.f);
j2(1,:) = 0; % barotropic mode is a mean?
Kh = sqrt(j2 + wvt.Kh.^2);
mid    = 0.5*(kp(1:end-1) + kp(2:end));
edges  = [-Inf, mid, +Inf];
bins_pm = discretize(Kh, edges); % 1..nK or NaN for out-of-range
end