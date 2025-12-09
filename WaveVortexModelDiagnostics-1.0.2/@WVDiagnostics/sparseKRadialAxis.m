function k = sparseKRadialAxis(self)
% Sparse KRadial Axis.
%
% sparseKRadialAxis is part of the WVDiagnostics toolbox. Update this description to explain its purpose, inputs, outputs, and how it is used in the overall diagnostics workflow.
%
% - Topic: Transformations — Axes
% - Declaration: k = sparseKRadialAxis(self)
% - Parameter self: WVDiagnostics object
% - Returns k: output value `k`
arguments
    self
end

wvt = self.wvt;

% Create kPseudoRadial *without* explicit de-aliasing
kRadial = reshape(wvt.kRadial,1,[]);    % 1×nK

% this is an expensive calculation so we need to be thoughtful about which
% bins we compute. The idea we use here is that we want to resolve three
% points, log-spaced, in between the damping bin and the maximum bin. So we
% compute that interval, then work backwards.
max_bin = length(kRadial);
mid    = 0.5*(kRadial(1:end-1) + kRadial(2:end));
edges  = [-Inf, mid, +Inf];
if wvt.hasForcingWithName("adaptive damping")
    svv = wvt.forcingWithName("adaptive damping");
    damp_bin = discretize(svv.k_damp,edges);
else
    damp_bin = round(2*max_bin/3);
end

dlog = (log10(max_bin)-log10(damp_bin))/4;
search_bins = cat(2,flip(round(10.^(log10(damp_bin):-dlog:log10(1)))), round(10.^((log10(damp_bin)+dlog):dlog:log10(max_bin))) );
search_bins(1:(find(diff(search_bins) == 0,1,'last')+1)) = [];
search_bins = cat(2,1:(search_bins(1)-1),search_bins);

% We will use this 'sparse' kPseudoRadial axes from here on out
k = reshape(kRadial(search_bins),[],1);

end