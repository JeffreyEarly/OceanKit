function [factor,isPositive] = positiveFitFactor(metric)
% Factor an explicitly positive physical fitting metric without changing it.
%
% When `isPositive` is true, `factor'*factor` represents the supplied
% Hermitian nonnegative metric. Zero fitting directions remain zero. Any
% negative eigenvalue selects signed physical-dual recovery instead: no
% eigenvalue floor or absolute value changes the physical pairing.
% The factor declares a fitting objective, not an observation/noise norm.
%
% - Topic: Developer topics — Construct fitting factors
% - Parameter metric: finite square real or complex physical metric
% - Returns factor: square fitting factor, or empty for an indefinite/non-Hermitian metric
% - Returns isPositive: whether this exact metric supports a nonnegative fit
% - Developer: true
arguments (Input)
    metric (:,:) double {mustBeFinite}
end
arguments (Output)
    factor (:,:) double
    isPositive (1,1) logical
end
if size(metric,1) ~= size(metric,2)
    error("IMProjection:InvalidMetricShape","A fitting metric must be square.");
end
factor = [];
isPositive = false;
if ~isequal(metric,metric')
    return
end
if isequal(metric,diag(diag(metric)))
    diagonal = diag(metric);
    isPositive = all(diagonal >= 0);
    if isPositive, factor = diag(sqrt(diagonal)); end
    return
end
[vectors,values] = eig(metric,"vector");
isPositive = all(values >= 0);
if isPositive, factor = diag(sqrt(values))*vectors'; end
end
