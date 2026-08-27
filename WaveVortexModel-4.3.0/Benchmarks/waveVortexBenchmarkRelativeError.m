function relativeError = waveVortexBenchmarkRelativeError(referenceOutputs,candidateOutputs)
% Return the maximum relative infinity error across operation outputs.
arguments
    referenceOutputs (1,:) cell
    candidateOutputs (1,:) cell
end
if numel(referenceOutputs) ~= numel(candidateOutputs)
    relativeError = Inf;
    return
end
relativeError = 0;
for iOutput = 1:numel(referenceOutputs)
    reference = referenceOutputs{iOutput};
    candidate = candidateOutputs{iOutput};
    if ~isequal(size(reference),size(candidate))
        relativeError = Inf;
        return
    end
    denominator = max(norm(reference(:),Inf),realmin("double"));
    relativeError = max(relativeError,norm(candidate(:)-reference(:),Inf)/denominator);
end
end
