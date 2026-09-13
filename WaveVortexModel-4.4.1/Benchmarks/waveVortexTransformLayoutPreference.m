function selection = waveVortexTransformLayoutPreference(strategyIds,medians,options)
% Select a strict winner while favoring the current path within a tie band.
arguments
    strategyIds (1,:) string
    medians (1,:) double
    options.currentStrategy (1,1) string = "wv-sorted-linear"
    options.tieTolerance (1,1) double {mustBeNonnegative} = 0.03
end
if numel(strategyIds) ~= numel(medians) || isempty(strategyIds)
    error("WaveVortexBenchmark:SelectionSizeMismatch","Strategy IDs and medians must be nonempty arrays with equal lengths.");
end
valid = isfinite(medians) & medians > 0;
if ~any(valid)
    error("WaveVortexBenchmark:NoValidLayoutStrategy","At least one strategy requires a finite positive median.");
end
validIndices = find(valid);
[strictMedian,relativeIndex] = min(medians(valid));
strictIndex = validIndices(relativeIndex);
currentIndex = find(strategyIds == options.currentStrategy,1);
preferredIndex = strictIndex;
currentRetained = false;
if ~isempty(currentIndex) && valid(currentIndex) && medians(currentIndex) <= (1+options.tieTolerance)*strictMedian
    preferredIndex = currentIndex;
    currentRetained = true;
end
selection = struct("strictFastestStrategy",strategyIds(strictIndex),"strictFastestMedianSeconds",strictMedian,"preferredStrategy",strategyIds(preferredIndex),"preferredMedianSeconds",medians(preferredIndex),"currentStrategy",options.currentStrategy,"currentMedianSeconds",NaN,"currentRelativeToFastest",NaN,"tieTolerance",options.tieTolerance,"currentRetained",currentRetained);
if ~isempty(currentIndex)
    selection.currentMedianSeconds = medians(currentIndex);
    selection.currentRelativeToFastest = medians(currentIndex)/strictMedian;
end
end
