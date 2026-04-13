function kappaEstimate = kappaEstimateFromFit(fit, totalDisplacementWeightsByTrajectory)
if nargin < 2
    totalDisplacementWeightsByTrajectory = {};
end

sampleData = fit.decompositionSampleData(fit.observedTrajectories);
nTrajectories = numel(sampleData.tCell);
if isempty(totalDisplacementWeightsByTrajectory)
    totalDisplacementWeightsByTrajectory = GriddedStreamfunctionBootstrap.totalDisplacementWeightsForSampleTimes(sampleData.tCell);
elseif numel(totalDisplacementWeightsByTrajectory) ~= nTrajectories
    error("GriddedStreamfunctionBootstrap:InvalidTotalDisplacementWeightCount", ...
        "Expected one total-displacement weight row for each sampled trajectory.");
end

[xTotalByTrajectory, yTotalByTrajectory, durationByTrajectory] = ...
    GriddedStreamfunctionBootstrap.submesoscaleTotalDisplacements(sampleData, totalDisplacementWeightsByTrajectory);
kappaEstimate = mean((xTotalByTrajectory.^2 + yTotalByTrajectory.^2) ./ (4 * durationByTrajectory));
end
