function [xTotalByTrajectory, yTotalByTrajectory, durationByTrajectory] = submesoscaleTotalDisplacements( ...
        sampleData, totalDisplacementWeightsByTrajectory)
nTrajectories = numel(sampleData.tCell);
xTotalByTrajectory = zeros(nTrajectories, 1);
yTotalByTrajectory = zeros(nTrajectories, 1);
durationByTrajectory = zeros(nTrajectories, 1);

sampleStartIndex = 1;
for iTrajectory = 1:nTrajectories
    ti = sampleData.tCell{iTrajectory};
    duration = ti(end) - ti(1);
    if duration <= 0
        error("GriddedStreamfunctionBootstrap:InvalidTrajectoryDuration", ...
            "Each submesoscale trajectory must span a strictly positive time interval.");
    end

    nSamples = numel(ti);
    sampleIndices = sampleStartIndex:(sampleStartIndex + nSamples - 1);
    totalDisplacementWeights = reshape(totalDisplacementWeightsByTrajectory{iTrajectory}, 1, []);
    xTotalByTrajectory(iTrajectory) = totalDisplacementWeights * reshape(sampleData.uSubmesoscaleObserved(sampleIndices), [], 1);
    yTotalByTrajectory(iTrajectory) = totalDisplacementWeights * reshape(sampleData.vSubmesoscaleObserved(sampleIndices), [], 1);
    durationByTrajectory(iTrajectory) = duration;
    sampleStartIndex = sampleStartIndex + nSamples;
end
end
