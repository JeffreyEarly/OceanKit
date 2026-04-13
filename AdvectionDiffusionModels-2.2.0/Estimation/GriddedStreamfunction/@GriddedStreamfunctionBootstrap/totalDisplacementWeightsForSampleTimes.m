function totalDisplacementWeightsByTrajectory = totalDisplacementWeightsForSampleTimes(tCell)
nTrajectories = numel(tCell);
totalDisplacementWeightsByTrajectory = cell(nTrajectories, 1);
cachedTimes = cell(0, 1);
cachedDegrees = zeros(0, 1);
cachedTotalDisplacementWeights = cell(0, 1);

for iTrajectory = 1:nTrajectories
    ti = reshape(tCell{iTrajectory}, [], 1);
    componentS = min(3, numel(ti) - 1);
    [totalDisplacementWeights, cachedTimes, cachedDegrees, cachedTotalDisplacementWeights] = lookupTotalDisplacementWeights( ...
        ti, componentS, cachedTimes, cachedDegrees, cachedTotalDisplacementWeights);
    totalDisplacementWeightsByTrajectory{iTrajectory} = totalDisplacementWeights;
end
end

function [totalDisplacementWeights, cachedTimes, cachedDegrees, cachedTotalDisplacementWeights] = lookupTotalDisplacementWeights( ...
        t, componentS, cachedTimes, cachedDegrees, cachedTotalDisplacementWeights)
t = reshape(t, [], 1);
for iCache = 1:numel(cachedTimes)
    if cachedDegrees(iCache) == componentS && isequal(cachedTimes{iCache}, t)
        totalDisplacementWeights = cachedTotalDisplacementWeights{iCache};
        return
    end
end

totalDisplacementWeights = computeTotalDisplacementWeights(t, componentS);
cachedTimes{end + 1, 1} = t;
cachedDegrees(end + 1, 1) = componentS;
cachedTotalDisplacementWeights{end + 1, 1} = totalDisplacementWeights;
end

function totalDisplacementWeights = computeTotalDisplacementWeights(t, componentS)
tKnot = BSpline.knotPointsForDataPoints(t, S=componentS);
totalDisplacementWeights = BSpline.integralMatrixForDataPoints(t, t(end), knotPoints=tKnot, S=componentS);
end
