function [queryTimes, scoreTimes, scoreIndices] = resolveBootstrapTimes( ...
        trajectories, queryTimesOption, scoreTimesOption, scoreStride)
tStart = zeros(numel(trajectories), 1);
tEnd = zeros(numel(trajectories), 1);
allTimes = cell(numel(trajectories), 1);

for iTrajectory = 1:numel(trajectories)
    ti = reshape(trajectories(iTrajectory).t, [], 1);
    tStart(iTrajectory) = ti(1);
    tEnd(iTrajectory) = ti(end);
    allTimes{iTrajectory} = ti;
end

commonStart = max(tStart);
commonEnd = min(tEnd);
if commonStart > commonEnd
    error("GriddedStreamfunctionBootstrap:EmptyCommonOverlap", ...
        "The original drifter trajectories do not share a common overlap interval.");
end

if isempty(queryTimesOption)
    pooledTimes = unique(vertcat(allTimes{:}), "sorted");
    queryTimes = pooledTimes(pooledTimes >= commonStart & pooledTimes <= commonEnd);
    if isempty(queryTimes)
        error("GriddedStreamfunctionBootstrap:EmptyCommonOverlap", ...
            "No original observation times lie inside the common overlap interval [%.15g, %.15g].", ...
            commonStart, commonEnd);
    end
else
    queryTimes = reshape(queryTimesOption, [], 1);
    if ~(isnumeric(queryTimes) && isreal(queryTimes) && all(isfinite(queryTimes)))
        error("GriddedStreamfunctionBootstrap:InvalidQueryTimes", ...
            "queryTimes must be a finite real vector.");
    end
    if any(diff(queryTimes) <= 0)
        error("GriddedStreamfunctionBootstrap:InvalidQueryTimes", ...
            "queryTimes must be strictly increasing.");
    end
    if any(queryTimes < commonStart | queryTimes > commonEnd)
        error("GriddedStreamfunctionBootstrap:InvalidQueryTimes", ...
            "queryTimes must lie inside the common overlap interval [%.15g, %.15g].", ...
            commonStart, commonEnd);
    end
end

if isempty(scoreTimesOption)
    scoreIndices = 1:scoreStride:numel(queryTimes);
    scoreTimes = queryTimes(scoreIndices);
else
    scoreTimes = reshape(scoreTimesOption, [], 1);
    if ~(isnumeric(scoreTimes) && isreal(scoreTimes) && all(isfinite(scoreTimes)))
        error("GriddedStreamfunctionBootstrap:InvalidScoreTimes", ...
            "scoreTimes must be a finite real vector.");
    end
    if any(diff(scoreTimes) <= 0)
        error("GriddedStreamfunctionBootstrap:InvalidScoreTimes", ...
            "scoreTimes must be strictly increasing.");
    end

    [isSubset, scoreIndices] = ismember(scoreTimes, queryTimes);
    if ~all(isSubset)
        error("GriddedStreamfunctionBootstrap:InvalidScoreTimes", ...
            "scoreTimes must be an exact subset of queryTimes.");
    end
end

queryTimes = reshape(queryTimes, [], 1);
scoreTimes = reshape(scoreTimes, [], 1);
scoreIndices = reshape(scoreIndices, 1, []);
end
