function dLevels = DensityLevelForCDF(gridVectors, density, pctTarget)
% DensityLevelForCDF Return contour levels that enclose target density mass.
%
% `gridVectors` should be the cell-array grid returned by
% `KernelDensityEstimate.densityOnGrid(...)`, with `density` stored in the
% same `ndgrid` ordering. `pctTarget` is a vector of target enclosed-mass
% fractions on `[0, 1]`.

arguments
    gridVectors (1,2) cell
    density (:,:) double {mustBeReal, mustBeFinite, mustBeNonempty}
    pctTarget (:,1) double {mustBeReal, mustBeFinite, mustBeGreaterThanOrEqual(pctTarget, 0), mustBeLessThanOrEqual(pctTarget, 1)}
end

xVector = reshape(gridVectors{1}, [], 1);
yVector = reshape(gridVectors{2}, [], 1);
if ~isnumeric(xVector) || ~isnumeric(yVector) || any(~isfinite(xVector)) || any(~isfinite(yVector))
    error("DensityLevelForCDF:InvalidGridVectors", ...
        "gridVectors must contain finite numeric coordinate vectors.");
end
if ~isequal(size(density), [numel(xVector) numel(yVector)])
    error("DensityLevelForCDF:SizeMismatch", ...
        "density must have size [numel(gridVectors{1}) numel(gridVectors{2})].");
end

nX = numel(xVector);
nY = numel(yVector);
if nX < 2 || nY < 2
    error("DensityLevelForCDF:InsufficientGrid", ...
        "gridVectors must each contain at least two points.");
end
if any(diff(xVector) <= 0) || any(diff(yVector) <= 0)
    error("DensityLevelForCDF:InvalidGridVectors", ...
        "gridVectors must be strictly increasing.");
end

xWeights = trapezoidalWeights(xVector);
yWeights = trapezoidalWeights(yVector);
pointMass = density .* (xWeights * yWeights.');
totalMass = sum(pointMass, "all");
if ~isfinite(totalMass) || totalMass <= 0
    dLevels = NaN(size(pctTarget));
    return
end

[sortedDensity, sortIndices] = sort(density(:), "descend");
sortedMass = pointMass(sortIndices);
cumulativeMass = cumsum(sortedMass) / totalMass;

[cumulativeMass, uniqueIndices] = unique(cumulativeMass, "stable");
sortedDensity = sortedDensity(uniqueIndices);

dLevels = NaN(size(pctTarget));
positiveTargets = pctTarget > 0;
if ~any(positiveTargets)
    return
end

positiveTargetValues = pctTarget(positiveTargets);
positiveTargetValues = min(positiveTargetValues, 1);
if numel(cumulativeMass) == 1
    dLevels(positiveTargets) = sortedDensity;
    return
end

dLevels(positiveTargets) = interp1(cumulativeMass, sortedDensity, positiveTargetValues, "linear", "extrap");
dLevels(pctTarget >= 1) = sortedDensity(end);
end

function weights = trapezoidalWeights(coordinates)
spacing = diff(coordinates);
weights = zeros(size(coordinates));
weights(1) = spacing(1)/2;
weights(end) = spacing(end)/2;
weights(2:end - 1) = (spacing(1:end - 1) + spacing(2:end))/2;
end
