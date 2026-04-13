function statistics = planarStatisticsFromData(data, options)
% Summarize a planar kernel density estimate from sampled data.
%
% `planarStatisticsFromData(...)` fits a two-dimensional Gaussian KDE,
% evaluates it on a regular grid, finds the KDE mode, and extracts one
% enclosed-mass contour around that mode. The returned summary stores only
% data-derived quantities, not rendering geometry.
%
% Relative to `originPoint = [x_0 y_0]`, the returned polar diagnostics use
%
% $$ r = \sqrt{(x-x_0)^2 + (y-y_0)^2}, $$
%
% and
%
% $$ \theta = \operatorname{atan2}(y-y_0, x-x_0). $$
%
% When explicit bounds are omitted, the default query box is the symmetric
% box centered on `originPoint` whose half-width is
% `boundPaddingFactor * max(abs(coordinates - originPoint), [], "all")`,
% clipped below by `minimumHalfWidth`.
%
% ```matlab
% data = [sigma_n(:), sigma_s(:)];
% statistics = KernelDensityEstimate.planarStatisticsFromData(data, referencePoint=[0 0]);
% ```
%
% - Topic: Evaluate the density estimate — Planar summary
% - Declaration: statistics = planarStatisticsFromData(data,referencePoint=...,originPoint=...,minimum=...,maximum=...,bandwidthGridSize=...,gridSize=...,summaryMass=...,boundPaddingFactor=...,minimumHalfWidth=...)
% - Parameter data: sampled planar coordinates as an `N-by-2` array
% - Parameter referencePoint: optional planar reference point stored in the returned summary
% - Parameter originPoint: planar origin used for radius and angle diagnostics
% - Parameter minimum: optional lower evaluation bound `[xmin ymin]`
% - Parameter maximum: optional upper evaluation bound `[xmax ymax]`
% - Parameter bandwidthGridSize: optional Botev-solver DCT grid size passed to `fromData(...)`
% - Parameter gridSize: optional evaluation grid size passed to `densityOnGrid(...)`
% - Parameter summaryMass: enclosed-mass fraction used to select the summary contour
% - Parameter boundPaddingFactor: multiplicative padding applied to the symmetric default bounds
% - Parameter minimumHalfWidth: lower bound for the symmetric default half-width
% - Returns statistics: planar KDE summary with density grid, mode point, selected contour, and radial/angular bounds
arguments
    data (:,2) double {mustBeReal, mustBeFinite, mustBeNonempty}
    options.referencePoint double {mustBePlanarVectorOrEmpty(options.referencePoint)} = []
    options.originPoint (1,2) double {mustBeReal, mustBeFinite} = [0 0]
    options.minimum double {mustBeBoundVectorOrEmpty(options.minimum)} = []
    options.maximum double {mustBeBoundVectorOrEmpty(options.maximum)} = []
    options.bandwidthGridSize = []
    options.gridSize = []
    options.summaryMass (1,1) double {mustBeReal, mustBeFinite, mustBeGreaterThanOrEqual(options.summaryMass, 0), mustBeLessThanOrEqual(options.summaryMass, 1)} = 0.68
    options.boundPaddingFactor (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 1.1
    options.minimumHalfWidth (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 1e-6
end

referencePoint = resolvedOptionalPoint(options.referencePoint);
originPoint = reshape(options.originPoint, 1, 2);
[minimum, maximum, halfWidth, usedExplicitBounds] = resolvedPlanarBounds( ...
    data, referencePoint, originPoint, options.minimum, options.maximum, options.boundPaddingFactor, options.minimumHalfWidth);

statisticsCore = evaluatePlanarStatisticsOnBounds( ...
    data, originPoint, minimum, maximum, options.bandwidthGridSize, options.gridSize, options.summaryMass);

if ~usedExplicitBounds
    iRetry = 0;
    while summaryContourNeedsExpandedBounds(statisticsCore.summaryContourVertices, statisticsCore.modePoint, minimum, maximum) && iRetry < 3
        halfWidth = 2 * halfWidth;
        minimum = originPoint - halfWidth;
        maximum = originPoint + halfWidth;
        statisticsCore = evaluatePlanarStatisticsOnBounds( ...
            data, originPoint, minimum, maximum, options.bandwidthGridSize, options.gridSize, options.summaryMass);
        iRetry = iRetry + 1;
    end
end

statistics = struct( ...
    "referencePoint", referencePoint, ...
    "originPoint", originPoint, ...
    "bounds", struct("minimum", minimum, "maximum", maximum, "halfWidth", halfWidth), ...
    "gridVectors", {statisticsCore.gridVectors}, ...
    "density", statisticsCore.density, ...
    "summaryMass", options.summaryMass, ...
    "summaryLevel", statisticsCore.summaryLevel, ...
    "summaryContourVertices", statisticsCore.summaryContourVertices, ...
    "modePoint", statisticsCore.modePoint, ...
    "modeRadius", statisticsCore.modeRadius, ...
    "modeAngle", statisticsCore.modeAngle, ...
    "radiusBounds", statisticsCore.radiusBounds, ...
    "angleBounds", statisticsCore.angleBounds);
end

function point = resolvedOptionalPoint(point)
if isempty(point)
    point = [NaN NaN];
    return
end

point = reshape(double(point), 1, 2);
end

function mustBePlanarVectorOrEmpty(value)
if isempty(value)
    return
end

if ~isnumeric(value) || ~isreal(value) || numel(value) ~= 2 || any(~isfinite(value))
    error("KernelDensityEstimate:InvalidPoint", ...
        "Planar coordinates must be finite two-element numeric vectors.");
end
end

function mustBeBoundVectorOrEmpty(value)
if isempty(value)
    return
end

if ~isnumeric(value) || ~isreal(value) || numel(value) ~= 2 || any(~isfinite(value))
    error("KernelDensityEstimate:InvalidBounds", ...
        "Bounds must be finite two-element numeric vectors.");
end
end

function [minimum, maximum, halfWidth, usedExplicitBounds] = resolvedPlanarBounds(data, referencePoint, originPoint, minimumOption, maximumOption, boundPaddingFactor, minimumHalfWidth)
if xor(isempty(minimumOption), isempty(maximumOption))
    error("KernelDensityEstimate:InvalidBounds", ...
        "minimum and maximum must be provided together when explicit bounds are used.");
end

if isempty(minimumOption)
    usedExplicitBounds = false;
    coordinates = data;
    if all(isfinite(referencePoint))
        coordinates = [coordinates; referencePoint]; %#ok<AGROW>
    end

    halfWidth = boundPaddingFactor * max(abs(coordinates - originPoint), [], "all");
    if ~isfinite(halfWidth) || halfWidth <= 0
        halfWidth = minimumHalfWidth;
    end
    halfWidth = max(halfWidth, minimumHalfWidth);
    minimum = originPoint - halfWidth;
    maximum = originPoint + halfWidth;
    return
end

usedExplicitBounds = true;
minimum = reshape(double(minimumOption), 1, []);
maximum = reshape(double(maximumOption), 1, []);
if numel(minimum) ~= 2 || numel(maximum) ~= 2 || any(~isfinite(minimum)) || any(~isfinite(maximum))
    error("KernelDensityEstimate:InvalidBounds", ...
        "minimum and maximum must be finite two-element numeric vectors.");
end
if any(maximum <= minimum)
    error("KernelDensityEstimate:InvalidBounds", ...
        "maximum must be strictly greater than minimum in every dimension.");
end

halfWidth = max(abs([minimum; maximum] - originPoint), [], "all");
halfWidth = max(halfWidth, minimumHalfWidth);
end

function statisticsCore = evaluatePlanarStatisticsOnBounds(data, originPoint, minimum, maximum, bandwidthGridSize, gridSize, summaryMass)
densityModel = KernelDensityEstimate.fromData(data, minimum=minimum, maximum=maximum, bandwidthGridSize=bandwidthGridSize);
[density, gridVectors] = densityModel.densityOnGrid(gridSize=gridSize);

[~, linearIndex] = max(density, [], "all");
[iXMode, iYMode] = ind2sub(size(density), linearIndex);
modePoint = [gridVectors{1}(iXMode), gridVectors{2}(iYMode)];
modeOffset = modePoint - originPoint;
modeRadius = hypot(modeOffset(1), modeOffset(2));
modeAngle = atan2(modeOffset(2), modeOffset(1));

summaryLevel = reshape(DensityLevelForCDF(gridVectors, density, summaryMass), 1, 1);
summaryContourVertices = zeros(0, 2);
radiusBounds = [NaN NaN];
angleBounds = [NaN NaN];

if isfinite(summaryLevel)
    contourVertices = contourPolygonsAtLevel(gridVectors, density, summaryLevel);
    if ~isempty(contourVertices)
        summaryContourVertices = selectedContourVertices(contourVertices, modePoint);
        if ~isempty(summaryContourVertices)
            contourOffsets = summaryContourVertices - originPoint;
            contourRadii = hypot(contourOffsets(:, 1), contourOffsets(:, 2));
            if any(isfinite(contourRadii))
                radiusBounds = [min(contourRadii), max(contourRadii)];
            end

            contourAngles = atan2(contourOffsets(:, 2), contourOffsets(:, 1));
            contourAngles = modeAngle + mod(contourAngles - modeAngle + pi, 2*pi) - pi;
            thetaMinimum = min(contourAngles);
            thetaMaximum = max(contourAngles);
            if isfinite(thetaMinimum) && isfinite(thetaMaximum) && thetaMaximum > thetaMinimum
                angleBounds = [thetaMinimum thetaMaximum];
            end
        end
    end
end

statisticsCore = struct( ...
    "gridVectors", {gridVectors}, ...
    "density", density, ...
    "summaryLevel", summaryLevel, ...
    "summaryContourVertices", summaryContourVertices, ...
    "modePoint", modePoint, ...
    "modeRadius", modeRadius, ...
    "modeAngle", modeAngle, ...
    "radiusBounds", radiusBounds, ...
    "angleBounds", angleBounds);
end

function tf = summaryContourNeedsExpandedBounds(vertices, modePoint, minimum, maximum)
if isempty(vertices)
    tf = true;
    return
end

if inpolygon(modePoint(1), modePoint(2), vertices(:, 1), vertices(:, 2))
    boundaryTolerance = 100 * eps(max(1, max(abs([minimum, maximum]), [], "all")));
    tf = any(abs(vertices(:, 1) - minimum(1)) <= boundaryTolerance | ...
        abs(vertices(:, 1) - maximum(1)) <= boundaryTolerance | ...
        abs(vertices(:, 2) - minimum(2)) <= boundaryTolerance | ...
        abs(vertices(:, 2) - maximum(2)) <= boundaryTolerance);
else
    tf = true;
end
end

function polygons = contourPolygonsAtLevel(gridVectors, density, contourLevel)
contourMatrix = contourc(gridVectors{1}, gridVectors{2}, density.', [contourLevel contourLevel]);
polygons = cell(0, 1);
iContour = 1;

while iContour < size(contourMatrix, 2)
    nVertices = contourMatrix(2, iContour);
    if nVertices >= 3
        vertices = contourMatrix(:, (iContour + 1):(iContour + nVertices)).';
        if ~isequal(vertices(1, :), vertices(end, :))
            vertices(end + 1, :) = vertices(1, :); %#ok<SAGROW>
        end
        polygons{end + 1, 1} = vertices; %#ok<SAGROW>
    end
    iContour = iContour + nVertices + 1;
end
end

function vertices = selectedContourVertices(polygons, modePoint)
containsMode = false(numel(polygons), 1);
distanceToMode = inf(numel(polygons), 1);

for iPolygon = 1:numel(polygons)
    vertices = polygons{iPolygon};
    containsMode(iPolygon) = inpolygon(modePoint(1), modePoint(2), vertices(:, 1), vertices(:, 2));
    distanceToMode(iPolygon) = min(hypot(vertices(:, 1) - modePoint(1), vertices(:, 2) - modePoint(2)));
end

if any(containsMode)
    polygonIndices = find(containsMode);
    [~, iSelected] = min(distanceToMode(polygonIndices));
    vertices = polygons{polygonIndices(iSelected)};
else
    [~, iSelected] = min(distanceToMode);
    vertices = polygons{iSelected};
end
end
