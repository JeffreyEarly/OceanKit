function plotPlanarStatistics(ax, statistics, options)
% Plot a planar KDE summary from precomputed statistics.
%
% `plotPlanarStatistics(...)` draws contours, a wedge summary, a scatter
% cloud, a reference point, and optional radial guides from the data
% stored in `statistics`. Rendering geometry such as contour levels and the
% wedge polygon is built inside this helper rather than stored in the
% statistical summary.
%
% ```matlab
% statistics = KernelDensityEstimate.planarStatisticsFromData(data, referencePoint=[0 0]);
% KernelDensityEstimate.plotPlanarStatistics(gca, statistics, samples=data, showWedge=true);
% ```
%
% - Topic: Evaluate the density estimate — Planar rendering
% - Declaration: plotPlanarStatistics(ax,statistics,samples=...,contourMasses=...,showContours=...,showScatter=...,showReferencePoint=...,showWedge=...,scatterColor=...,scatterAlpha=...,scatterSize=...,referenceFaceColor=...,referenceEdgeColor=...,referenceSize=...,wedgeColor=...,wedgeAlpha=...,wedgeResolution=...,ringRadii=...,ringColor=...,ringLineStyle=...,ringLineWidth=...,applyEqualAxes=...,applyStoredBounds=...)
% - Parameter ax: target axes
% - Parameter statistics: planar KDE statistics returned by `planarStatisticsFromData(...)`
% - Parameter samples: optional `N-by-2` sample cloud to draw as a scatter layer
% - Parameter contourMasses: enclosed-mass targets used to draw filled contours
% - Parameter showContours: toggle filled-density contours
% - Parameter showScatter: toggle the sample cloud
% - Parameter showReferencePoint: toggle the stored reference point marker
% - Parameter showWedge: toggle the wedge defined by the stored mode radius and angle bounds
% - Parameter scatterColor: scatter-marker face color
% - Parameter scatterAlpha: scatter-marker face alpha
% - Parameter scatterSize: scatter-marker size
% - Parameter referenceFaceColor: reference-point face color
% - Parameter referenceEdgeColor: reference-point edge color
% - Parameter referenceSize: reference-point marker size
% - Parameter wedgeColor: wedge fill color
% - Parameter wedgeAlpha: wedge fill alpha
% - Parameter wedgeResolution: number of arc samples used to draw the wedge
% - Parameter ringRadii: optional radii for dotted reference rings centered on `originPoint`
% - Parameter ringColor: reference-ring color
% - Parameter ringLineStyle: reference-ring line style
% - Parameter ringLineWidth: reference-ring line width
% - Parameter applyEqualAxes: toggle `axis equal`
% - Parameter applyStoredBounds: toggle use of `statistics.bounds`
arguments
    ax (1,1) matlab.graphics.axis.Axes
    statistics (1,1) struct {mustBeRenderablePlanarStatistics(statistics)}
    options.samples (:,2) double {mustBeReal, mustBeFinite} = zeros(0, 2)
    options.contourMasses (:,1) double {mustBeReal, mustBeFinite, mustBeGreaterThanOrEqual(options.contourMasses, 0), mustBeLessThanOrEqual(options.contourMasses, 1)} = (0.8:-0.1:0.1).'
    options.showContours (1,1) logical = true
    options.showScatter (1,1) logical = true
    options.showReferencePoint (1,1) logical = true
    options.showWedge (1,1) logical = false
    options.scatterColor (1,3) double {mustBeReal, mustBeFinite} = [0 0.4470 0.7410]
    options.scatterAlpha (1,1) double {mustBeReal, mustBeFinite, mustBeGreaterThanOrEqual(options.scatterAlpha, 0), mustBeLessThanOrEqual(options.scatterAlpha, 1)} = 0.18
    options.scatterSize (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 12
    options.referenceFaceColor (1,3) double {mustBeReal, mustBeFinite} = [1 1 1]
    options.referenceEdgeColor (1,3) double {mustBeReal, mustBeFinite} = [0 0 0]
    options.referenceSize (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 42
    options.wedgeColor (1,3) double {mustBeReal, mustBeFinite} = 0.7 * [1 1 1]
    options.wedgeAlpha (1,1) double {mustBeReal, mustBeFinite, mustBeGreaterThanOrEqual(options.wedgeAlpha, 0), mustBeLessThanOrEqual(options.wedgeAlpha, 1)} = 0.35
    options.wedgeResolution (1,1) double {mustBeReal, mustBeFinite, mustBeInteger, mustBePositive} = 100
    options.ringRadii (:,1) double {mustBeReal, mustBeFinite} = zeros(0, 1)
    options.ringColor (1,3) double {mustBeReal, mustBeFinite} = 0.75 * [1 1 1]
    options.ringLineStyle {mustBeTextScalar} = ":"
    options.ringLineWidth (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 0.8
    options.applyEqualAxes (1,1) logical = true
    options.applyStoredBounds (1,1) logical = true
end

samples = options.samples;
gridVectors = statistics.gridVectors;
density = statistics.density;
referencePoint = reshape(statistics.referencePoint, 1, 2);
originPoint = reshape(statistics.originPoint, 1, 2);
bounds = statistics.bounds;

hold(ax, "on");

ringRadii = reshape(options.ringRadii, [], 1);
if any(ringRadii <= 0)
    error("KernelDensityEstimate:InvalidRingRadii", ...
        "ringRadii must contain only positive values.");
end
for iRadius = 1:numel(ringRadii)
    rectangle(ax, Position=[originPoint(1) - ringRadii(iRadius), originPoint(2) - ringRadii(iRadius), 2 * ringRadii(iRadius), 2 * ringRadii(iRadius)], ...
        Curvature=[1 1], EdgeColor=options.ringColor, LineStyle=options.ringLineStyle, LineWidth=options.ringLineWidth);
end

if options.showScatter && ~isempty(samples)
    scatter(ax, samples(:, 1), samples(:, 2), options.scatterSize, options.scatterColor, "filled", ...
        MarkerEdgeColor="none", MarkerFaceAlpha=options.scatterAlpha);
end

if options.showContours
    contourLevels = finiteContourLevels(gridVectors, density, options.contourMasses);
    contourf(ax, gridVectors{1}, gridVectors{2}, density.', contourLevels, LineStyle="none");
end

if options.showWedge
    wedgeVertices = wedgeVerticesFromStatistics(statistics, options.wedgeResolution);
    if ~isempty(wedgeVertices)
        patch(ax, wedgeVertices(:, 1), wedgeVertices(:, 2), options.wedgeColor, EdgeColor="none", FaceAlpha=options.wedgeAlpha);
    end
end

if options.showReferencePoint && all(isfinite(referencePoint))
    scatter(ax, referencePoint(1), referencePoint(2), options.referenceSize, options.referenceFaceColor, "filled", ...
        MarkerEdgeColor=options.referenceEdgeColor, LineWidth=1.1);
end

if options.applyEqualAxes
    axis(ax, "equal");
end
if options.applyStoredBounds
    xlim(ax, [bounds.minimum(1), bounds.maximum(1)]);
    ylim(ax, [bounds.minimum(2), bounds.maximum(2)]);
end
box(ax, "on");
end

function contourLevels = finiteContourLevels(gridVectors, density, contourMasses)
contourLevels = DensityLevelForCDF(gridVectors, density, contourMasses);
contourLevels = reshape(contourLevels, [], 1);
contourLevels = contourLevels(isfinite(contourLevels));
contourLevels = unique(contourLevels, "sorted");
if numel(contourLevels) >= 2
    return
end

densityMinimum = min(density, [], "all");
densityMaximum = max(density, [], "all");
if densityMaximum <= densityMinimum
    contourLevels = [densityMinimum; densityMinimum + eps(densityMinimum + 1)];
else
    contourLevels = linspace(densityMinimum, densityMaximum, max(2, numel(contourMasses))).';
end
end

function wedgeVertices = wedgeVerticesFromStatistics(statistics, wedgeResolution)
originPoint = reshape(statistics.originPoint, 1, 2);
modeRadius = statistics.modeRadius;
angleBounds = reshape(statistics.angleBounds, 1, 2);

wedgeVertices = zeros(0, 2);
if ~isfinite(modeRadius) || modeRadius <= 0 || any(~isfinite(angleBounds)) || angleBounds(2) <= angleBounds(1)
    return
end

wedgeAngles = linspace(angleBounds(1), angleBounds(2), wedgeResolution);
wedgeVertices = [originPoint; originPoint + modeRadius * [cos(wedgeAngles(:)), sin(wedgeAngles(:))]; originPoint];
end

function mustBeRenderablePlanarStatistics(statistics)
requiredFields = ["gridVectors", "density", "referencePoint", "originPoint", "bounds", "modeRadius", "angleBounds"];
missingFields = requiredFields(~isfield(statistics, requiredFields));
if ~isempty(missingFields)
    error("KernelDensityEstimate:MissingPlanarStatisticsField", ...
        "statistics is missing required field(s): %s.", strjoin(cellstr(missingFields), ", "));
end

gridVectors = statistics.gridVectors;
if ~iscell(gridVectors) || numel(gridVectors) ~= 2
    error("KernelDensityEstimate:InvalidPlanarStatistics", ...
        "statistics.gridVectors must be a two-element cell array.");
end
if ~isnumeric(statistics.density) || ~isreal(statistics.density) || ndims(statistics.density) ~= 2 || ...
        ~isequal(size(statistics.density), [numel(gridVectors{1}), numel(gridVectors{2})]) || any(~isfinite(statistics.density), "all")
    error("KernelDensityEstimate:InvalidPlanarStatistics", ...
        "statistics.density must be a finite matrix with size [numel(gridVectors{1}) numel(gridVectors{2})].");
end

mustBePlanarPoint(statistics.referencePoint, "statistics.referencePoint", true);
mustBePlanarPoint(statistics.originPoint, "statistics.originPoint", false);

bounds = statistics.bounds;
if ~isstruct(bounds) || ~isfield(bounds, "minimum") || ~isfield(bounds, "maximum")
    error("KernelDensityEstimate:InvalidPlanarStatistics", ...
        "statistics.bounds must contain minimum and maximum fields.");
end
mustBeFinitePlanarVector(bounds.minimum, "statistics.bounds.minimum");
mustBeFinitePlanarVector(bounds.maximum, "statistics.bounds.maximum");
if any(reshape(bounds.maximum, 1, 2) <= reshape(bounds.minimum, 1, 2))
    error("KernelDensityEstimate:InvalidPlanarStatistics", ...
        "statistics.bounds.maximum must be strictly greater than statistics.bounds.minimum.");
end

if ~isscalar(statistics.modeRadius) || ~isnumeric(statistics.modeRadius) || ~isreal(statistics.modeRadius)
    error("KernelDensityEstimate:InvalidPlanarStatistics", ...
        "statistics.modeRadius must be a real scalar.");
end

angleBounds = statistics.angleBounds;
if ~isnumeric(angleBounds) || ~isreal(angleBounds) || numel(angleBounds) ~= 2
    error("KernelDensityEstimate:InvalidPlanarStatistics", ...
        "statistics.angleBounds must be a two-element real vector.");
end
angleBounds = reshape(angleBounds, 1, 2);
if any(isfinite(angleBounds)) && (~all(isfinite(angleBounds)) || angleBounds(2) <= angleBounds(1))
    error("KernelDensityEstimate:InvalidPlanarStatistics", ...
        "statistics.angleBounds must contain either two finite increasing values or two NaNs.");
end
end

function mustBePlanarPoint(value, fieldName, allowNaNPair)
if ~isnumeric(value) || ~isreal(value) || numel(value) ~= 2
    error("KernelDensityEstimate:InvalidPlanarStatistics", ...
        "%s must be a two-element real vector.", fieldName);
end

value = reshape(value, 1, 2);
if all(isfinite(value))
    return
end
if allowNaNPair && all(isnan(value))
    return
end

error("KernelDensityEstimate:InvalidPlanarStatistics", ...
    "%s must contain either two finite values or two NaNs.", fieldName);
end

function mustBeFinitePlanarVector(value, fieldName)
if ~isnumeric(value) || ~isreal(value) || numel(value) ~= 2 || any(~isfinite(value))
    error("KernelDensityEstimate:InvalidPlanarStatistics", ...
        "%s must be a finite two-element real vector.", fieldName);
end
end
