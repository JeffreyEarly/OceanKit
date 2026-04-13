function polarSummary = polarSummaryFromPlanarStatistics(statistics, options)
% Reduce planar KDE statistics to radius-angle uncertainty summaries.
%
% `polarSummaryFromPlanarStatistics(...)` converts the planar mode and the
% selected enclosed-mass contour into one-dimensional radius and angle
% summaries relative to `statistics.originPoint`. The returned angular
% quantities are the raw plane angle
% $$\alpha = \operatorname{atan2}(y-y_0, x-x_0),$$
% not any problem-specific half-angle reduction.
%
% If `angleReference` is supplied, the returned angular quantities are
% shifted by integer multiples of $$2\pi$$ so the mode angle lies on the
% branch nearest that reference.
%
% If the selected contour encloses `statistics.originPoint`, the radius
% summary remains mode-centered but the lower radial bound is set to zero
% and the angular bounds are reported as undefined. This origin-inclusive
% condition is exposed through `containsOrigin`.
%
% ```matlab
% statistics = KernelDensityEstimate.planarStatisticsFromData(data);
% polarSummary = KernelDensityEstimate.polarSummaryFromPlanarStatistics(statistics);
% ```
%
% - Topic: Evaluate the density estimate — Polar reduction
% - Declaration: polarSummary = polarSummaryFromPlanarStatistics(statistics,angleReference=...)
% - Parameter statistics: planar KDE statistics returned by `planarStatisticsFromData(...)`
% - Parameter angleReference: optional `2*pi`-periodic plane-angle reference in radians used to choose the returned angle branch
% - Returns polarSummary: radius-plane-angle summary with mode values, lower/upper uncertainty bounds, and an origin-inclusion flag
arguments
    statistics (1,1) struct {mustBePolarReduciblePlanarStatistics(statistics)}
    options.angleReference (1,1) double {mustBeReal} = NaN
end

modeRadius = statistics.modeRadius;
modeAngle = statistics.modeAngle;
radiusBounds = reshape(statistics.radiusBounds, 1, 2);
angleBounds = reshape(statistics.angleBounds, 1, 2);
referencePoint = reshape(statistics.referencePoint, 1, 2);
originPoint = reshape(statistics.originPoint, 1, 2);
containsOrigin = contourContainsPoint(statistics.summaryContourVertices, originPoint);

if containsOrigin
    contourOffsets = statistics.summaryContourVertices - originPoint;
    contourRadii = hypot(contourOffsets(:, 1), contourOffsets(:, 2));
    contourRadii = contourRadii(isfinite(contourRadii));
    if isempty(contourRadii)
        radiusBounds = [0 modeRadius];
    else
        radiusBounds = [0 max([modeRadius; contourRadii])];
    end
    angleBounds = [NaN NaN];
end

angleReference = options.angleReference;
if isfinite(angleReference) && isfinite(modeAngle)
    wrappedModeAngle = wrapNearReference(modeAngle, angleReference);
    angleShift = wrappedModeAngle - modeAngle;
    modeAngle = wrappedModeAngle;
    if all(isfinite(angleBounds))
        angleBounds = angleBounds + angleShift;
    end
else
    angleShift = 0;
end

if all(isfinite(referencePoint))
    referenceOffset = referencePoint - originPoint;
    referenceRadius = hypot(referenceOffset(1), referenceOffset(2));
    referenceAngle = atan2(referenceOffset(2), referenceOffset(1));
    if isfinite(angleReference)
        referenceAngle = wrapNearReference(referenceAngle, angleReference);
    end
else
    referenceRadius = NaN;
    referenceAngle = NaN;
end

if all(isfinite(radiusBounds))
    radiusErrors = [modeRadius - radiusBounds(1), radiusBounds(2) - modeRadius];
else
    radiusErrors = [NaN NaN];
end

if all(isfinite(angleBounds))
    angleErrors = [modeAngle - angleBounds(1), angleBounds(2) - modeAngle];
else
    angleErrors = [NaN NaN];
end

polarSummary = struct( ...
    "modeRadius", modeRadius, ...
    "modeAngle", modeAngle, ...
    "radiusBounds", radiusBounds, ...
    "radiusErrors", radiusErrors, ...
    "angleBounds", angleBounds, ...
    "angleErrors", angleErrors, ...
    "referenceRadius", referenceRadius, ...
    "referenceAngle", referenceAngle, ...
    "containsOrigin", containsOrigin);
end

function angle = wrapNearReference(angle, angleReference)
angle = angleReference + mod(angle - angleReference + pi, 2*pi) - pi;
end

function tf = contourContainsPoint(vertices, point)
tf = false;
if isempty(vertices)
    return
end

[isInterior, isOnBoundary] = inpolygon(point(1), point(2), vertices(:, 1), vertices(:, 2));
tf = isInterior | isOnBoundary;
end

function mustBePolarReduciblePlanarStatistics(statistics)
requiredFields = ["modeRadius", "modeAngle", "radiusBounds", "angleBounds", "referencePoint", "originPoint", "summaryContourVertices"];
missingFields = requiredFields(~isfield(statistics, requiredFields));
if ~isempty(missingFields)
    error("KernelDensityEstimate:MissingPlanarStatisticsField", ...
        "statistics is missing required field(s): %s.", strjoin(cellstr(missingFields), ", "));
end

mustBeRealScalar(statistics.modeRadius, "statistics.modeRadius");
mustBeRealScalar(statistics.modeAngle, "statistics.modeAngle");
mustBeBoundPair(statistics.radiusBounds, "statistics.radiusBounds");
mustBeBoundPair(statistics.angleBounds, "statistics.angleBounds");
mustBePlanarPoint(statistics.referencePoint, "statistics.referencePoint", true);
mustBePlanarPoint(statistics.originPoint, "statistics.originPoint", false);
mustBeContourVertices(statistics.summaryContourVertices);
end

function mustBeRealScalar(value, fieldName)
if ~isscalar(value) || ~isnumeric(value) || ~isreal(value)
    error("KernelDensityEstimate:InvalidPlanarStatistics", ...
        "%s must be a real scalar.", fieldName);
end
end

function mustBeBoundPair(value, fieldName)
if ~isnumeric(value) || ~isreal(value) || numel(value) ~= 2
    error("KernelDensityEstimate:InvalidPlanarStatistics", ...
        "%s must be a two-element real vector.", fieldName);
end

value = reshape(value, 1, 2);
if all(isnan(value))
    return
end
if all(isfinite(value)) && value(2) >= value(1)
    return
end

error("KernelDensityEstimate:InvalidPlanarStatistics", ...
    "%s must contain either two finite ordered values or two NaNs.", fieldName);
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

function mustBeContourVertices(value)
if ~(isnumeric(value) && isreal(value) && size(value, 2) == 2 && all(isfinite(value), "all"))
    error("KernelDensityEstimate:InvalidPlanarStatistics", ...
        "statistics.summaryContourVertices must be a finite N-by-2 real array.");
end
end
