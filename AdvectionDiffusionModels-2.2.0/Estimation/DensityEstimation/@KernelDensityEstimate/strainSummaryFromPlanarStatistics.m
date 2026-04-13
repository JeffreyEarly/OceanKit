function strainSummary = strainSummaryFromPlanarStatistics(statistics, options)
% Reduce planar KDE statistics to strain magnitude-angle summaries.
%
% `strainSummaryFromPlanarStatistics(...)` converts the planar KDE mode and
% selected contour in the $$(\sigma_n,\sigma_s)$$ plane into the physical
% strain variables
%
% $$ \sigma = \sqrt{\sigma_n^2 + \sigma_s^2}, \qquad
% \theta = \tfrac{1}{2}\operatorname{atan2}(\sigma_s,\sigma_n). $$
%
% If `thetaReference` is supplied, the returned angular quantities are
% shifted by integer multiples of $$\pi$$ so the mode angle lies on the
% branch nearest that reference.
%
% If the selected contour encloses the origin, the returned magnitude
% summary remains mode-centered but is made zero-compatible by forcing the
% lower magnitude bound to zero. In that same case the angular bounds are
% reported as undefined through `thetaBounds = [NaN NaN]`.
%
% ```matlab
% statistics = KernelDensityEstimate.planarStatisticsFromData(data);
% strainSummary = KernelDensityEstimate.strainSummaryFromPlanarStatistics(statistics);
% ```
%
% - Topic: Evaluate the density estimate — Strain reduction
% - Declaration: strainSummary = strainSummaryFromPlanarStatistics(statistics,thetaReference=...)
% - Parameter statistics: planar KDE statistics returned by `planarStatisticsFromData(...)`
% - Parameter thetaReference: optional `pi`-periodic extensional-axis reference angle in radians used to choose the returned angle branch
% - Returns strainSummary: strain magnitude-angle summary with mode values, lower/upper uncertainty bounds, and a zero-compatibility flag
arguments
    statistics (1,1) struct
    options.thetaReference (1,1) double {mustBeReal} = NaN
end

if isfinite(options.thetaReference)
    polarSummary = KernelDensityEstimate.polarSummaryFromPlanarStatistics( ...
        statistics, ...
        angleReference=2 * options.thetaReference);
else
    polarSummary = KernelDensityEstimate.polarSummaryFromPlanarStatistics(statistics);
end

modeSigma = polarSummary.modeRadius;
sigmaBounds = polarSummary.radiusBounds;
sigmaErrors = polarSummary.radiusErrors;
modeTheta = polarSummary.modeAngle / 2;
thetaBounds = polarSummary.angleBounds / 2;
thetaErrors = polarSummary.angleErrors / 2;
referenceSigma = polarSummary.referenceRadius;
referenceTheta = polarSummary.referenceAngle / 2;
containsZero = polarSummary.containsOrigin;

if containsZero
    sigmaBounds(1) = 0;
    sigmaErrors = [modeSigma - sigmaBounds(1), sigmaBounds(2) - modeSigma];
    thetaBounds = [NaN NaN];
    thetaErrors = [NaN NaN];
end

strainSummary = struct( ...
    "modeSigma", modeSigma, ...
    "modeTheta", modeTheta, ...
    "sigmaBounds", sigmaBounds, ...
    "sigmaErrors", sigmaErrors, ...
    "thetaBounds", thetaBounds, ...
    "thetaErrors", thetaErrors, ...
    "referenceSigma", referenceSigma, ...
    "referenceTheta", referenceTheta, ...
    "containsZero", containsZero);
end
