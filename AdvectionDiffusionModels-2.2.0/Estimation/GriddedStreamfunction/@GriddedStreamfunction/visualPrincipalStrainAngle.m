function theta = visualPrincipalStrainAngle(sigma_n, sigma_s, options)
% Evaluate a jump-free principal strain angle for visualization.
%
% The principal strain angle is defined by
% $$\theta_p = \tfrac{1}{2}\operatorname{atan2}(\sigma_s,\sigma_n),$$
% with a principal representative in
% $$[-45^\circ,45^\circ)$$ or $$[-\pi/4,\pi/4).$$ For plotting, this
% helper chooses the nearest equivalent branch at each step by adding
% integer multiples of $$90^\circ$$ or $$\pi/2$$ along the first
% dimension, so the plotted angle stays continuous when it crosses the
% principal-range boundary.
%
% Vectors are treated as one time series and preserve their input
% orientation. Matrices and higher-dimensional arrays are processed
% independently down the first dimension.
%
% Use this helper after evaluating the derived strain diagnostics when you
% want a smoothly plotted angle rather than the wrapped principal value.
%
% ```matlab
% tFit = fit.fitSupportTimes;
% xCom = fit.centerOfMassTrajectory.x(tFit);
% yCom = fit.centerOfMassTrajectory.y(tFit);
% theta = GriddedStreamfunction.visualPrincipalStrainAngle( ...
%     fit.sigma_n(tFit, xCom, yCom), fit.sigma_s(tFit, xCom, yCom));
% plot(tFit, theta)
% xlabel("t (s)")
% ylabel("\theta_p (deg)")
% ```
%
% - Topic: Evaluate derived fields — Visualization helper
% - nav_order: 1
% - Declaration: theta = visualPrincipalStrainAngle(sigma_n,sigma_s,units=...)
% - Parameter sigma_n: normal strain array
% - Parameter sigma_s: shear strain array with the same size as `sigma_n`
% - Parameter units: optional output units `"degrees"` or `"radians"`, default `"degrees"`
% - Returns theta: continuity-preserving principal strain angle for visualization
arguments
    sigma_n {mustBeNumeric, mustBeReal}
    sigma_s {mustBeNumeric, mustBeReal}
    options.units {mustBeTextScalar, mustBeMember(options.units, ["degrees", "radians"])} = "degrees"
end

if ~isequal(size(sigma_n), size(sigma_s))
    error("GriddedStreamfunction:StrainAngleSizeMismatch", ...
        "sigma_n and sigma_s must have the same size.");
end

principalTheta = mod(atan2(sigma_s, sigma_n)/2 + pi/4, pi/2) - pi/4;
principalTheta(~(isfinite(sigma_n) & isfinite(sigma_s))) = NaN;

theta = principalTheta;
if isvector(principalTheta)
    theta = reshape(unwrapPrincipalBranch(principalTheta(:), pi/2), size(principalTheta));
else
    nTime = size(principalTheta, 1);
    theta = reshape(principalTheta, nTime, []);
    for iSeries = 1:size(theta, 2)
        theta(:, iSeries) = unwrapPrincipalBranch(theta(:, iSeries), pi/2);
    end
    theta = reshape(theta, size(principalTheta));
end

if string(options.units) == "degrees"
    theta = theta * 180/pi;
end
end

function theta = unwrapPrincipalBranch(thetaPrincipal, fullRange)
theta = thetaPrincipal;
previousTheta = NaN;

for iTime = 1:numel(thetaPrincipal)
    if isnan(thetaPrincipal(iTime))
        previousTheta = NaN;
        continue
    end

    if isnan(previousTheta)
        theta(iTime) = thetaPrincipal(iTime);
    else
        branchShift = round((previousTheta - thetaPrincipal(iTime))/fullRange);
        theta(iTime) = thetaPrincipal(iTime) + branchShift * fullRange;
    end

    previousTheta = theta(iTime);
end
end
