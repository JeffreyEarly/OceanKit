function [tEval, q, r] = centeredCoordinates(self, t, x, y)
% Convert fixed-frame coordinates to COM-frame coordinates.
%
% This method evaluates the fitted center-of-mass trajectory and returns
% the centered coordinates
% $$\tilde{x} = x - m_x(t)$$ and $$\tilde{y} = y - m_y(t)$$ with the same
% shape rules used by the mesoscale evaluators.
%
% `centeredCoordinates` is a derived coordinate transform built from the
% solved `centerOfMassTrajectory`; it does not refit the estimator or
% introduce additional state beyond the stored COM trajectory.
%
% ```matlab
% trajectory = fit.observedTrajectories(1);
% ti = trajectory.t;
% [tEval, q, r] = fit.centeredCoordinates(ti, trajectory.x(ti), trajectory.y(ti));
% plot(q, r)
% axis equal
% xlabel("\tilde{x} (m)")
% ylabel("\tilde{y} (m)")
% ```
%
% - Topic: Evaluate derived fields — Coordinate transform
% - nav_order: 1
% - Declaration: [tEval,q,r] = centeredCoordinates(self,t,x,y)
% - Parameter t: scalar time, vector matching the first dimension of `x` and `y`, or array matching `x` and `y`
% - Parameter x: x-coordinate array in meters
% - Parameter y: y-coordinate array in meters
% - Returns tEval: expanded evaluation-time array matching `q` and `r`
% - Returns q: centered x-coordinate array $$\tilde{x}$$
% - Returns r: centered y-coordinate array $$\tilde{y}$$
arguments
    self
    t {mustBeNumeric, mustBeReal}
    x {mustBeNumeric, mustBeReal}
    y {mustBeNumeric, mustBeReal}
end

if ~isequal(size(x), size(y))
    error("GriddedStreamfunction:EvaluationSizeMismatch", ...
        "x and y must have the same size.");
end

if isscalar(t)
    tEval = repmat(t, size(x));
elseif isvector(t) && ismatrix(x) && size(x, 1) == numel(t)
    tEval = repmat(reshape(t, [], 1), 1, size(x, 2));
else
    if ~isequal(size(t), size(x))
        error("GriddedStreamfunction:EvaluationTimeSizeMismatch", ...
            "t must be scalar or have the same size as x and y.");
    end
    tEval = t;
end

q = x - self.centerOfMassTrajectory.x(tEval);
r = y - self.centerOfMassTrajectory.y(tEval);
end
