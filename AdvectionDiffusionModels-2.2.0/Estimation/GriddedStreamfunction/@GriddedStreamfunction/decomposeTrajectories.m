function decomposition = decomposeTrajectories(self, trajectories)
% Apply the fitted decomposition to a supplied drifter ensemble.
%
% Use this method when a fitted `GriddedStreamfunction` should be applied
% to a different collection of `TrajectorySpline` drifters without
% refitting the mesoscale streamfunction or the center-of-mass trajectory.
% The returned struct has the same field layout as `self.decomposition`,
% but all trajectory vectors are aligned with the supplied
% `trajectories` instead of `observedTrajectories`.
%
% This method applies the solved COM trajectory, mesoscale field, and
% common anchored background path to a new drifter ensemble; it does not
% change the fitted state stored on `self`.
%
% The returned decomposition keeps the fitted
% $$m_x(t), m_y(t), u^{\mathrm{meso}}, v^{\mathrm{meso}}$$ fixed and
% recomputes the common background and submesoscale residuals on the
% supplied trajectories. For the supplied drifters, the method enforces
% the same fixed-frame and centered-frame relations used by the fitted
% object,
%
% $$
% \dot{x} = u^{\mathrm{meso}} + u^{\mathrm{bg}} + u^{\mathrm{sm}},
% \qquad
% \dot{y} = v^{\mathrm{meso}} + v^{\mathrm{bg}} + v^{\mathrm{sm}},
% $$
%
% together with the centered-frame decomposition
%
% $$
% \dot{\tilde{x}} = u^{\mathrm{meso}}_{\mathrm{rel}} + u^{\mathrm{sm}}_{\mathrm{rel}},
% \qquad
% \dot{\tilde{y}} = v^{\mathrm{meso}}_{\mathrm{rel}} + v^{\mathrm{sm}}_{\mathrm{rel}}.
% $$
%
% Each `decomposition.fixedFrame.background(k)` is the common fitted
% background path re-anchored to zero at the first sample time of
% `trajectories(k)`. The fixed-frame mesoscale path carries the supplied
% drifter's observed initial position, while the fixed-frame and
% centered-frame submesoscale paths are zero-anchored.
%
% ```matlab
% decomposition = fit.decomposeTrajectories(otherTrajectories);
% iDrifter = 1;
% trajectory = otherTrajectories(iDrifter);
% ti = trajectory.t;
%
% background = decomposition.fixedFrame.background(iDrifter);
% mesoscale = decomposition.fixedFrame.mesoscale(iDrifter);
% submesoscale = decomposition.fixedFrame.submesoscale(iDrifter);
%
% plot(trajectory.x(ti), trajectory.y(ti), "k")
% hold on
% plot(background.x(ti), background.y(ti))
% plot(mesoscale.x(ti), mesoscale.y(ti))
% plot(submesoscale.x(ti), submesoscale.y(ti))
% axis equal
% ```
%
% - Topic: Apply fitted decomposition
% - nav_order: 1
% - Declaration: decomposition = decomposeTrajectories(self,trajectories)
% - Parameter trajectories: nonempty vector of `TrajectorySpline` drifters evaluated against the current fit
% - Returns decomposition: struct of fixed-frame and centered-frame `TrajectorySpline` vectors aligned with `trajectories`
arguments (Input)
    self (1,1) GriddedStreamfunction
    trajectories {mustBeA(trajectories, "TrajectorySpline"), mustBeVector, mustBeNonempty}
end
arguments (Output)
    decomposition struct
end

trajectories = reshape(trajectories, [], 1);
[~, decomposition] = decomposeTrajectorySet(self, trajectories);
end
