function tKnotDidChange(self)
% Clear cached piecewise-polynomial data after knot updates.
%
% - Topic: Maintain cached state
% - Developer: true
% - Declaration: tKnotDidChange(self)
% - Parameter self: BSpline instance
self.Xtpp = [];
self.C = [];
self.t_pp = [];
self.xi_ = [];
