function S = smoothingMatrix(self)
% Return the smoothing matrix that maps observations to fitted values.
%
% Use this to inspect the linear action of the final weighted
% fit on the observed data.
%
% For an unconstrained fit with diagonal weights $$W$$, this matrix is the
% familiar linear smoother
%
% $$
% S = \mathbf{B}(\mathbf{B}^{T} W \mathbf{B})^{-1}\mathbf{B}^{T}W.
% $$
%
% When the fit uses a full observation covariance, the implementation uses
% the same linear map but applies the covariance through the stored solve
% object rather than explicitly forming an inverse.
%
% ```matlab
% S = spline.smoothingMatrix();
% valuesFit = S * spline.dataValues;
% ```
%
% - Topic: Analyze the fit
% - Declaration: S = smoothingMatrix(self)
% - Parameter self: ConstrainedSpline instance
% - Returns S: smoothing matrix
if ~isempty(self.Aeq) || ~isempty(self.Aineq)
    error('ConstrainedSpline:UnavailableSmoothingMatrix',  'smoothingMatrix is only available for unconstrained tensor fits.');
end

if size(self.W,1) == length(self.dataValues) && size(self.W,2) == 1
    S = (self.X*ConstrainedSpline.leftSolve(self.CmInv, self.X.')).*(self.W.');
elseif isa(self.W, 'decomposition')
    S = self.X*ConstrainedSpline.leftSolve(self.CmInv, self.X.');
    S = (self.W \ S.').';
else
    S = (self.X*ConstrainedSpline.leftSolve(self.CmInv, self.X.'))*self.W;
end
