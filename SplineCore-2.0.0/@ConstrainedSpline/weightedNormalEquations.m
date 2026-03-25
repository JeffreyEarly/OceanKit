function [normalMatrix, rhs] = weightedNormalEquations(X, x, W)
% Assemble weighted normal equations.
%
% For the weighted least-squares problem
%
% $$
% \min_{\xi}\ (x - X\xi)^T W (x - X\xi),
% $$
%
% this helper returns the normal matrix and right-hand side
%
% $$
% N = X^T W X, \qquad r = X^T W x.
% $$
%
% `W` may be a dense matrix, a diagonal weight vector, a scalar, or a
% MATLAB `decomposition` object representing a correlated observation
% covariance solve.
%
% - Topic: Solve fit systems
% - Developer: true
% - Declaration: [normalMatrix,rhs] = weightedNormalEquations(X,x,W)
% - Parameter X: design matrix
% - Parameter x: observed values
% - Parameter W: weight matrix, vector, scalar, or decomposition
% - Returns normalMatrix: weighted normal-equation matrix
% - Returns rhs: weighted normal-equation right-hand side
XT = X';
if isa(W, 'decomposition')
    weightedX = W \ X;
    normalMatrix = XT*weightedX;
    rhs = XT*(W \ x);
elseif size(W,1) == length(x) && size(W,2) == length(x)
    normalMatrix = XT*W*X;
    rhs = XT*W*x;
elseif isscalar(W)
    normalMatrix = (XT*X)*W;
    rhs = XT*W*x;
elseif isvector(W) && numel(W) == length(x)
    normalMatrix = XT*(W.*X);
    rhs = XT*(W.*x);
else
    error('W must have the same length as x and X.');
end
