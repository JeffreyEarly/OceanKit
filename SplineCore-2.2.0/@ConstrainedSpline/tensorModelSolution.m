function [xi,CmInv,W] = tensorModelSolution(values, designMatrix, distribution, rho_X, Aeq, beq, Aineq, bineq)
% Solve the tensor noisy-data model with iteratively reweighted least squares.
%
% This is the core fitting loop behind
% [`ConstrainedSpline`](/spline-core/classes/constrainedspline/constrainedspline.html).
% At each iteration it solves
%
% $$
% \min_{\xi}\ (y - \mathbf{B}\xi)^{T} W^{(n)} (y - \mathbf{B}\xi)
% $$
%
% subject to the supplied equality and inequality constraints, then updates
% the per-observation variances from the current residuals
% $$r^{(n)} = y - \mathbf{B}\xi^{(n)}$$ through the distribution model.
%
% When `rho_X` is supplied, the observation covariance is modeled as
%
% $$
% \Sigma^{(n)}_{ij} = \sigma_i^{(n)} \rho_{ij} \sigma_j^{(n)},
% $$
%
% and `W` represents the corresponding weighted solve rather than an
% explicitly formed inverse matrix.
%
% - Topic: Solve fit systems
% - Developer: true
% - Declaration: [xi,CmInv,W] = tensorModelSolution(values,designMatrix,distribution,rho_X,Aeq,beq,Aineq,bineq)
% - Parameter values: observation values as an N-by-1 vector
% - Parameter designMatrix: splines on the observation grid, N-by-M
% - Parameter distribution: distribution describing the errors
% - Parameter rho_X: optional observation correlation matrix
% - Parameter Aeq: optional equality-constraint matrix
% - Parameter beq: optional equality-constraint values
% - Parameter Aineq: optional inequality-constraint matrix
% - Parameter bineq: optional inequality-constraint values
% - Returns xi: fitted tensor spline coefficients
% - Returns CmInv: inverse coefficient covariance or system matrix
% - Returns W: final weight matrix or weights
arguments
    values (:,1) double
    designMatrix (:,:) double
    distribution
    rho_X = []
    Aeq = []
    beq = []
    Aineq = []
    bineq = []
end

sigma2_previous = (distribution.sigma0)^2 * ones(size(values));
W = ConstrainedSpline.weightMatrixFromSigma2(sigma2_previous, rho_X);

rel_error = 1.0;
repeats = 1;
while rel_error > 0.01 && repeats < 250
    [normalMatrix, rhs] = ConstrainedSpline.weightedNormalEquations(designMatrix, values, W);
    [xi, CmInv] = ConstrainedSpline.constrainedWeightedSolution(normalMatrix, rhs, Aeq, beq, Aineq, bineq);

    sigma2 = distribution.w(values - designMatrix*xi);
    rel_error = max(abs((sigma2-sigma2_previous)./sigma2), [], 'all');
    sigma2_previous = sigma2;
    W = ConstrainedSpline.weightMatrixFromSigma2(sigma2, rho_X);
    repeats = repeats + 1;
end
