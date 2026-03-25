function W = weightMatrixFromSigma2(sigma2, rho_X)
% Build the observation-weight matrix from per-observation variances.
%
% For independent errors this returns the diagonal weights
%
% $$
% W = \mathrm{diag}(\sigma_1^{-2},\ldots,\sigma_N^{-2}),
% $$
%
% represented here by the vector `1./sigma2`.
%
% With correlated errors, it builds the observation covariance
%
% $$
% \Sigma_{ij} = \sigma_i \rho_{ij} \sigma_j,
% $$
%
% and returns a MATLAB `decomposition` object for $$\Sigma$$, which is then
% used in left-division solves by
% [`weightedNormalEquations`](/spline-core/classes/constrainedspline/weightednormalequations.html).
%
% - Topic: Solve fit systems
% - Developer: true
% - Declaration: W = weightMatrixFromSigma2(sigma2,rho_X)
% - Parameter sigma2: per-observation variances
% - Parameter rho_X: optional observation correlation matrix
% - Returns W: weight representation for the fit
if ~isempty(rho_X)
    Sigma2 = (sqrt(sigma2) * sqrt(sigma2).') .* rho_X;
    W = decomposition(Sigma2);
else
    W = 1./sigma2;
end
