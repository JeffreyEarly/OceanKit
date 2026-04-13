function [Aineq, bineq] = compileGlobalConstraints(globalConstraints, tKnot, K)
% Compile global constraints into coefficient inequalities.
%
% This helper turns global shape requests into coefficient-space
% inequalities of the form
%
% $$
% A_{\mathrm{ineq}}\xi \le 0.
% $$
%
% Positivity uses $$\xi_j \ge 0$$, which is sufficient because the terminated
% B-spline basis is nonnegative. Monotonicity uses first differences of
% adjacent coefficients along a selected tensor dimension, via
% [`monotonicDifferenceMatrix`](/spline-core/classes/constrainedspline/monotonicdifferencematrix.html).
%
% - Topic: Compile constraints
% - Developer: true
% - Declaration: [Aineq,bineq] = compileGlobalConstraints(globalConstraints,tKnot,K)
% - Parameter globalConstraints: GlobalConstraint array
% - Parameter tKnot: per-dimension knot vectors
% - Parameter K: spline order per dimension
% - Returns Aineq: inequality-constraint matrix
% - Returns bineq: inequality-constraint right-hand side
basisSize = reshape(cellfun(@numel, tKnot), 1, []) - reshape(K, 1, []);
numCoefficients = prod(basisSize);
Aineq = sparse([], [], [], 0, numCoefficients);
bineq = zeros(0,1);

for iConstraint = 1:numel(globalConstraints)
    constraint = globalConstraints(iConstraint);
    switch constraint.shape
        case GlobalConstraint.positiveShape
            constraintMatrix = -speye(numCoefficients);
        case GlobalConstraint.monotonicIncreasingShape
            constraintMatrix = ConstrainedSpline.monotonicDifferenceMatrix(basisSize, constraint.dimension, "increasing");
        case GlobalConstraint.monotonicDecreasingShape
            constraintMatrix = ConstrainedSpline.monotonicDifferenceMatrix(basisSize, constraint.dimension, "decreasing");
        otherwise
            error('ConstrainedSpline:UnsupportedGlobalConstraint',  'Unsupported global constraint shape.');
    end

    Aineq = [Aineq; constraintMatrix];
    bineq = [bineq; zeros(size(constraintMatrix,1),1)];
end
