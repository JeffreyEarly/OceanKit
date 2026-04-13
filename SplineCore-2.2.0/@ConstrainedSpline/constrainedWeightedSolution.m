function [xi, systemMatrix] = constrainedWeightedSolution(normalMatrix, rhs, Aeq, beq, Aineq, bineq)
% Solve weighted least squares with optional linear constraints.
%
% If only equality constraints are present, this forms and solves the KKT
% system
%
% $$
% \begin{bmatrix}
% N & A_{\mathrm{eq}}^T \\
% A_{\mathrm{eq}} & 0
% \end{bmatrix}
% \begin{bmatrix}
% \xi \\ \lambda
% \end{bmatrix}
% =
% \begin{bmatrix}
% r \\ b_{\mathrm{eq}}
% \end{bmatrix}.
% $$
%
% If inequality constraints are also present, the method switches to
% `quadprog` and solves the equivalent convex quadratic program with
% Hessian `N`.
%
% - Topic: Solve fit systems
% - Developer: true
% - Declaration: [xi,systemMatrix] = constrainedWeightedSolution(normalMatrix,rhs,Aeq,beq,Aineq,bineq)
% - Parameter normalMatrix: weighted normal-equation matrix
% - Parameter rhs: weighted right-hand side
% - Parameter Aeq: equality-constraint matrix
% - Parameter beq: equality-constraint right-hand side
% - Parameter Aineq: inequality-constraint matrix
% - Parameter bineq: inequality-constraint right-hand side
% - Returns xi: fitted coefficient vector
% - Returns systemMatrix: solved system matrix or quadratic Hessian proxy
numCoefficients = size(normalMatrix, 1);

if isempty(Aeq)
    Aeq = zeros(0, numCoefficients);
    beq = zeros(0,1);
end

if isempty(Aineq)
    Aineq = zeros(0, numCoefficients);
    bineq = zeros(0,1);
end

if isempty(Aineq)
    if isempty(Aeq)
        xi = ConstrainedSpline.leftSolve(normalMatrix, rhs);
        systemMatrix = normalMatrix;
    else
        systemMatrix = [normalMatrix, Aeq'; Aeq, zeros(size(Aeq,1))];
        solution = ConstrainedSpline.leftSolve(systemMatrix, [rhs; beq]);
        xi = solution(1:numCoefficients);
    end
    return;
end

H = (normalMatrix + normalMatrix')*0.5;
options = optimoptions('quadprog', 'Display', 'off', 'Algorithm', 'interior-point-convex');
[xi, ~, exitflag] = quadprog(2*H, -2*rhs, Aineq, bineq, Aeq, beq, [], [], [], options);
if exitflag <= 0
    error('ConstrainedSpline:OptimizationFailed',  'The constrained tensor-spline fit failed to converge.');
end

systemMatrix = normalMatrix;
