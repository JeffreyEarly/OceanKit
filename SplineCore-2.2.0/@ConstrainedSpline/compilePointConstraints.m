function [Aeq, beq, Aineq, bineq] = compilePointConstraints(pointConstraints, tKnot, K)
% Compile point constraints into equality and inequality systems.
%
% For each point constraint, this function evaluates the tensor-product
% basis or derivative row at the constrained points and turns the requested
% relation into linear constraints on the coefficient vector `xi`.
%
% If a row of the evaluated basis is denoted by $$b_m^T$$, then the three
% supported relations are compiled as
%
% $$
% b_m^{T}\xi = v_m, \qquad
% b_m^{T}\xi \le v_m, \qquad
% b_m^{T}\xi \ge v_m,
% $$
%
% with lower bounds rewritten as $$-b_m^T \xi \le -v_m$$.
%
% Rows sharing the same derivative order are batched together so the basis
% matrix only has to be evaluated once per derivative multi-index.
%
% - Topic: Compile constraints
% - Developer: true
% - Declaration: [Aeq,beq,Aineq,bineq] = compilePointConstraints(pointConstraints,tKnot,K)
% - Parameter pointConstraints: PointConstraint array
% - Parameter tKnot: per-dimension knot vectors
% - Parameter K: spline order per dimension
% - Returns Aeq: equality-constraint matrix
% - Returns beq: equality-constraint right-hand side
% - Returns Aineq: inequality-constraint matrix
% - Returns bineq: inequality-constraint right-hand side
basisSize = reshape(cellfun(@numel, tKnot), 1, []) - reshape(K, 1, []);
numCoefficients = prod(basisSize);
Aeq = sparse([], [], [], 0, numCoefficients);
beq = zeros(0,1);
Aineq = sparse([], [], [], 0, numCoefficients);
bineq = zeros(0,1);

for iConstraint = 1:numel(pointConstraints)
    constraint = pointConstraints(iConstraint);
    [groupOrders, ~, groupIndex] = unique(constraint.D, 'rows', 'stable');
    for iGroup = 1:size(groupOrders, 1)
        isGroup = groupIndex == iGroup;
        B = sparse(TensorSpline.matrixForPointMatrix(constraint.points(isGroup,:), knotPoints=tKnot, S=K - 1, D=groupOrders(iGroup,:)));
        values = constraint.value(isGroup);

        switch constraint.relation
            case PointConstraint.equalRelation
                Aeq = [Aeq; B];
                beq = [beq; values];
            case PointConstraint.lowerBoundRelation
                Aineq = [Aineq; -B];
                bineq = [bineq; -values];
            case PointConstraint.upperBoundRelation
                Aineq = [Aineq; B];
                bineq = [bineq; values];
            otherwise
                error('ConstrainedSpline:InvalidConstraintRelation',  'Unsupported point-constraint relation.');
        end
    end
end
