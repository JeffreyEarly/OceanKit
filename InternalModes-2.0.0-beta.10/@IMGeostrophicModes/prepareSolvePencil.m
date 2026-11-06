function [A,B,map] = prepareSolvePencil(self,solver,~,~)
% Balance inverse-depth units and expand the all-neutral finite subspace exactly.
if isa(solver,"IMSolverFiniteDifference")
    [A,B,map] = IMGeostrophicFiniteDifferenceTools.preparePencil(self,solver);
    return
end
n = numel(solver.zNative);
[~,~,IN] = self.endpointClassification(solver);
if IN==0
    [r,V] = self.fullyNeutralMode(solver);
    rScale = 1/diff(self.zDomain);
    A = r/rScale;
    B = 1;
    fieldSizes = 1;
    columnScale = 1;
    expand = @(coefficients) V*coefficients;
else
    [A,B] = self.inverseDepthPencil(solver);
    rScale = self.g/(IN*diff(self.zDomain));
    B = rScale*B;
    fieldSizes = [n n];
    columnScale = [IN/self.g*ones(1,n) ones(1,n)];
    expand = @(coefficients) coefficients;
end
map = struct("nativeFieldSizes",fieldSizes,"hasInverseDepth",true, ...
    "nativeColumnScale",columnScale,"restoreEigenvalues",@(lambda) rScale*lambda, ...
    "expandEigenvectors",expand,"publicEigenvalues",@(V,r) self.publicEigenvalues(solver,V,r));
if IN>0 && isa(solver,"IMSolverSpectral") && self.name=="geostrophicAPVModes"
    staticRows = all(B==0,2);
    if self.parameters.g0==0
        % F vanishes at the surface, so the surface G derivative is
        % another eigenvalue-independent constraint.
        staticRows(n+solver.boundaryIndex("surface")) = true;
    end
    map.reduceEquilibratedPencil = @(Aeq,Beq) reduceAPVConstraints(Aeq,Beq,staticRows);
end
end

function [A,B,lift] = reduceAPVConstraints(A,B,staticRows)
% Eliminate the static equations before QZ so infinite roots cannot
% enter the finite signed-mode inventory through roundoff.
lift = [];
dynamicRows = ~staticRows;
candidateLift = null(A(staticRows,:));
if size(candidateLift,2) ~= nnz(dynamicRows)
    return
end
finiteB = B(dynamicRows,:)*candidateLift;
% Degenerate stratification or endpoint laws can leave further static
% constraints. Keep the original descriptor solve for those cases.
if rank(finiteB) < size(finiteB,1)
    return
end
A = A(dynamicRows,:)*candidateLift;
B = finiteB;
lift = candidateLift;
end
