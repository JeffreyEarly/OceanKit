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
end
