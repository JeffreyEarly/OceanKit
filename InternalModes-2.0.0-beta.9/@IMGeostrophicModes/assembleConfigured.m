function [A,B,samples] = assembleConfigured(self,solver)
% Keep public eigenvalue units while assembling the regular coupled fields.
if isa(solver,"IMSolverFiniteDifference")
    [A,B,samples] = IMGeostrophicFiniteDifferenceTools.assemble(self,solver);
    return
end
[A,B] = self.inverseDepthPencil(solver);
n = numel(solver.zNative);
samples = struct("N2",self.stratificationValues(solver.zNative));
if self.name == "geostrophicGeneralizedPotentialEnstrophyModes"
    A(n+1:end,:) = (self.f0^2/self.g)*A(n+1:end,:)+self.parameters.k^2*B(n+1:end,:);
    V = solver.physicalDerivativeMatrix(0);
    s = solver.boundaryIndex("surface");
    b = solver.boundaryIndex("bottom");
    % Write physical endpoint rows directly, avoiding cancellation in A.
    A(s,:) = [-V(s,:) V(s,:)];
    B(s,:) = [-self.g/self.parameters.alpha0*V(s,:) zeros(1,n)];
    A(n+b,:) = [zeros(1,n) V(b,:)];
    B(n+b,:) = [self.g/self.parameters.alphaD*V(b,:) zeros(1,n)];
end
end
