function [V,r] = prepareGeneralizedEigenpairs(self,solver,A,B,V,r)
% Qualify raw pairs in solve units before taking their real parts.
if isa(solver,"IMSolverFiniteDifference")
    [V,r] = IMGeostrophicFiniteDifferenceTools.qualify(self,solver,A,B,V,r);
    return
end
finite = isfinite(real(r)) & isfinite(imag(r));
V = V(:,finite);
r = r(finite);
if isempty(r), return; end
[~,zeroPresent,IN] = self.endpointClassification(solver);
% The solve eigenvalue is r/rScale, and physical field units are balanced
% before measuring pencil-null directions. Weak N2 must not make F null.
realPair = abs(imag(r)) <= 1e-10*max(abs(real(r)),1);
for j = 1:numel(r)
    [~,i] = max(abs(V(:,j)));
    V(:,j) = V(:,j)*exp(-1i*angle(V(i,j)));
end
realPair = realPair & (vecnorm(imag(V),2,1).' <= 1e-8*vecnorm(real(V),2,1).');
V = real(V(:,realPair));
r = real(r(realPair));
if isempty(r), return; end
fieldScale = ones(size(V,1),1);
if size(V,1)>1, fieldScale(1:numel(solver.zNative)) = IN/self.g; end
balancedV = V./fieldScale;
balancedA = A.*fieldScale.';
balancedB = B.*fieldScale.';
metricAction = vecnorm(B*V,2,1)./max(norm(balancedB,1)*vecnorm(balancedV,2,1),realmin);
valid = metricAction>1e3*eps;
if size(V,1)>1
    n = numel(solver.zNative);
    F = solver.physicalDerivativeMatrix(0)*V(1:n,:);
    valid = valid & vecnorm(F,2,1)/(IN/self.g)>1e3*eps*vecnorm(balancedV,2,1);
end
V = V(:,valid);
r = r(valid);
if zeroPresent && ~isempty(r)
    [~,j] = min(abs(r));
    zeroResidual = norm(A*V(:,j))/(max(norm(balancedA,1)*norm(V(:,j)./fieldScale),realmin));
    if zeroResidual<1e-9 && abs(r(j))<1e-7
        r(j) = 0;
    end
end
[V,r] = self.refineSmallInverseDepthModes(solver,V,r,IN);
end
