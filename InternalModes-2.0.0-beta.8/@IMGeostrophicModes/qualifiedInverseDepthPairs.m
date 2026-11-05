function valid = qualifiedInverseDepthPairs(self,solver,V,r)
% Assess the regular pencil in balanced field and inverse-depth units.
[A,B] = self.inverseDepthPencil(solver);
[~,~,IN] = self.endpointClassification(solver);
n = numel(solver.zNative);
if IN>0
    fieldScale = [IN/self.g*ones(n,1);ones(n,1)];
    rScale = self.g/(IN*diff(self.zDomain));
else
    fieldScale = ones(2*n,1);
    rScale = 1/diff(self.zDomain);
end
B = rScale*B;
rowScale = max(max(abs(A.*fieldScale.'),[],2),max(abs(B.*fieldScale.'),[],2));
rowScale(rowScale==0) = 1;
A = A./rowScale;
B = B./rowScale;
balancedV = V./fieldScale;
lambda = r/rScale;
residual = vecnorm(A*V-(B*V).*lambda.',2,1) ...
    ./max((norm(A.*fieldScale.',1)+abs(lambda.')*norm(B.*fieldScale.',1)).*vecnorm(balancedV,2,1),realmin);
valid = residual.'<=1e-8;
end
