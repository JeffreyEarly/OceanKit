function mapped = publicEigenvalues(self,solver,V,r)
% Compute positive physical eigenvalues without cancelling k2 against r.
if self.name=="geostrophicAPVModes"
    mapped = struct("values",r,"qualified",self.qualifiedInverseDepthPairs(solver,V,r));
    return
end
n = numel(solver.zNative);
T = solver.physicalDerivativeMatrix(0);
F = T*V(1:n,:);
G = T*V(n+1:end,:);
w = solver.innerProductWeights(solver.zNative,self.zDomain);
N2 = self.stratificationValues(solver.zNative);
Fs = F(solver.boundaryIndex("surface"),:);
Fb = F(solver.boundaryIndex("bottom"),:);
volumeF = sum(w.*F.^2,1);
numerator = self.parameters.k^2*volumeF+self.f0^2/self.g^2*sum(w.*N2.*G.^2,1)+self.f0^2/self.g*Fs.^2;
denominator = volumeF+self.f0^2/self.parameters.alpha0*Fs.^2+self.f0^2/self.parameters.alphaD*Fb.^2;
lambda = (numerator./denominator).';
lambda(r==0) = self.parameters.k^2;
[A,B] = self.assembleConfigured(solver);
residual = vecnorm(A*V-(B*V).*lambda.',2,1) ...
    ./max((norm(A,1)+abs(lambda.')*norm(B,1)).*vecnorm(V,2,1),realmin);
affine = self.parameters.k^2+self.f0^2/self.g*r;
scale = self.parameters.k^2+self.f0^2/self.g*abs(r);
% Preserve failed candidates in their physical ordering. Selection must
% fail if one lies in the requested band, never relabel a later mode.
valid = isfinite(lambda) & lambda>0 & residual.'<=1e-7 & abs(lambda-affine)<=1e-7*scale ...
    & self.qualifiedInverseDepthPairs(solver,V,r);
% A failed candidate below the positive physical spectrum belongs before
% every qualified mode; it must not disappear through sign classification.
lambda(~valid) = max(affine(~valid),realmin);
mapped = struct("values",lambda,"qualified",valid);
end
