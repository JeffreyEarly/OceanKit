function [V,lambda] = refineSmallInverseDepthModes(self,solver,V,lambda,IN)
% Refine near-zero inverse depths and small companion fields regularly.
% For weak stratification the finite external branch can be much smaller
% than the interior spectrum. Endpoint Rayleigh identities retain it
% without subtracting physical eigenvalues. Integrating G'=rF removes
% derivative-amplified roundoff when G is small relative to F.
if IN==0 || isempty(lambda), return; end
small = find(abs(lambda)<1e-6);
if isempty(small), return; end
n = numel(solver.zNative);
T = solver.physicalDerivativeMatrix(0);
F = T*V(1:n,small);
G = T*V(n+1:end,small);
s = solver.boundaryIndex("surface");
b = solver.boundaryIndex("bottom");
Fs = F(s,:);
Fb = F(b,:);
N2 = self.stratificationValues(solver.zNative);
w = solver.innerProductWeights(solver.zNative,self.zDomain);
volumeF = sum(w.*F.^2,1);
volumeG = sum(w.*N2.*G.^2,1)/self.g;
if self.name=="geostrophicAPVModes"
    sigma = double(self.parameters.surfaceBoundary=="freeSurface");
    surface = zeros(size(Fs));
    bottom = zeros(size(Fb));
    if self.parameters.g0~=0, surface = (sigma+self.g/self.parameters.g0)*Fs.^2; end
    if self.parameters.gd~=0, bottom = self.g/self.parameters.gd*Fb.^2; end
    r = (volumeG+surface+bottom)./volumeF;
else
    cS = 1-self.g*self.parameters.k^2/self.parameters.alpha0;
    cB = self.g*self.parameters.k^2/self.parameters.alphaD;
    r = (volumeG+cS*Fs.^2-cB*Fb.^2) ...
        ./(volumeF+self.f0^2/self.parameters.alpha0*Fs.^2+self.f0^2/self.parameters.alphaD*Fb.^2);
end
r(lambda(small)==0) = 0;
integralF = solver.integrateGridValuesFromSurface(F);
integralValues = solver.evaluateNativeModes(integralF,solver.zNative);
if self.name=="geostrophicAPVModes" && self.parameters.gd==0
    Gs = (sigma+self.g/self.parameters.g0)*Fs;
    G = Gs-r.*integralValues;
else
    if self.name=="geostrophicAPVModes"
        Gb = -self.g/self.parameters.gd*Fb;
    else
        Gb = (self.g*self.parameters.k^2+self.f0^2*r)/self.parameters.alphaD.*Fb;
    end
    G = Gb+r.*(integralValues(b,:)-integralValues);
end
V(n+1:end,small) = T\G;
lambda(small) = (r*(IN*diff(self.zDomain)/self.g)).';
end
