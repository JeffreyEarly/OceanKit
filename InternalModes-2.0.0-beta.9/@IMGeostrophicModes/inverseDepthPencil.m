function [A,B] = inverseDepthPencil(self,solver)
% Assemble the regular pair directly in inverse-depth units.
z = solver.zNative(:);
n = numel(z);
N2 = self.stratificationValues(z);
V = solver.physicalDerivativeMatrix(0);
D = solver.physicalDerivativeMatrix(1);
A = [D diag(N2/self.g)*V; zeros(n) D];
B = [zeros(n,2*n); V zeros(n)];
s = solver.boundaryIndex("surface");
b = solver.boundaryIndex("bottom");
B([s n+b],:) = 0;
if self.name == "geostrophicAPVModes"
    sigma = double(self.parameters.surfaceBoundary=="freeSurface");
    g0 = self.parameters.g0;
    gd = self.parameters.gd;
    if isinf(g0)
        A(s,:) = [-sigma*V(s,:) V(s,:)];
    else
        A(s,:) = [-(self.g+sigma*g0)*V(s,:) g0*V(s,:)];
    end
    if isinf(gd)
        A(n+b,:) = [zeros(1,n) V(b,:)];
    else
        A(n+b,:) = [self.g*V(b,:) gd*V(b,:)];
    end
else
    k2 = self.parameters.k^2;
    a0 = self.parameters.alpha0;
    aD = self.parameters.alphaD;
    A(s,:) = [(self.g*k2/a0-1)*V(s,:) V(s,:)];
    B(s,:) = [-self.f0^2/a0*V(s,:) zeros(1,n)];
    A(n+b,:) = [-self.g*k2/aD*V(b,:) V(b,:)];
    B(n+b,:) = [self.f0^2/aD*V(b,:) zeros(1,n)];
end
end
