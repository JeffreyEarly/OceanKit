function [r,coefficients] = fullyNeutralMode(self,solver)
% Deflate the singular neutral pencil to its constant-F, affine-G family.
D = diff(self.zDomain);
z = solver.zNative(:);
F = ones(size(z));
if self.name == "geostrophicAPVModes"
    if self.parameters.g0==0 || self.parameters.gd==0
        error("IMSolver:NoValidEigenvalues","The fully neutral APV problem has no F-normalizable mode when an endpoint imposes F=0.");
    end
    betaS = double(self.parameters.surfaceBoundary=="freeSurface")/self.g+1/self.parameters.g0;
    betaB = 1/self.parameters.gd;
    r = self.g*(betaS+betaB)/D;
    G = -self.g*betaB+r*(z-self.zDomain(1));
else
    a = 1/self.parameters.alpha0+1/self.parameters.alphaD;
    cS = 1-self.g*self.parameters.k^2/self.parameters.alpha0;
    cB = self.g*self.parameters.k^2/self.parameters.alphaD;
    if isinf(self.parameters.alphaD) && (self.parameters.alpha0==self.g*self.parameters.k^2 || self.parameters.alpha0==(self.g*self.parameters.k)*self.parameters.k)
        cS = 0;
    end
    r = (cS-cB)/(D+self.f0^2*a);
    G = (self.g*self.parameters.k^2+self.f0^2*r)/self.parameters.alphaD+r*(z-self.zDomain(1));
end
V = solver.physicalDerivativeMatrix(0);
coefficients = [V\F;V\G];
end
