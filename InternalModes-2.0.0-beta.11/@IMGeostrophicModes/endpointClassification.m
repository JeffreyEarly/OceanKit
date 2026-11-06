function [negativeCount,zeroPresent,IN] = endpointClassification(self,solver)
% Use endpoint inertia and the regular zero-inverse-depth compatibility.
IN = integral(@(z) reshape(self.stratificationValues(z),size(z)),self.zDomain(1),self.zDomain(2),AbsTol=1e-13,RelTol=1e-12);
if IN==0
    if self.name=="geostrophicAPVModes" && (self.parameters.g0==0 || self.parameters.gd==0)
        negativeCount = 0;
        zeroPresent = false;
    else
        r = self.fullyNeutralMode(solver);
        negativeCount = double(self.name=="geostrophicAPVModes" && r<0);
        zeroPresent = r==0;
    end
elseif self.name=="geostrophicAPVModes"
    betaS = double(self.parameters.surfaceBoundary=="freeSurface")/self.g+1/self.parameters.g0;
    betaB = 1/self.parameters.gd;
    [negativeCount,zeroPresent] = IMGeostrophicAPVEndpointInertia(IN,betaS,betaB,self.parameters.g0~=0,self.parameters.gd~=0);
else
    cS = 1-self.g*self.parameters.k^2/self.parameters.alpha0;
    cB = self.g*self.parameters.k^2/self.parameters.alphaD;
    if isinf(self.parameters.alphaD) && (self.parameters.alpha0==self.g*self.parameters.k^2 || self.parameters.alpha0==(self.g*self.parameters.k)*self.parameters.k)
        cS = 0;
    end
    difference = cS-cB;
    product = cS*cB*IN/self.g;
    zeroPresent = abs(difference-product)<=1e3*eps*(abs(difference)+abs(product));
    negativeCount = 0;
end
end
