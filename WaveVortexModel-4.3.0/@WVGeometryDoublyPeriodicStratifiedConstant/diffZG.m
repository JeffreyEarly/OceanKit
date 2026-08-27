function du = diffZG(self,u,options)
% Differentiate a G-grid field with respect to z.
%
% `u` must use the gridded layout `[Nx Ny Nz]`. Orders 1 through 4 are
% supported. Odd orders return an F-representation and even orders return
% a G-representation.
%
% - Topic: Operations — Calculus
% - Declaration: du = diffZG(u,n=n)
% - Parameter u: G-grid field with dimensions `[Nx Ny Nz]`
% - Parameter n: derivative order from 1 through 4 (default 1)
% - Returns du: vertical derivative in the alternating G/F representation
arguments
    self        WVTransform
    u           double
    options.n (1,1) double {mustBeMember(options.n,1:4)} = 1
end
mustBeMember(ndims(u),3);
mustBeMember(size(u,1),self.Nx);
mustBeMember(size(u,2),self.Ny);
mustBeMember(size(u,3),self.Nz);
n = options.n;
u = permute(u,[3 1 2]); % keep adjacent in memory
u = reshape(u,self.Nz,[]);

% sine goes to, [1,-1,-1,1] for numDerivs = [1,2,3,4]
thesign = [1,-1,-1,1];
m = reshape(pi*self.j/self.Lz,[],1);
du_bar = thesign(mod(n-1,4)+1)*(m.^n) .* (self.DST*u);

if mod(n,2) == 0
    du = self.iDST*du_bar;
else
    du = self.iDCT*du_bar;
end
du = reshape(du,self.Nz,self.Nx,self.Ny);
du = permute(du,[2 3 1]);

end
