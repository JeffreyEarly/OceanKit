function du = diffZF(self,u,options)
% Differentiate an F-grid field with respect to z.
%
% `u` must use the gridded layout `[Nx Ny Nz]`. Orders 1 through 4 are
% supported. Odd orders return a G-representation and even orders return
% an F-representation.
%
% - Topic: Operations — Calculus
% - Declaration: du = diffZF(u,n=n)
% - Parameter u: F-grid field with dimensions `[Nx Ny Nz]`
% - Parameter n: derivative order from 1 through 4 (default 1)
% - Returns du: vertical derivative in the alternating F/G representation
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

if n == 1
    DzF = -(self.N2/self.g) .* self.QG0inv*(squeeze(self.Q0 ./ self.P0).*self.PF0);
    D = DzF;
elseif n == 2
    DzF = -(self.N2/self.g) .* self.QG0inv*(squeeze(self.Q0 ./ self.P0).*self.PF0);
    DzG = self.PF0inv*(squeeze(self.P0./(self.Q0 .* self.h_0)).*self.QG0);
    D = DzG*DzF;
elseif n == 3
    DzF = -(self.N2/self.g) .* self.QG0inv*(squeeze(self.Q0 ./ self.P0).*self.PF0);
    DzzG = -(self.N2/self.g) .* self.QG0inv*(squeeze(1./self.h_0).*self.QG0);
    D = DzzG*DzF;
elseif n == 4
    DzF = -(self.N2/self.g) .* self.QG0inv*(squeeze(self.Q0 ./ self.P0).*self.PF0);
    DzzG = -(self.N2/self.g) .* self.QG0inv*(squeeze(1./self.h_0).*self.QG0);
    DzG = self.PF0inv*(squeeze(self.P0./(self.Q0 .* self.h_0)).*self.QG0);
    D = DzG * DzzG * DzF;
end
du =  D*u;

du = reshape(du,self.Nz,self.Nx,self.Ny);
du = permute(du,[2 3 1]);

end
