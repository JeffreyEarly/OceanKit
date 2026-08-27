function U = intZF(self,u,options)
% Return the first antiderivative of an F-representation.
%
% The result is the antiderivative representable in G space. The
% barotropic F mode is removed because it has no corresponding G mode, so
% the result vanishes at both vertical boundaries. `u` may use the gridded
% layout `[Nx Ny Nz]` or a vertical-first matrix `[Nz N]`; the returned
% array preserves that layout.
%
% - Topic: Operations — Calculus
% - Declaration: U = intZF(u,n=1)
% - Parameter u: F-representation in `[Nx Ny Nz]` or `[Nz N]` layout
% - Parameter n: antiderivative order; only 1 is supported (default 1)
% - Returns U: G-representation antiderivative in the input layout
arguments
    self        WVTransform
    u           double
    options.n (1,1) double = 1
end
mustBeMember(options.n,1);
didShift = false;
if ndims(u) == 3
    mustBeMember(size(u,1),self.Nx);
    mustBeMember(size(u,2),self.Ny);
    mustBeMember(size(u,3),self.Nz);
    u = permute(u,[3 1 2]); % keep adjacent in memory
    u = reshape(u,self.Nz,[]);
    didShift = true;
else
    mustBeMatrix(u);
    mustBeMember(size(u,1),self.Nz);
end
IntZF = self.QG0inv*(squeeze(self.h_0 .* self.Q0 ./ self.P0).*self.PF0);
U = IntZF*u;

if didShift
    U = reshape(U,self.Nz,self.Nx,self.Ny);
    U = permute(U,[2 3 1]);
end

end
