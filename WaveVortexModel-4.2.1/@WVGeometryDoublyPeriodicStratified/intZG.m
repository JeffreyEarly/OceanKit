function W = intZG(self,w,options)
% Return the bottom-zero first antiderivative of a G-representation.
%
% A G-to-F antiderivative is defined up to an additive constant. This
% method selects the representative that vanishes at the bottom boundary.
% `w` may use the gridded layout `[Nx Ny Nz]` or a vertical-first matrix
% `[Nz N]`; the returned array preserves that layout.
%
% - Topic: Operations — Calculus
% - Declaration: W = intZG(w,n=1)
% - Parameter w: G-representation in `[Nx Ny Nz]` or `[Nz N]` layout
% - Parameter n: antiderivative order; only 1 is supported (default 1)
% - Returns W: bottom-zero F-representation antiderivative in the input layout
arguments
    self        WVTransform
    w           double
    options.n (1,1) double = 1
end
mustBeMember(options.n,1);
didShift = false;
if ndims(w) == 3
    mustBeMember(size(w,1),self.Nx);
    mustBeMember(size(w,2),self.Ny);
    mustBeMember(size(w,3),self.Nz);
    w = permute(w,[3 1 2]); % keep adjacent in memory
    w = reshape(w,self.Nz,[]);
    didShift = true;
else
    mustBeMatrix(w);
    mustBeMember(size(w,1),self.Nz);
end

PF0inv = self.PF0inv - self.PF0inv(1,:);
IntZG = - self.g*PF0inv*(squeeze(self.P0./self.Q0).*(self.QG0 .* shiftdim(1./self.N2,-1)));
W = IntZG*w;

if didShift
    W = reshape(W,self.Nz,self.Nx,self.Ny);
    W = permute(W,[2 3 1]);
end

end
