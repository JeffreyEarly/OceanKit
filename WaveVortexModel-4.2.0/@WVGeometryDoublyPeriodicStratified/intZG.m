function W = intZG(self,w,options)
arguments
    self        WVTransform
    w           double
    options.n (1,1)     double = 1
end
n = options.n;
didShift = false;
if size(w,3) == self.Nz
    w = permute(w,[3 1 2]); % keep adjacent in memory
    w = reshape(w,self.Nz,[]);
    didShift = true;
elseif size(w,1) ~= self.Nz
    error("incompatible size for operation");
end

if n == 1
    PF0inv = self.PF0inv - self.PF0inv(1,:);
    IntZG = - self.g*PF0inv*(squeeze(self.P0./self.Q0).*(self.QG0 .* shiftdim(1./self.N2,-1)));
    Int = IntZG;
else
    error('Not yet implemented')
end
W =  Int*w;

if didShift
    W = reshape(W,self.Nz,self.Nx,self.Ny);
    W = permute(W,[2 3 1]);
end

end