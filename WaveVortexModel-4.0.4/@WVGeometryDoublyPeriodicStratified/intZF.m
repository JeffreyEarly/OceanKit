function U = intZF(self,u,options)
arguments
    self        WVTransform
    u (:,:,:)   double
    options.n (1,1)     double = 1
end
n = options.n;
didShift = false;
if size(u,3) == self.Nz
    u = permute(u,[3 1 2]); % keep adjacent in memory
    u = reshape(u,self.Nz,[]);
    didShift = true;
elseif size(u,1) ~= self.Nz
    error("incompatible size for operation");
end
if n == 1
    IntZF = self.QG0inv*(squeeze(self.h_0 .* self.Q0 ./ self.P0).*self.PF0);
    Int = IntZF;
else
    error('Not yet implemented')
end
U =  Int*u;

if didShift
    U = reshape(U,self.Nz,self.Nx,self.Ny);
    U = permute(U,[2 3 1]);
end

end