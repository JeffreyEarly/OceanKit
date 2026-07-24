function W = intZG(self,w,options)
arguments
    self        WVTransform
    w           double
    options.n (1,1)     double = 1
end
error("Not yet implemented for constant stratification")
% n = options.n;
% if size(w,3) == self.Nz
%     w = permute(w,[3 1 2]); % keep adjacent in memory
%     w = reshape(w,self.Nz,[]);
% elseif size(w,1) ~= self.Nz
%     error("incompatible size for operation");
% end
% 
% if n == 1
%     PF0inv = self.PF0inv - self.PF0inv(1,:);
%     IntZG = - self.g*PF0inv*(squeeze(self.P0./self.Q0).*(wvt.QG0 .* shiftdim(1./wvt.N2,-1)));
%     Int = IntZG;
% else
%     error('Not yet implemented')
% end
% W =  Int*w;
% 
% W = reshape(W,self.Nz,self.Nx,self.Ny);
% W = permute(W,[2 3 1]);

end