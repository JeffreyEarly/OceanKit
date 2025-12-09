function [varargout] = transformToPseudoRadialWavenumberApm(self,varargin)     
% transforms Ap/Am modes in the from (j,kRadial) to kPseudoRadial.
%
% transforms Ap/Am modes in the from (j,kRadial) to kPseudoRadial
% Sums all the variance/energy in radial bins `kPseudoRadial`.
%
% - Topic: Transformations — Axes
% - Declaration: [varargout] = transformToRadialWavenumber(varargin)
% - Parameter self: WVDiagnostics object
% - Parameter varargin: variables with dimensions $$(j,kl)$$
% - Returns varargout: variables with dimensions $$(kRadial)$$ or $$(kRadial,j)$$

wvt = self.wvt;
if size(wvt.h_pm,2) > 1
    kj2 = 1./self.Lr2_pm;
    kj2(1,1) = 0;
    kr2 = repmat(reshape(self.kRadial.^2,1,[]),[length(self.j) 1]);
    Kh = sqrt(kj2 + kr2);
else
    [kj,kr] = ndgrid(self.jWavenumber,self.kRadial);
    kj(1,1) = 0;
    Kh = sqrt(kj.^2 + kr.^2);
end

k = self.kPseudoRadial;
edges  = [-Inf; 0.5*(k(1:end-1) + k(2:end)); +Inf];
bin = discretize(Kh, edges);
valid = ~isnan(bin);
S = sparse(find(valid), bin(valid), 1, numel(Kh), numel(k), nnz(valid));


varargout = cell(size(varargin));
spectralMatrixSize = [length(self.j) length(self.kRadial)];
for iVar=1:length(varargin)
    if size(varargin{iVar},2) ~= spectralMatrixSize(2)
        error('The input matrix must be of size [Nj NkRadial]');
    end
    
    varargout{iVar} = reshape((varargin{iVar}(:)).' * S,[],1);
end

end