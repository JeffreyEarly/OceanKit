function [varargout] = transformToPseudoRadialWavenumberA0(self,varargin)     
% transforms in the from (j,kRadial) to kPseudoRadial.
%
% transforms in the from (j,kRadial) to kPseudoRadial
% Sums all the variance/energy in radial bins `kPseudoRadial`.
%
% - Topic: Transformations — Axes
% - Declaration: [varargout] = transformToRadialWavenumber(varargin)
% - Parameter self: WVDiagnostics object
% - Parameter varargin: variables with dimensions $$(j,kl)$$
% - Returns varargout: variables with dimensions $$(kRadial)$$ or $$(kRadial,j)$$

% Thi is the final output axis for wavenumber

jWavenumber = 1./sqrt(self.Lr2);
jWavenumber(1) = 0; % barotropic mode is a mean?
[kj,kr] = ndgrid(jWavenumber,self.kRadial);
Kh = sqrt(kj.^2 + kr.^2);

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