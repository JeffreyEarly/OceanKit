function [varargout] = transformToKePeAxis(self,varargin)
% Transform To Ke Pe Axis.
%
% transformToKePeAxis is part of the WVDiagnostics toolbox. Update this description to explain its purpose, inputs, outputs, and how it is used in the overall diagnostics workflow.
%
% - Topic: Transformations — Axes
% - Declaration: [varargout] = transformToKePeAxis(self,varargin)
% - Parameter self: WVDiagnostics object
% - Parameter varargin: input argument `varargin`
% - Returns varargout: output value `varargout`
    kePeAxis = reshape(self.kePeAxis,1,[]);
    kePeFraction = self.geo_hke_jk./(self.geo_hke_jk+self.geo_pe_jk);

    mid    = 0.5*(kePeAxis(1:end-1) + kePeAxis(2:end));
    edges  = [-Inf, mid, +Inf];

    bins = discretize(kePeFraction, edges);

    valid = ~isnan(bins);
    S = sparse(find(valid), bins(valid), 1, numel(kePeFraction), numel(kePeAxis), nnz(valid));

    varargout = cell(size(varargin));
    for iVar=1:length(varargin)
        varargout{iVar} =  reshape(varargin{iVar}(:).' * S,[],1);
    end

end