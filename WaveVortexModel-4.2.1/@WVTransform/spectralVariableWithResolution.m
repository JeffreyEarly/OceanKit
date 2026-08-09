function [varargout] = spectralVariableWithResolution(self,wvtX2,varargin)
% create a new variable with different resolution
%
% Given a variable with dimensions `[Nj Nkl]`, this returns a new variable
% with dimensions matching `wvtX2`. Coefficients are matched by their
% integer `(kMode,lMode,jMode)` identities. Modes absent from the target are
% discarded and target modes absent from the source are initialized to zero.
%
% - Topic: Utility function
% - Declaration: varX2 = spectralVariableWithResolution(wvtX2,var)
% - Parameter var: a variable with dimensions [Nj Nkl]
% - Parameter wvtX2: a WVTransform of different size.
% - Returns varX2: matrix the size Nklj

if ~isequal(self.Lx,wvtX2.Lx) || ~isequal(self.Ly,wvtX2.Ly) || ~isequal(self.Lz,wvtX2.Lz)
    error('These transforms are not compatible.')
end

sourceK = repmat(self.kMode_wv.',self.Nj,1);
sourceL = repmat(self.lMode_wv.',self.Nj,1);
sourceJ = repmat(self.j,1,self.Nkl);
targetK = repmat(wvtX2.kMode_wv.',wvtX2.Nj,1);
targetL = repmat(wvtX2.lMode_wv.',wvtX2.Nj,1);
targetJ = repmat(wvtX2.j,1,wvtX2.Nkl);
[isCommon,sourceIndex] = ismember([targetK(:),targetL(:),targetJ(:)],[sourceK(:),sourceL(:),sourceJ(:)],'rows');

varargout = cell(size(varargin));
for iVar=1:length(varargin)
    if isempty(varargin{iVar})
        varargout{iVar} = [];
    else
        varargout{iVar} = zeros(wvtX2.spectralMatrixSize);
        varargout{iVar}(isCommon) = varargin{iVar}(sourceIndex(isCommon));
    end
end

end
