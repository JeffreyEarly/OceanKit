function [F,G,h,varargout] = ModesFromGEP(self,A,B,varargin,options)
% Return modes from matrices defining a generalized eigenvalue problem.
%
% - Topic: Developer topics
% - Developer: true
% - Declaration: [F,G,h,varargout] = ModesFromGEP(self,A,B,varargin,options)
% - Parameter self: InternalModesSpectral instance
% - Parameter A: left generalized-eigenproblem matrix
% - Parameter B: right generalized-eigenproblem matrix
% - Parameter varargin: additional normalization or integral diagnostics
% - Parameter options.negativeEigenvalues: number of negative-depth modes to retain first
% - Returns F: horizontal-velocity mode matrix on `zOut`
% - Returns G: vertical-velocity mode matrix on `zOut`
% - Returns h: equivalent-depth row vector
% - Returns varargout: requested diagnostics
arguments
    self InternalModesSpectral
    A (:,:) double
    B (:,:) double
end
arguments (Repeating)
    varargin
end
arguments
    options.negativeEigenvalues (1,1) double = 0
end

[V_cheb,h] = self.solveGEP(A, B, negativeEigenvalues=options.negativeEigenvalues);

if self.nModes == 0
    maxModes = ceil(find(h > 0, 1, 'last')/2);
    if isempty(maxModes) || maxModes == 0
        error('InternalModesSpectral:NoUsableModes', 'No usable modes found. Try with higher resolution.');
    end
else
    maxModes = self.nModes;
end

if isempty(varargin)
    varargout = {};
    [F,G,h] = self.transformModesToSpatialDomain(V_cheb, h, maxModes);
else
    varargout = cell(size(varargin));
    [F,G,h,varargout{:}] = self.transformModesToSpatialDomain(V_cheb, h, maxModes, varargin{:});
end
end
