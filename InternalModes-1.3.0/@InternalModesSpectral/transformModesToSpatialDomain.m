function [F,G,h,varargout] = transformModesToSpatialDomain(self,V_cheb,h,maxModes,varargin)
% Transform Chebyshev eigenvectors onto the active output grid.
%
% - Topic: Developer topics
% - Developer: true
% - Declaration: [F,G,h,varargout] = transformModesToSpatialDomain(self,V_cheb,h,maxModes,varargin)
% - Parameter self: InternalModesSpectral instance
% - Parameter V_cheb: eigenvectors in Chebyshev coefficient space
% - Parameter h: equivalent-depth vector
% - Parameter maxModes: number of modes to transform
% - Parameter varargin: additional normalization or integral diagnostics
% - Returns F: horizontal-velocity mode matrix on `zOut`
% - Returns G: vertical-velocity mode matrix on `zOut`
% - Returns h: equivalent-depth row vector
% - Returns varargout: requested diagnostics

F = zeros(length(self.z), maxModes);
G = zeros(length(self.z), maxModes);
h = reshape(h(1:maxModes), 1, []);

varargout = cell(size(varargin));
for iArg = 1:length(varargin)
    varargout{iArg} = zeros(1, maxModes);
end

maxIndexZ = find(self.N2_xLobatto - self.gridFrequency*self.gridFrequency > 0, 1, 'first');
if maxIndexZ > 1
    maxIndexZ = maxIndexZ - 1;
elseif isempty(maxIndexZ)
    maxIndexZ = 1;
end

for j = 1:maxModes
    Fj = self.FFromVCheb(V_cheb(:,j), h(j));
    Gj = self.GFromVCheb(V_cheb(:,j), h(j));
    switch self.normalization
        case Normalization.uMax
            A = max(abs(Fj));
        case Normalization.wMax
            A = max(abs(Gj));
        case Normalization.kConstant
            A = sqrt(self.GNorm(Gj));
        case Normalization.omegaConstant
            A = sqrt(self.FNorm(Fj));
        case Normalization.geostrophic
            A = sqrt(self.GeostrophicNorm(Gj));
    end
    if Fj(maxIndexZ) < 0
        A = -A;
    end

    G(:,j) = self.GOutFromVCheb(V_cheb(:,j), h(j))/A;
    F(:,j) = self.FOutFromVCheb(V_cheb(:,j), h(j))/A;

    for iArg = 1:length(varargin)
        switch string(varargin{iArg})
            case "F2"
                varargout{iArg}(j) = self.Lz*self.FNorm(Fj/A);
            case "G2"
                varargout{iArg}(j) = self.Lz*self.FNorm(Gj/A);
            case "N2G2"
                varargout{iArg}(j) = self.g*(self.GNorm(Gj/A)-Gj(1)*Gj(1)) + self.f0*self.f0*self.Lz*self.FNorm(Gj/A);
            case "uMax"
                B = max(abs(Fj));
                varargout{iArg}(j) = abs(A/B);
            case "wMax"
                B = max(abs(Gj));
                varargout{iArg}(j) = abs(A/B);
            case "kConstant"
                B = sqrt(self.GNorm(Gj));
                varargout{iArg}(j) = abs(A/B);
            case "omegaConstant"
                B = sqrt(self.FNorm(Fj));
                varargout{iArg}(j) = abs(A/B);
            case "geostrophicNorm"
                B = sqrt(self.GeostrophicNorm(Gj));
                varargout{iArg}(j) = abs(A/B);
            case "int_N2_G_dz/g"
                integrand = (1/self.g) * self.N2_xLobatto .* (Gj/A);
                varargout{iArg}(j) = sum(self.Int_xCheb .* InternalModesSpectral.fct(integrand));
            case "int_F_dz"
                varargout{iArg}(j) = sum(self.Int_xCheb .* InternalModesSpectral.fct(Fj/A));
            case "int_G_dz"
                varargout{iArg}(j) = sum(self.Int_xCheb .* InternalModesSpectral.fct(Gj/A));
            otherwise
                error('InternalModesSpectral:InvalidDiagnostic', ...
                    'Invalid diagnostic request "%s".', char(string(varargin{iArg})));
        end
    end
end
end
