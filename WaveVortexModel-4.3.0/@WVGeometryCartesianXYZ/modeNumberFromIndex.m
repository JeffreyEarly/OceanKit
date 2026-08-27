function [kMode,lMode,jMode] = modeNumberFromIndex(self,linearIndex)
% Return mode numbers for spectral linear indices.
%
% Scalar and column-vector inputs preserve their shape and ordering. Each
% index must lie within the current spectral matrix.
%
% - Topic: Index Gymnastics
% - Declaration: [kMode,lMode,jMode] = modeNumberFromIndex(linearIndex)
% - Parameter linearIndex: positive integer scalar or column vector
% - Returns kMode: integer scalar or column vector
% - Returns lMode: integer scalar or column vector
% - Returns jMode: integer scalar or column vector
arguments (Input)
    self WVTransform {mustBeNonempty}
    linearIndex (:,1) double {mustBeInteger,mustBePositive}
end
arguments (Output)
    kMode (:,1) double {mustBeInteger}
    lMode (:,1) double {mustBeInteger}
    jMode (:,1) double {mustBeInteger,mustBeNonnegative}
end
mustBeLessThanOrEqual(linearIndex,prod(self.spectralMatrixSize))
[jIndex,klIndex] = ind2sub(self.spectralMatrixSize,linearIndex);
[kMode,lMode] = self.klModeNumberFromIndex(klIndex);
jMode = self.j(jIndex);
end
