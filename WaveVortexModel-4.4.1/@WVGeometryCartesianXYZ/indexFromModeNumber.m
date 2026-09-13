function index = indexFromModeNumber(self,kMode,lMode,jMode)
% return the linear index into a spectral matrix given (k,l,j)
%
% This function will return the linear index in a spectral
% matrix given a mode number. Scalar and column-vector inputs preserve
% their shape and ordering; conjugate mode numbers map to their primary
% coefficient.
%
% - Topic: Index Gymnastics
% - Declaration: index = indexFromModeNumber(kMode,lMode,jMode)
% - Parameter kMode: integer
% - Parameter lMode: integer
% - Parameter jMode: integer vertical mode number present in `self.j`
% - Returns index: a non-negative integer number
arguments (Input)
    self WVTransform {mustBeNonempty}
    kMode (:,1) double {mustBeInteger}
    lMode (:,1) double {mustBeInteger}
    jMode (:,1) double {mustBeInteger,mustBeNonnegative}
end
arguments (Output)
    index (:,1) double {mustBeInteger,mustBePositive}
end
if ~all(self.isValidModeNumber(kMode,lMode,jMode))
    error('Invalid WV mode number!');
end
[kMode,lMode] = self.primaryKLModeNumberFromKLModeNumber(kMode,lMode);
klIndex = self.indexFromKLModeNumber(kMode,lMode);
[isPresent,jIndex] = ismember(jMode,self.j);
if ~all(isPresent)
    error('Invalid WV vertical mode number!');
end
index = sub2ind(self.spectralMatrixSize,jIndex,klIndex);
end
