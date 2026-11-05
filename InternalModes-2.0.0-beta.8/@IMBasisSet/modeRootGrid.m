function [z, gridDesign] = modeRootGrid(self, options)
% Design a mode-root grid and report how it was generated.
%
% Specify either an exact physical point count or the represented leading
% mode count. The returned `gridDesign` records the source EVP and family,
% generating variable and mode, and the interpretation of the points for
% internal-mode `G` structures.
%
% ```matlab
% [z,gridDesign] = basisSet.modeRootGrid(nPoints=128);
% gridDesign.generatingVariable
% gridDesign.interpretationForG
% ```
%
% - Topic: Build discrete transforms
% - Declaration: [z,gridDesign] = modeRootGrid(basisSet,options)
% - Parameter options.nPoints: exact requested physical point count
% - Parameter options.modeCount: leading mode count represented by the generating mode
% - Returns z: increasing physical grid
% - Returns gridDesign: scalar provenance struct
arguments
    self IMBasisSet
    options.nPoints double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.modeCount double {mustBeReal, mustBeFinite} = zeros(0,1)
end

IMDiscreteTransformTools.validateOptionalPositiveInteger(options.nPoints,"nPoints");
IMDiscreteTransformTools.validateOptionalPositiveInteger(options.modeCount,"modeCount");
if isempty(options.nPoints) == isempty(options.modeCount)
    error("IMBasisSet:InvalidModeRootGridSpecification", "Specify exactly one of nPoints or modeCount.");
end
if ~isempty(options.nPoints)
    [z,~,gridDesign] = IMDiscreteTransformTools.pointsForExactCount(self,options.nPoints,size(self.nativeModes,2));
else
    [z,gridDesign] = self.pointsFromModeRoots(nModes=options.modeCount);
end
end
