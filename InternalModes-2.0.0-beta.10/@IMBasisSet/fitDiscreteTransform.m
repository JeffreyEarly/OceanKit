function [transform, assessment] = fitDiscreteTransform(self, options)
% Fit and certify one exact modal band on caller-selected points.
%
% This is the strict, diagnostic API. `modeCount` is the number of leading
% family modes whose quadrature weights are fitted. The same exact band is
% then required to pass every enabled policy. Construction throws rather
% than silently returning a shorter prefix.
%
% Use `certifiedDiscreteTransform` when the retained count should be chosen
% automatically. Use `modeRootGrid` when the physical points should be
% designed from a particular modal family before fitting one or more
% families independently on that shared grid.
%
% - Topic: Build discrete transforms
% - Declaration: [transform,assessment] = fitDiscreteTransform(basisSet,options)
% - Parameter options.z: increasing physical sample points
% - Parameter options.modeCount: exact number of leading modes to fit and certify
% - Parameter options.gridDesign: optional provenance returned by `modeRootGrid`
% - Parameter options.gramTolerance: normalized-Gram operator-error tolerance
% - Parameter options.leakageTolerance: optional rejected-mode leakage tolerance
% - Parameter options.quadraticAliasingTolerance: optional quadratic-aliasing tolerance
% - Parameter options.nCheckModes: optional rejected-mode check count
% - Returns transform: exact fitted transform
% - Returns assessment: exact-band diagnostics and grid provenance
arguments
    self IMBasisSet
    options.z (:,1) double {mustBeReal, mustBeFinite}
    options.modeCount (1,1) double {mustBeInteger, mustBePositive}
    options.gridDesign struct = struct.empty
    options.gramTolerance (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative} = 1e-2
    options.leakageTolerance double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.quadraticAliasingTolerance double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.nCheckModes double {mustBeReal, mustBeFinite} = zeros(0,1)
end

gridDesign = IMDiscreteTransformTools.validatedGridDesign(self,options.z,options.gridDesign);
[transform,assessment] = self.discreteTransform(z=options.z,nModes=options.modeCount, ...
    gramTolerance=options.gramTolerance,leakageTolerance=options.leakageTolerance, ...
    quadraticAliasingTolerance=options.quadraticAliasingTolerance,nCheckModes=options.nCheckModes);
assessment = assessment.withCertificationMetadata(gridDesign,table());
end
