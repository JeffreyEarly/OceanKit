function [transform, assessment] = fitDiscreteTransform(self, options)
% Fit and certify one exact aligned F/G family band.
%
% `modeCount` is explicit and strict: the weights are fitted to exactly
% that many aligned family columns, and all requested variables and enabled
% policies must accept the complete band. Use `certifiedDiscreteTransform`
% to select a retained count automatically.
%
% - Topic: Build discrete transforms
% - Declaration: [transform,assessment] = fitDiscreteTransform(basisSet,options)
% - Parameter options.z: increasing physical sample points
% - Parameter options.modeCount: exact number of leading family modes
% - Parameter options.variables: requested direct channels, F and/or G
% - Parameter options.gridDesign: optional provenance returned by `modeRootGrid`
% - Parameter options.gramTolerance: per-channel normalized-Gram tolerance
% - Parameter options.leakageTolerance: optional rejected-mode leakage tolerance
% - Parameter options.quadraticAliasingTolerance: optional coupled-product tolerance
% - Parameter options.nCheckModes: optional rejected-mode check count
% - Returns transform: exact fitted aligned transform
% - Returns assessment: exact-band diagnostics and grid provenance
arguments
    self IMInternalModesBasis
    options.z (:,1) double {mustBeReal, mustBeFinite}
    options.modeCount (1,1) double {mustBeInteger, mustBePositive}
    options.variables (1,:) string = strings(1,0)
    options.gridDesign struct = struct.empty
    options.gramTolerance (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative} = 1e-2
    options.leakageTolerance double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.quadraticAliasingTolerance double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.nCheckModes double {mustBeReal, mustBeFinite} = zeros(0,1)
end

gridDesign = IMDiscreteTransformTools.validatedGridDesign(self,options.z,options.gridDesign);
[transform,assessment] = self.discreteTransform(z=options.z,nModes=options.modeCount,variables=options.variables, ...
    gramTolerance=options.gramTolerance,leakageTolerance=options.leakageTolerance, ...
    quadraticAliasingTolerance=options.quadraticAliasingTolerance,nCheckModes=options.nCheckModes);
assessment = assessment.withCertificationMetadata(gridDesign,table());
end
