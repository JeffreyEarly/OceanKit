function recipe = projectionRecipe(self,options)
% Report the unavailable scalar F/G projection for boundary-response bases.
%
% Endpoint response and energy matrices are coefficient-space forms. They
% do not define projection of arbitrary scalar F or G samples. Endpoint
% labels remain intact; no modal-prefix truncation is inferred.
% - Topic: Evaluate geostrophic zero-APV modes
% - Parameter options.variable: requested scalar variable, "F" or "G"
% - Returns recipe: unavailable recipe with an explicit reason
arguments (Input)
    self IMAnalyticalGeostrophicZeroAPVModesBasis
    options.variable (1,1) string {mustBeMember(options.variable,["F","G"])} = "F"
end
arguments (Output)
    recipe (1,1) struct
end
reason = "Zero-APV coefficient-space response and energy forms do not define a scalar F/G sample projection. Paired endpoint-observation projection is not implemented.";
provenance = struct("representation","analytical","solveAccuracy","unverified","requestedKappa",self.k);
recipe = IMProjection.createRecipe(isAvailable=false,reason=reason,variable=options.variable,columnLabels=self.endpoints,normalization=string(self.normalizationConvention),provenance=provenance,zDomain=self.zDomain);
end
