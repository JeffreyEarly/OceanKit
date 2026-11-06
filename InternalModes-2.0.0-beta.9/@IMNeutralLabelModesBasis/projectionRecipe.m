function recipe = projectionRecipe(self,options)
% Describe the prescribed continuous APV projector and its exact reference.
%
% This metadata-only recipe is unavailable to the generic Galerkin
% constructor. Collection construction delegates to projectionOnGrid.
% Fixed endpoint labels give strictly positive vertical wavenumbers. APV
% has a positive continuous norm for every retained label when f0 is nonzero.
% - Topic: Project APV
% - Parameter options.variable: APV only
% - Returns recipe: family-owned capability and reference metadata
arguments (Input)
    self IMNeutralLabelModesBasis
    options.variable (1,1) string {mustBeMember(options.variable,"APV")} = "APV"
end
arguments (Output)
    recipe (1,1) struct
end
provenance = struct("representation","analyticalNeutralLabels","projection","continuousAPVDual","coefficientSystem","identity after dividing each APV pairing by its exact norm","continuousAPVGramMatrix",self.gramMatrix(variable="APV"),"referenceIntegration","exact","boundaryConditions","Dirichlet at both interval endpoints","g0",self.g0,"gd",self.gd,"endpointDerivativeTrace","inside neutral interval");
recipe = IMProjection.createRecipe(isAvailable=false,reason="Neutral labels require their prescribed continuous APV dual. Use basis.projectionOnGrid or collection.projectionOnGrid.",variable=options.variable,columnLabels=self.columnLabels,normalization=self.normalization,targetGramMatrix=eye(numel(self.m_n)),majorantGramMatrix=self.gramMatrix(variable="APV"),supportsLeakage=false,supportsQuadratic=false,provenance=provenance,zDomain=self.zDomain);
end
