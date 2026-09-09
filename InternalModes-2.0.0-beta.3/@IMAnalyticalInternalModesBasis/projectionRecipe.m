function recipe = projectionRecipe(self,options)
% Bind an analytical aligned basis without requiring a numerical solver.
%
% Gram targets use the existing 1024-point trapezoidal reference rule;
% analytical evaluation does not make this quadrature convergence-qualified.
% - Topic: Analyze Gram matrices
% - Parameter options.variable: "F" or "G", the solved formulation by default
% - Returns recipe: continuous basis-owned projection recipe
arguments (Input)
    self IMAnalyticalInternalModesBasis
    options.variable (1,1) string {mustBeMember(options.variable,["F","G"])} = string(self.evp.formulation)
end
arguments (Output)
    recipe (1,1) struct
end
spec = self.evp.innerProduct(options.variable);
provenance = struct("representation","analytical","referenceIntegration","1024-point trapezoidal rule","referencePointCount",1024,"solveAccuracy","unverified","referenceConvergence","unverified");
if ~spec.hasInnerProduct
    recipe = IMProjection.createRecipe(isAvailable=false,reason=string(spec.reason),variable=options.variable,columnLabels=string(self.modeNumber),normalization=self.normalizationName(self.normalization),provenance=provenance,zDomain=self.zDomain);
    return;
end
context = self.context();
factors = self.normalizationFactors(self.normalization);
target = self.gramMatrix(variable=options.variable);
majorant = self.majorantGramMatrix(variable=options.variable);
active = abs(diag(target)) > 1e3*eps(max(1,diag(majorant)));
[~,notPositive] = chol(target(active,active));
recipe = IMProjection.createRecipe(variable=options.variable,columnLabels=string(self.modeNumber),normalization=self.normalizationName(self.normalization),targetGramMatrix=target,majorantGramMatrix=majorant, ...
    supportsLeakage=notPositive==0,supportsQuadratic=true,provenance=provenance,zDomain=self.zDomain, ...
    evaluateFunction=@(z) self.rawVariable(options.variable,z)./factors,weightFunction=@(z) IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight,z,context),spec=spec);
end
