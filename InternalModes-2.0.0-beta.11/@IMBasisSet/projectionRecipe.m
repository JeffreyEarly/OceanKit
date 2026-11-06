function recipe = projectionRecipe(self,options)
% Bind the scalar signed pairing and continuous evaluation to a recipe.
%
% Full targets retain the actual continuous pairing of the represented columns.
% A signed scalar pairing supports Gram measurements but not leakage or quadratic policies.
% Reference integration does not independently qualify the eigenproblem.
% - Topic: Build discrete transforms
% - Parameter options.variable: solved scalar variable, "u"
% - Returns recipe: continuous basis-owned projection recipe
arguments (Input)
    self IMBasisSet
    options.variable (1,1) string {mustBeMember(options.variable,"u")} = "u"
end
arguments (Output)
    recipe (1,1) struct
end
spec = self.evp.innerProduct();
context = self.evp.contextForSolver(self.solver);
weightFunction = @(z) IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight,z,context);
[target,majorant,isPositiveMetric] = self.scalarGramMatrices(self.zDomain);
provenance = struct("representation","numerical","referenceIntegration","basis reference quadrature","solverClass",string(class(self.solver)),"solveAccuracy","unverified","referenceConvergence","unverified");
recipe = IMProjection.createRecipe(variable=options.variable,columnLabels=self.columnLabels,normalization=self.normalizationName(self.normalization),targetGramMatrix=target,majorantGramMatrix=majorant, ...
    supportsLeakage=isPositiveMetric,supportsQuadratic=isPositiveMetric,provenance=provenance,zDomain=self.zDomain,evaluateFunction=@(z) self.u(z),weightFunction=weightFunction,spec=spec);
end
