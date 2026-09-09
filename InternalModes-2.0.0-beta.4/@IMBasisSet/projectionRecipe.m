function recipe = projectionRecipe(self,options)
% Bind the scalar signed pairing and continuous evaluation to a recipe.
%
% Scalar targets retain the existing diagonal convention. A signed scalar
% pairing supports Gram measurements but not leakage or quadratic policies.
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
gram = self.gramMatrix();
target = diag(diag(gram));
% Build the scalar majorant on the same reference rule, replacing signed
% interior and endpoint contributions individually before summing.
z = self.solver.innerProductGrid(self.zDomain);
values = self.u(z);
weight = weightFunction(z);
if isscalar(weight)
    weight = repmat(weight,size(z));
end
majorant = zeros(size(gram));
for i = 1:size(values,2)
    for j = i:size(values,2)
        majorant(i,j) = self.solver.integrateInnerProduct(z,abs(weight(:)).*values(:,i).*values(:,j),self.zDomain);
        majorant(j,i) = majorant(i,j);
    end
end
terms = self.endpointGramTerms();
isPositiveMetric = all(weight(:) >= 0);
for i = 1:numel(terms)
    endpointValues = terms(i).values(:);
    majorant = majorant+abs(terms(i).coefficient)*(endpointValues*endpointValues.');
    isPositiveMetric = isPositiveMetric && terms(i).coefficient >= 0;
end
provenance = struct("representation","numerical","referenceIntegration","solver inner-product rule","solverClass",string(class(self.solver)),"solveAccuracy","unverified","referenceConvergence","unverified");
recipe = IMProjection.createRecipe(variable=options.variable,columnLabels=string(self.modeNumber),normalization=self.normalizationName(self.normalization),targetGramMatrix=target,majorantGramMatrix=majorant, ...
    supportsLeakage=isPositiveMetric,supportsQuadratic=isPositiveMetric,provenance=provenance,zDomain=self.zDomain,evaluateFunction=@(z) self.u(z),weightFunction=weightFunction,spec=spec);
end
