function recipe = projectionRecipe(self,options)
% Bind an analytical aligned basis without requiring a numerical solver.
%
% Gram targets use the basis's physical Gauss-Legendre reference rule;
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
referencePoints = self.integrationGrid(self.zDomain);
provenance = struct("representation","analytical","referenceIntegration","Gauss-Legendre rule","referencePointCount",numel(referencePoints),"solveAccuracy","unverified","referenceConvergence","unverified");
if ~spec.hasInnerProduct
    recipe = IMProjection.createRecipe(isAvailable=false,reason=string(spec.reason),variable=options.variable,columnLabels=self.columnLabels,normalization=self.normalizationName(self.normalization),provenance=provenance,zDomain=self.zDomain);
    return;
end
context = self.context();
factors = self.normalizationFactors(self.normalization);
target = self.gramMatrix(variable=options.variable);
majorant = self.majorantGramMatrix(variable=options.variable);
active = true(1,numel(self.modeNumber));
surface = self.evp.surfaceBoundary;
bottom = self.evp.bottomBoundary;
ordinaryNeumann = surface.a == 0 && surface.c == 0 && surface.d == 0 ...
    && bottom.a == 0 && bottom.c == 0 && bottom.d == 0;
if options.variable == "G" && self.evp.modeFamily == "hydrostatic" ...
        && self.evp.formulation == "F" && ordinaryNeumann
    % The exact zero-eigenvalue F-constant sector has G=0. This includes
    % rigid hydrostatic modes and literal APV Neumann/Neumann limits.
    % No sampled values or magnitude cutoff decide scientific activity.
    active = ~isinf(self.h);
end
[~,notPositive] = chol(target(active,active));
recipe = IMProjection.createRecipe(variable=options.variable,columnLabels=self.columnLabels,normalization=self.normalizationName(self.normalization),targetGramMatrix=target,majorantGramMatrix=majorant,activeColumnMask=active, ...
    supportsLeakage=notPositive==0,supportsQuadratic=true,provenance=provenance,zDomain=self.zDomain, ...
    evaluateFunction=@(z) self.rawVariable(options.variable,z)./factors,weightFunction=@(z) IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight,z,context),spec=spec);
end
