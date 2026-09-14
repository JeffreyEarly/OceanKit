function [projection,recipe] = prepareProjection(self,z,weights,options)
% Construct a fixed projection for one requested page and variable.
%
% This operation preserves the supplied samples, weights, and explicit
% columns. It performs no quadrature fitting or acceptance policy. The
% continuous basis supplies the signed pairing and positive error metric.
%
% - Topic: Construct projections
% - Declaration: [projection,recipe] = prepareProjection(collection,z,weights,options)
% - Parameter z: physical sample column
% - Parameter weights: fixed quadrature weights aligned with z
% - Parameter options.page: requested-page index, not distinct-solve index
% - Parameter options.variable: scalar u or aligned F/G variable
% - Parameter options.columns: explicit column indices; empty means all
% - Returns projection: IMProjection with the requested columns
% - Returns recipe: basis-owned metric, capabilities, and provenance
arguments (Input)
    self (1,1) IMBasisCollection
    z (:,1) double {mustBeReal,mustBeFinite}
    weights (:,1) double {mustBeReal,mustBeFinite}
    options.page (1,1) double {mustBeInteger,mustBePositive} = 1
    options.variable (1,1) string = ""
    options.columns (1,:) double {mustBeInteger,mustBePositive} = []
end
if options.page > numel(self.basisIndex)
    error("IMBasisCollection:InvalidPage","page must select an existing requested page.");
end
basis = self.bases{self.basisIndex(options.page)};
if options.variable == ""
    recipe = basis.projectionRecipe();
else
    recipe = basis.projectionRecipe(variable=options.variable);
end
projection = IMProjection.fromRecipe(recipe,z,weights,columns=options.columns);
end
