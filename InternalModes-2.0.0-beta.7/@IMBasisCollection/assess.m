function assessment = assess(self,z,weights,options)
% Measure a fixed requested page through the common result contract.
%
% Products and leakage are explicit prepared recipes as documented by
% IMProjection.assess. This operation neither changes scientific modes nor
% creates an independent reference solve. Apply tolerances afterward with
% checkBasisAssessment(assessment,...).
%
% - Topic: Assess collections
% - Declaration: assessment = assess(collection,z,weights,options)
arguments (Input)
    self (1,1) IMBasisCollection
    z (:,1) double {mustBeReal,mustBeFinite}
    weights (:,1) double {mustBeReal,mustBeFinite}
    options.page (1,1) double {mustBeInteger,mustBePositive} = 1
    options.variable (1,1) string = ""
    options.columns (1,:) double {mustBeInteger,mustBePositive} = []
    options.prefixColumnCounts (1,:) double {mustBeInteger,mustBePositive} = []
    options.leakage struct = struct.empty
    options.products struct = struct.empty
    options.coverage (1,1) struct = struct()
end
arguments (Output)
    assessment (1,1) struct
end
started = tic;
[p,recipe] = self.prepareProjection(z,weights,page=options.page,variable=options.variable,columns=options.columns);
constructionSeconds = toc(started);
if ~isempty(options.leakage) && ~recipe.supportsLeakage
    error("IMBasisCollection:UnsupportedAssessment","This recipe does not support leakage with its continuous target metric.");
end
if ~isempty(options.products) && ~recipe.supportsQuadratic
    error("IMBasisCollection:UnsupportedAssessment","This recipe does not support quadratic products; consult the family's projection capabilities.");
end
metadata = self.metadata(self.basisIndex(options.page));
identity = struct(family=metadata.family,variable=recipe.variable,page=options.page,basisIndex=self.basisIndex(options.page),sourcePage=self.sourcePage(options.page),columnLabels=p.columnLabels,normalization=recipe.normalization,projectionProvenance=recipe.provenance);
if ~isempty(self.kappa)
    identity.kappa = self.kappa(options.page);
end
assessment = p.assess(identity=identity,columnKind=metadata.columnKind,prefixColumnCounts=options.prefixColumnCounts,leakage=options.leakage,products=options.products,coverage=options.coverage,constructionSeconds=constructionSeconds);
end
