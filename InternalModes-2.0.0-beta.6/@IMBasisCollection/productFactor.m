function factor = productFactor(collection,options)
% Bind one collection page, variable, derivative, and coefficient function.
%
% Independent reference collections must be supplied and matched explicitly
% by the caller through a custom factor evaluator. This factory evaluates
% the same scientific basis on every requested quadrature grid.
%
% - Topic: Assess sampled bases
% - Parameter collection: continuous collection snapshot
% - Parameter options.id: distinct factor identity
% - Parameter options.page: requested collection page
% - Parameter options.variable: available continuous variable
% - Parameter options.derivativeOrder: supported vertical derivative order
% - Parameter options.countRole: retained or fixed; endpoints require fixed
% - Parameter options.coefficient: function of z multiplying each column
% - Returns factor: immutable metadata and deferred continuous evaluator
arguments (Input)
    collection (1,1) IMBasisCollection
    options.id (1,1) string
    options.page (1,1) double {mustBeInteger,mustBePositive} = 1
    options.variable (1,1) string
    options.derivativeOrder (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    options.countRole (1,1) string {mustBeMember(options.countRole,["retained","fixed"])} = "retained"
    options.coefficient (1,1) function_handle = @(z) ones(size(z))
end
if options.page > numel(collection.basisIndex)
    error("IMBasisCollection:InvalidPage","The requested page is absent from the collection.");
end
metadata = collection.metadata(collection.basisIndex(options.page));
variableIndex = find(metadata.variables==options.variable,1);
if isempty(variableIndex) || ~ismember(options.derivativeOrder,metadata.derivativeOrders{variableIndex})
    error("IMBasisCollection:UnsupportedFactor","The collection does not support the requested variable/derivative.");
end
if metadata.columnKind=="endpoint" && options.countRole~="fixed"
    error("IMBasisCollection:EndpointPrefix","Boundary endpoint factors must preserve the complete fixed coordinate set.");
end
factor = struct(id=options.id,family=metadata.family,labels=metadata.columnLabels,ordinals=1:numel(metadata.columnLabels),frequencySigns=zeros(size(metadata.columnLabels)),countRole=options.countRole,columnKind=metadata.columnKind,evaluate=@evaluate,provenance=struct(collectionMetadata=metadata,page=options.page,kappa=collection.kappa,variable=options.variable,derivativeOrder=options.derivativeOrder,independentSolve=false));
    function values = evaluate(z,columns,referenceId) %#ok<INUSD>
        values = collection.evaluate(z,pages=options.page,columns=columns,variable=options.variable,derivativeOrder=options.derivativeOrder);
        coefficient = options.coefficient(z);
        if ~isnumeric(coefficient) || ~(isscalar(coefficient) || isequal(size(coefficient),size(z))) || any(~isfinite(coefficient),"all")
            error("IMBasisCollection:InvalidCoefficient","A coefficient function must return one finite value or one per grid point.");
        end
        values = values.*coefficient;
    end
end
