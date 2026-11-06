classdef IMProductAssessmentPlan
    % Reserve individual product identities before invoking expensive callbacks.
    %
    % The stored inventory is immutable. Execution caches are local to assess;
    % a later execution cannot reuse stale evaluated fields or reference data.
    %
    % ```matlab
    % plan = IMProductAssessmentPlan(factors,products,outputs,retainedCounts=1:8,productBudget=10000);
    % result = plan.assess(grids,chunkSize=128);
    % ```
    %
    % - Topic: Inspect reservations
    % - Topic: Execute assessments
    % - Declaration: classdef IMProductAssessmentPlan
    properties (SetAccess = private)
        % Deferred continuous factor metadata and evaluators.
        % - Topic: Inspect reservations
        factors
        % Ordered declared interaction and channel identities.
        % - Topic: Inspect reservations
        products
        % Deferred output projection metadata and preparation callbacks.
        % - Topic: Inspect reservations
        outputs
        % Selection rule, fixed sparse or bounded all-products validation.
        % - Topic: Inspect reservations
        selection
        % Examined retained-count values, unchanged by assessment decisions.
        % - Topic: Inspect reservations
        retainedCounts
        % Product-family rows selected in original inventory order.
        % - Topic: Inspect reservations
        rows
        % Two-by-product input column indices for each selected row.
        % - Topic: Inspect reservations
        pairs
        % Earliest retained count at which each selected product is examined.
        % - Topic: Inspect reservations
        firstCounts
        % Exact reserved product count, including future structural zeros.
        % - Topic: Inspect reservations
        reservedProducts
        % Explicit maximum authorized product count.
        % - Topic: Inspect reservations
        productBudget
        % Inventory construction/planning elapsed seconds.
        % - Topic: Inspect reservations
        planningSeconds
    end
    methods
        function self = IMProductAssessmentPlan(factors,products,outputs,options)
            % Validate declared products and reserve all work before callbacks.
            %
            % Factors declare id, family, labels, ordinals, frequencySigns,
            % countRole, evaluate(z,columns,referenceId), and provenance.
            % Outputs declare id, family, labels, nondecreasing ordinals,
            % countRole, prepare(grids), and provenance. The product table
            % requires interactionId, channel, factorA, factorB, and output;
            % an optional complex coefficient multiplies each entire product.
            %
            % - Topic: Inspect reservations
            % - Parameter factors: cell row of deferred continuous factors
            % - Parameter products: ordered valid interaction/channel table
            % - Parameter outputs: cell row of deferred physical output projections
            % - Parameter options.interactionIds: supplied interaction IDs; default all
            % - Parameter options.retainedCounts: strictly increasing retained ordinals to assess
            % - Parameter options.productBudget: maximum reservation, including exact zeros
            % - Parameter options.selection: fixedSparse stresses or bounded allProducts control
            arguments (Input)
                factors (1,:) cell
                products table
                outputs (1,:) cell
                options.interactionIds (1,:) string = string.empty(1,0)
                options.retainedCounts (1,:) double {mustBeInteger,mustBePositive,mustBeFinite}
                options.productBudget (1,1) double {mustBeInteger,mustBePositive,mustBeFinite}
                options.selection (1,1) string {mustBeMember(options.selection,["fixedSparse","allProducts"])} = "fixedSparse"
            end
            timer = tic;
            if isempty(factors) || isempty(outputs) || isempty(products)
                error("IMProductAssessmentPlan:EmptyInventory","Supply factors, outputs, and at least one valid product-family row.");
            end
            required = ["interactionId","channel","factorA","factorB","output"];
            if ~all(ismember(required,string(products.Properties.VariableNames)))
                error("IMProductAssessmentPlan:InvalidProducts","Product table requires interactionId, channel, factorA, factorB, and output.");
            end
            for kind = ["factor","output"]
                if kind == "factor", entries = factors; else, entries = outputs; end
                ids = strings(1,numel(entries));
                for i = 1:numel(entries)
                    entry = entries{i};
                    fields = ["id","family","labels","ordinals","countRole","provenance"];
                    if kind == "factor", fields = [fields,"frequencySigns","evaluate"]; else, fields = [fields,"prepare"]; end %#ok<AGROW>
                    if ~isstruct(entry) || ~isscalar(entry) || ~all(isfield(entry,fields))
                        error("IMProductAssessmentPlan:InvalidEntry","Each %s needs %s.",kind,join(fields,", "));
                    end
                    entry.id = string(entry.id); entry.family = string(entry.family);
                    entry.labels = reshape(string(entry.labels),1,[]);
                    entry.ordinals = reshape(entry.ordinals,1,[]);
                    entry.countRole = string(entry.countRole);
                    if ~isscalar(entry.id) || ismissing(entry.id) || strlength(entry.id)==0 || ~isscalar(entry.family) || ismissing(entry.family) || strlength(entry.family)==0 || isempty(entry.labels) || any(ismissing(entry.labels) | strlength(entry.labels)==0) || numel(unique(entry.labels)) ~= numel(entry.labels)
                        error("IMProductAssessmentPlan:InvalidIdentity","Entry identities and distinct column labels must be nonempty.");
                    end
                    if ~isnumeric(entry.ordinals) || ~isreal(entry.ordinals) || numel(entry.ordinals) ~= numel(entry.labels) || any(~isfinite(entry.ordinals) | entry.ordinals < 1 | fix(entry.ordinals) ~= entry.ordinals) || ~isscalar(entry.countRole) || ~ismember(entry.countRole,["retained","fixed"])
                        error("IMProductAssessmentPlan:InvalidOrdinals","Column ordinals must be positive integers; countRole must be retained or fixed.");
                    end
                    if ~isstruct(entry.provenance) || ~isscalar(entry.provenance) || isempty(fieldnames(entry.provenance))
                        error("IMProductAssessmentPlan:MissingProvenance","Every factor/output requires explicit provenance.");
                    end
                    if isfield(entry,"columnKind") && string(entry.columnKind)=="endpoint" && entry.countRole~="fixed"
                        error("IMProductAssessmentPlan:EndpointPrefix","Boundary endpoint coordinates require fixed countRole.");
                    end
                    if kind == "factor"
                        entry.frequencySigns = reshape(entry.frequencySigns,1,[]);
                        if ~isa(entry.evaluate,"function_handle") || numel(entry.frequencySigns) ~= numel(entry.labels) || ~isreal(entry.frequencySigns) || any(~ismember(entry.frequencySigns,[-1 0 1]))
                            error("IMProductAssessmentPlan:InvalidFactor","Factors require an evaluator and one frequency sign per column.");
                        end
                        factors{i} = entry;
                    else
                        if ~isa(entry.prepare,"function_handle") || any(diff(entry.ordinals)<0)
                            error("IMProductAssessmentPlan:InvalidOutput","Outputs require a preparation callback and nondecreasing ordinals.");
                        end
                        outputs{i} = entry;
                    end
                    ids(i) = entry.id;
                end
                if numel(unique(ids))~=numel(ids)
                    error("IMProductAssessmentPlan:DuplicateIdentity","%s IDs must be distinct.",kind);
                end
            end
            products.interactionId = string(products.interactionId);
            products.channel = string(products.channel);
            if size(products.interactionId,2)~=1 || size(products.channel,2)~=1 || any(ismissing(products.interactionId) | strlength(products.interactionId)==0 | ismissing(products.channel) | strlength(products.channel)==0)
                error("IMProductAssessmentPlan:InvalidProducts","Interaction and channel identities must be nonempty scalar strings per row.");
            end
            for field = ["factorA","factorB","output"]
                values = products.(field);
                bound = numel(factors); if field=="output", bound=numel(outputs); end
                if ~isnumeric(values) || ~isreal(values) || size(values,2)~=1 || any(~isfinite(values) | values<1 | values>bound | fix(values)~=values)
                    error("IMProductAssessmentPlan:InvalidProducts","%s must index declared entries.",field);
                end
            end
            if ~ismember("coefficient",string(products.Properties.VariableNames))
                products.coefficient = ones(height(products),1);
            elseif ~isnumeric(products.coefficient) || size(products.coefficient,2)~=1 || any(~isfinite(products.coefficient))
                error("IMProductAssessmentPlan:InvalidProducts","Product coefficients must be finite scalars per row.");
            end
            retainedCounts = options.retainedCounts; productBudget = options.productBudget; selection = options.selection;
            interactionIds = options.interactionIds;
            if isempty(interactionIds), interactionIds = unique(products.interactionId,"stable").'; end
            if isempty(retainedCounts) || any(diff(retainedCounts)<=0) || isempty(interactionIds) || numel(unique(interactionIds))~=numel(interactionIds) || any(~ismember(interactionIds,products.interactionId))
                error("IMProductAssessmentPlan:InvalidSelection","Counts must strictly increase; interaction IDs must be distinct supplied identities.");
            end
            rows = find(ismember(products.interactionId,interactionIds)).';
            pairs = cell(size(rows)); firstCounts = cell(size(rows)); reserved = 0;
            for j=1:numel(rows)
                row = products(rows(j),:);
                a = factors{row.factorA}; b = factors{row.factorB};
                for factor={a,b}
                    if factor{1}.countRole=="retained" && ~all(ismember(1:retainedCounts(end),factor{1}.ordinals))
                        error("IMProductAssessmentPlan:InsufficientInputBand","Retained factors must contain every explicitly requested input ordinal; no count is silently reduced.");
                    end
                end
                [pairs{j},firstCounts{j},number] = reservePairs(a,b,retainedCounts(end),selection,productBudget-reserved);
                reserved=reserved+number;
                output = outputs{row.output};
                if output.countRole=="retained" && ~all(ismember(1:retainedCounts(end),output.ordinals))
                    error("IMProductAssessmentPlan:InvalidOutputCounts","Every ordinal through the requested count must exist in each retained output; sparse count lists do not permit missing output modes.");
                end
            end
            if reserved==0
                error("IMProductAssessmentPlan:EmptySelection","The selected counts contain no input products.");
            end
            self.factors=factors; self.products=products; self.outputs=outputs; self.selection=selection; self.retainedCounts=retainedCounts;
            self.rows=rows; self.pairs=pairs; self.firstCounts=firstCounts;
            self.reservedProducts=reserved; self.productBudget=productBudget; self.planningSeconds=toc(timer);
        end
    end
end

function [pairs,firstCounts,number] = reservePairs(a,b,maximum,selection,budget)
% Enumerate small ordinal stresses first; expand sign/column multiplicities
% only after their reservation fits. Large unrequested basis tails cost no
% Cartesian product storage.
if selection=="allProducts"
    ia=1:numel(a.labels); ib=1:numel(b.labels);
    if a.countRole=="retained", ia=ia(a.ordinals<=maximum); end
    if b.countRole=="retained", ib=ib(b.ordinals<=maximum); end
    number=numel(ia)*numel(ib); requireBudget(number,budget);
    [i,j]=ndgrid(ia,ib); pairs=[i(:).';j(:).']; firstCounts=ones(1,number);
    if a.countRole=="retained", firstCounts=max(firstCounts,a.ordinals(pairs(1,:))); end
    if b.countRole=="retained", firstCounts=max(firstCounts,b.ordinals(pairs(2,:))); end
    return
end
ordinalPairs=zeros(0,2); first=zeros(0,1); multiplicities=zeros(0,1); number=0;
for count=1:maximum
    anchors=unique([1 2 ceil(count/2) max(1,count-2):count]); anchors=anchors(anchors<=count);
    left=anchors; right=anchors;
    if a.countRole=="fixed", left=0; end
    if b.countRole=="fixed", right=0; end
    [i,j]=ndgrid(left,right); candidates=[i(:),j(:)];
    additions=candidates(~ismember(candidates,ordinalPairs,"rows"),:);
    for q=1:size(additions,1)
        na=nnz(a.ordinals==additions(q,1)); nb=nnz(b.ordinals==additions(q,2));
        if additions(q,1)==0, na=numel(a.labels); end
        if additions(q,2)==0, nb=numel(b.labels); end
        amount=na*nb; number=number+amount; requireBudget(number,budget);
        ordinalPairs(end+1,:)=additions(q,:); first(end+1,1)=count; multiplicities(end+1,1)=amount; %#ok<AGROW>
    end
end
pairs=zeros(2,number); firstCounts=zeros(1,number); position=0;
for q=1:size(ordinalPairs,1)
    ia=find(a.ordinals==ordinalPairs(q,1)); ib=find(b.ordinals==ordinalPairs(q,2));
    if ordinalPairs(q,1)==0, ia=1:numel(a.labels); end
    if ordinalPairs(q,2)==0, ib=1:numel(b.labels); end
    [i,j]=ndgrid(ia,ib); batch=position+(1:multiplicities(q));
    pairs(:,batch)=[i(:).';j(:).']; firstCounts(batch)=first(q); position=position+multiplicities(q);
end
[~,order]=sortrows(pairs.',[2 1]); pairs=pairs(:,order); firstCounts=firstCounts(order);
end

function requireBudget(number,budget)
if number>budget
    error("IMProductAssessmentPlan:ProductBudgetExceeded","The metadata reservation exceeds the remaining product budget of %d (at least %d products, including zeros). No evaluator has run.",budget,number);
end
end
