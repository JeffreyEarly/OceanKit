function projection = fromRecipe(self,z,weights,options)
    % Construct a projection from a hidden family descriptor on an unchanged grid.
    % - Topic: Evaluate projection recipes
    % - Parameter z: physical sample coordinates
    % - Parameter weights: fixed quadrature weights
    % - Parameter options.columns: explicit column selection, all by default
    % - Parameter options.prepareRecovery: also prepare numerical recovery, true by default
    % - Parameter options.failureMode: diagnostic or error on recovery failure
    % - Parameter options.rankTolerance: relative cutoff for numerical recovery
    % - Returns projection: physical pairing and optional numerical recovery
    arguments (Input)
        self (1,1) struct
        z (:,1) double {mustBeReal,mustBeFinite}
        weights (:,1) double {mustBeReal,mustBeFinite}
        options.columns (1,:) double {mustBeInteger,mustBePositive} = 1:numel(self.columnLabels)
        options.prepareRecovery (1,1) logical = true
        options.failureMode (1,1) string {mustBeMember(options.failureMode,["diagnostic","error"])} = "diagnostic"
        options.rankTolerance (1,1) double {mustBeReal,mustBePositive,mustBeFinite} = 1e-12
    end
    arguments (Output)
        projection (1,1) IMProjection
    end
    columns = selectedColumns(self,options.columns);
    sampled = evaluate(self,z,columns=columns);
    matrix = metric(self,z,weights);
    target = self.targetGramMatrix(columns,columns);
    majorant = self.majorantGramMatrix(columns,columns);
    % A vanishing sampled column does not change scientific activity. In
    % particular, an active signed-null physical norm must fail explicitly.
    projection = IMProjection(sampled,matrix,target,majorantGramMatrix=majorant,activeColumnMask=self.activeColumnMask(columns),columnLabels=self.columnLabels(columns),provenance=self.provenance, ...
        prepareRecovery=options.prepareRecovery,failureMode=options.failureMode,rankTolerance=options.rankTolerance);
end

function values = evaluate(self,z,options)
    % Evaluate normalized continuous columns at physical points.
    % - Topic: Evaluate projection recipes
    % - Parameter z: column vector of physical coordinates
    % - Parameter options.columns: selected column indices, all by default
    % - Returns values: nZ by nSelectedColumns values
    arguments (Input)
        self (1,1) struct
        z (:,1) double {mustBeReal,mustBeFinite}
        options.columns (1,:) double {mustBeInteger,mustBePositive} = 1:numel(self.columnLabels)
    end
    arguments (Output)
        values (:,:) double
    end
    assertAvailable(self);
    validatePoints(self,z);
    columns = selectedColumns(self,options.columns);
    values = self.evaluateFunction(z);
    values = values(:,columns);
end
function matrix = metric(self,z,weights,options)
    % Construct the signed or majorant pairing on fixed samples.
    %
    % Endpoint traces must be represented by supplied variable
    % samples. Derivative and companion-variable traces are explicit
    % unsupported capabilities, never inferred from the basis span.
    % Quadrature weights are unchanged; negative quadrature weights
    % can make the sampled majorant indefinite.
    % - Topic: Evaluate projection recipes
    % - Parameter z: strictly increasing physical sample coordinates
    % - Parameter weights: fixed quadrature weights
    % - Parameter options.majorant: use absolute endpoint coefficients
    % - Returns matrix: nZ by nZ sample metric
    arguments (Input)
        self (1,1) struct
        z (:,1) double {mustBeReal,mustBeFinite}
        weights (:,1) double {mustBeReal,mustBeFinite}
        options.majorant (1,1) logical = false
    end
    arguments (Output)
        matrix (:,:) double
    end
    assertAvailable(self);
    validatePoints(self,z);
    if numel(z) < 2 || any(diff(z) <= 0)
        error("IMProjection:InvalidGrid","z must contain at least two strictly increasing points.");
    end
    if numel(weights) ~= numel(z) || ~any(weights ~= 0)
        error("IMProjection:InvalidWeights","Supply one weight per sample and at least one nonzero weight.");
    end
    weight = self.weightFunction(z);
    if isscalar(weight)
        weight = repmat(weight,size(z));
    end
    if numel(weight) ~= numel(z) || ~isreal(weight) || any(~isfinite(weight(:)))
        error("IMProjection:InvalidInteriorWeight","The interior weight must return one finite real value per sample.");
    end
    if options.majorant
        weight = abs(weight);
    end
    matrix = diag(weights.*weight(:));
    endpoints = [self.spec.surfaceWeights;self.spec.bottomWeights];
    for i = 1:numel(endpoints)
        term = endpoints(i);
        if term.coefficient == 0
            continue
        end
        if term.d ~= 0
            error("IMProjection:UnsupportedDerivativeTrace","The %s endpoint requires a derivative trace; value-only samples cannot represent this projection.",term.location);
        end
        if term.c ~= 0
            matrix = addEndpoint(self,matrix,z,term.location,term.coefficient*term.c^2,options.majorant);
        end
    end
    if isfield(self.spec,"endpointInnerProductTerms")
        for i = 1:numel(self.spec.endpointInnerProductTerms)
            term = self.spec.endpointInnerProductTerms(i);
            if term.coefficient == 0
                continue
            end
            if string(term.variable) ~= self.variable
                error("IMProjection:UnsupportedCompanionTrace","The %s endpoint requires %s values in the %s pairing; paired-state projection is unsupported.",term.location,term.variable,self.variable);
            end
            matrix = addEndpoint(self,matrix,z,term.location,term.coefficient,options.majorant);
        end
    end
end
function assertAvailable(self)
    if ~self.isAvailable
        error("IMProjection:UnavailableProjection","%s",self.reason);
    end
end
function columns = selectedColumns(self,columns)
    if any(columns > numel(self.columnLabels)) || numel(unique(columns)) ~= numel(columns)
        error("IMProjection:InvalidColumns","Select distinct existing column indices; scientific labels are not array indices.");
    end
end
function validatePoints(self,z)
    tolerance = 100*eps(max(1,max(abs(self.zDomain))));
    if any(z < self.zDomain(1)-tolerance | z > self.zDomain(2)+tolerance)
        error("IMProjection:InvalidGrid","All points must lie inside the basis domain.");
    end
end
function matrix = addEndpoint(self,matrix,z,location,coefficient,majorant)
    if ~isfinite(coefficient)
        error("IMProjection:InvalidEndpointCoefficient","The %s endpoint coefficient must be finite.",location);
    end
    if string(location) == "surface"
        endpoint = self.zDomain(2);
    else
        endpoint = self.zDomain(1);
    end
    index = find(abs(z-endpoint) <= 100*eps(max(1,max(abs(self.zDomain)))),1);
    if isempty(index)
        error("IMProjection:MissingEndpointSample","Include the %s endpoint at z=%g to represent its value-only pairing.",location,endpoint);
    end
    if majorant
        coefficient = abs(coefficient);
    end
    matrix(index,index) = matrix(index,index)+coefficient;
end
