classdef IMInternalModesDiscreteTransform
    % Store aligned discrete transforms for internal-mode F/G variables.
    %
    % A single instance stores one point-and-weight rule and the complete
    % aligned internal-mode family.  The sampled `F` and `G` columns share
    % physical mode labels, equivalent depths, and normalization, while
    % each directly representable variable owns its metric, Gram matrix,
    % active-column projector, and Galerkin forward matrix.
    %
    % Use variable-qualified accessors because the two channels generally
    % have different continuous inner products:
    %
    % Let $$n_z$$ be the sample count and $$n_m$$ the aligned family count.
    % `inverseMatrix(variable=V)` is $$n_z\times n_m$$, while
    % `forwardMatrix(variable=V)` is $$n_m\times n_z$$ and has zero rows
    % for inactive variable columns. For the active projector $$Q_V$$,
    %
    % $$A_\mathrm{f}^{V}A_\mathrm{i}^{V}=Q_V.$$
    %
    % `modeProjectionFunctional` returns
    % $$(A_\mathrm{i}^{V})^\mathsf{T}W_VX$$ before the active Gram solve;
    % `transformForward` returns coefficients after that solve. Both use
    % the signed Pontryagin pairing. `targetMajorantGramMatrix` returns the
    % induced positive Hilbert majorant used for error magnitudes; it does
    % not replace the signed projection metric.
    % When `variable` is omitted, accessors use `primaryVariable`, which is
    % the solved formulation when that channel was requested and otherwise
    % the first directly representable requested channel.
    %
    % ```matlab
    % [transform,assessment] = basisSet.discreteTransform(nPoints=24,variables=["F","G"]);
    % aG = transform.transformForward(G,variable="G");
    % F = transform.transformBack(aG,variable="F");
    % ```
    %
    % - Topic: Create discrete transforms
    % - Topic: Apply discrete transforms
    % - Topic: Inspect samples and modes
    % - Topic: Assess transform quality
    % - Declaration: classdef IMInternalModesDiscreteTransform

    properties (SetAccess = private)
        % Shared physical sample locations.
        z
        % Shared quadrature weights.
        weights
        % Physical labels for the aligned modal family.
        modeNumber
        % Equivalent depths aligned with `modeNumber`.
        h
        % Normalization captured when the transform was built.
        normalization
        % Directly representable forward channels in canonical order.
        availableVariables
        % Default F/G channel used when `variable` is omitted.
        %
        % `primaryVariable` is the solved formulation when it is among
        % `availableVariables`; otherwise it is the first available requested
        % channel in canonical order.
        %
        % - Topic: Inspect samples and modes
        primaryVariable
        % Physical vertical domain `[zBottom zSurface]`.
        zDomain
        % Full physical depth.
        depth
        % Gravitational acceleration.
        g
        % Internal-mode family name.
        modeFamily
        % Buoyancy frequency squared sampled at `z`.
        N2Values
        % Canonical endpoint identities, `["surface";"bottom"]`.
        endpointLocations
        % Immutable basis metadata snapshot.
        problemMetadata
        % Whether any shared quadrature weight is negative.
        hasNegativeWeights
    end

    properties (Access = private)
        inverseMatrices
        endpointMatrices
        channelData
    end

    methods
        function self = IMInternalModesDiscreteTransform(options)
            % Create an aligned internal-mode transform from prepared data.
            %
            % `channelData.F` and `channelData.G` are scalar structs with
            % fields `available`, `reason`, `activeModeMask`,
            % `metricMatrix`, `targetGramMatrix`, and
            % `targetMajorantGramMatrix`. Diagnostic inverse
            % matrices and endpoint traces are supplied for both variables,
            % including variables that have no direct sampled projection.
            %
            % - Topic: Developer topics — Construct transforms
            % - Declaration: transform = IMInternalModesDiscreteTransform(options)
            % - Developer: true
            arguments
                options.z (:,1) double {mustBeReal, mustBeFinite}
                options.weights (:,1) double {mustBeReal, mustBeFinite}
                options.modeNumber (1,:) double {mustBeInteger}
                options.h (1,:) double {mustBeReal}
                options.normalization {mustBeTextScalar}
                options.inverseF (:,:) double {mustBeReal, mustBeFinite}
                options.inverseG (:,:) double {mustBeReal, mustBeFinite}
                options.endpointF (2,:) double {mustBeReal, mustBeFinite}
                options.endpointG (2,:) double {mustBeReal, mustBeFinite}
                options.channelData (1,1) struct
                options.primaryVariable {mustBeTextScalar, mustBeMember(options.primaryVariable,["","F","G"])} = ""
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive}
                options.modeFamily {mustBeTextScalar}
                options.N2Values (:,1) double {mustBeReal, mustBeFinite}
                options.problemMetadata (1,1) struct = struct()
            end

            z = options.z(:);
            weights = options.weights(:);
            nSamples = length(z);
            nModes = length(options.modeNumber);
            if length(weights) ~= nSamples || size(options.inverseF,1) ~= nSamples || size(options.inverseG,1) ~= nSamples
                error("IMInternalModesDiscreteTransform:InvalidShape", "z, weights, and both inverse matrices must describe the same sample count.");
            end
            if size(options.inverseF,2) ~= nModes || size(options.inverseG,2) ~= nModes || length(options.h) ~= nModes
                error("IMInternalModesDiscreteTransform:InvalidShape", "F, G, h, and modeNumber must describe the same aligned modal family.");
            end
            if size(options.endpointF,2) ~= nModes || size(options.endpointG,2) ~= nModes
                error("IMInternalModesDiscreteTransform:InvalidShape", "Endpoint traces must have two rows and one column per mode.");
            end
            if length(options.N2Values) ~= nSamples || options.zDomain(1) >= options.zDomain(2)
                error("IMInternalModesDiscreteTransform:InvalidPhysicalSnapshot", "N2Values must align with z and zDomain must be increasing.");
            end

            inverseMatrices = struct(F=options.inverseF,G=options.inverseG);
            endpointMatrices = struct(F=options.endpointF,G=options.endpointG);
            channelData = options.channelData;
            variables = ["F","G"];
            availableVariables = strings(1,0);
            for variable = variables
                field = char(variable);
                if ~isfield(channelData,field)
                    error("IMInternalModesDiscreteTransform:InvalidChannelData", "channelData must contain F and G records.");
                end
                data = channelData.(field);
                requiredFields = ["available","reason","activeModeMask","metricMatrix","targetGramMatrix","targetMajorantGramMatrix"];
                if ~all(isfield(data,requiredFields))
                    error("IMInternalModesDiscreteTransform:InvalidChannelData", "Each channel record must contain availability, reason, active mask, metric, and target Gram data.");
                end
                data.available = logical(data.available);
                data.reason = string(data.reason);
                data.activeModeMask = reshape(logical(data.activeModeMask),1,[]);
                if length(data.activeModeMask) ~= nModes
                    error("IMInternalModesDiscreteTransform:InvalidChannelData", "Each activeModeMask must contain one value per aligned mode.");
                end
                if data.available
                    if ~isequal(size(data.metricMatrix),[nSamples nSamples]) || ~isequal(size(data.targetGramMatrix),[nModes nModes]) ...
                            || ~isequal(size(data.targetMajorantGramMatrix),[nModes nModes])
                        error("IMInternalModesDiscreteTransform:InvalidChannelData", "Available channel matrices have incompatible dimensions.");
                    end
                    metricTolerance = 100*eps(max(1,norm(data.metricMatrix,2)));
                    targetTolerance = 100*eps(max(1,norm(data.targetGramMatrix,2)));
                    majorantTolerance = 100*eps(max(1,norm(data.targetMajorantGramMatrix,2)));
                    if norm(data.metricMatrix-data.metricMatrix.',2) > metricTolerance ...
                            || norm(data.targetGramMatrix-data.targetGramMatrix.',2) > targetTolerance ...
                            || norm(data.targetMajorantGramMatrix-data.targetMajorantGramMatrix.',2) > majorantTolerance
                        error("IMInternalModesDiscreteTransform:NonSymmetricMatrix", "Metric, target Gram, and target majorant Gram matrices must be symmetric.");
                    end
                    data.metricMatrix = 0.5*(data.metricMatrix+data.metricMatrix.');
                    data.targetGramMatrix = 0.5*(data.targetGramMatrix+data.targetGramMatrix.');
                    data.targetMajorantGramMatrix = 0.5*(data.targetMajorantGramMatrix+data.targetMajorantGramMatrix.');
                    inverse = inverseMatrices.(field);
                    data.gramMatrix = 0.5*(inverse.'*data.metricMatrix*inverse + inverse.'*data.metricMatrix.'*inverse);
                    active = find(data.activeModeMask);
                    data.forwardMatrix = zeros(nModes,nSamples);
                    if isempty(active)
                        data.relativeGramOperatorError = 0;
                        data.roundTripError = 0;
                        data.inverseMatrixConditionNumber = NaN;
                        data.gramConditionNumber = NaN;
                        data.targetGramIsPositiveDefinite = true;
                        data.sampledGramRank = 0;
                    else
                        gramActive = data.gramMatrix(active,active);
                        singularValues = svd(gramActive);
                        rankTolerance = max(size(gramActive))*eps(max(1,norm(gramActive,2)));
                        data.sampledGramRank = sum(singularValues > rankTolerance);
                        if data.sampledGramRank < length(active)
                            data.forwardMatrix(active,:) = pinv(gramActive,rankTolerance)*(inverse(:,active).'*data.metricMatrix);
                        else
                            data.forwardMatrix(active,:) = gramActive \ (inverse(:,active).'*data.metricMatrix);
                        end
                        targetActive = data.targetGramMatrix(active,active);
                        targetNorms = diag(targetActive);
                        scale = 1./sqrt(abs(targetNorms));
                        scaledDifference = scale.*(gramActive-targetActive).*scale.';
                        projector = diag(double(data.activeModeMask));
                        data.relativeGramOperatorError = norm(scaledDifference,2);
                        if data.sampledGramRank < length(active)
                            data.relativeGramOperatorError = Inf;
                        end
                        data.roundTripError = norm(data.forwardMatrix*inverse-projector,2);
                        data.inverseMatrixConditionNumber = cond(inverse(:,active));
                        data.gramConditionNumber = cond(gramActive);
                        targetTolerance = 100*eps(max(1,norm(targetActive,2)));
                        data.targetGramIsPositiveDefinite = min(eig(targetActive)) > targetTolerance;
                    end
                    availableVariables(end+1) = variable; %#ok<AGROW>
                else
                    data.metricMatrix = zeros(0,0);
                    data.targetGramMatrix = zeros(0,0);
                    data.targetMajorantGramMatrix = zeros(0,0);
                    data.gramMatrix = zeros(0,0);
                    data.forwardMatrix = zeros(0,nSamples);
                    data.relativeGramOperatorError = NaN;
                    data.roundTripError = NaN;
                    data.inverseMatrixConditionNumber = NaN;
                    data.gramConditionNumber = NaN;
                    data.targetGramIsPositiveDefinite = false;
                    data.sampledGramRank = 0;
                end
                channelData.(field) = data;
            end

            primaryVariable = string(options.primaryVariable);
            if primaryVariable == ""
                if isempty(availableVariables)
                    error("IMInternalModesDiscreteTransform:NoAvailableVariable", "At least one directly representable F or G channel is required.");
                end
                primaryVariable = availableVariables(1);
            elseif ~ismember(primaryVariable,availableVariables)
                error("IMInternalModesDiscreteTransform:InvalidPrimaryVariable", "primaryVariable must name one of the directly representable requested channels: %s.",join(availableVariables,", "));
            end

            self.z = z;
            self.weights = weights;
            self.modeNumber = reshape(options.modeNumber,1,[]);
            self.h = reshape(options.h,1,[]);
            self.normalization = string(options.normalization);
            self.availableVariables = availableVariables;
            self.primaryVariable = primaryVariable;
            self.zDomain = options.zDomain;
            self.depth = diff(options.zDomain);
            self.g = options.g;
            self.modeFamily = string(options.modeFamily);
            self.N2Values = options.N2Values(:);
            self.endpointLocations = ["surface";"bottom"];
            self.problemMetadata = options.problemMetadata;
            self.hasNegativeWeights = any(weights < 0);
            self.inverseMatrices = inverseMatrices;
            self.endpointMatrices = endpointMatrices;
            self.channelData = channelData;
        end

        function tf = hasForwardTransform(self, options)
            % Return whether a direct sampled projection exists.
            arguments
                self IMInternalModesDiscreteTransform
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = self.primaryVariable
            end
            tf = self.channelData.(char(string(options.variable))).available;
        end

        function reason = forwardTransformReason(self, options)
            % Return the reason a direct projection is unavailable.
            arguments
                self IMInternalModesDiscreteTransform
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = self.primaryVariable
            end
            data = self.channelData.(char(string(options.variable)));
            if data.available
                reason = "";
            else
                reason = data.reason;
            end
        end

        function mask = activeModeMask(self, options)
            % Return the full-family active-column projector mask.
            arguments
                self IMInternalModesDiscreteTransform
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = self.primaryVariable
            end
            data = self.requireChannel(string(options.variable));
            mask = data.activeModeMask;
        end

        function matrix = inverseMatrix(self, options)
            % Return the sampled modal synthesis matrix for F or G.
            arguments
                self IMInternalModesDiscreteTransform
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = self.primaryVariable
            end
            matrix = self.inverseMatrices.(char(string(options.variable)));
        end

        function matrix = forwardMatrix(self, options)
            % Return the variable-qualified Galerkin forward matrix.
            arguments
                self IMInternalModesDiscreteTransform
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = self.primaryVariable
            end
            data = self.requireChannel(string(options.variable));
            matrix = data.forwardMatrix;
        end

        function matrix = metricMatrix(self, options)
            % Return the variable-qualified sampled signed metric.
            arguments
                self IMInternalModesDiscreteTransform
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = self.primaryVariable
            end
            data = self.requireChannel(string(options.variable));
            matrix = data.metricMatrix;
        end

        function matrix = gramMatrix(self, options)
            % Return the sampled signed full-family Gram matrix.
            arguments
                self IMInternalModesDiscreteTransform
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = self.primaryVariable
            end
            data = self.requireChannel(string(options.variable));
            matrix = data.gramMatrix;
        end

        function matrix = targetGramMatrix(self, options)
            % Return the continuous signed full-family target Gram matrix.
            arguments
                self IMInternalModesDiscreteTransform
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = self.primaryVariable
            end
            data = self.requireChannel(string(options.variable));
            matrix = data.targetGramMatrix;
        end

        function matrix = targetMajorantGramMatrix(self, options)
            % Return the continuous positive Hilbert-majorant Gram matrix.
            %
            % `targetGramMatrix` remains the signed Pontryagin pairing used
            % for projection. This accessor returns the positive matrix used
            % for error magnitudes and coupled quadratic certification.
            %
            % - Topic: Assess transform quality
            % - Declaration: matrix = targetMajorantGramMatrix(transform,options)
            % - Parameter options.variable: `"F"` or `"G"`
            % - Returns matrix: continuous positive majorant Gram matrix
            arguments
                self IMInternalModesDiscreteTransform
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = self.primaryVariable
            end
            data = self.requireChannel(string(options.variable));
            matrix = data.targetMajorantGramMatrix;
        end

        function values = endpointValues(self, options)
            % Return normalized endpoint traces, surface then bottom.
            arguments
                self IMInternalModesDiscreteTransform
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = self.primaryVariable
            end
            values = self.endpointMatrices.(char(string(options.variable)));
        end

        function pairings = modeProjectionFunctional(self, values, options)
            % Apply `(A_i^V)' W_V` without solving the modal Gram system.
            arguments
                self IMInternalModesDiscreteTransform
                values double
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = self.primaryVariable
            end
            variable = string(options.variable);
            data = self.requireChannel(variable);
            if size(values,1) ~= length(self.z)
                error("IMInternalModesDiscreteTransform:InvalidSampleCount", "values must have one row for each sample point in z.");
            end
            pairings = self.inverseMatrices.(char(variable)).'*data.metricMatrix*values;
            pairings(~data.activeModeMask,:) = 0;
        end

        function coefficients = transformForward(self, values, options)
            % Project sampled values into aligned family coefficients.
            arguments
                self IMInternalModesDiscreteTransform
                values double
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = self.primaryVariable
            end
            variable = string(options.variable);
            data = self.requireChannel(variable);
            if size(values,1) ~= length(self.z)
                error("IMInternalModesDiscreteTransform:InvalidSampleCount", "values must have one row for each sample point in z.");
            end
            coefficients = data.forwardMatrix*values;
        end

        function values = transformBack(self, coefficients, options)
            % Synthesize F or G values from full-family coefficients.
            arguments
                self IMInternalModesDiscreteTransform
                coefficients double
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = self.primaryVariable
            end
            if size(coefficients,1) ~= length(self.modeNumber)
                error("IMInternalModesDiscreteTransform:InvalidCoefficientCount", "coefficients must have one row for each aligned family mode.");
            end
            values = self.inverseMatrices.(char(string(options.variable)))*coefficients;
        end

        function diagnostics = channelDiagnostics(self, options)
            % Return scalar quality diagnostics for one built channel.
            % - Developer: true
            arguments
                self IMInternalModesDiscreteTransform
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = self.primaryVariable
            end
            data = self.requireChannel(string(options.variable));
            diagnostics = rmfield(data,["metricMatrix","targetGramMatrix","targetMajorantGramMatrix","gramMatrix","forwardMatrix"]);
        end

        function value = relativeGramOperatorError(self, options)
            % Return the normalized Gram operator error for one channel.
            arguments
                self IMInternalModesDiscreteTransform
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = self.primaryVariable
            end
            data = self.requireChannel(string(options.variable));
            value = data.relativeGramOperatorError;
        end

        function value = roundTripError(self, options)
            % Return the active-projector round-trip error for one channel.
            arguments
                self IMInternalModesDiscreteTransform
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = self.primaryVariable
            end
            data = self.requireChannel(string(options.variable));
            value = data.roundTripError;
        end

        function value = inverseMatrixConditionNumber(self, options)
            % Return the active sampled-basis condition number.
            arguments
                self IMInternalModesDiscreteTransform
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = self.primaryVariable
            end
            data = self.requireChannel(string(options.variable));
            value = data.inverseMatrixConditionNumber;
        end

        function value = gramConditionNumber(self, options)
            % Return the active sampled-Gram condition number.
            arguments
                self IMInternalModesDiscreteTransform
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = self.primaryVariable
            end
            data = self.requireChannel(string(options.variable));
            value = data.gramConditionNumber;
        end

        function tf = targetGramIsPositiveDefinite(self, options)
            % Return whether the active target channel defines a norm.
            arguments
                self IMInternalModesDiscreteTransform
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = self.primaryVariable
            end
            data = self.requireChannel(string(options.variable));
            tf = data.targetGramIsPositiveDefinite;
        end
    end

    methods (Hidden)
        function transform = prefixTransform(self, nModes)
            % Derive a leading family prefix without reevaluating the basis.
            arguments
                self IMInternalModesDiscreteTransform
                nModes (1,1) double {mustBeInteger, mustBePositive}
            end
            if nModes > length(self.modeNumber)
                error("IMInternalModesDiscreteTransform:InvalidModeCount", "The transform contains %d modes, but a %d-mode prefix was requested.",length(self.modeNumber),nModes);
            end

            preparedChannels = self.channelData;
            for variable = ["F","G"]
                field = char(variable);
                data = preparedChannels.(field);
                data.activeModeMask = data.activeModeMask(1:nModes);
                if data.available
                    data.targetGramMatrix = data.targetGramMatrix(1:nModes,1:nModes);
                    data.targetMajorantGramMatrix = data.targetMajorantGramMatrix(1:nModes,1:nModes);
                end
                preparedChannels.(field) = data;
            end
            transform = IMInternalModesDiscreteTransform(z=self.z,weights=self.weights,modeNumber=self.modeNumber(1:nModes),h=self.h(1:nModes), ...
                normalization=self.normalization,inverseF=self.inverseMatrices.F(:,1:nModes),inverseG=self.inverseMatrices.G(:,1:nModes), ...
                endpointF=self.endpointMatrices.F(:,1:nModes),endpointG=self.endpointMatrices.G(:,1:nModes),channelData=preparedChannels, ...
                primaryVariable=self.primaryVariable,zDomain=self.zDomain,g=self.g,modeFamily=self.modeFamily,N2Values=self.N2Values,problemMetadata=self.problemMetadata);
        end

        function transform = transformWithWeights(self, weights)
            % Replace shared weights while reusing prepared modal data.
            arguments
                self IMInternalModesDiscreteTransform
                weights (:,1) double {mustBeReal, mustBeFinite}
            end
            weights = weights(:);
            if length(weights) ~= length(self.z) || ~any(weights ~= 0)
                error("IMInternalModesDiscreteTransform:InvalidWeights", "weights must contain one value per sample and at least one nonzero value.");
            end

            preparedChannels = self.channelData;
            for variable = self.availableVariables
                field = char(variable);
                data = preparedChannels.(field);
                if ~isfield(data,"interiorWeight") || ~isfield(data,"endpointMetricMatrix")
                    error("IMInternalModesDiscreteTransform:MissingPreparedMetric", "The %s channel does not contain reusable sampled-metric components.",variable);
                end
                data.metricMatrix = diag(data.interiorWeight.*weights)+data.endpointMetricMatrix;
                preparedChannels.(field) = data;
            end
            transform = IMInternalModesDiscreteTransform(z=self.z,weights=weights,modeNumber=self.modeNumber,h=self.h,normalization=self.normalization, ...
                inverseF=self.inverseMatrices.F,inverseG=self.inverseMatrices.G,endpointF=self.endpointMatrices.F,endpointG=self.endpointMatrices.G, ...
                channelData=preparedChannels,primaryVariable=self.primaryVariable,zDomain=self.zDomain,g=self.g,modeFamily=self.modeFamily,N2Values=self.N2Values,problemMetadata=self.problemMetadata);
        end

        function context = quadratureFitContext(self, options)
            % Return the prepared normalized-Gram system for one channel.
            arguments
                self IMInternalModesDiscreteTransform
                options.variable {mustBeTextScalar, mustBeMember(options.variable,["F","G"])} = self.primaryVariable
            end
            variable = string(options.variable);
            data = self.requireChannel(variable);
            if ~isfield(data,"interiorWeight") || ~isfield(data,"endpointMetricMatrix")
                error("IMInternalModesDiscreteTransform:MissingPreparedMetric", "The %s channel does not contain reusable sampled-metric components.",variable);
            end
            inverseMatrix = self.inverseMatrices.(char(variable));
            endpointGramMatrix = inverseMatrix.'*data.endpointMetricMatrix*inverseMatrix;
            endpointGramMatrix = 0.5*(endpointGramMatrix+endpointGramMatrix.');
            [A,b,pairs] = normalizedGramSystem(inverseMatrix,data.interiorWeight,data.targetGramMatrix,endpointGramMatrix,data.activeModeMask);
            context = struct(variable=variable,z=self.z,modeNumber=self.modeNumber,normalization=self.normalization,inverseMatrix=inverseMatrix, ...
                interiorWeight=data.interiorWeight,targetGramMatrix=data.targetGramMatrix,endpointGramMatrix=endpointGramMatrix, ...
                geometricWeights=self.weights,activeModeMask=data.activeModeMask,normalizedGramA=A,normalizedGramB=b,normalizedGramModePairs=pairs);
        end
    end

    methods (Access = private)
        function data = requireChannel(self, variable)
            data = self.channelData.(char(variable));
            if ~data.available
                error("IMInternalModesDiscreteTransform:UnavailableForwardTransform", "%s", data.reason);
            end
        end
    end
end

function [A,b,pairs] = normalizedGramSystem(inverseMatrix,interiorWeight,targetGramMatrix,endpointGramMatrix,activeMask)
active = find(activeMask);
nRows = length(active)*(length(active)+1)/2;
A = zeros(nRows,size(inverseMatrix,1));
b = zeros(nRows,1);
pairs = zeros(nRows,2);
iRow = 0;
targetNorms = diag(targetGramMatrix);
for iActive = 1:length(active)
    iMode = active(iActive);
    for jActive = iActive:length(active)
        jMode = active(jActive);
        iRow = iRow+1;
        rowFactor = 1;
        if iMode ~= jMode
            rowFactor = sqrt(2);
        end
        scale = sqrt(abs(targetNorms(iMode)*targetNorms(jMode)));
        A(iRow,:) = rowFactor*(interiorWeight.*inverseMatrix(:,iMode).*inverseMatrix(:,jMode)).'/scale;
        b(iRow) = rowFactor*(targetGramMatrix(iMode,jMode)-endpointGramMatrix(iMode,jMode))/scale;
        pairs(iRow,:) = [iMode jMode];
    end
end
end
