classdef IMNeutralLabelModesBasis
    % Evaluate independent zero-energy labels on exact neutral intervals.
    %
    % The caller specifies maximal intervals where N2 is exactly zero; no
    % inference from sampled density or stratification is performed. With
    % D = diff(zDomain), the default depth-normalized profiles satisfy
    % $$G_{zz} = -m_{\mathrm n}^2 G,\quad D^{-1}\int G_i G_j\,dz = \delta_{ij}.$$
    % Each profile is zero outside its interval and at both interval ends,
    % including physical surface and bottom boundaries. Endpoint labels
    % are held fixed; independent boundary anomalies belong to other mode
    % families. All m_n are positive, so no constant label is included.
    % Generalized-energy coefficients g0 and gd are provenance only and do
    % not change these boundary conditions. These modes have no equivalent
    % depth or physical-energy norm.
    % Amplitudes have units of length and require particle-label information
    % in addition to velocity and density. At finite amplitude, deleting a
    % component in this fixed linear basis need not preserve physical density.
    %
    % ```matlab
    % basis = IMNeutralLabelModesBasis([-200 -100],zDomain=[-1000 0],nModesPerInterval=4);
    % G = basis.G(linspace(-1000,0,1001).');
    % projection = basis.projectionOnGrid(z,weights,variable="APV");
    % amplitudes = projection.projectPhysical(APV);
    % ```
    %
    % - Topic: Create neutral labels
    % - Topic: Inspect neutral labels
    % - Topic: Evaluate neutral labels
    % - Topic: Project APV
    % - Declaration: classdef IMNeutralLabelModesBasis

    properties (SetAccess = private)
        % Increasing physical bottom and surface coordinates, in meters.
        % - Topic: Inspect neutral labels
        zDomain (1,2) double
        % Exact maximal neutral intervals, ordered from bottom to surface.
        % - Topic: Inspect neutral labels
        neutralIntervals (:,2) double
        % Nonnegative retained mode count for each neutral interval.
        % - Topic: Inspect neutral labels
        nModesPerInterval (1,:) double
        % Interval index for each retained column.
        % - Topic: Inspect neutral labels
        intervalIndex (1,:) double
        % Positive local mode number within each interval.
        % - Topic: Inspect neutral labels
        localModeNumber (1,:) double
        % Distinct interval/local-mode identities in retained-column order.
        % - Topic: Inspect neutral labels
        columnLabels (1,:) string
        % Positive vertical wavenumbers, in inverse meters.
        % - Topic: Inspect neutral labels
        m_n (1,:) double
        % Fixed lower/upper Dirichlet boundary types for each interval.
        % - Topic: Inspect neutral labels
        boundaryTypes (:,2) string
        % Planetary rotation rate in radians per second.
        % - Topic: Inspect neutral labels
        rotationRate (1,1) double
        % Central latitude in degrees north.
        % - Topic: Inspect neutral labels
        latitude (1,1) double
        % Surface generalized-energy provenance in meters per second squared; does not alter labels.
        % - Topic: Inspect neutral labels
        g0 (1,1) double
        % Bottom generalized-energy provenance in meters per second squared; does not alter labels.
        % - Topic: Inspect neutral labels
        gd (1,1) double
    end

    properties
        % Active profile normalization; per-evaluation overrides do not solve again.
        % - Topic: Inspect neutral labels
        normalization = "rmsG"
    end

    properties (Access = private)
        normalizationRules = configureDictionary("string","cell")
    end

    properties (Hidden, SetAccess = protected)
        % Immutable construction identities; derived normalization stays separate.
        basisContractState struct = struct()
    end

    properties (Dependent, SetAccess = private)
        % Supported continuous operations and conditional preparation requirements.
        % - Topic: Inspect basis sets
        % Requested/resolved scientific identity, without evaluating normalization.
        % - Topic: Inspect basis sets
        scientificSpecification
        % Original backend request and current selection, with unverified accuracy.
        % - Topic: Inspect basis sets
        constructionRecord
        % Fixed, per-column or independent horizontal identity.
        % - Topic: Inspect basis sets
        horizontalKind
        k
        kSquared
        % Signs of the active normalization form.
        % - Topic: Normalize modes
        normSignatures
        % Exact number of retained label columns.
        % - Topic: Inspect neutral labels
        nColumns
        % Scientific family name.
        % - Topic: Inspect neutral labels
        family
        % Retained interval and local-mode identities.
        % - Topic: Inspect neutral labels
        modeInfo
        % Zero physical frequency of particle-label coordinates.
        % - Topic: Inspect neutral labels
        frequency
        % Supported profile and derivative channels.
        % - Topic: Inspect neutral labels
        capabilities
        % Derived Coriolis frequency, f0 = 2 rotationRate sin(latitude).
        % - Topic: Inspect neutral labels
        f0
    end

    methods
        function value = get.scientificSpecification(self)
            value = IMBasisContracts.scientificSpecification(self);
        end
        function value = get.constructionRecord(self)
            value = IMBasisContracts.constructionRecord(self);
        end
        function value = get.horizontalKind(self)
            description = IMBasisContracts.horizontalIdentity(self); value = description.kind;
        end
        function value = get.k(self)
            description = IMBasisContracts.horizontalIdentity(self); value = description.k;
        end
        function value = get.kSquared(self)
            description = IMBasisContracts.horizontalIdentity(self); value = description.kSquared;
        end
        function value = get.normSignatures(self)
            value = IMBasisContracts.normalizationSignatures(self);
        end

        function self = IMNeutralLabelModesBasis(neutralIntervals,options)
            % Create analytical neutral labels on explicitly specified intervals.
            %
            % Use IMNeutralLabelModesBasis(background,neutralIntervals=...,nModesPerInterval=...)
            % to take domain and rotation from an authoritative background and
            % reject interval declarations contradicting its N2 profile.
            %
            % Intervals must be nonempty, strictly separated, increasing and
            % contained in zDomain. Merge touching intervals before calling.
            % Counts retain positive modes 1:nModesPerInterval per interval, in
            % interval order. Zero counts return no columns; there is no
            % additional constant sector, including in a fully neutral column.
            %
            % - Topic: Create neutral labels
            % - Parameter neutralIntervals: exact neutral intervals, one [lower upper] row each
            % - Parameter options.zDomain: physical [bottom surface] domain in meters
            % - Parameter options.nModesPerInterval: nonnegative positive-mode count, scalar or one per interval
            % - Parameter options.rotationRate: nonnegative planetary rotation rate
            % - Parameter options.latitude: central latitude in degrees north
            % - Parameter options.g0: optional surface generalized-energy provenance, zero by default
            % - Parameter options.gd: optional bottom generalized-energy provenance, zero by default
            % - Returns self: continuous neutral-label value basis
            arguments (Input)
                neutralIntervals
                options.neutralIntervals (:,2) double {mustBeReal,mustBeFinite} = zeros(0,2)
                options.zDomain double {mustBeReal,mustBeFinite} = []
                options.nModesPerInterval (1,:) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = []
                options.rotationRate double {mustBeReal,mustBeFinite,mustBeNonnegative} = []
                options.latitude double {mustBeReal,mustBeFinite} = []
                options.normalization = "rmsG"
                options.g0 (1,1) double {mustBeReal} = 0
                options.gd (1,1) double {mustBeReal} = 0
            end
            arguments (Output)
                self (1,1) IMNeutralLabelModesBasis
            end
            if isempty(options.nModesPerInterval)
                error("IMNeutralLabelModesBasis:MissingModeCount","Supply nModesPerInterval as an explicit scalar or one nonnegative count per declared interval; zero is valid.");
            end
            if isa(neutralIntervals,"IMBackground")
                background = neutralIntervals;
                if ~isempty(options.zDomain) || ~isempty(options.rotationRate) || ~isempty(options.latitude)
                    error("IMNeutralLabelModesBasis:ConflictingBackground","Physical inputs come from the background. Do not repeat zDomain, rotationRate or latitude.");
                end
                neutralIntervals = options.neutralIntervals;
                IMDeclaredIntervals.validate(background,neutralIntervals,0,"neutral","IMNeutralLabelModesBasis");
                options.zDomain = background.zDomain;
                options.rotationRate = background.rotationRate;
                options.latitude = background.latitude;
            elseif ~isempty(options.neutralIntervals)
                error("IMNeutralLabelModesBasis:ConflictingIntervals","Supply intervals positionally, or supply a background with neutralIntervals=.");
            end
            validateattributes(neutralIntervals,{'double'},{'real','finite','2d','ncols',2});
            validateattributes(options.zDomain,{'double'},{'real','finite','size',[1 2]});
            if isempty(options.rotationRate), options.rotationRate = 7.2921e-5; end
            if isempty(options.latitude), options.latitude = 33; end
            validateattributes(options.rotationRate,{'double'},{'real','finite','scalar','nonnegative'});
            validateattributes(options.latitude,{'double'},{'real','finite','scalar','>=',-90,'<=',90});
            if options.zDomain(1) >= options.zDomain(2)
                error("IMNeutralLabelModesBasis:InvalidDomain","zDomain must be strictly increasing.");
            end
            if isempty(neutralIntervals) || any(neutralIntervals(:,1) >= neutralIntervals(:,2)) || any(neutralIntervals(:,1) < options.zDomain(1)) || any(neutralIntervals(:,2) > options.zDomain(2))
                error("IMNeutralLabelModesBasis:InvalidIntervals","Supply nonempty increasing neutral intervals contained in zDomain.");
            end
            if any(neutralIntervals(2:end,1) <= neutralIntervals(1:end-1,2))
                error("IMNeutralLabelModesBasis:InvalidIntervals","Neutral intervals must be ordered, disjoint and strictly separated. Merge touching intervals.");
            end
            if isnan(options.g0) || isnan(options.gd)
                error("IMNeutralLabelModesBasis:InvalidEndpointAcceleration","g0 and gd must be real numbers or infinities, not NaN.");
            end
            counts = options.nModesPerInterval;
            if isscalar(counts)
                counts = repmat(counts,1,size(neutralIntervals,1));
            elseif numel(counts) ~= size(neutralIntervals,1)
                error("IMNeutralLabelModesBasis:InvalidModeCounts","nModesPerInterval must be scalar or contain one nonnegative count per neutral interval.");
            end
            if isempty(counts)
                error("IMNeutralLabelModesBasis:InvalidModeCounts","Supply at least one nonnegative positive-mode count.");
            end
            self.zDomain = options.zDomain;
            self.neutralIntervals = neutralIntervals;
            self.nModesPerInterval = counts;
            self.rotationRate = options.rotationRate;
            self.latitude = options.latitude;
            self.g0 = options.g0;
            self.gd = options.gd;
            self.normalization = options.normalization;
            self.boundaryTypes = repmat("D",size(neutralIntervals));
            self.intervalIndex = repelem(1:numel(counts),counts);
            self.localModeNumber = zeros(1,sum(counts));
            self.m_n = zeros(1,sum(counts));
            for iInterval = 1:numel(counts)
                columns = self.intervalIndex == iInterval;
                local = 1:counts(iInterval);
                self.localModeNumber(columns) = local;
                self.m_n(columns) = pi*local/diff(neutralIntervals(iInterval,:));
            end
            self.columnLabels = "interval" + string(self.intervalIndex) + ":mode" + string(self.localModeNumber);
            self.basisContractState = IMBasisContracts.initialState(self);
            self.basisContractState.requestedScientific = struct(neutralIntervals=neutralIntervals,nModesPerInterval=options.nModesPerInterval,energyBoundary=struct(g0=options.g0,gd=options.gd));
            self.basisContractState.request = struct(nModesPerInterval=options.nModesPerInterval);
            if exist("background","var"), self.basisContractState.physicalBackground = IMBasisContracts.physicalBackground(background); end
        end

        function value = get.f0(self)
            value = 2*self.rotationRate*sind(self.latitude);
        end

        function value = get.nColumns(self)
            value = numel(self.m_n);
        end

        function value = get.family(~)
            value = "neutralLabel";
        end

        function value = get.frequency(self)
            value = zeros(1,self.nColumns);
        end

        function value = get.capabilities(self)
            value = IMBasisContracts.capabilities(self);
        end

        function info = get.modeInfo(self)
            branch = repmat("neutralLabel",self.nColumns,1);
            info = table(self.columnLabels(:),branch,self.localModeNumber(:),self.intervalIndex(:),self.localModeNumber(:),VariableNames=["columnLabel","branch","branchIndex","intervalIndex","localModeNumber"]);
        end

        function self = selectColumns(self,columns)
            % Select aligned labels without reconstructing a basis.
            % - Topic: Inspect neutral labels
            arguments
                self IMNeutralLabelModesBasis
                columns (1,:) double {mustBeInteger,mustBePositive}
            end
            if any(columns > self.nColumns) || numel(unique(columns)) ~= numel(columns)
                error("IMNeutralLabelModesBasis:InvalidColumns","Select distinct existing label columns.");
            end
            source = self;
            self.basisContractState = IMBasisContracts.selectState(self.basisContractState,columns);
            for name = reshape(string(keys(self.normalizationRules)),1,[])
                self.normalizationRules{name} = @(~,j) IMInternalModesNormalizationTools.selectedNormFactor(source,name,columns(j));
            end
            self.intervalIndex = self.intervalIndex(columns);
            self.localModeNumber = self.localModeNumber(columns);
            self.columnLabels = self.columnLabels(columns);
            self.m_n = self.m_n(columns);
            self.nModesPerInterval = arrayfun(@(j) nnz(self.intervalIndex == j),1:size(self.neutralIntervals,1));
        end

        function self = addNormalization(self,name,rule)
            % Install a custom positive factor rule, factor = rule(basis,j).
            % - Topic: Evaluate neutral labels
            arguments
                self IMNeutralLabelModesBasis
                name
                rule (1,1) function_handle
            end
            name = IMInternalModesNormalizationTools.customName(name,"IMNeutralLabelModesBasis");
            self.normalizationRules{name} = rule;
        end

        function self = set.normalization(self,name)
            name = IMInternalModesNormalizationTools.name(name,"IMNeutralLabelModesBasis");
            IMBasisContracts.validateNormalization(self,name);
            self.normalization = name;
        end

        function self = withNormalization(self,name)
            % Select a structurally available normalization without another solve.
            % - Topic: Evaluate neutral labels
            arguments
                self IMNeutralLabelModesBasis
                name
            end
            name = IMInternalModesNormalizationTools.name(name,"IMNeutralLabelModesBasis");
            if ~ismember(name,self.normalizationNames())
                error("IMNeutralLabelModesBasis:UnsupportedNormalization","This basis does not define normalization %s.",name);
            end
            self.normalization = name;
        end

        function names = normalizationNames(self)
            % List exact and custom profile normalization rules.
            % - Topic: Evaluate neutral labels
            names = unique(["rmsG";"maxAbsG";reshape(string(keys(self.normalizationRules)),[],1)]);
        end

        function factors = normalizationFactors(self,normalization)
            % Return factors dividing both G and its APV diagnostic.
            % - Topic: Evaluate neutral labels
            arguments
                self IMNeutralLabelModesBasis
                normalization = self.normalization
            end
            normalization = IMInternalModesNormalizationTools.name(normalization,"IMNeutralLabelModesBasis");
            if isKey(self.normalizationRules,normalization)
                rule = self.normalizationRules{normalization};
                factors = zeros(1,self.nColumns);
                for j = 1:self.nColumns
                    value = rule(self,j);
                    IMInternalModesNormalizationTools.validateCustomFactor(value,"IMNeutralLabelModesBasis");
                    factors(j) = value;
                end
            else
            switch normalization
                case "rmsG"
                    factors = ones(1,self.nColumns);
                case "maxAbsG"
                    lengths = diff(self.neutralIntervals,1,2);
                    factors = sqrt(2*diff(self.zDomain)./reshape(lengths(self.intervalIndex),1,[]));
                otherwise
                    error("IMNeutralLabelModesBasis:UnsupportedNormalization","The basis does not define normalization %s.",normalization);
            end
            end
            if ~isreal(factors) || any(~isfinite(factors) | factors <= 0)
                error("IMNeutralLabelModesBasis:NullNormalization","Normalization factors must be positive and finite.");
            end
        end

        function values = evaluate(self,z,options)
            % Evaluate a supported label channel using one normalization.
            % - Topic: Evaluate neutral labels
            % - Declaration: values = evaluate(self,z,options)
            % - Parameter z: physical coordinate column within zDomain
            % - Parameter options.variable: required channel G or APV
            % - Parameter options.derivativeOrder: supported nonnegative vertical derivative order, zero by default
            % - Parameter options.normalization: normalization rule for this call
            % - Returns values: sample-by-mode profiles or derivatives
            arguments (Input)
                self IMNeutralLabelModesBasis
                z (:,1) double {mustBeReal,mustBeFinite}
                options.variable (1,1) string = ""
                options.derivativeOrder (1,1) double {mustBeInteger,mustBeNonnegative} = 0
                options.normalization = self.normalization
            end
            arguments (Output)
                values double
            end
            variable = options.variable;
            derivativeOrder = options.derivativeOrder;
            if variable == ""
                error("IMNeutralLabelModesBasis:MissingVariable","Specify the channel with variable= followed by a supported channel name.");
            end
            if ~ismember(variable,["G","APV"])
                error("IMNeutralLabelModesBasis:UnsupportedVariable","Unsupported channel %s. Inspect capabilities for supported channels and derivatives.",variable);
            end
            IMBasisContracts.validateEvaluation(self,variable,derivativeOrder);
            IMBasisContracts.validateNormalization(self,options.normalization);
            if variable == "G" && derivativeOrder <= 2
                values = self.evaluateDerivative(z,derivativeOrder);
            elseif variable == "APV" && derivativeOrder == 0
                values = -self.f0*self.evaluateDerivative(z,1);
            else
                error("IMNeutralLabelModesBasis:UnsupportedDerivative","G supports orders zero through two; APV supports order zero.");
            end
            values = values./self.normalizationFactors(options.normalization);
        end

        function values = G(self,z,options)
            % Evaluate dimensionless label profiles, zero outside their interval.
            % - Topic: Evaluate neutral labels
            % - Parameter z: physical coordinate column within zDomain
            % - Returns values: sample-by-mode label profiles
            arguments (Input)
                self IMNeutralLabelModesBasis
                z (:,1) double {mustBeReal,mustBeFinite}
                options.normalization = self.normalization
            end
            arguments (Output)
                values (:,:) double
            end
            values = self.evaluate(z,variable="G",derivativeOrder=0,normalization=options.normalization);
        end

        function values = Gz(self,z,options)
            % Evaluate interval-interior first derivatives, in inverse meters.
            %
            % At interval endpoints use the one-sided trace from inside the
            % neutral interval; derivatives may jump to zero outside it.
            % Quadrature must resolve intervals and use these traces consistently.
            %
            % - Topic: Evaluate neutral labels
            % - Parameter z: physical coordinate column within zDomain
            % - Returns values: sample-by-mode first derivatives
            arguments (Input)
                self IMNeutralLabelModesBasis
                z (:,1) double {mustBeReal,mustBeFinite}
                options.normalization = self.normalization
            end
            arguments (Output)
                values (:,:) double
            end
            values = self.evaluate(z,variable="G",derivativeOrder=1,normalization=options.normalization);
        end

        function values = Gzz(self,z,options)
            % Evaluate classical interval-interior second derivatives.
            %
            % Values use the inside trace at interval endpoints. This is not
            % the distributional derivative of the zero-extended Gz.
            %
            % - Topic: Evaluate neutral labels
            % - Parameter z: physical coordinate column within zDomain
            % - Returns values: sample-by-mode second derivatives in inverse meters squared
            arguments (Input)
                self IMNeutralLabelModesBasis
                z (:,1) double {mustBeReal,mustBeFinite}
                options.normalization = self.normalization
            end
            arguments (Output)
                values (:,:) double
            end
            values = self.evaluate(z,variable="G",derivativeOrder=2,normalization=options.normalization);
        end

        function values = APV(self,z,options)
            % Evaluate APV per unit displacement amplitude, q = -f0 Gz.
            % - Topic: Evaluate neutral labels
            % - Parameter z: physical coordinate column within zDomain
            % - Returns values: sample-by-mode APV profiles in inverse meter seconds
            arguments (Input)
                self IMNeutralLabelModesBasis
                z (:,1) double {mustBeReal,mustBeFinite}
                options.normalization = self.normalization
            end
            arguments (Output)
                values (:,:) double
            end
            values = self.evaluate(z,variable="APV",derivativeOrder=0,normalization=options.normalization);
        end

        function gram = gramMatrix(self,options)
            % Return exact continuous label or potential-enstrophy pairings.
            %
            % G uses integral G_i G_j / D = delta_ij. APV uses the
            % unscaled integral q_i q_j = f0^2 D m_n^2 delta_ij, with no
            % factor of one half. APV is degenerate on every column when
            % f0 is zero.
            %
            % - Topic: Inspect neutral labels
            % - Parameter options.variable: G (default) or APV
            % - Returns gram: exact continuous Gram matrix
            arguments (Input)
                self IMNeutralLabelModesBasis
                options.variable (1,1) string {mustBeMember(options.variable,["G","APV"])} = "G"
                options.normalization = self.normalization
            end
            arguments (Output)
                gram (:,:) double
            end
            if options.variable == "G"
                gram = eye(numel(self.m_n));
            else
                gram = diag(self.f0^2*diff(self.zDomain)*self.m_n.^2);
            end
            factors = self.normalizationFactors(options.normalization);
            gram = gram./(factors.'*factors);
        end

        projection = projectionOnGrid(self,z,weights,options)
    end

    methods
        transform = scalarTransform(self,z,weights,options)
    end

    methods (Hidden)
        function self = withConstructionRequest(self,request,solver,selection)
            % Attach original creation evidence without evaluating the modes.
            arguments
                self
                request (1,1) struct
                solver
                selection (1,1) struct = struct()
            end
            state = self.basisContractState;
            state.request = request;
            state.requestedSolver = IMBasisContracts.solverRecord(solver,false);
            if isfield(selection,"sortIndex"), state.sourceEigenpairIndices = reshape(selection.sortIndex,1,[]); end
            if isfield(selection,"modeSelectionDiagnostics"), state.selectionDiagnostics = selection.modeSelectionDiagnostics; end
            self.basisContractState = state;
        end
        recipe = projectionRecipe(self,options)
    end

    methods (Access = private)
        function values = evaluateDerivative(self,z,order)
            if any(z < self.zDomain(1) | z > self.zDomain(2))
                error("IMNeutralLabelModesBasis:OutsideDomain","Evaluation coordinates must lie inside zDomain.");
            end
            values = zeros(numel(z),numel(self.m_n));
            for iInterval = 1:size(self.neutralIntervals,1)
                bounds = self.neutralIntervals(iInterval,:);
                rows = z >= bounds(1) & z <= bounds(2);
                if ~any(rows)
                    continue
                end
                columns = self.intervalIndex == iInterval;
                m = self.m_n(columns);
                phase = (z(rows)-bounds(1))*m;
                amplitude = sqrt(2*diff(self.zDomain)/diff(bounds));
                shape = sin(phase + order*pi/2);
                values(rows,columns) = amplitude*(m.^order).*shape;
                if order == 0 || order == 2
                    endpointRows = z == bounds(1) | z == bounds(2);
                    values(endpointRows,columns) = 0;
                end
            end
        end
    end
end
