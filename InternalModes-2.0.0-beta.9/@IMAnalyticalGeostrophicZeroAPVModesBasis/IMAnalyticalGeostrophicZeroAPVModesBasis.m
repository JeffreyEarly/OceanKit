classdef IMAnalyticalGeostrophicZeroAPVModesBasis
    % Evaluate canonical and generalized-energy orthogonal boundary-buoyancy modes.
    %
    % Canonical columns have unit surface or bottom buoyancy response. At a
    % free surface that response is G-F; the bottom response is G. The
    % profiles obey Fz=-N2 G/g and Gz=-g k^2 F/f0^2, including neutral layers.
    % Select endpoints in the problem, solve once, and optionally construct
    % the manuscript's generalized-energy orthogonal basis before selecting
    % columns. Profiles have shape nZ-by-nColumns-by-nWavenumbers.
    %
    % ```matlab
    % boundary = solver.solveGeostrophicZeroAPVModes(problem);
    % orthogonal = boundary.orthogonalizeGeneralizedEnergy(g0=-.035,gd=.01);
    % hbb = orthogonal.boundaryDepth;
    % ```
    %
    % - Topic: Evaluate boundary-buoyancy modes
    % - Topic: Inspect boundary-buoyancy modes
    % - Topic: Choose boundary-buoyancy coordinates
    % - Topic: Form boundary-buoyancy quadratic forms
    % - Declaration: classdef IMAnalyticalGeostrophicZeroAPVModesBasis

    properties (SetAccess = private)
        % Canonical physical endpoint problem.
        % - Topic: Inspect boundary-buoyancy modes
        problem
        % Analytical solution that supplied the canonical profiles.
        % - Topic: Inspect boundary-buoyancy modes
        solution
        % Positive horizontal-wavenumber pages.
        % - Topic: Inspect boundary-buoyancy modes
        k
        % Canonical endpoint-row identities, independent of retained columns.
        % - Topic: Inspect boundary-buoyancy modes
        endpoints
        % Signed Coriolis frequency in inverse seconds.
        % - Topic: Inspect boundary-buoyancy modes
        f0
        % Gravitational acceleration.
        % - Topic: Inspect boundary-buoyancy modes
        g
        % Physical free-surface or rigid-lid convention.
        % - Topic: Inspect boundary-buoyancy modes
        surfaceBoundary
        % Physical domain [bottom surface].
        % - Topic: Inspect boundary-buoyancy modes
        zDomain
        % Shape-preserving buoyancy-frequency-squared profile.
        % - Topic: Inspect boundary-buoyancy modes
        N2
        % Construction and orthogonality evidence.
        % - Topic: Inspect boundary-buoyancy modes
        metadata
        % Signed manuscript boundary depths; empty before orthogonalization.
        % - Topic: Inspect boundary-buoyancy modes
        boundaryDepth
    end

    properties (Hidden, SetAccess = private)
        basisContractState = struct()
    end

    properties (Dependent)
        % Active normalization; assigning it changes only this value object.
        % - Topic: Evaluate boundary-buoyancy modes
        normalization
    end

    properties (Dependent, SetAccess = private)
        % Exact retained column count per page.
        % - Topic: Inspect boundary-buoyancy modes
        nColumns
        % Scientific family identity.
        % - Topic: Inspect boundary-buoyancy modes
        family
        % Current scientific column labels.
        % - Topic: Inspect boundary-buoyancy modes
        columnLabels
        % Current column identities and branches.
        % - Topic: Inspect boundary-buoyancy modes
        modeInfo
        % Zero geostrophic frequencies.
        % - Topic: Inspect boundary-buoyancy modes
        frequency
        % Canonical source endpoint identities.
        % - Topic: Inspect boundary-buoyancy modes
        canonicalEndpointLabels
        % Canonical endpoint responses per current coefficient, pagewise.
        % - Topic: Inspect boundary-buoyancy modes
        boundaryResponseMatrix
        % Physical-energy form in current normalized coordinates.
        % - Topic: Inspect boundary-buoyancy modes
        energyMatrix
        % Surface-buoyancy form in current normalized coordinates.
        % - Topic: Inspect boundary-buoyancy modes
        surfaceBuoyancyMatrix
        % Bottom-buoyancy form in current normalized coordinates.
        % - Topic: Inspect boundary-buoyancy modes
        bottomBuoyancyMatrix
        % Total endpoint-response form in current coordinates.
        % - Topic: Inspect boundary-buoyancy modes
        endpointResponseMetric
        % Signs of the active normalization form.
        % - Topic: Inspect boundary-buoyancy modes
        normSignatures
        % Structured supported-variable and operation contract.
        % - Topic: Inspect boundary-buoyancy modes
        capabilities
        % Requested and resolved scientific identity.
        % - Topic: Inspect boundary-buoyancy modes
        scientificSpecification
        % Original construction and later selection evidence.
        % - Topic: Inspect boundary-buoyancy modes
        constructionRecord
        % Fixed-wavenumber horizontal identity.
        % - Topic: Inspect boundary-buoyancy modes
        horizontalKind
        % Empty because every page has one fixed wavenumber.
        % - Topic: Inspect boundary-buoyancy modes
        kSquared
    end

    properties (Access = private)
        canonicalFFunction
        canonicalGFunction
        canonicalEndpointResponseMetric
        canonicalEnergyMatrix
        canonicalSurfaceBuoyancyMatrix
        canonicalBottomBuoyancyMatrix
        endpointDirections_
        columnLabels_
        columnOrdinals_
        normalization_ = "endpointResponse"
        normalizationRules_ = configureDictionary("string","cell")
        normalizationCache_
        isOrthogonal_ = false
    end

    methods
        function self = IMAnalyticalGeostrophicZeroAPVModesBasis(options)
            % Create an exact canonical boundary-normalized basis.
            %
            % Concrete analytical solution families supply exact `F` and
            % `G` evaluators. Users normally construct this basis through
            % `geostrophicZeroAPVModesAtWavenumber` on one of those families.
            %
            % - Topic: Developer topics
            % - Declaration: basisSet = IMAnalyticalGeostrophicZeroAPVModesBasis(options)
            % - Parameter options.solution: analytical solution family
            % - Parameter options.problem: canonical zero-APV problem
            % - Parameter options.FFunction: exact canonical `F` evaluator, `(z,pages)` for bounded page evaluation
            % - Parameter options.GFunction: exact canonical `G` evaluator, `(z,pages)` for bounded page evaluation
            % - Parameter options.metadata: creation metadata
            % - Returns basisSet: exact canonical basis
            % - Developer: true
            arguments
                options.solution IMAnalyticalSolution
                options.problem IMGeostrophicZeroAPVModes
                options.FFunction (1,1) function_handle
                options.GFunction (1,1) function_handle
                options.metadata struct = struct()
            end

            self.solution = options.solution;
            self.problem = options.problem;
            self.k = options.problem.k;
            self.endpoints = options.problem.endpoints;
            self.f0 = options.problem.f0;
            self.g = options.problem.g;
            self.surfaceBoundary = options.problem.surfaceBoundary;
            self.zDomain = options.problem.zDomain;
            self.N2 = options.problem.N2;
            self.canonicalFFunction = options.FFunction;
            self.canonicalGFunction = options.GFunction;

            nEndpoints = numel(self.endpoints);
            nK = numel(self.k);
            [self.canonicalEndpointResponseMetric,self.canonicalEnergyMatrix,self.canonicalSurfaceBuoyancyMatrix,self.canonicalBottomBuoyancyMatrix,responseResidual] = self.formCanonicalQuadraticMatrices();
            self.metadata = options.metadata;
            self.metadata.canonicalEndpointResponseResidual = responseResidual;
            self.endpointDirections_ = repmat(eye(nEndpoints),1,1,nK);
            self.columnLabels_ = reshape(string(self.endpoints),1,[]);
            self.columnOrdinals_ = 1:nEndpoints;
            self.normalizationCache_ = IMInternalModesNormalizationTools.normalizationCache();
            self.metadata.basisChoice = "boundaryNormalized";
            self.basisContractState = IMBasisContracts.initialState(self);
        end
        function value = get.nColumns(self)
            value = size(self.endpointDirections_,2);
        end
        function value = get.family(~)
            value = "geostrophicZeroAPV";
        end
        function value = get.columnLabels(self)
            value = self.columnLabels_;
        end
        function value = get.canonicalEndpointLabels(self)
            value = self.endpoints;
        end
        function value = get.modeInfo(self)
            branch = repmat("endpointResponse",self.nColumns,1);
            if self.isOrthogonal_, branch(:) = "generalizedEnergyOrthogonal"; end
            value = table(self.columnLabels(:),branch,self.columnOrdinals_(:),VariableNames=["columnLabel","branch","branchIndex"]);
        end
        function value = get.frequency(self)
            value = zeros(1,self.nColumns);
        end
        function value = get.normalization(self)
            value = self.normalization_;
        end
        function self = set.normalization(self,name)
            name = IMInternalModesNormalizationTools.name(name,"IMAnalyticalGeostrophicZeroAPVModesBasis");
            if ~isscalar(name) || ~ismember(name,self.normalizationNames())
                error("IMAnalyticalGeostrophicZeroAPVModesBasis:UnsupportedNormalization","Choose a rule from normalizationNames().");
            end
            self.normalization_ = name;
        end
        function value = get.normSignatures(self)
            value = ones(self.nColumns,numel(self.k));
        end
        function value = get.boundaryResponseMatrix(self)
            factors = reshape(self.normalizationFactors(),1,self.nColumns,numel(self.k));
            value = self.endpointDirections_./factors;
        end
        function value = get.energyMatrix(self)
            value = IMGeostrophicZeroAPVFormTools.transformForm(self.canonicalEnergyMatrix,self.boundaryResponseMatrix);
        end
        function value = get.surfaceBuoyancyMatrix(self)
            value = IMGeostrophicZeroAPVFormTools.transformForm(self.canonicalSurfaceBuoyancyMatrix,self.boundaryResponseMatrix);
        end
        function value = get.bottomBuoyancyMatrix(self)
            value = IMGeostrophicZeroAPVFormTools.transformForm(self.canonicalBottomBuoyancyMatrix,self.boundaryResponseMatrix);
        end
        function value = get.endpointResponseMetric(self)
            value = IMGeostrophicZeroAPVFormTools.transformForm(self.canonicalEndpointResponseMetric,self.boundaryResponseMatrix);
        end
        function value = get.capabilities(self)
            value = IMBasisContracts.capabilities(self);
        end
        function value = get.scientificSpecification(self)
            value = IMBasisContracts.scientificSpecification(self);
        end
        function value = get.constructionRecord(self)
            value = IMBasisContracts.constructionRecord(self);
        end
        function value = get.horizontalKind(~)
            value = "fixedWavenumber";
        end
        function value = get.kSquared(~)
            value = zeros(1,0);
        end

        function names = normalizationNames(self)
            % List admissible profile normalizations and installed rules.
            % - Topic: Evaluate boundary-buoyancy modes
            names = ["endpointResponse","rmsF","rmsG","maxAbsF","maxAbsG"];
            names = unique([names,reshape(string(keys(self.normalizationRules_)),1,[])],"stable");
        end

        function factors = normalizationFactors(self,name,options)
            % Return positive divisors for selected columns and source pages.
            % - Topic: Evaluate boundary-buoyancy modes
            arguments
                self IMAnalyticalGeostrophicZeroAPVModesBasis
                name = self.normalization
                options.pages (1,:) double {mustBeInteger,mustBePositive} = 1:numel(self.k)
            end
            name = IMInternalModesNormalizationTools.name(name,"IMAnalyticalGeostrophicZeroAPVModesBasis");
            if ~ismember(name,self.normalizationNames())
                error("IMAnalyticalGeostrophicZeroAPVModesBasis:UnsupportedNormalization","Choose a rule from normalizationNames().");
            end
            if any(options.pages>numel(self.k))
                error("IMAnalyticalGeostrophicZeroAPVModesBasis:InvalidPage","pages must select existing source wavenumber pages.");
            end
            factors = zeros(self.nColumns,numel(options.pages));
            for iPage = 1:numel(options.pages)
                page = options.pages(iPage);
                key = name+"@"+string(page);
                factor = self.normalizationCache_(key);
                if isempty(factor) && self.nColumns>0
                    if isKey(self.normalizationRules_,name)
                        factor = zeros(self.nColumns,1);
                        rule = self.normalizationRules_{name};
                        for j = 1:self.nColumns
                            value = rule(self,j);
                            if ~isnumeric(value) || ~isreal(value) || ~(isscalar(value) || numel(value)==numel(self.k))
                                error("IMAnalyticalGeostrophicZeroAPVModesBasis:InvalidNormalization","A custom rule must return one positive scalar or one value per wavenumber page.");
                            end
                            if isscalar(value), factor(j) = value; else, factor(j) = value(page); end
                        end
                    else
                        factor = IMGeostrophicZeroAPVFormTools.normalizationFactors(self,name,page);
                    end
                    if any(~isfinite(factor) | factor<=0)
                        error("IMAnalyticalGeostrophicZeroAPVModesBasis:NullNormalization","Normalization divisors must be positive and finite for every retained column.");
                    end
                    self.normalizationCache_(key,factor);
                end
                if self.nColumns>0, factors(:,iPage) = factor; end
            end
        end

        function basis = withNormalization(self,name)
            % Choose a copied basis normalization without another solve.
            % - Topic: Choose boundary-buoyancy coordinates
            arguments
                self IMAnalyticalGeostrophicZeroAPVModesBasis
                name
            end
            basis = self;
            basis.normalization = name;
        end

        function basis = addNormalization(self,name,rule)
            % Install or replace a positive divisor rule, factor=rule(basis,j).
            % Built-in names retain their scientific definitions on every
            % family; use a distinct name for a custom alternative.
            % - Topic: Choose boundary-buoyancy coordinates
            arguments
                self IMAnalyticalGeostrophicZeroAPVModesBasis
                name
                rule (1,1) function_handle
            end
            name = IMInternalModesNormalizationTools.customName(name,"IMAnalyticalGeostrophicZeroAPVModesBasis");
            basis = self;
            basis.normalizationRules_{name} = rule;
            basis.normalizationCache_ = IMInternalModesNormalizationTools.normalizationCache();
        end

        function basis = selectColumns(self,columns)
            % Retain distinct existing columns in caller order without a solve.
            % - Topic: Choose boundary-buoyancy coordinates
            arguments
                self IMAnalyticalGeostrophicZeroAPVModesBasis
                columns (1,:) double {mustBeInteger,mustBePositive}
            end
            if any(columns > self.nColumns) || numel(unique(columns)) ~= numel(columns)
                error("IMAnalyticalGeostrophicZeroAPVModesBasis:InvalidColumns","Select distinct existing columns.");
            end
            basis = self;
            for name = reshape(string(keys(self.normalizationRules_)),1,[])
                rule = self.normalizationRules_{name};
                basis.normalizationRules_{name} = @(~,j) rule(self,columns(j));
            end
            basis.endpointDirections_ = self.endpointDirections_(:,columns,:);
            basis.columnLabels_ = self.columnLabels_(columns);
            basis.columnOrdinals_ = self.columnOrdinals_(columns);
            if self.isOrthogonal_, basis.boundaryDepth = self.boundaryDepth(columns,:); end
            basis.normalizationCache_ = IMInternalModesNormalizationTools.normalizationCache();
            basis.basisContractState = IMBasisContracts.selectState(self.basisContractState,columns);
        end

        function matrix = generalizedEnergyMatrix(self,options)
            % Form signed generalized energy in the current coefficient coordinates.
            % - Topic: Form boundary-buoyancy quadratic forms
            % - Parameter options.g0: signed surface coefficient; infinity only for an absent endpoint
            % - Parameter options.gd: signed bottom coefficient; infinity only for an absent endpoint
            % - Returns matrix: current-coordinate H + g0 B0 + gd Bd pages
            arguments
                self IMAnalyticalGeostrophicZeroAPVModesBasis
                options.g0 (1,1) double {mustBeReal} = 0
                options.gd (1,1) double {mustBeReal} = 0
            end
            matrix = IMGeostrophicZeroAPVFormTools.generalizedEnergy(self.energyMatrix,self.surfaceBuoyancyMatrix,self.bottomBuoyancyMatrix,self.endpoints,options.g0,options.gd,"IMAnalyticalGeostrophicZeroAPVModesBasis");
        end

        function basis = orthogonalizeGeneralizedEnergy(self,options)
            % Construct the manuscript's generalized-energy orthogonal modes.
            %
            % Solve Hg C=C diag(boundaryDepth/(2 k^2)), C' C=I, in the
            % complete originally chosen endpoint coordinates. Negative,
            % zero and repeated depths remain valid. A proper retained
            % selection must be made after this operation.
            %
            % - Topic: Choose boundary-buoyancy coordinates
            % - Parameter options.g0: required signed surface coefficient
            % - Parameter options.gd: required signed bottom coefficient
            % - Returns basis: orthogonal modes with unit endpoint-response geometry
            arguments
                self IMAnalyticalGeostrophicZeroAPVModesBasis
                options.g0 (1,1) double {mustBeReal}
                options.gd (1,1) double {mustBeReal}
            end
            if self.nColumns ~= numel(self.endpoints)
                error("IMAnalyticalGeostrophicZeroAPVModesBasis:IncompleteEndpointBasis","Orthogonalize the complete chosen endpoint family before selecting columns.");
            end
            matrix = IMGeostrophicZeroAPVFormTools.generalizedEnergy(self.canonicalEnergyMatrix,self.canonicalSurfaceBuoyancyMatrix,self.canonicalBottomBuoyancyMatrix,self.endpoints,options.g0,options.gd,"IMAnalyticalGeostrophicZeroAPVModesBasis");
            [map,depth,evidence] = IMGeostrophicZeroAPVFormTools.orthogonalizeEnergy(matrix,self.k,"IMAnalyticalGeostrophicZeroAPVModesBasis");
            basis = self;
            basis.endpointDirections_ = map;
            basis.boundaryDepth = depth;
            basis.columnLabels_ = "generalizedEnergy"+string(1:self.nColumns);
            basis.columnOrdinals_ = 1:self.nColumns;
            basis.normalization_ = "endpointResponse";
            basis.isOrthogonal_ = true;
            basis.normalizationCache_ = IMInternalModesNormalizationTools.normalizationCache();
            basis.metadata.basisChoice = "generalizedEnergyOrthogonal";
            basis.metadata.generalizedEnergyCoefficients = struct(g0=options.g0,gd=options.gd);
            basis.metadata.orthogonality = evidence;
            basis.basisContractState.initialColumnIndices = zeros(1,0);
            basis.basisContractState.activeSourceColumns = zeros(1,0);
            basis.basisContractState.sourceDirections = map;
        end

        function projection = projectionOnGrid(self,z,weights,options)
            % Request an explicit scalar physical projection capability.
            % - Topic: Evaluate boundary-buoyancy modes
            arguments
                self IMAnalyticalGeostrophicZeroAPVModesBasis
                z (:,1) double {mustBeReal,mustBeFinite}
                weights (:,1) double {mustBeReal,mustBeFinite}
                options.variable (1,1) string
            end
            projection = IMScalarProjectionTools.projectionOnGrid(self,z,weights,variable=options.variable);
        end

        function transform = scalarTransform(self,z,weights,options)
            % Request an explicit scalar physical transform capability.
            % - Topic: Evaluate boundary-buoyancy modes
            arguments
                self IMAnalyticalGeostrophicZeroAPVModesBasis
                z (:,1) double {mustBeReal,mustBeFinite}
                weights (:,1) double {mustBeReal,mustBeFinite}
                options.variable (1,1) string
            end
            transform = IMScalarProjectionTools.scalarTransform(self,z,weights,variable=options.variable);
        end

        values = evaluate(self,z,options)
        values = F(self,z,options)
        values = G(self,z,options)
        values = Fz(self,z,options)
        values = Gz(self,z,options)
        values = Gzz(self,z,options)

        function summarize(self)
            % Print the scientific family, physical domain and retained count.
            % - Topic: Inspect boundary-buoyancy modes
            fprintf("%s: %d columns on [%g,%g], %d wavenumber pages\n",class(self),self.nColumns,self.zDomain(1),self.zDomain(2),numel(self.k));
        end
    end

    methods (Hidden)
        recipe = projectionRecipe(self,options)

        function values = endpointUnscaledProfiles(self,variable,z,pages)
            % Evaluate the retained coordinate directions before profile divisors.
            if isempty(z) || self.nColumns==0 || isempty(pages)
                values = zeros(numel(z),self.nColumns,numel(pages));
                return
            end
            if variable == "F", evaluator = self.canonicalFFunction; else, evaluator = self.canonicalGFunction; end
            canonical = self.evaluateCanonical(evaluator,z,pages);
            values = IMGeostrophicZeroAPVFormTools.applyMap(canonical,self.endpointDirections_(:,:,pages));
        end

        function map = endpointUnscaledResponseMatrix(self)
            map = self.endpointDirections_;
        end

        function [z,weights] = endpointNormalizationQuadrature(self)
            [z,weights] = legpts(256,self.zDomain);
            z = z(:); weights = weights(:);
        end

        function self = withProblem(self,problem)
            % Retain the authoritative descriptor after exact source dispatch.
            arguments
                self IMAnalyticalGeostrophicZeroAPVModesBasis
                problem (1,1) IMGeostrophicZeroAPVModes
            end
            if ~isequal(problem.zDomain,self.zDomain) || ~isequal(problem.k,self.k) || ~isequal(problem.endpoints,self.endpoints) || problem.g ~= self.g || problem.f0 ~= self.f0 || problem.surfaceBoundary ~= self.surfaceBoundary
                error("IMAnalyticalGeostrophicZeroAPVModesBasis:IncompatibleProblem","The descriptor must represent the same exact physical endpoint problem.");
            end
            self.problem = problem;
            authoritative = IMBasisContracts.initialState(self);
            self.basisContractState.requestedScientific = authoritative.requestedScientific;
            if isfield(authoritative,"physicalBackground")
                self.basisContractState.physicalBackground = authoritative.physicalBackground;
            end
        end
        function tf = sharesCanonicalBasisWith(self,other)
            tf = isa(other,"IMAnalyticalGeostrophicZeroAPVModesBasis") && isequal(self.canonicalFFunction,other.canonicalFFunction) && isequal(self.canonicalGFunction,other.canonicalGFunction) && isequal(self.problem,other.problem);
        end
    end

    methods (Access = private)
        function [responseMetric,physicalEnergy,surfaceBuoyancy,bottomBuoyancy,responseResidual] = formCanonicalQuadraticMatrices(self)
            zSurface = self.zDomain(2);
            zBottom = self.zDomain(1);
            FSurface = self.evaluateCanonical(self.canonicalFFunction,zSurface);
            FBottom = self.evaluateCanonical(self.canonicalFFunction,zBottom);
            GSurface = self.evaluateCanonical(self.canonicalGFunction,zSurface);
            GBottom = self.evaluateCanonical(self.canonicalGFunction,zBottom);
            surfaceResponse = GSurface;
            if self.surfaceBoundary == "freeSurface"
                surfaceResponse = surfaceResponse-FSurface;
            end
            bottomResponse = GBottom;
            fullResponse = cat(1,surfaceResponse,bottomResponse);
            endpointRows = zeros(1,numel(self.endpoints));
            endpointRows(self.endpoints == "surface") = 1;
            endpointRows(self.endpoints == "bottom") = 2;
            selectedResponse = fullResponse(endpointRows,:,:);
            responseResidual = zeros(1,numel(self.k));
            for iK = 1:numel(self.k)
                responseResidual(iK) = norm(selectedResponse(:,:,iK)-eye(numel(self.endpoints)),"fro")/max(1,sqrt(numel(self.endpoints)));
            end
            if any(responseResidual > 1e-10)
                error("IMAnalyticalGeostrophicZeroAPVModesBasis:EndpointNormalizationFailure", "Exact canonical modes did not satisfy unit endpoint responses.");
            end
            [responseMetric,physicalEnergy,surfaceBuoyancy,bottomBuoyancy] = IMGeostrophicZeroAPVFormTools.formMatrices(FSurface,FBottom,surfaceResponse,bottomResponse,selectedResponse,self.endpoints,self.k,self.f0,self.g);
        end
        function values = evaluateCanonical(self,evaluator,z,pages)
            if nargin < 4
                pages = 1:numel(self.k);
            end
            if nargin(evaluator) == 1
                values = evaluator(z(:));
                values = values(:,:,pages);
            else
                values = evaluator(z(:),pages);
            end
            nEndpoints = numel(self.endpoints);
            nK = numel(pages);
            if size(values,1) ~= numel(z) || size(values,2) ~= nEndpoints || size(values,3) ~= nK
                error("IMAnalyticalGeostrophicZeroAPVModesBasis:InvalidEvaluatorShape", "Exact evaluators must return nZ x nEndpoints x nK arrays.");
            end
            if ~isreal(values) || any(~isfinite(values(:)))
                error("IMAnalyticalGeostrophicZeroAPVModesBasis:InvalidEvaluatorValues", "Exact evaluators must return finite real values.");
            end
        end
    end
end
