classdef IMBasisCollection
    % Preserve continuous bases and their requested wavenumber-page mapping.
    %
    % A collection keeps the original scientific columns and normalization.
    % It performs no solve, truncation, acceptance decision, or approximate
    % sharing. `basisIndex` maps requested pages to stored basis objects;
    % `sourcePage` selects a page within a multi-wavenumber boundary basis.
    % Supported bases are value objects, so subsequent changes to the input
    % basis normalization do not change this collection.
    % Evaluation can override normalization for one call without changing
    % these stored defaults or the collection's scientific column identity.
    %
    % ```matlab
    % collection = IMBasisCollection({basis1,basis2},kappa=[k2 k1 k2],basisIndex=[2 1 2]);
    % G = collection.evaluate(z,variable="G",pages=[3 1]);
    % ```
    %
    % - Topic: Create collections
    % - Topic: Inspect collections
    % - Topic: Evaluate collections
    % - Declaration: classdef IMBasisCollection

    properties (SetAccess = private)
        % Continuous value bases, in stored-basis order.
        % - Topic: Inspect collections
        bases (1,:) cell
        % Requested horizontal wavenumbers, in requested-page order.
        %
        % Empty for collections without a wavenumber identity. Values are
        % never rounded or replaced by representative wavenumbers.
        % - Topic: Inspect collections
        kappa (1,:) double
        % Requested-page to stored-basis mapping.
        % - Topic: Inspect collections
        basisIndex (1,:) double
        % Page within each mapped basis, normally one for scalar solves.
        % - Topic: Inspect collections
        sourcePage (1,:) double
        % Scientific column identity and evaluation capabilities per basis.
        %
        % `derivativeOrders` is a cell row aligned with `variables`. These
        % describe evaluation only; no projection capability is implied.
        % Mode labels, endpoint labels, and wavenumber pages remain separate.
        % Frequency signs are not introduced by a vertical basis collection.
        % - Topic: Inspect collections
        metadata (1,:) struct
    end

    methods (Access = private)
        [projection,recipe] = prepareProjection(self,z,weights,options)
        basis = selectProjectionBasis(self,page,columns)
        columns = allColumnIndices(self,page)
    end

    methods
        function self = IMBasisCollection(bases,options)
            % Construct a collection without solving or resampling modes.
            %
            % Multi-page boundary bases require explicit `sourcePage` when
            % mapping multiple requested pages. When all mapped bases have
            % known wavenumbers, omitted `kappa` is inferred from them.
            %
            % - Topic: Create collections
            % - Parameter bases: nonempty cell row of continuous value bases
            % - Parameter options.kappa: requested nonnegative wavenumber row
            % - Parameter options.basisIndex: mapping to stored bases
            % - Parameter options.sourcePage: page within each mapped basis
            % - Returns self: immutable collection of continuous bases
            arguments (Input)
                bases (1,:) cell
                options.kappa (1,:) double {mustBeReal,mustBeFinite,mustBeNonnegative} = []
                options.basisIndex (1,:) double {mustBeInteger,mustBePositive} = []
                options.sourcePage (1,:) double {mustBeInteger,mustBePositive} = []
            end
            if isempty(bases)
                error("IMBasisCollection:EmptyCollection","Supply at least one continuous basis.");
            end
            for iBasis = 1:numel(bases)
                basis = bases{iBasis};
                supported = isa(basis,"IMBasisSet") || isa(basis,"IMAnalyticalInternalModesBasis") || isa(basis,"IMGeostrophicZeroAPVModesBasis") || isa(basis,"IMAnalyticalGeostrophicZeroAPVModesBasis") || isa(basis,"IMNeutralLabelModesBasis") || isa(basis,"IMZeroDepthWaveModesBasis");
                if ~supported || ~isscalar(basis) || isa(basis,"handle")
                    error("IMBasisCollection:UnsupportedBasis","Entry %d must be a supported scalar continuous value basis.",iBasis);
                end
                entry = IMBasisCollection.describeBasis(basis);
                if iBasis == 1
                    self.metadata = entry;
                else
                    self.metadata(iBasis) = entry;
                end
            end
            self.bases = bases;
            self.basisIndex = options.basisIndex;
            if isempty(self.basisIndex)
                self.basisIndex = 1:numel(bases);
            end
            if any(self.basisIndex > numel(bases))
                error("IMBasisCollection:InvalidBasisIndex","basisIndex must select existing entries in bases.");
            end
            nPages = numel(self.basisIndex);
            self.sourcePage = options.sourcePage;
            if isempty(self.sourcePage)
                self.sourcePage = ones(1,nPages);
            elseif numel(self.sourcePage) ~= nPages
                error("IMBasisCollection:InvalidSourcePage","sourcePage must contain one index per requested page.");
            end
            if ~isempty(options.kappa) && numel(options.kappa) ~= nPages
                error("IMBasisCollection:InvalidKappaCount","kappa must contain one value per requested page.");
            end
            knownKappa = nan(1,nPages);
            for iPage = 1:nPages
                basis = bases{self.basisIndex(iPage)};
                if isprop(basis,"horizontalKind") && ismember(basis.horizontalKind,["perColumn","unspecified"]) && ~isempty(options.kappa)
                    error("IMBasisCollection:UnsupportedHorizontalPage","A per-column or unspecified horizontal identity cannot be assigned one scalar kappa. Preserve the declared dispersion or use a fixed-wavenumber family.");
                end
                basisKappa = IMBasisCollection.basisWavenumbers(basis);
                if self.sourcePage(iPage) > max(1,numel(basisKappa))
                    error("IMBasisCollection:InvalidSourcePage","Requested page %d selects a nonexistent source page.",iPage);
                end
                if ~isempty(basisKappa)
                    knownKappa(iPage) = basisKappa(self.sourcePage(iPage));
                end
            end
            self.kappa = options.kappa;
            if isempty(self.kappa)
                if all(isfinite(knownKappa))
                    self.kappa = knownKappa;
                end
            elseif numel(self.kappa) ~= nPages
                error("IMBasisCollection:InvalidKappaCount","kappa must contain one value per requested page.");
            elseif any(self.kappa(isfinite(knownKappa)) ~= knownKappa(isfinite(knownKappa)))
                error("IMBasisCollection:KappaMismatch","Requested kappa must equal the underlying solved wavenumber. Approximate sharing is not supported.");
            end
        end

        values = evaluate(self,z,options)
        evaluateChunks(self,z,consumer,options)
    end

    methods (Static, Access = private)
        function entry = describeBasis(basis)
            description = basis.capabilities;
            entry = struct("family",string(class(basis)),"columnKind","mode","columnLabels",basis.columnLabels,"modeNumber",[],"endpointLabels",strings(1,0),"normalization","","zDomain",basis.zDomain,"variables",description.variables,"derivativeOrders",{description.derivativeOrders},"capabilities",description,"horizontalIdentity",struct(),"scientificSpecification",basis.scientificSpecification,"constructionRecord",basis.constructionRecord,"boundaryResponseMatrixStatus","notApplicable");
            if isprop(basis,"normalization"), entry.normalization = string(basis.normalization); end
            if isprop(basis,"modeNumber"), entry.modeNumber = reshape(basis.modeNumber,1,[]); end
            if isprop(basis,"evp"), entry.family = string(basis.evp.name); end
            if isprop(basis,"family"), entry.family = string(basis.family); end
            if isa(basis,"IMGeostrophicZeroAPVModesBasis") || isa(basis,"IMAnalyticalGeostrophicZeroAPVModesBasis")
                entry.family = "geostrophicZeroAPVModes";
                entry.columnKind = "endpoint";
                entry.endpointLabels = reshape(string(basis.endpoints),1,[]);
                entry.boundaryResponseMatrixStatus = "availableOnBasis";
                if entry.normalization == "", entry.normalization = string(basis.normalizationConvention); end
            end
            if isprop(basis,"horizontalKind")
                entry.horizontalIdentity = struct(kind=basis.horizontalKind,k=basis.k,kSquared=basis.kSquared);
            else
                entry.horizontalIdentity = struct(kind="unspecified",k=[],kSquared=[]);
            end
            if numel(entry.columnLabels) ~= basis.nColumns || numel(unique(entry.columnLabels)) ~= numel(entry.columnLabels)
                error("IMBasisCollection:InvalidColumnLabels","Every retained column requires a distinct scientific label; empty inventories remain empty.");
            end
        end

        function kappa = basisWavenumbers(basis)
            kappa = [];
            if isprop(basis,"horizontalKind")
                if ismember(basis.horizontalKind,["fixedWavenumber","fixedWavenumberPages"])
                    kappa = reshape(basis.k,1,[]);
                end
            elseif isprop(basis,"k")
                kappa = reshape(basis.k,1,[]);
            elseif isprop(basis,"evp") && isfield(basis.evp.parameters,"k")
                kappa = reshape(basis.evp.parameters.k,1,[]);
            end
        end
    end
end
