classdef (Hidden, Sealed) IMModeSelection
    % Select scientific ladders separately from finite exceptional sectors.
    methods (Static)
        function validateRequest(evp,nModes,branch)
            explicit = ~(isscalar(branch) && ismember(branch,["","all"]));
            if (~explicit && ~isscalar(nModes)) || (explicit && numel(branch) ~= numel(nModes)) || numel(unique(branch)) ~= numel(branch)
                error("IMSolver:InvalidModeRequest","Use one total nModes count, or one count for each distinct explicitly requested branch.");
            end
            if ~explicit, return; end
            known = ["negative","zeroEigenvalue","positive"];
            if isa(evp,"IMInternalModes")
                if evp.name == "waveModesAtWavenumber"
                    known = ["surface","internal","subinertial"];
                elseif evp.name == "waveModesAtFrequency"
                    known = ["surface","propagating","evanescent"];
                elseif evp.name == "hydrostaticFModes"
                    known = ["negative","barotropic","positive"];
                elseif evp.modeFamily == "meanDensityAnomaly"
                    known = ["negative","constant","positive"];
                end
            end
            unknown = branch(~ismember(branch,known));
            if ~isempty(unknown)
                error("IMSolver:UnknownBranch","Unknown branch %s for %s. Supported branches are %s.",join(unknown,", "),evp.name,join(known,", "));
            end
        end

        function plan = wavePagePlan(kappa,nModes)
            if isempty(kappa)
                error("IMSolver:EmptyWavenumberCollection","Supply at least one requested wavenumber.");
            end
            if isempty(nModes)
                error("IMSolver:MissingModeCount","Supply nModes for every requested page, including kappa=0. Use a scalar or one count per entry of kappa; zero is valid.");
            end
            if isscalar(nModes)
                counts = repmat(nModes,size(kappa));
            elseif numel(nModes) == numel(kappa)
                counts = nModes;
            else
                error("IMSolver:InvalidWaveModeCounts","nModes must be scalar or contain one count per entry of kappa, including zero and repeated wavenumbers.");
            end
            distinct = unique(kappa,"stable");
            solveCounts = arrayfun(@(k) max(counts(kappa==k)),distinct);
            [requests,~,mapping] = unique([kappa(:),counts(:)],"rows","stable");
            plan = struct(kappa=kappa,pageCounts=counts,solvedKappa=distinct,solveCounts=solveCounts,basisKappa=requests(:,1).',basisCounts=requests(:,2).',basisIndex=reshape(mapping,1,[]));
        end

        function collection = wavePrefixCollection(sources,plan,solver)
            bases = cell(size(plan.basisCounts));
            for j = 1:numel(bases)
                sourceIndex = find(plan.solvedKappa==plan.basisKappa(j),1);
                source = sources{sourceIndex};
                count = plan.basisCounts(j);
                basis = source;
                if count ~= source.nColumns, basis = source.selectColumns(1:count); end
                sharedSolve = struct(kappa=plan.basisKappa(j),nModes=plan.solveCounts(sourceIndex),request=source.constructionRecord.selection.request);
                bases{j} = basis.withConstructionRequest(struct(nModes=count,branch="",sharedSolve=sharedSolve),solver);
            end
            collection = IMBasisCollection(bases,kappa=plan.kappa,basisIndex=plan.basisIndex);
        end

        function selection = selectCounted(evp,lambda,nModes,branch,options)
            % Select an exact total inventory, or explicit counted branches.
            arguments (Input)
                evp IMEigenvalueProblem
                lambda (:,1) double {mustBeReal,mustBeFinite}
                nModes (1,:) double {mustBeInteger,mustBeNonnegative,mustBeFinite}
                branch (1,:) string = ""
                options.solver = []
                options.A double = []
                options.diagnostics (1,1) struct = struct()
                options.zeroIndices (1,:) double = []
            end
            IMModeSelection.validateRequest(evp,nModes,branch);
            explicit = ~(isscalar(branch) && ismember(branch,["","all"]));
            diagnostics = options.diagnostics;
            wave = isa(evp,"IMInternalModes") && ismember(string(evp.name),["waveModesAtWavenumber","waveModesAtFrequency"]);
            if wave
                positive = find(lambda > 0);
                [~,order] = sort(lambda(positive)); positive = positive(order);
                negative = find(lambda < 0);
                [~,order] = sort(lambda(negative),"descend"); negative = negative(order);
                surface = zeros(0,1);
                if evp.surfaceBoundary.isEigenvalueDependent() && ~isempty(positive)
                    surface = positive(1); positive(1) = [];
                end
                if evp.name == "waveModesAtWavenumber"
                    k = evp.parameters.k;
                    if k > 0
                        negative = negative(evp.f0^2+evp.g*k^2./lambda(negative) > 0);
                    end
                    index = [surface;positive;negative];
                    branches = [repmat("surface",numel(surface),1);repmat("internal",numel(positive),1);repmat("subinertial",numel(negative),1)];
                    numbers = [ones(numel(surface),1);(1:numel(positive))';(1:numel(negative))'];
                else
                    k2 = (evp.parameters.omega^2-evp.f0^2)*lambda/evp.g;
                    propagating = find(k2 > 0 & ~ismember((1:numel(lambda))',surface));
                    evanescent = find(k2 < 0 & ~ismember((1:numel(lambda))',surface));
                    [~,order] = sort(abs(lambda(propagating))); propagating = propagating(order);
                    [~,order] = sort(abs(lambda(evanescent))); evanescent = evanescent(order);
                    index = [surface;propagating;evanescent];
                    branches = [repmat("surface",numel(surface),1);repmat("propagating",numel(propagating),1);repmat("evanescent",numel(evanescent),1)];
                    numbers = [ones(numel(surface),1);(1:numel(propagating))';(1:numel(evanescent))'];
                end
                labels = 1:numel(index);
                labels(lambda(index)<0) = -1:-1:-nnz(lambda(index)<0);
            else
                if ~isempty(options.solver) && ~isempty(lambda)
                    base = evp.selectModes(lambda,max(1,numel(lambda)),options.solver,options.A,diagnostics=diagnostics);
                    index = base.sortIndex;
                    labels = base.modeNumber;
                    diagnostics = base.modeSelectionDiagnostics;
                else
                    negative = find(lambda < 0); [~,order] = sort(lambda(negative)); negative = negative(order);
                    zero = options.zeroIndices(:);
                    positive = find(lambda > 0 & ~ismember((1:numel(lambda))',zero));
                    [~,order] = sort(lambda(positive)); positive = positive(order);
                    index = [negative;zero;positive];
                    labels = [-1:-1:-numel(negative),zeros(1,numel(zero)),1:numel(positive)];
                end
                branches = repmat("positive",numel(index),1);
                branches(labels < 0) = "negative";
                branches(labels == 0) = "zeroEigenvalue";
                if isa(evp,"IMInternalModes") && evp.modeFamily == "meanDensityAnomaly"
                    branches(labels == 0) = "constant";
                end
                if isa(evp,"IMInternalModes") && evp.name == "hydrostaticFModes"
                    branches(labels == 0) = "barotropic";
                end
                numbers = abs(labels(:));
            end
            inventory = table(branches+":"+string(numbers),branches,numbers,lambda(index),VariableNames=["columnLabel","branch","branchIndex","eigenvalue"]);
            if explicit
                positions = zeros(0,1);
                for j = 1:numel(branch)
                    available = find(inventory.branch == branch(j));
                    if numel(available) < nModes(j)
                        error("IMSolver:InsufficientModes","Requested %d %s columns, but only %d are available. Increase resolution or revise the request.",nModes(j),branch(j),numel(available));
                    end
                    positions = [positions;available(1:nModes(j))]; %#ok<AGROW>
                end
            else
                if numel(index) < nModes
                    error("IMSolver:InsufficientModes","Requested %d total columns, but only %d are available. Increase resolution or revise the request.",nModes,numel(index));
                end
                positions = (1:nModes)';
            end
            selection = struct("sortIndex",index(positions),"modeNumber",reshape(labels(positions),1,[]),"modeSelectionDiagnostics",diagnostics,"modeInfo",inventory(positions,:));
        end

    end
end
