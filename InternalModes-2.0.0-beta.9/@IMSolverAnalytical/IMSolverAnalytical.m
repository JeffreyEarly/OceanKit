classdef IMSolverAnalytical
    % Solve supported family problems with their declared exact background.
    %
    % This backend obtains the physical profile and source parameters from
    % the problem's authoritative IMBackground. It returns the same
    % capable analytical family basis as direct analytical construction.
    % Counts include every returned column; explicit branch requests never
    % append an unrequested surface, negative, zero or barotropic sector.
    %
    % ```matlab
    % background = IMBackground.constantStratification(N0=1e-2,zDomain=[-1000 0],surfaceCondition="freeSurface");
    % problem = IMInternalModes.waveModesAtWavenumber(background,k=1e-3);
    % solver = IMSolverAnalytical();
    % modes = solver.solveEVP(problem,nModes=4);
    % ```
    %
    % - Topic: Create an analytical solver
    % - Topic: Solve analytical problems
    % - Declaration: classdef IMSolverAnalytical

    methods
        function self = IMSolverAnalytical()
            % Create an exact family backend without duplicating physical state.
            % - Topic: Create an analytical solver
            % - Returns self: analytical family backend
        end

        function basis = solveEVP(self,problem,options)
            % Solve a supported scalar family with an exact counted inventory.
            %
            % Supply an explicit nModes total count. With branch specified,
            % supply one count for each requested branch. Zero is valid.
            % Analytical construction may resolve additional candidates, but
            % they are never appended to the retained output inventory.
            %
            % - Topic: Solve analytical problems
            % - Parameter problem: background-aware internal-mode descriptor
            % - Parameter options.nModes: required total count, or counts matching branch
            % - Parameter options.branch: optional explicit branch names
            % - Parameter options.normalization: optional active profile normalization
            % - Returns basis: capable analytical basis with aligned metadata
            arguments
                self (1,1) IMSolverAnalytical
                problem (1,1) IMInternalModes
                options.nModes (1,:) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = []
                options.branch (1,:) string = ""
                options.normalization = []
            end
            if isempty(options.nModes)
                error("IMSolver:MissingModeCount","Supply an explicit nModes count. Use nModes=0 to request an empty basis.");
            end
            IMModeSelection.validateRequest(problem,options.nModes,options.branch);
            if ismember(string(problem.name),["meanDensityAnomalyModes","geostrophicGeneralizedPotentialEnstrophyModes"])
                error("IMSolverAnalytical:UnsupportedFamily","The analytical backend does not implement %s. Choose a numerical solver.",problem.name);
            end
            if isempty(problem.background)
                error("IMSolverAnalytical:MissingBackground","Construct the problem from a source-specific IMBackground to use the analytical backend.");
            end
            prefixRequest = isscalar(options.branch) && ismember(options.branch,["","all"]) && isscalar(options.nModes);
            onlyBarotropic = numel(options.branch) == numel(options.nModes) && all(options.nModes(options.branch ~= "barotropic") == 0) && sum(options.nModes(options.branch == "barotropic")) <= 1;
            exactConstant = problem.name == "hydrostaticFModes" && problem.background.surfaceCondition == "rigidLid" && ((prefixRequest && options.nModes <= 1) || onlyBarotropic);
            if exactConstant
                source = IMAnalyticalInternalModesBasis(solution=[],evp=problem,h=Inf,modeNumber=0,N2=problem.N2,rawVariableFunction=@(~,variable,z) double(variable=="F")*ones(numel(z),1),rawUzFunction=@(~,z) zeros(numel(z),1),metadata=struct(analyticalFamily="barotropic"));
                selection = IMModeSelection.selectCounted(problem,source.eigenvalues(:),options.nModes,options.branch,zeroIndices=1);
                basis = source.selectColumns(reshape(selection.sortIndex,1,[]),modeInfo=selection.modeInfo);
                basis = basis.withConstructionRequest(struct(nModes=options.nModes,branch=options.branch),self,selection);
                if ~isempty(options.normalization), basis.normalization = options.normalization; end
                return
            end
            solution = problem.background.analyticalSolution();
            % Existing exact catalogs count total columns. Three additional
            % candidates cover finite endpoint/zero sectors when selecting a
            % regular branch, without changing the returned count contract.
            count = max(1,sum(options.nModes)+3);
            if problem.name == "hydrostaticFModes" && problem.background.surfaceCondition == "freeSurface"
                pairedProblem = IMInternalModes.geostrophicAPVModes(problem.background,g0=Inf,gd=Inf);
                paired = solution.internalModes(pairedProblem,nModes=count);
                source = IMAnalyticalInternalModesBasis(solution=solution,evp=problem,h=paired.h,modeNumber=paired.modeNumber,N2=problem.N2,rawVariableFunction=@(~,variable,z) IMSolverAnalytical.pairedVariable(paired,variable,z),rawUzFunction=@(~,z) paired.rawUz(z),metadata=paired.metadata);
            else
                source = solution.internalModes(problem,nModes=count);
            end
            selection = IMModeSelection.selectCounted(problem,source.eigenvalues(:),options.nModes,options.branch,zeroIndices=find(source.modeNumber==0));
            basis = source.selectColumns(reshape(selection.sortIndex,1,[]),modeInfo=selection.modeInfo);
            basis = basis.withConstructionRequest(struct(nModes=options.nModes,branch=options.branch),self,selection);
            if ~isempty(options.normalization), basis.normalization = options.normalization; end
        end

        [collection,diagnostics] = solveWaveModesAtWavenumbers(self,background,kappa,options)

        function basis = solveGeostrophicZeroAPVModes(self,problem)
            % Solve explicitly requested canonical endpoint responses exactly.
            % - Topic: Solve analytical problems
            % - Parameter problem: background-aware canonical endpoint problem
            % - Returns basis: capable canonical analytical endpoint basis
            arguments
                self (1,1) IMSolverAnalytical
                problem (1,1) IMGeostrophicZeroAPVModes
            end
            if isempty(problem.background)
                error("IMSolverAnalytical:MissingBackground","Construct the endpoint problem from a source-specific IMBackground to use the analytical backend.");
            end
            solution = problem.background.analyticalSolution();
            basis = solution.geostrophicZeroAPVModesAtWavenumber(problem.k,endpoints=problem.endpoints,surfaceBoundary=problem.surfaceBoundary,metadata=problem.metadata);
            basis = basis.withProblem(problem);
        end
    end

    methods (Static, Access = private)
        function values = pairedVariable(paired,variable,z)
            values = paired.rawVariable(variable,z);
        end
    end
end
