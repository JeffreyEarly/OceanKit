classdef IMGeostrophicModesBasis < IMInternalModesBasis
    % Store aligned geostrophic pressure and displacement modes.
    %
    % The regular numerical EVP solves both fields,
    % $$F_\xi=-N^2G/g,\qquad G_\xi=rF,\qquad r=1/h.$$
    % `nativeModes` remains the native `F` coefficient array. The companion
    % `G` coefficients are stored separately, so evaluation remains regular
    % wherever `N2` vanishes. Both fields share normalization and orientation.
    % The inverse depth is retained directly from the internal solve rather
    % than recovered by subtracting public physical eigenvalues.
    %
    % ```matlab
    % evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=[-1000 0],g0=0,gd=-0.01,surfaceBoundary="freeSurface");
    % basis = IMSolverSpectral(nEVP=128).solveEVP(evp,nModes=10);
    % G = basis.G(z);
    % ```
    %
    % - Topic: Create internal-mode bases
    % - Topic: Evaluate modes
    % - Topic: Inspect basis sets
    % - Declaration: classdef IMGeostrophicModesBasis < IMInternalModesBasis

    properties (SetAccess = private)
        % Solved inverse equivalent depths, in inverse meters.
        %
        % Exact zero values represent `h=Inf`, including generalized-
        % enstrophy modes whose public eigenvalue is nonzero.
        %
        % - Topic: Inspect basis sets
        inverseDepth
    end

    properties (Access = private)
        nativeGModes
        gSolver
    end

    methods
        function self = IMGeostrophicModesBasis(options)
            % Create a basis from aligned native `F` and `G` columns.
            %
            % - Topic: Create internal-mode bases
            % - Declaration: basisSet = IMGeostrophicModesBasis(options)
            % - Parameter options.solver: configured solver
            % - Parameter options.evp: geostrophic EVP descriptor
            % - Parameter options.nativeModes: native solved `F` columns
            % - Parameter options.nativeGModes: aligned native solved `G` columns
            % - Parameter options.companionSolver: optional G discretization for staggered finite differences
            % - Parameter options.eigenvalues: retained public eigenvalues
            % - Parameter options.inverseDepth: retained solved inverse depths
            % - Parameter options.modeNumber: physical mode labels
            % - Parameter options.modeSelectionDiagnostics: selection diagnostics
            % - Parameter options.normalization: normalization rule name or enum value
            % - Parameter options.metadata: captured problem metadata
            % - Returns basisSet: aligned geostrophic basis set
            arguments
                options.solver IMSolver
                options.evp IMInternalModes
                options.nativeModes (:,:) double {mustBeReal, mustBeFinite}
                options.nativeGModes (:,:) double {mustBeReal, mustBeFinite}
                options.companionSolver = []
                options.eigenvalues (1,:) double {mustBeReal, mustBeFinite}
                options.inverseDepth (1,:) double {mustBeReal, mustBeFinite}
                options.modeNumber (1,:) double {mustBeInteger}
                options.modeSelectionDiagnostics struct = struct()
                options.normalization = []
                options.metadata struct = struct()
            end

            if size(options.nativeModes,2)~=size(options.nativeGModes,2)
                error("IMGeostrophicModesBasis:InconsistentModeArrays", "Native F and G arrays must have the same number of columns.");
            end
            if numel(options.inverseDepth) ~= size(options.nativeModes,2)
                error("IMGeostrophicModesBasis:InvalidInverseDepthCount", "inverseDepth must contain one solved value per F/G mode pair.");
            end
            self@IMInternalModesBasis(solver=options.solver,evp=options.evp,nativeModes=options.nativeModes, ...
                eigenvalues=options.eigenvalues,modeNumber=options.modeNumber, ...
                modeSelectionDiagnostics=options.modeSelectionDiagnostics,normalization=options.normalization,metadata=options.metadata);
            self.nativeGModes = options.nativeGModes;
            if isempty(options.companionSolver)
                self.gSolver = options.solver;
            elseif isa(options.companionSolver,"IMSolver")
                self.gSolver = options.companionSolver;
            else
                error("IMGeostrophicModesBasis:InvalidCompanionSolver","companionSolver must be an IMSolver.");
            end
            self.inverseDepth = options.inverseDepth;
            self.h = 1./self.inverseDepth;
        end

        function Gz = Gz(self,z,options)
            % Evaluate derivatives of the stored companion `G` modes.
            %
            % This differentiates the native representation independently of
            % the first-order differential equations.
            %
            % - Topic: Evaluate modes
            % - Declaration: Gz = Gz(basisSet,z,options)
            % - Parameter z: physical coordinate
            % - Parameter options.normalization: normalization rule name or enum value
            % - Returns Gz: evaluated `G` derivatives
            arguments
                self IMGeostrophicModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
                options.normalization = self.normalization
            end
            Gz = self.gSolver.evaluatePhysicalDerivative(self.nativeGModes,z,1) ./ self.normalizationFactors(options.normalization);
        end

        function [self,signs] = orientModeSigns(self)
            % Apply one deterministic sign to each aligned `F`/`G` pair.
            %
            % Resolved surface traces retain the shallow-interior convention.
            % A trapped signed mode with unresolved surface traces uses its
            % dominant stored G value; metadata records this phase fallback.
            %
            % - Topic: Developer topics — Diagnostic variables
            % - Declaration: basisSet = orientModeSigns(basisSet)
            % - Returns basisSet: basis set with aligned signs oriented
            % - Returns signs: applied sign for each pair
            % - Developer: true
            signs = ones(1,size(self.nativeModes,2));
            if isempty(self.nativeModes)
                return;
            end
            z = self.solver.zNative;
            surfaceIndex = self.solver.boundaryIndex("surface");
            GValues = self.rawVariable("G",z);
            GzSurface = self.gSolver.evaluatePhysicalDerivative(self.nativeGModes,self.zDomain(2),1);
            FValues = self.rawVariable("F",z);
            fallback = false(size(signs));
            try
                signs = IMModeOrientationTools.shallowInteriorGPositive( ...
                    GValues=GValues,GzSurface=GzSurface,FValues=FValues, ...
                    depth=diff(self.zDomain),surfaceIndex=surfaceIndex,allowFFallback=isinf(self.h));
            catch exception
                if string(exception.identifier)~="IMModeOrientationTools:IndeterminateOrientation", rethrow(exception); end
                for j = 1:numel(signs)
                    try
                        signs(j) = IMModeOrientationTools.shallowInteriorGPositive( ...
                            GValues=GValues(:,j),GzSurface=GzSurface(j),FValues=FValues(:,j), ...
                            depth=diff(self.zDomain),surfaceIndex=surfaceIndex,allowFFallback=isinf(self.h(j)));
                    catch modeException
                        [amplitude,index] = max(abs(GValues(:,j)));
                        if string(modeException.identifier)~="IMModeOrientationTools:IndeterminateOrientation" ...
                                || self.inverseDepth(j)>=0 || amplitude==0
                            rethrow(modeException)
                        end
                        signs(j) = sign(GValues(index,j));
                        fallback(j) = true;
                    end
                end
            end
            self.nativeModes = self.nativeModes.*signs;
            self.nativeGModes = self.nativeGModes.*signs;
            self.metadata.modeOrientation = IMModeOrientationTools.convention;
            self.metadata.dominantGOrientationFallback = fallback;
            self = self.resetNormalizationCache();
        end

        function values = rawVariable(self,variable,z)
            % Evaluate stored unnormalized `F` or `G` modes.
            %
            % - Topic: Developer topics — Diagnostic variables
            % - Declaration: values = rawVariable(basisSet,variable,z)
            % - Parameter variable: `"F"` or `"G"`
            % - Parameter z: physical coordinate
            % - Returns values: unnormalized aligned mode values
            % - Developer: true
            arguments
                self IMGeostrophicModesBasis
                variable {mustBeTextScalar, mustBeMember(variable,["F","G"])}
                z (:,1) double {mustBeReal, mustBeFinite}
            end
            if string(variable) == "F"
                values = self.solver.evaluateNativeModes(self.nativeModes,z);
            else
                values = self.gSolver.evaluateNativeModes(self.nativeGModes,z);
            end
        end
    end
    methods (Access = protected)
        function points = referenceBreakpoints(self,zBounds)
            % Stored F and G can have different staggered reconstruction knots.
            points = referenceBreakpoints@IMBasisSet(self,zBounds);
            companion = self.gSolver.zNative(:);
            points = unique([points;companion(companion > zBounds(1) & companion < zBounds(2))]);
        end

        function values = rawVariableDerivative(self,variable,z,order)
            if variable == "G" && ismember(order,[1 2])
                values = self.gSolver.evaluatePhysicalDerivative(self.nativeGModes,z,order);
            else
                values = rawVariableDerivative@IMInternalModesBasis(self,variable,z,order);
            end
        end

        function self = selectAlignedColumns(self,columns)
            self = selectAlignedColumns@IMInternalModesBasis(self,columns);
            self.nativeGModes = self.nativeGModes(:,columns);
            self.inverseDepth = self.inverseDepth(columns);
        end
    end

end
