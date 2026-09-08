classdef IMMeanDensityAnomalyModesBasis < IMInternalModesBasis
    % Store aligned mean-density-anomaly `F` and `G` modes.
    %
    % The EVP solves `G`. This basis stores the companion pressure shapes
    % obtained by surface-referenced integration of the native spectral or
    % finite-difference representation:
    %
    % $$
    % F_j(z)=\frac{1}{g}\int_z^{z_s}N^2(z')G_j(z')\,dz'.
    % $$
    %
    % The default `unity` normalization divides each aligned `F`/`G` pair
    % by the positive magnitude of its signed generalized-energy norm.
    % `signatures(j)` records the remaining sign, so the normalized Gram
    % matrix is `diag(signatures)`.
    %
    % Discrete transforms project `G` and synthesize both `G` and the
    % diagnostic `F`. Omitted transform variables therefore select `G`.
    %
    % - Topic: Create internal-mode bases
    % - Topic: Evaluate modes
    % - Topic: Analyze modes
    % - Topic: Build discrete transforms
    % - Declaration: classdef IMMeanDensityAnomalyModesBasis < IMInternalModesBasis

    properties (SetAccess = private)
        % Surface generalized-energy acceleration.
        %
        % - Topic: Inspect basis sets
        g0

        % Bottom generalized-energy acceleration.
        %
        % - Topic: Inspect basis sets
        gd

        % Active endpoint identities in canonical surface-bottom order.
        %
        % - Topic: Inspect basis sets
        activeEndpoints

        % Signs of the normalized generalized-energy mode norms.
        %
        % Each value is `-1` or `+1`. A zero signed norm is rejected when
        % the basis is constructed.
        %
        % - Topic: Analyze modes
        signatures
    end

    properties (Access = private)
        nativeFModes
    end

    methods
        function self = IMMeanDensityAnomalyModesBasis(options)
            % Create an aligned mean-density-anomaly basis set.
            %
            % - Topic: Create internal-mode bases
            % - Declaration: basisSet = IMMeanDensityAnomalyModesBasis(options)
            % - Parameter options.solver: configured solver
            % - Parameter options.evp: mean-density-anomaly EVP descriptor
            % - Parameter options.nativeModes: native solved `G` columns
            % - Parameter options.eigenvalues: retained eigenvalues
            % - Parameter options.modeNumber: physical mode labels
            % - Parameter options.modeSelectionDiagnostics: selection diagnostics
            % - Parameter options.metadata: captured problem metadata
            % - Returns basisSet: aligned mean-density-anomaly basis set
            arguments
                options.solver IMSolver
                options.evp IMMeanDensityAnomalyModes
                options.nativeModes (:,:) double
                options.eigenvalues (1,:) double {mustBeReal, mustBeFinite}
                options.modeNumber (1,:) double {mustBeInteger}
                options.modeSelectionDiagnostics struct = struct()
                options.metadata struct = struct()
            end

            self@IMInternalModesBasis(solver=options.solver,evp=options.evp,nativeModes=options.nativeModes, ...
                eigenvalues=options.eigenvalues,modeNumber=options.modeNumber, ...
                modeSelectionDiagnostics=options.modeSelectionDiagnostics,metadata=options.metadata);
            self.g0 = options.evp.g0;
            self.gd = options.evp.gd;
            self.activeEndpoints = options.evp.activeEndpoints;

            zNative = options.solver.zNative;
            N2Values = options.evp.N2(zNative);
            N2Values = N2Values(:);
            if numel(N2Values) ~= numel(zNative) || any(~isfinite(N2Values)) || any(N2Values <= 0)
                error("IMMeanDensityAnomalyModesBasis:InvalidStratification", "N2 must return one finite positive value for each solver grid point.");
            end
            GNativeValues = options.solver.evaluateNativeModes(options.nativeModes,zNative);
            self.nativeFModes = options.solver.integrateGridValuesFromSurface((N2Values/options.evp.g).*GNativeValues);

            rawGram = self.variableGramMatrix("G",self.zDomain,false);
            rawNorms = real(diag(rawGram)).';
            normScale = max(1,norm(rawGram,2));
            zeroNorm = abs(rawNorms) <= 1e3*eps(normScale);
            if any(zeroNorm)
                labels = self.modeNumber(zeroNorm);
                error("IMMeanDensityAnomalyModesBasis:ZeroNormMode", ...
                    "Mean-density-anomaly mode label(s) %s have zero signed generalized-energy norm. Adjust g0 or gd away from the singular endpoint balance or change the retained band.", ...
                    join(string(labels),", "));
            end
            self.signatures = sign(rawNorms);
            self.normalization = "unity";
        end

        function Gz = Gz(self,z,options)
            % Evaluate vertical derivatives of the solved `G` modes.
            %
            % - Topic: Evaluate modes
            % - Declaration: Gz = Gz(basisSet,z,options)
            % - Parameter z: physical coordinate
            % - Parameter options.normalization: normalization rule name
            % - Returns Gz: evaluated `G` derivatives
            arguments
                self IMMeanDensityAnomalyModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
                options.normalization = self.normalization
            end
            Gz = self.uz(z,normalization=options.normalization);
        end

        function self = orientModeSigns(self)
            % Apply one deterministic sign to each aligned `F`/`G` pair.
            %
            % - Topic: Developer topics — Diagnostic variables
            % - Declaration: basisSet = orientModeSigns(basisSet)
            % - Returns basisSet: basis set with aligned signs oriented
            % - Developer: true
            [self,signs] = orientModeSigns@IMInternalModesBasis(self);
            self.nativeFModes = self.nativeFModes.*signs;
        end
    end

    methods
        function values = rawVariable(self,variable,z)
            % Evaluate unnormalized aligned `F` or `G` modes.
            %
            % - Topic: Developer topics — Diagnostic variables
            % - Declaration: values = rawVariable(basisSet,variable,z)
            % - Parameter variable: `"F"` or `"G"`
            % - Parameter z: physical coordinate
            % - Returns values: unnormalized aligned mode values
            % - Developer: true
            arguments
                self IMMeanDensityAnomalyModesBasis
                variable {mustBeTextScalar, mustBeMember(variable,["F","G"])}
                z (:,1) double {mustBeReal, mustBeFinite}
            end
            if string(variable) == "F"
                values = self.solver.evaluateNativeModes(self.nativeFModes,z);
            else
                values = self.solver.evaluateNativeModes(self.nativeModes,z);
            end
        end
    end
end
