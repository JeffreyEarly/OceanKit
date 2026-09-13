classdef (Sealed) WVCompiledTransformBackend < handle
    % Own one compiled variable-evaluation transform session.
    %
    % This developer-facing adapter provides the in-memory bridge for all
    % six built-in v4 configurations. It owns only the native handle
    % and immutable shape/module identity; the native session owns modal data,
    % prepared operators, and evaluation workspaces.
    %
    % - Developer: true
    % - Topic: Compiled transform internals
    % - Declaration: classdef (Sealed) WVCompiledTransformBackend

    properties (GetAccess=public, SetAccess=private)
        capabilities
    end

    properties (Access=private)
        moduleName (1,1) string = ""
        transformHandle = []
        transformClass (1,1) string = ""
        Nx (1,1) double = 0
        Ny (1,1) double = 0
        Nz (1,1) double = 0
        Nj (1,1) double = 0
        Nkl (1,1) double = 0
        compiledSourceIdentity = []
        hasSourceIdentityLease (1,1) logical = false
        geometryValues = []
        evaluationToken (1,1) uint64 = uint64(0)
        leasedEvaluationToken (1,1) uint64 = uint64(0)
        evaluationGeneration (1,1) uint64 = uint64(0)
        supportedVariableNames = strings(0,1)
    end

    methods (Static)
        function backend = create(wvt)
            % Create a compiled transform session.
            %
            % - Developer: true
            % - Topic: Compiled transform internals
            % - Declaration: backend = WVCompiledTransformBackend.create(wvt)
            % - Parameter wvt: transform supplying its geometry and solved modal source
            % - Returns backend: owning compiled transform session
            arguments
                wvt (1,1) WVTransform
            end
            if ~ismember(string(class(wvt)),["WVTransformHydrostatic" "WVTransformBoussinesq" "WVTransformStratifiedQG" "WVTransformBarotropicQG" "WVTransformConstantStratification"])
                error("WaveVortexModel:CompiledTransformUnsupportedFamily","Compiled MATLAB support covers the five built-in transform families, including both constant-stratification configurations.")
            end
            capabilities = WVCompiledBackend.capabilities();
            WVCompiledTransformBackend.validateCapabilities(capabilities);
            backend = WVCompiledTransformBackend(wvt,capabilities);
        end
    end

    methods
        function cleanup = scopedEvaluation(self,wvt)
            % Open or join one immutable state evaluation.
            %
            % Keep the returned onCleanup object alive for all queries. A new
            % native-owned coefficient snapshot is copied once per outer scope.
            % Mutating the MATLAB state invalidates further queries in that scope.
            %
            % - Developer: true
            % - Topic: Compiled transform internals
            arguments
                self (1,1) WVCompiledTransformBackend
                wvt (1,1) WVTransform
            end
            self.assertActive(); self.assertMatchingTransform(wvt);
            self.assertScopeUsable(wvt);
            if self.evaluationToken ~= 0
                self.assertEvaluationGeneration(wvt);
                cleanup = onCleanup(@()[]);
                return
            end
            [Ap,Am,A0] = WVCompiledTransformBackend.stateArrays(wvt,self.Nj,self.Nkl);
            token = feval(char(self.moduleName),'transformBeginEvaluation',self.transformHandle,Ap,Am,A0,wvt.t,wvt.t0);
            self.evaluationToken = token;
            self.leasedEvaluationToken = token;
            self.evaluationGeneration = wvt.compiledStateGeneration;
            cleanup = onCleanup(@()WVCompiledTransformBackend.finishScopeIfValid(self,token));
        end

        function result = evaluateVariables(self,wvt,variableNames)
            % Query dependencies in an explicit scope, or use a fresh scope.
            %
            % - Developer: true
            % - Topic: Compiled transform internals
            arguments
                self (1,1) WVCompiledTransformBackend
                wvt (1,1) WVTransform
                variableNames
            end
            cleanup = self.scopedEvaluation(wvt); %#ok<NASGU>
            variableNames = WVCompiledTransformBackend.normalizeVariableNames(variableNames);
            actualDensity = ~isprop(wvt,'shouldUseTrueNoMotionProfile') || wvt.shouldUseTrueNoMotionProfile;
            result = feval(char(self.moduleName),'transformQueryEvaluation',self.transformHandle,self.evaluationToken,variableNames,logical(actualDensity));
        end

        function flag = supportsVariables(self,names)
            % - Developer: true
            % - Topic: Compiled transform internals
            flag = all(ismember(string(names),self.supportedVariableNames));
        end

        function result = nonlinearCoefficients(self,wvt)
            % Return coefficients while retaining shared native dependencies.
            %
            % - Developer: true
            % - Topic: Compiled transform internals
            arguments
                self (1,1) WVCompiledTransformBackend
                wvt (1,1) WVTransform
            end
            cleanup = self.scopedEvaluation(wvt); %#ok<NASGU>
            result = feval(char(self.moduleName),'transformNonlinearCoefficients',self.transformHandle,self.evaluationToken);
        end

        function invalidateEvaluation(self)
            % Close an event before its owning MATLAB state is changed.
            %
            % - Developer: true
            % - Topic: Compiled transform internals
            self.endEvaluation(self.evaluationToken);
        end

        function flag = canBeginCoefficientOnlyEvaluation(self)
            % Preserve active and invalidated external evaluation leases.
            %
            % - Developer: true
            % - Topic: Compiled transform internals
            flag = self.evaluationToken == 0 && self.leasedEvaluationToken == 0;
        end

        function result = coefficientOnlyRightHandSide(self,wvt)
            % Execute a sealed native RHS while borrowing inputs for this call.
            %
            % - Developer: true
            % - Topic: Compiled transform internals
            arguments
                self (1,1) WVCompiledTransformBackend
                wvt (1,1) WVTransform
            end
            [used,result] = self.tryCoefficientOnlyRightHandSide(wvt);
            if ~used
                error("WaveVortexModel:CompiledCoefficientOnlyWorkload","A coefficient-only RHS requires the default nonlinear forcing and no active evaluation lease.");
            end
        end

        function [used,result] = tryCoefficientOnlyRightHandSide(self,wvt)
            % Execute the sealed native RHS when the receiver and workload allow it.
            %
            % - Developer: true
            % - Topic: Compiled transform internals
            arguments
                self (1,1) WVCompiledTransformBackend
                wvt (1,1) WVTransform
            end
            self.assertActive(); self.assertMatchingTransform(wvt);
            self.assertScopeUsable(wvt);
            used = false;
            result = [];
            if ~self.canBeginCoefficientOnlyEvaluation() || ~wvt.canUseCompiledCoefficientOnlyRightHandSide()
                return
            end
            [Ap,Am,A0] = WVCompiledTransformBackend.stateArrays(wvt,self.Nj,self.Nkl);
            result = feval(char(self.moduleName),'transformCoefficientOnlyRightHandSide',self.transformHandle,Ap,Am,A0,wvt.t,wvt.t0);
            used = true;
        end

        function prepare(self,variableNames)
            % Prepare one immutable variable evaluation plan.
            %
            % - Developer: true
            % - Topic: Compiled transform internals
            % - Declaration: prepare(variableNames)
            % - Parameter variableNames: names accepted by the portable variable catalog
            self.assertActive();
            variableNames = WVCompiledTransformBackend.normalizeVariableNames(variableNames);
            feval(char(self.moduleName),'transformPrepare',self.transformHandle,variableNames);
        end

        function result = evaluate(self,wvt,variableNames,options)
            % Evaluate registered variables and optionally the nonlinear flux.
            %
            % - Developer: true
            % - Topic: Compiled transform internals
            % - Declaration: result = evaluate(wvt,variableNames,options)
            % - Parameter wvt: source transform whose state is evaluated
            % - Parameter variableNames: prepared variable names
            % - Parameter options.shouldEvaluateFlux: also return three wave flux cells or one QG flux cell; default false
            % - Returns result: scalar native result struct with values, flux, and metrics
            arguments
                self (1,1) WVCompiledTransformBackend
                wvt (1,1) WVTransform
                variableNames
                options.shouldEvaluateFlux (1,1) logical = false
            end
            self.assertActive();
            self.assertMatchingTransform(wvt);
            variableNames = WVCompiledTransformBackend.normalizeVariableNames(variableNames);
            [Ap,Am,A0] = WVCompiledTransformBackend.stateArrays(wvt,self.Nj,self.Nkl);
            result = feval(char(self.moduleName),'transformEvaluate',self.transformHandle,Ap,Am,A0,wvt.t,wvt.t0,variableNames,options.shouldEvaluateFlux);
        end

        function result = operation(self,wvt,operationName,inputs,options)
            % Invoke one explicitly supported compiled numerical primitive.
            %
            % - Developer: true
            % - Topic: Compiled transform internals
            % - Declaration: result = operation(wvt,operationName,inputs,options)
            % - Parameter wvt: source transform identifying this session
            % - Parameter operationName: native primitive identifier
            % - Parameter inputs: ordered input arrays borrowed for this call
            % - Parameter options: primitive-specific options
            % - Returns result: ordered values and native execution metrics
            arguments
                self (1,1) WVCompiledTransformBackend
                wvt (1,1) WVTransform
                operationName (1,1) string
                inputs (1,:) cell = cell(1,0)
                options (1,1) struct = struct()
            end
            self.assertActive();
            self.assertMatchingTransform(wvt);
            if ~ismember(string(class(wvt)),["WVTransformHydrostatic" "WVTransformBoussinesq" "WVTransformStratifiedQG" "WVTransformBarotropicQG" "WVTransformConstantStratification"])
                error("WaveVortexModel:CompiledTransformUnsupportedFamily","Compiled operation family is not supported.")
            end
            self.assertScopeUsable(wvt);
            result = feval(char(self.moduleName),'transformOperation',self.transformHandle,char(operationName),inputs,options);
        end

        function value = metadata(self)
            % Return native transform-session metadata.
            %
            % - Developer: true
            % - Topic: Compiled transform internals
            % - Declaration: value = metadata()
            % - Returns value: scalar native metadata struct, or deleted status
            if isempty(self.transformHandle)
                value = struct("status","deleted");
            else
                value = feval(char(self.moduleName),'transformMetadata',self.transformHandle);
            end
        end

        function delete(self)
            if isempty(self.transformHandle) && ~self.hasSourceIdentityLease
                return
            end
            handle = self.transformHandle;
            self.transformHandle = [];
            self.evaluationToken = uint64(0);
            self.leasedEvaluationToken = uint64(0);
            if ~isempty(handle)
                try
                    feval(char(self.moduleName),'transformDelete',handle);
                catch exception
                    warning("WaveVortexModel:CompiledTransformCleanup","Unable to delete the compiled transform cleanly: %s",exception.message)
                end
            end
            self.releaseSourceIdentityLease();
        end
    end

    methods (Access=private)
        function assertScopeUsable(self,wvt)
            if self.leasedEvaluationToken ~= 0 && self.evaluationToken == 0
                error("WaveVortexModel:CompiledTransformStateChanged","The enclosing evaluation was invalidated. Release its scope before starting a new evaluation.")
            end
            if self.evaluationToken ~= 0
                self.assertEvaluationGeneration(wvt);
            end
        end

        function finishEvaluationScope(self,token)
            if token ~= self.leasedEvaluationToken, return, end
            self.endEvaluation(token);
            self.leasedEvaluationToken = uint64(0);
        end

        function assertEvaluationGeneration(self,wvt)
            if self.evaluationGeneration ~= wvt.compiledStateGeneration
                self.invalidateEvaluation();
                error("WaveVortexModel:CompiledTransformStateChanged","The MATLAB state changed during the compiled evaluation. Start a new scope.")
            end
        end

        function endEvaluation(self,token)
            if token == 0 || isempty(self.transformHandle) || token ~= self.evaluationToken
                return
            end
            feval(char(self.moduleName),'transformEndEvaluation',self.transformHandle,token);
            self.evaluationToken = uint64(0);
        end

        function self = WVCompiledTransformBackend(wvt,capabilities)
            self.capabilities = capabilities;
            self.moduleName = WVCompiledBackend.activateModule(capabilities);
            self.transformClass = string(class(wvt));
            self.Nx = wvt.Nx; self.Ny = wvt.Ny; self.Nj = wvt.Nj; self.Nkl = wvt.Nkl;
            self.geometryValues = WVCompiledTransformBackend.sourceGeometry(wvt);
            self.Nz = self.geometryValues(3);
            self.compiledSourceIdentity = wvt.compiledSourceIdentity;
            configuration = wvCompiledStratifiedModalConfiguration(wvt);
            handle = [];
            try
                handle = feval(char(self.moduleName),'transformCreate',configuration,capabilities.contract.threadCount);
                self.transformHandle = handle;
                self.supportedVariableNames = string(feval(char(self.moduleName),'transformVariables',handle));
                self.compiledSourceIdentity.acquireLease();
                self.hasSourceIdentityLease = true;
            catch exception
                self.transformHandle = [];
                if ~isempty(handle)
                    try
                        feval(char(self.moduleName),'transformDelete',handle);
                    catch
                    end
                end
                rethrow(exception)
            end
        end

        function assertActive(self)
            if isempty(self.transformHandle)
                error("WaveVortexModel:CompiledTransformDeleted","The compiled transform session has already been deleted.")
            end
        end

        function releaseSourceIdentityLease(self)
            if ~self.hasSourceIdentityLease, return, end
            self.hasSourceIdentityLease = false;
            identity = self.compiledSourceIdentity;
            if ~isempty(identity) && isvalid(identity)
                identity.releaseLease();
            end
        end

        function assertMatchingTransform(self,wvt)
            if string(class(wvt)) ~= self.transformClass || ~isvalid(self.compiledSourceIdentity) || ~(wvt.compiledSourceIdentity == self.compiledSourceIdentity) || wvt.Nx ~= self.Nx || wvt.Ny ~= self.Ny || ~isequal(WVCompiledTransformBackend.sourceGeometry(wvt),self.geometryValues) || wvt.Nj ~= self.Nj || wvt.Nkl ~= self.Nkl
                error("WaveVortexModel:CompiledTransformMismatch","The evaluated transform does not match the compiled session configuration.")
            end
        end
    end

    methods (Static, Access=private)
        function finishScopeIfValid(backend,token)
            if isvalid(backend)
                backend.finishEvaluationScope(token);
            end
        end

        function values = sourceGeometry(wvt)
            if isa(wvt,"WVTransformBarotropicQG"), nz = 1; else, nz = wvt.Nz; end
            values = [wvt.Nx wvt.Ny nz wvt.Nj wvt.Nkl wvt.Lx wvt.Ly wvt.Lz wvt.g wvt.latitude wvt.rotationRate wvt.planetaryRadius double(wvt.shouldAntialias)];
        end

        function validateCapabilities(capabilities)
            if ~isstruct(capabilities) || ~isfield(capabilities,"isAvailable") || ~capabilities.isAvailable
                error("WaveVortexModel:CompiledTransformUnavailable","The compiled transform backend is unavailable. Build native support explicitly before creating a session.")
            end
            if ~isfield(capabilities,"module") || ~isfield(capabilities.module,"name") || strlength(string(capabilities.module.name)) == 0 || ...
                    ~isfield(capabilities.module,"identityValidated") || ~capabilities.module.identityValidated
                error("WaveVortexModel:CompiledTransformCapabilityMismatch","The installed compiled module does not have a validated identity.")
            end
            if ~isfield(capabilities.module,"matlabTransformBridgeVersion") || capabilities.module.matlabTransformBridgeVersion < 9
                error("WaveVortexModel:CompiledTransformCapabilityMismatch","The installed module predates the MATLAB transform operation bridge. Rebuild with WVCompiledBackend.build().")
            end
            if ~isfield(capabilities,"contract") || ~isfield(capabilities.contract,"threadCount") || ...
                    ~isnumeric(capabilities.contract.threadCount) || ~isscalar(capabilities.contract.threadCount) || ...
                    ~isfinite(capabilities.contract.threadCount) || capabilities.contract.threadCount < 1 || ...
                    capabilities.contract.threadCount ~= floor(capabilities.contract.threadCount)
                error("WaveVortexModel:CompiledTransformCapabilityMismatch","The compiled capability record has no validated positive FFT thread count.")
            end
        end

        function variableNames = normalizeVariableNames(variableNames)
            if isstring(variableNames)
                variableNames = cellstr(variableNames(:));
            elseif ischar(variableNames)
                variableNames = {variableNames};
            elseif iscell(variableNames) && all(cellfun(@(name) ischar(name) && isrow(name),variableNames(:)))
                variableNames = variableNames(:);
            else
                error("WaveVortexModel:CompiledTransformVariables","Variable names must be a string array, character vector, or cell array of character vectors.")
            end
            if isempty(variableNames)
                error("WaveVortexModel:CompiledTransformVariables","At least one variable name is required.")
            end
        end

        function [Ap,Am,A0] = stateArrays(wvt,Nj,Nkl)
            empty = complex(zeros(0,0));
            if ismember(string(class(wvt)),["WVTransformHydrostatic" "WVTransformBoussinesq" "WVTransformConstantStratification"])
                Ap = complex(wvt.Ap); Am = complex(wvt.Am); A0 = complex(wvt.A0);
            else
                Ap = empty; Am = empty; A0 = complex(wvt.A0);
            end
            if ~isequal(size(A0),[Nj Nkl]) || (~isempty(Ap) && ~isequal(size(Ap),[Nj Nkl])) || (~isempty(Am) && ~isequal(size(Am),[Nj Nkl]))
                error("WaveVortexModel:CompiledTransformShape","The evaluated transform coefficient arrays do not match the compiled session shape.")
            end
        end
    end
end
