classdef (Sealed) WVCompiledConstantStratificationBackend < handle
    % Own one compiled constant-stratification nonlinear-flux kernel.
    %
    % This developer-facing adapter is the MATLAB ownership boundary for
    % the source-only compiled preview. It validates an already-built native
    % provider, creates exactly one MEX kernel, and returns MATLAB-owned
    % `[Nj,Nkl]` flux arrays. It never downloads, builds, or falls back.
    %
    % Construct this object through `create`. Callers must delete it or its
    % owning transform to release all FFT plans and the MEX module lock.
    %
    % - Developer: true
    % - Topic: Compiled preview internals
    % - Declaration: classdef (Sealed) WVCompiledConstantStratificationBackend < handle

    properties (GetAccess=public, SetAccess=private)
        % Capability record validated during construction.
        %
        % - Developer: true
        % - Topic: Compiled preview internals
        capabilities
    end

    properties (Access=private)
        moduleName (1,1) string = ""
        kernelHandle = []
        configuration
        storageEstimate
    end

    methods (Static)
        function backend = create(wvt)
            % Create a compiled adapter for one constant-stratification transform.
            %
            % The native provider must already have been built explicitly by
            % calling `WVCompiledBackend.build()`.
            %
            % - Developer: true
            % - Topic: Compiled preview internals
            % - Declaration: backend = WVCompiledConstantStratificationBackend.create(wvt)
            % - Parameter wvt: constant-stratification transform defining the immutable kernel configuration
            % - Returns backend: validated compiled-kernel owner
            arguments
                wvt (1,1) WVTransformConstantStratification
            end
            capabilities = WVCompiledBackend.capabilities();
            WVCompiledConstantStratificationBackend.validateCapabilities(capabilities,wvt.isHydrostatic);
            backend = WVCompiledConstantStratificationBackend(wvt,capabilities);
        end
    end

    methods
        function [Fp,Fm,F0] = nonlinearFlux(self,Ap,Am,A0,t,t0)
            % Evaluate ordinary nonlinear advection in the compiled kernel.
            %
            % - Developer: true
            % - Topic: Compiled preview internals
            % - Declaration: [Fp,Fm,F0] = nonlinearFlux(Ap,Am,A0,t,t0)
            % - Parameter Ap: positive-frequency `[Nj,Nkl]` coefficients
            % - Parameter Am: negative-frequency `[Nj,Nkl]` coefficients
            % - Parameter A0: zero-frequency `[Nj,Nkl]` coefficients
            % - Parameter t: evaluation time in seconds
            % - Parameter t0: coefficient reference time in seconds
            % - Returns Fp: positive-frequency `[Nj,Nkl]` flux
            % - Returns Fm: negative-frequency `[Nj,Nkl]` flux
            % - Returns F0: zero-frequency `[Nj,Nkl]` flux
            if isempty(self.kernelHandle)
                error("WaveVortexModel:CompiledBackendDeleted","The compiled backend has already been deleted.")
            end
            [Fp,Fm,F0] = feval(char(self.moduleName),'nonlinearFlux',self.kernelHandle,Ap,Am,A0,t,t0);
        end

        function value = metadata(self)
            % Return JSON-safe identity, storage, and runtime metadata.
            %
            % - Developer: true
            % - Topic: Compiled preview internals
            % - Declaration: value = metadata()
            % - Returns value: immutable provider identity plus current kernel metrics
            if isempty(self.kernelHandle)
                runtimeMetrics = struct("status","deleted");
            else
                runtimeMetrics = feval(char(self.moduleName),'metrics',self.kernelHandle);
                runtimeMetrics.status = "active";
            end
            value = struct( ...
                "schemaVersion","1.0.0", ...
                "requestedBackend","compiled", ...
                "activeBackend","compiled", ...
                "scope","constant stratification with exactly the default WVNonlinearAdvection forcing", ...
                "provider",self.capabilities.provider, ...
                "libraries",self.capabilities.libraries, ...
                "module",self.capabilities.module, ...
                "contract",self.capabilities.contract, ...
                "storage",self.storageEstimate, ...
                "runtimeMetrics",runtimeMetrics);
        end

        function delete(self)
            if isempty(self.kernelHandle)
                return
            end
            handle = self.kernelHandle;
            self.kernelHandle = [];
            try
                feval(char(self.moduleName),'delete',handle);
            catch exception
                warning("WaveVortexModel:CompiledBackendCleanup","Unable to delete the compiled backend cleanly: %s",exception.message)
            end
        end
    end

    methods (Access=private)
        function self = WVCompiledConstantStratificationBackend(wvt,capabilities)
            self.capabilities = capabilities;
            self.moduleName = string(capabilities.module.name);
            self.configuration = WVCompiledConstantStratificationBackend.configurationForTransform(wvt);
            self.storageEstimate = feval(char(self.moduleName),'estimate',self.configuration);
            handle = [];
            try
                handle = feval(char(self.moduleName),'create',self.configuration,capabilities.contract.threadCount);
                self.kernelHandle = handle;
            catch exception
                if ~isempty(handle)
                    try
                        feval(char(self.moduleName),'delete',handle);
                    catch
                    end
                end
                rethrow(exception)
            end
        end
    end

    methods (Static, Access=private)
        function validateCapabilities(capabilities,isHydrostatic)
            if ~isstruct(capabilities) || ~isfield(capabilities,"schemaVersion") || string(capabilities.schemaVersion) ~= "1.0.0"
                error("WaveVortexModel:CompiledBackendCapabilitySchema","The compiled preview requires WVCompiledBackend capability schema 1.0.0.")
            end
            if ~capabilities.isAvailable
                message = "The compiled preview is unavailable. Call WVCompiledBackend.capabilities() for details and WVCompiledBackend.build() to build native support explicitly.";
                if isfield(capabilities,"failure") && isfield(capabilities.failure,"message") && strlength(string(capabilities.failure.message)) > 0
                    message = message + " " + string(capabilities.failure.message);
                end
                error("WaveVortexModel:CompiledBackendUnavailable","%s",message)
            end
            required = isfield(capabilities,"provider") && string(capabilities.provider.id) == "native-neon-pthreads" && ...
                isfield(capabilities,"module") && capabilities.module.identityValidated && ...
                isfield(capabilities,"contract") && capabilities.contract.version == 4 && capabilities.contract.planCount == 17 && ...
                isfield(capabilities,"libraries") && ~capabilities.libraries.openmp.detected && ...
                isfield(capabilities,"featureValidation") && string(capabilities.featureValidation.status) == "passed" && ...
                capabilities.featureValidation.maximumRelativeError <= 1e-12;
            if ~required
                error("WaveVortexModel:CompiledBackendCapabilityMismatch","The installed compiled backend does not satisfy the native provider, library identity, contract, plan-count, OpenMP, and numerical-validation requirements.")
            end
            feature = capabilities.featureValidation.nonhydrostatic;
            if isHydrostatic
                feature = capabilities.featureValidation.hydrostatic;
            end
            if string(feature.status) ~= "passed" || feature.maximumRelativeError > 1e-12
                error("WaveVortexModel:CompiledBackendSelfTest","The requested physical configuration did not pass the compiled-backend self-test.")
            end
        end

        function configuration = configurationForTransform(wvt)
            configuration = struct( ...
                "Nx",wvt.Nx,"Ny",wvt.Ny,"Nz",wvt.Nz,"Nj",wvt.Nj, ...
                "Lx",wvt.Lx,"Ly",wvt.Ly,"Lz",wvt.Lz, ...
                "N0",wvt.N0,"rho0",wvt.rho0,"g",wvt.g, ...
                "planetaryRadius",wvt.planetaryRadius,"rotationRate",wvt.rotationRate,"latitude",wvt.latitude, ...
                "isHydrostatic",wvt.isHydrostatic,"shouldAntialias",wvt.shouldAntialias);
        end
    end
end
