classdef IMConstantStratificationSolution < IMAnalyticalSolution
    % Analytical solution family for constant stratification.
    %
    % `IMConstantStratificationSolution` owns the closed-form formulas for
    % $$N^2(z)=N_0^2$$. It can create exact internal-mode bases for the
    % supported canonical internal-mode EVPs and exact geostrophic zero-APV
    % boundary-response modes.
    %
    % The public `IMInternalModes.geostrophicAPVModes` descriptor selects
    % an exact generalized-energy APV catalog. Positive eigendepths use
    % trigonometric modes, a zero eigenvalue uses the affine solution, and
    % negative eigendepths use hyperbolic modes. Endpoint inertia determines
    % whether zero, one, or two negative modes precede the optional exact
    % zero mode and positive branch. All finite, zero, and positive-infinite
    % `g0` and `gd` limits are supported under both surface conventions.
    % APV bases use the volume-only `depth` normalization by default and
    % expose roots, residuals, branch labels, endpoint inertia, `g0`, `gd`,
    % and `surfaceBoundary` in their metadata.
    %
    % Exact geostrophic zero-APV modes use scaled hyperbolic functions with
    % $$m=kN_0/|f_0|$$. A two-column response solve produces unit surface
    % and bottom coordinates under either the free-surface response `G-F`
    % or the rigid-lid response `G`. The older one-boundary profile formula
    % used unit `f0*Fz` at its active boundary; the public exact basis instead
    % uses the canonical response normalization shared with the numerical
    % zero-APV basis.
    %
    % ```matlab
    % solution = IMConstantStratificationSolution(N0=5.2e-3,zDomain=[-5000 0]);
    % evp = IMInternalModes.hydrostaticGModes(N2=@(z) solution.N2(z), zDomain=solution.zDomain);
    % basisSet = solution.internalModes(evp,nModes=4);
    % exactModes = solution.geostrophicZeroAPVModesAtWavenumber(1e-4);
    % ```
    %
    % - Topic: Create analytical solutions
    % - Topic: Compute internal modes
    % - Topic: Compute geostrophic zero-APV modes
    % - Topic: Inspect analytical solutions
    % - Topic: Developer topics
    % - Declaration: classdef IMConstantStratificationSolution < IMAnalyticalSolution

    properties (SetAccess = private)
        % Constant buoyancy frequency $$N_0$$ in radians per second.
        %
        % - Topic: Inspect analytical solutions
        N0
    end

    methods
        function self = IMConstantStratificationSolution(options)
            % Create a constant-stratification analytical solution family.
            %
            % - Topic: Create analytical solutions
            % - Declaration: solution = IMConstantStratificationSolution(options)
            % - Parameter options.N0: constant buoyancy frequency
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Returns solution: analytical solution family
            arguments
                options.N0 (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 5.2e-3
                options.zDomain (1,2) double {mustBeReal, mustBeFinite} = [-1 0]
                options.f0 (1,1) double {mustBeReal, mustBeFinite} = 0
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
            end

            self@IMAnalyticalSolution(zDomain=options.zDomain, f0=options.f0, g=options.g);
            self.N0 = options.N0;
        end

        function values = N2(self, z)
            % Evaluate $$N^2(z)=N_0^2$$.
            %
            % - Topic: Inspect analytical solutions
            % - Declaration: values = N2(solution,z)
            % - Parameter z: physical coordinate
            % - Returns values: buoyancy frequency squared
            arguments
                self IMConstantStratificationSolution
                z double {mustBeReal, mustBeFinite}
            end

            values = self.N0*self.N0*ones(size(z));
        end

        function basisSet = internalModes(self, evp, options)
            % Create an exact internal-mode basis.
            %
            % Generalized-energy APV modes are returned in ascending
            % eigenvalue order $$1/h$$. Hyperbolic negative branches come
            % first, followed by an optional affine zero branch with
            % `h=Inf`, then trigonometric positive branches. The same
            % positive depth-normalization factor scales exact `F` and `G`.
            %
            % - Topic: Compute internal modes
            % - Declaration: basisSet = internalModes(solution,evp,options)
            % - Parameter evp: internal-mode EVP
            % - Parameter options.nModes: number of retained modes
            % - Parameter options.normalization: active normalization
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: exact analytical internal-mode basis
            arguments
                self IMConstantStratificationSolution
                evp = []
                options.nModes (1,1) double {mustBeInteger, mustBePositive} = 64
                options.normalization = []
                options.metadata struct = struct()
            end

            evp = self.resolveEVP(evp);
            if IMConstantGeostrophicAPVCatalog.isTargetEVP(evp)
                [modeData, h, modeNumber, analyticalMetadata, diagnostics] = IMConstantGeostrophicAPVCatalog.solve(evp, self.N0, self.zDomain, options.nModes);
                metadata = options.metadata;
                metadataFields = fieldnames(analyticalMetadata);
                for iField = 1:length(metadataFields)
                    fieldName = metadataFields{iField};
                    metadata.(fieldName) = analyticalMetadata.(fieldName);
                end
                metadata.analyticalSolution = "constantStratification";
                metadata.analyticalFamily = "generalizedEnergyAPV";
                basisSet = IMAnalyticalInternalModesBasis(solution=self, evp=evp, h=h, modeNumber=modeNumber, N2=@(z) self.N2(z), rawVariableFunction=@(~,variable,z) IMConstantGeostrophicAPVCatalog.rawVariable(modeData, self.N0, evp.g, variable, z), rawUzFunction=@(~,z) IMConstantGeostrophicAPVCatalog.rawFz(modeData, self.N0, evp.g, z), normalization=options.normalization, metadata=metadata, modeSelectionDiagnostics=diagnostics);
                return;
            end
            f0 = evp.f0;
            gValue = evp.g;
            [h, verticalWavenumbers, solutionTypes, isBoundaryMode, baroclinicNumbers, modeNumber] = IMConstantStratificationSolution.solveSpectrum(evp, self.N0, self.zDomain, options.nModes, f0, gValue);
            modeData = struct();
            modeData.h = h;
            modeData.verticalWavenumbers = verticalWavenumbers;
            modeData.solutionTypes = solutionTypes;
            modeData.isBoundaryMode = isBoundaryMode;
            modeData.baroclinicNumbers = baroclinicNumbers;

            metadata = options.metadata;
            metadata.analyticalSolution = "constantStratification";
            metadata.verticalWavenumbers = verticalWavenumbers;
            metadata.solutionTypes = solutionTypes;
            metadata.isBoundaryMode = isBoundaryMode;
            metadata.baroclinicNumbers = baroclinicNumbers;

            basisSet = IMAnalyticalInternalModesBasis(solution=self, evp=evp, h=h, modeNumber=modeNumber, N2=@(z) self.N2(z), rawVariableFunction=@(basisSet,variable,z) self.rawVariable(modeData, evp, variable, z), rawUzFunction=@(basisSet,z) self.rawUz(modeData, evp, z), normalization=options.normalization, metadata=metadata);
        end

        function basisSet = geostrophicZeroAPVModesAtWavenumber(self,k,options)
            % Create exact canonical geostrophic zero-APV modes.
            %
            % The exact columns solve
            %
            % $$
            % \frac{f_0^2}{N_0^2}F_{zz}-k^2F=0,
            % \qquad G=-\frac{g}{N_0^2}F_z,
            % $$
            %
            % and have unit response at each requested endpoint and zero
            % response at the other endpoint.
            %
            % - Topic: Compute geostrophic zero-APV modes
            % - Declaration: exactModes = geostrophicZeroAPVModesAtWavenumber(solution,k,options)
            % - Parameter k: positive horizontal wavenumbers
            % - Parameter options.endpoints: requested surface and bottom coordinates
            % - Parameter options.surfaceBoundary: `"freeSurface"` or `"rigidLid"`
            % - Parameter options.metadata: additional metadata
            % - Returns exactModes: exact boundary-normalized basis
            arguments
                self IMConstantStratificationSolution
                k double {mustBeReal, mustBeFinite, mustBePositive}
                options.endpoints {mustBeText} = ["surface", "bottom"]
                options.surfaceBoundary {mustBeTextScalar} = "freeSurface"
                options.metadata struct = struct()
            end

            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=@(z) self.N2(z),zDomain=self.zDomain,f0=self.f0,g=self.g,k=k,endpoints=options.endpoints,surfaceBoundary=options.surfaceBoundary,metadata=options.metadata);
            modeData = self.zeroAPVModeData(problem);
            metadata = options.metadata;
            metadata.analyticalSolution = "constantStratification";
            metadata.analyticalFamily = "geostrophicZeroAPV";
            metadata.endpoints = problem.endpoints;
            metadata.surfaceBoundary = problem.surfaceBoundary;
            metadata.responseReciprocalCondition = modeData.responseReciprocalCondition;
            metadata.formulaFamily = "scaledHyperbolic";
            basisSet = IMAnalyticalGeostrophicZeroAPVModesBasis(solution=self,problem=problem,FFunction=@(z) self.zeroAPVVariable(modeData,problem,"F",z),GFunction=@(z) self.zeroAPVVariable(modeData,problem,"G",z),metadata=metadata);
        end

        function summarize(self)
            % Print a readable solution-family summary.
            %
            % - Topic: Inspect analytical solutions
            % - Declaration: summarize(solution)
            arguments
                self IMConstantStratificationSolution
            end

            summarize@IMAnalyticalSolution(self);
            fprintf("  N0: %g\n", self.N0);
            fprintf("  exact operations: internalModes, geostrophicZeroAPVModesAtWavenumber\n");
        end
    end

    methods (Access = private)
        function evp = resolveEVP(self, evp)
            if isempty(evp)
                evp = IMInternalModes.hydrostaticGModes(N2=@(z) self.N2(z), zDomain=self.zDomain, f0=self.f0, g=self.g);
                return;
            end
            if ~isa(evp, "IMInternalModes")
                error("IMConstantStratificationSolution:InvalidEVP", "The EVP must be an IMInternalModes instance.");
            end
            IMConstantStratificationSolution.validateEVPDomain(evp, self.zDomain);
        end

        function values = rawVariable(self, modeData, evp, variable, z)
            variable = string(variable);
            if variable ~= "G" && variable ~= "F"
                error("IMConstantStratificationSolution:UnsupportedVariable", "Variable must be ""F"" or ""G"".");
            end

            z = z(:);
            s = z - self.zDomain(1);
            values = zeros(length(z), length(modeData.h));
            for iMode = 1:length(modeData.h)
                k_z = modeData.verticalWavenumbers(iMode);
                hMode = modeData.h(iMode);
                if modeData.isBoundaryMode(iMode)
                    [G, F] = self.rawBoundaryMode(modeData.solutionTypes(iMode), k_z, hMode, s);
                elseif evp.formulation == "F"
                    [G, F] = self.rawHydrostaticFMode(k_z, modeData.baroclinicNumbers(iMode), s, evp.g);
                else
                    signValue = (-1)^modeData.baroclinicNumbers(iMode);
                    G = signValue*sin(k_z*s);
                    F = signValue*hMode*k_z*cos(k_z*s);
                end

                if variable == "G"
                    values(:,iMode) = G;
                else
                    values(:,iMode) = F;
                end
            end
        end

        function values = rawUz(self, modeData, evp, z)
            z = z(:);
            s = z - self.zDomain(1);
            values = zeros(length(z), length(modeData.h));
            for iMode = 1:length(modeData.h)
                k_z = modeData.verticalWavenumbers(iMode);
                hMode = modeData.h(iMode);
                if modeData.isBoundaryMode(iMode)
                    [dGdz, dFdz] = self.rawBoundaryModeDerivatives(modeData.solutionTypes(iMode), k_z, hMode, s);
                elseif evp.formulation == "F"
                    [dGdz, dFdz] = self.rawHydrostaticFModeDerivatives(k_z, modeData.baroclinicNumbers(iMode), s, evp.g);
                else
                    signValue = (-1)^modeData.baroclinicNumbers(iMode);
                    dGdz = signValue*k_z*cos(k_z*s);
                    dFdz = -signValue*hMode*k_z*k_z*sin(k_z*s);
                end

                if evp.formulation == "G"
                    values(:,iMode) = dGdz;
                else
                    values(:,iMode) = dFdz;
                end
            end
        end

        function [G, F] = rawHydrostaticFMode(self, k_z, modeNumber, s, gValue)
            if abs(k_z) <= eps
                F = ones(size(s));
                G = zeros(size(s));
                return;
            end

            signValue = (-1)^modeNumber;
            F = signValue*cos(k_z*s);
            G = signValue*(gValue*k_z/(self.N0*self.N0))*sin(k_z*s);
        end

        function [dGdz, dFdz] = rawHydrostaticFModeDerivatives(self, k_z, modeNumber, s, gValue)
            if abs(k_z) <= eps
                dFdz = zeros(size(s));
                dGdz = zeros(size(s));
                return;
            end

            signValue = (-1)^modeNumber;
            dFdz = -signValue*k_z*sin(k_z*s);
            dGdz = signValue*(gValue*k_z*k_z/(self.N0*self.N0))*cos(k_z*s);
        end

        function [G, F] = rawBoundaryMode(self, solutionType, k_z, hMode, s)
            depth = diff(self.zDomain);
            switch string(solutionType)
                case "linear"
                    G = s;
                    F = depth*ones(size(s));
                case "hyperbolic"
                    G = sinh(k_z*s);
                    F = hMode*k_z*cosh(k_z*s);
                case "trig"
                    G = sin(k_z*s);
                    F = hMode*k_z*cos(k_z*s);
                otherwise
                    error("IMConstantStratificationSolution:InvalidSolutionType", "Unknown solution type ""%s"".", string(solutionType));
            end
        end

        function [dGdz, dFdz] = rawBoundaryModeDerivatives(~, solutionType, k_z, hMode, s)
            switch string(solutionType)
                case "linear"
                    dGdz = ones(size(s));
                    dFdz = zeros(size(s));
                case "hyperbolic"
                    dGdz = k_z*cosh(k_z*s);
                    dFdz = hMode*k_z*k_z*sinh(k_z*s);
                case "trig"
                    dGdz = k_z*cos(k_z*s);
                    dFdz = -hMode*k_z*k_z*sin(k_z*s);
                otherwise
                    error("IMConstantStratificationSolution:InvalidSolutionType", "Unknown solution type ""%s"".", string(solutionType));
            end
        end

        function modeData = zeroAPVModeData(self,problem)
            nEndpoints = numel(problem.endpoints);
            nK = numel(problem.k);
            coefficients = zeros(2,nEndpoints,nK);
            m = problem.k*self.N0/abs(self.f0);
            responseReciprocalCondition = zeros(1,nK);
            for iK = 1:nK
                [FSurface,FzSurface] = self.zeroAPVFundamental(m(iK),problem.zDomain(2));
                [~,FzBottom] = self.zeroAPVFundamental(m(iK),problem.zDomain(1));
                GSurface = -(self.g/(self.N0*self.N0))*FzSurface;
                GBottom = -(self.g/(self.N0*self.N0))*FzBottom;
                surfaceResponse = GSurface;
                if problem.surfaceBoundary == "freeSurface"
                    surfaceResponse = surfaceResponse-FSurface;
                end
                responseMatrix = [surfaceResponse; GBottom];
                responseReciprocalCondition(iK) = rcond(responseMatrix);
                if responseReciprocalCondition(iK) <= 100*eps
                    error("IMConstantStratificationSolution:DegenerateEndpointResponse", "The exact endpoint-response matrix is numerically singular at k=%g.",problem.k(iK));
                end
                canonicalCoefficients = responseMatrix\eye(2);
                columnIndices = 1 + (problem.endpoints == "bottom");
                coefficients(:,:,iK) = canonicalCoefficients(:,columnIndices);
            end
            modeData = struct(m=m,coefficients=coefficients,responseReciprocalCondition=responseReciprocalCondition);
        end

        function values = zeroAPVVariable(self,modeData,problem,variable,z)
            z = z(:);
            values = zeros(numel(z),numel(problem.endpoints),numel(problem.k));
            for iK = 1:numel(problem.k)
                [fundamental,FzFundamental] = self.zeroAPVFundamental(modeData.m(iK),z);
                if string(variable) == "F"
                    raw = fundamental;
                else
                    raw = -(self.g/(self.N0*self.N0))*FzFundamental;
                end
                values(:,:,iK) = raw*modeData.coefficients(:,:,iK);
            end
        end

        function [F,Fz] = zeroAPVFundamental(self,m,z)
            midpoint = mean(self.zDomain);
            halfDepth = diff(self.zDomain)/2;
            y = z(:)-midpoint;
            if m*halfDepth < 50
                scale = exp(-m*halfDepth);
                first = scale*cosh(m*y);
                second = scale*sinh(m*y)/m;
            else
                upperExponential = exp(m*(y-halfDepth));
                lowerExponential = exp(-m*(y+halfDepth));
                first = 0.5*(upperExponential+lowerExponential);
                second = 0.5*(upperExponential-lowerExponential)/m;
            end
            F = [first second];
            Fz = [m*m*second first];
        end
    end

    methods (Static, Access = private)
        function validateEVPDomain(evp, zDomain)
            tolerance = 100*eps(max([1 abs(evp.zDomain) abs(zDomain)]));
            if max(abs(evp.zDomain - zDomain)) > tolerance
                error("IMConstantStratificationSolution:DomainMismatch", "The analytical solution zDomain must match evp.zDomain.");
            end
        end

        function [h, verticalWavenumbers, solutionTypes, isBoundaryMode, baroclinicNumbers, modeNumber] = solveSpectrum(evp, N0, zDomain, nModes, f0, gValue)
            [surfaceBoundary, bottomBoundary] = IMConstantStratificationSolution.validateEVP(evp);
            depth = diff(zDomain);
            evpName = string(evp.name);

            switch evpName
                case "waveModesAtWavenumber"
                    if ~isfield(evp.parameters, "k")
                        error("IMConstantStratificationSolution:UnsupportedEVP", "A fixed-wavenumber EVP must include parameters.k.");
                    end
                    k = evp.parameters.k;
                    [hBaroclinic, k_zBaroclinic] = IMConstantStratificationSolution.baroclinicAtWavenumber(k, N0, depth, nModes, f0, gValue, surfaceBoundary);
                    if surfaceBoundary == "free"
                        [h0, k_z0, solutionType0] = IMConstantStratificationSolution.surfaceBoundaryAtWavenumber(k, N0, depth, f0, gValue);
                        h = [h0 hBaroclinic(1:end-1)];
                        verticalWavenumbers = [k_z0 k_zBaroclinic(1:end-1)];
                        solutionTypes = IMConstantStratificationSolution.freeSolutionTypes(solutionType0, nModes);
                        isBoundaryMode = [true false(1,nModes-1)];
                        baroclinicNumbers = [0 1:(nModes-1)];
                    else
                        h = hBaroclinic;
                        verticalWavenumbers = k_zBaroclinic;
                        solutionTypes = repmat("baroclinic",1,nModes);
                        isBoundaryMode = false(1,nModes);
                        baroclinicNumbers = 1:nModes;
                    end
                case "waveModesAtFrequency"
                    if ~isfield(evp.parameters, "omega")
                        error("IMConstantStratificationSolution:UnsupportedEVP", "A fixed-frequency EVP must include parameters.omega.");
                    end
                    omega = evp.parameters.omega;
                    if surfaceBoundary == "free"
                        [h0, k_z0, solutionType0] = IMConstantStratificationSolution.surfaceBoundaryAtFrequency(omega, N0, depth, gValue);
                        if omega < N0
                            [hBaroclinic, k_zBaroclinic] = IMConstantStratificationSolution.freeSurfaceBaroclinicAtFrequency(omega, N0, depth, max(nModes - 1,0), gValue);
                        else
                            hBaroclinic = zeros(1,0);
                            k_zBaroclinic = zeros(1,0);
                        end
                        h = [h0 hBaroclinic];
                        verticalWavenumbers = [k_z0 k_zBaroclinic];
                        solutionTypes = IMConstantStratificationSolution.freeSolutionTypes(solutionType0, length(h));
                        isBoundaryMode = [true false(1,length(h)-1)];
                        baroclinicNumbers = [0 1:(length(h)-1)];
                    else
                        if omega >= N0
                            error("IMConstantStratificationSolution:UnsupportedFrequency", "Rigid-surface fixed-frequency constant-stratification modes require omega < N0.");
                        end
                        [hBaroclinic, k_zBaroclinic] = IMConstantStratificationSolution.baroclinicAtFrequency(omega, N0, depth, nModes, gValue);
                        h = hBaroclinic;
                        verticalWavenumbers = k_zBaroclinic;
                        solutionTypes = repmat("baroclinic",1,nModes);
                        isBoundaryMode = false(1,nModes);
                        baroclinicNumbers = 1:nModes;
                    end
                case "hydrostaticGModes"
                    if surfaceBoundary == "free"
                        [h0, k_z0, solutionType0] = IMConstantStratificationSolution.surfaceBoundaryAtFrequency(0, N0, depth, gValue);
                        [hBaroclinic, k_zBaroclinic] = IMConstantStratificationSolution.freeSurfaceBaroclinicAtFrequency(0, N0, depth, max(nModes - 1,0), gValue);
                        h = [h0 hBaroclinic];
                        verticalWavenumbers = [k_z0 k_zBaroclinic];
                        solutionTypes = IMConstantStratificationSolution.freeSolutionTypes(solutionType0, length(h));
                        isBoundaryMode = [true false(1,length(h)-1)];
                        baroclinicNumbers = [0 1:(length(h)-1)];
                    else
                        [h, verticalWavenumbers] = IMConstantStratificationSolution.baroclinicAtFrequency(0, N0, depth, nModes, gValue);
                        solutionTypes = repmat("baroclinic",1,nModes);
                        isBoundaryMode = false(1,nModes);
                        baroclinicNumbers = 1:nModes;
                    end
                case "hydrostaticFModes"
                    if surfaceBoundary ~= "rigid" || bottomBoundary ~= "rigid"
                        error("IMConstantStratificationSolution:UnsupportedBoundary", "Hydrostatic F constant-stratification modes currently require rigid surface and rigid bottom boundaries.");
                    end
                    interiorNumbers = 1:(nModes-1);
                    verticalWavenumbers = [0 interiorNumbers*pi/depth];
                    h = [Inf (N0*N0)./(gValue*verticalWavenumbers(2:end).*verticalWavenumbers(2:end))];
                    solutionTypes = ["null" repmat("baroclinic",1,nModes-1)];
                    isBoundaryMode = false(1,nModes);
                    baroclinicNumbers = 0:(nModes-1);
                otherwise
                    error("IMConstantStratificationSolution:UnsupportedEVP", "Constant stratification supports fixed-wavenumber, fixed-frequency, hydrostatic G, and hydrostatic F EVPs.");
            end
            modeNumber = baroclinicNumbers;
            modeNumber(isBoundaryMode) = -1;
        end

        function [h, k_z] = baroclinicAtWavenumber(k, N0, depth, nModes, f0, gValue, surfaceBoundary)
            k_z = (1:nModes)*pi/depth;
            if surfaceBoundary == "free"
                for iMode = 1:nModes
                    residual = @(xi) (xi + iMode*pi)*(N0*N0 - f0*f0)*depth - gValue*(k*k*depth*depth + (xi + iMode*pi)*(xi + iMode*pi))*tan(xi);
                    k_z(iMode) = k_z(iMode) + fzero(residual, 0)/depth;
                end
            end
            h = (N0*N0 - f0*f0)./(gValue*(k*k + k_z.*k_z));
        end

        function solutionTypes = freeSolutionTypes(boundaryType, nModes)
            solutionTypes = repmat("baroclinic",1,nModes);
            solutionTypes(1) = boundaryType;
        end

        function [h, k_z] = baroclinicAtFrequency(omega, N0, depth, nModes, gValue)
            k_z = (1:nModes)*pi/depth;
            if omega >= N0
                error("IMConstantStratificationSolution:UnsupportedFrequency", "Interior fixed-frequency constant-stratification modes require omega < N0.");
            end
            h = (N0*N0 - omega*omega)./(gValue*k_z.*k_z);
        end

        function [h, k_z] = freeSurfaceBaroclinicAtFrequency(omega, N0, depth, nModes, gValue)
            if omega >= N0
                error("IMConstantStratificationSolution:UnsupportedFrequency", "Free-surface interior fixed-frequency constant-stratification modes require omega < N0.");
            end
            k_z = (1:nModes)*pi/depth;
            for iMode = 1:nModes
                residual = @(xi) gValue*tan(xi)/(xi + iMode*pi) - (N0*N0 - omega*omega)*depth/((xi + iMode*pi)^2);
                k_z(iMode) = k_z(iMode) + fzero(residual, 0)/depth;
            end
            h = (N0*N0 - omega*omega)./(gValue*k_z.*k_z);
        end

        function [h0, k_z, solutionType] = surfaceBoundaryAtWavenumber(k, N0, depth, f0, gValue)
            kStar = sqrt((N0*N0 - f0*f0)/(gValue*depth));
            if abs(k - kStar)/kStar < 1e-6
                solutionType = "linear";
                h0 = depth;
                k_z = 0;
            elseif k > kStar
                solutionType = "hyperbolic";
                residual = @(q) depth*(N0*N0 - f0*f0) - (1./q).*(gValue*(k*k*depth*depth - q.*q)).*tanh(q);
                kInitial = sqrt(k*k*depth*depth - depth*(N0*N0 - f0*f0)/gValue);
                k_z = fzero(residual, kInitial)/depth;
                h0 = (N0*N0 - f0*f0)/(gValue*(k*k - k_z*k_z));
            else
                solutionType = "trig";
                residual = @(q) depth*(N0*N0 - f0*f0) - (1./q).*(gValue*(k*k*depth*depth + q.*q)).*tan(q);
                kInitial = sqrt(-k*k*depth*depth + depth*(N0*N0 - f0*f0)/gValue);
                k_z = fzero(residual, kInitial)/depth;
                h0 = (N0*N0 - f0*f0)/(gValue*(k*k + k_z*k_z));
            end
        end

        function [h0, k_z, solutionType] = surfaceBoundaryAtFrequency(omega, N0, depth, gValue)
            if abs(omega - N0)/N0 < 1e-6
                solutionType = "linear";
                h0 = depth;
                k_z = 0;
            elseif omega > N0
                solutionType = "hyperbolic";
                residual = @(q) depth*(omega*omega - N0*N0) - gValue*q.*tanh(q);
                k_z = fzero(residual, sqrt(depth*(omega*omega - N0*N0)/gValue))/depth;
                h0 = (omega*omega - N0*N0)/(gValue*k_z*k_z);
            else
                solutionType = "trig";
                residual = @(q) depth*(N0*N0 - omega*omega) - gValue*q.*tan(q);
                k_z = fzero(residual, sqrt(depth*(N0*N0 - omega*omega)/gValue))/depth;
                h0 = (N0*N0 - omega*omega)/(gValue*k_z*k_z);
            end
        end

        function [surfaceBoundary, bottomBoundary] = validateEVP(evp)
            surfaceBoundary = IMConstantStratificationSolution.boundaryKind(evp.surfaceBoundary, evp.formulation);
            bottomBoundary = IMConstantStratificationSolution.boundaryKind(evp.bottomBoundary, evp.formulation);
            if evp.formulation == "F"
                if evp.name ~= "hydrostaticFModes"
                    error("IMConstantStratificationSolution:UnsupportedEVP", "Constant stratification supports F-formulation EVPs only for hydrostatic F modes.");
                end
                if surfaceBoundary ~= "rigid" || bottomBoundary ~= "rigid"
                    error("IMConstantStratificationSolution:UnsupportedBoundary", "Hydrostatic F constant-stratification modes currently require rigid surface and rigid bottom boundaries.");
                end
                return;
            end
            if evp.formulation ~= "G"
                error("IMConstantStratificationSolution:UnsupportedEVP", "Constant stratification supports only F or G formulations.");
            end
            if bottomBoundary ~= "rigid"
                error("IMConstantStratificationSolution:UnsupportedBoundary", "Unsupported constant-stratification bottom boundary ""%s"".", bottomBoundary);
            end
            if ~ismember(surfaceBoundary, ["rigid", "free"])
                error("IMConstantStratificationSolution:UnsupportedBoundary", "Unsupported constant-stratification surface boundary ""%s"".", surfaceBoundary);
            end
        end

        function kind = boundaryKind(boundary, formulation)
            formulation = string(formulation);
            if formulation == "G"
                if IMConstantStratificationSolution.isCondition(boundary, [1 0 0 0])
                    kind = "rigid";
                elseif IMConstantStratificationSolution.isCondition(boundary, [0 1 1 0])
                    kind = "free";
                elseif IMConstantStratificationSolution.isCondition(boundary, [0 1 0 0])
                    kind = "noSlip";
                else
                    kind = "custom";
                end
            else
                if IMConstantStratificationSolution.isCondition(boundary, [0 1 0 0])
                    kind = "rigid";
                elseif IMConstantStratificationSolution.isCondition(boundary, [1 0 0 0])
                    kind = "noSlip";
                else
                    kind = "custom";
                end
            end
        end

        function tf = isCondition(boundary, coefficients)
            actual = [boundary.a boundary.b boundary.c boundary.d];
            tolerance = 100*eps(max(1,max(abs([actual coefficients]))));
            tf = max(abs(actual - coefficients)) <= tolerance;
        end
    end
end
