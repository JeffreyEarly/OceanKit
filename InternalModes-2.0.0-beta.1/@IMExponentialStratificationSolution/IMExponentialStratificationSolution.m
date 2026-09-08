classdef IMExponentialStratificationSolution < IMAnalyticalSolution
    % Analytical solution family for exponential stratification.
    %
    % `IMExponentialStratificationSolution` owns the closed-form formulas for
    % $$N^2(z)=N_0^2\exp(2z/b)$$ on domains with the surface at $$z=0$$. It
    % can create exact internal-mode bases for supported rigid-bottom
    % internal-mode EVPs, generalized-energy geostrophic APV EVPs, and
    % exact geostrophic zero-APV boundary-response modes.
    %
    % The generalized-energy APV branch recognizes a hydrostatic `F` EVP
    % named `"geostrophicAPVModes"` with canonical coefficients
    % $$p=1/N^2$$, $$q=0$$, and $$r=1/g$$. Its parameter struct must contain
    % the signed endpoint accelerations `g0` and `gd`, together with
    % `surfaceBoundary="freeSurface"` or `"rigidLid"`. Finite, zero, and
    % positive-infinite endpoint values are supported. The returned exact
    % basis is ordered by $$1/h$$ and may contain negative modes, an exact
    % zero-eigenvalue mode represented by `h=Inf`, and positive modes.
    % The public factory supplies the endpoint descriptor and direct APV
    % metadata. Exact APV bases use the volume-only `depth` normalization
    % by default; the Bessel formulas and endpoint inner products remain
    % available under every other supported normalization.
    %
    % Exact geostrophic zero-APV modes use scaled modified-Bessel functions
    % and exact derivatives. A two-column response solve supplies unit
    % surface and bottom coordinates for free-surface and rigid-lid models.
    % The older one-boundary profile formula used unit `f0*Fz` at its active
    % boundary; the public exact basis uses the canonical response
    % normalization shared with the numerical zero-APV basis.
    %
    % ```matlab
    % solution = IMExponentialStratificationSolution(N0=5.2e-3,b=1300,zDomain=[-5000 0]);
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
    % - Declaration: classdef IMExponentialStratificationSolution < IMAnalyticalSolution

    properties (SetAccess = private)
        % Surface buoyancy frequency $$N_0$$ in radians per second.
        %
        % - Topic: Inspect analytical solutions
        N0

        % Exponential e-folding depth $$b$$ in meters.
        %
        % - Topic: Inspect analytical solutions
        b
    end

    methods
        function self = IMExponentialStratificationSolution(options)
            % Create an exponential-stratification analytical solution family.
            %
            % - Topic: Create analytical solutions
            % - Declaration: solution = IMExponentialStratificationSolution(options)
            % - Parameter options.N0: surface buoyancy frequency
            % - Parameter options.b: exponential e-folding depth
            % - Parameter options.zDomain: physical vertical domain with surface at zero
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Returns solution: analytical solution family
            arguments
                options.N0 (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 5.2e-3
                options.b (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 1300
                options.zDomain (1,2) double {mustBeReal, mustBeFinite} = [-1 0]
                options.f0 (1,1) double {mustBeReal, mustBeFinite} = 0
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
            end

            zDomain = sort(options.zDomain);
            IMExponentialStratificationSolution.validateDomain(zDomain);
            self@IMAnalyticalSolution(zDomain=zDomain, f0=options.f0, g=options.g);
            self.N0 = options.N0;
            self.b = options.b;
        end

        function values = N2(self, z)
            % Evaluate $$N^2(z)=N_0^2\exp(2z/b)$$.
            %
            % - Topic: Inspect analytical solutions
            % - Declaration: values = N2(solution,z)
            % - Parameter z: physical coordinate
            % - Returns values: buoyancy frequency squared
            arguments
                self IMExponentialStratificationSolution
                z double {mustBeReal, mustBeFinite}
            end

            values = self.N0*self.N0*exp(2*z/self.b);
        end

        function basisSet = internalModes(self, evp, options)
            % Create an exact internal-mode basis.
            %
            % Generalized-energy APV modes use ordinary Bessel functions on
            % positive-$$h$$ branches, modified Bessel functions on
            % negative-$$h$$ branches, and the exact integrated solution on
            % a zero branch. Endpoint inertia determines whether zero, one,
            % or two negative modes precede the zero and positive modes.
            % Both `geostrophic` and volume-only `depth` normalization
            % factors remain real and positive for every retained branch;
            % `depth` is the APV default.
            %
            % - Topic: Compute internal modes
            % - Declaration: basisSet = internalModes(solution,evp,options)
            % - Parameter evp: internal-mode EVP
            % - Parameter options.nModes: number of retained modes
            % - Parameter options.normalization: active normalization
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: exact analytical internal-mode basis
            arguments
                self IMExponentialStratificationSolution
                evp = []
                options.nModes (1,1) double {mustBeInteger, mustBePositive} = 64
                options.normalization = []
                options.metadata struct = struct()
            end

            evp = self.resolveEVP(evp);
            if IMExponentialGeostrophicAPVCatalog.isTargetEVP(evp)
                [modeData, h, modeNumber, analyticalMetadata, diagnostics] = ...
                    IMExponentialGeostrophicAPVCatalog.solve( ...
                    evp, self.N0, self.b, self.zDomain, options.nModes);
                metadata = options.metadata;
                metadataFields = fieldnames(analyticalMetadata);
                for iField = 1:length(metadataFields)
                    fieldName = metadataFields{iField};
                    metadata.(fieldName) = analyticalMetadata.(fieldName);
                end
                metadata.analyticalSolution = "exponentialStratification";
                metadata.analyticalFamily = "generalizedEnergyAPV";
                basisSet = IMAnalyticalInternalModesBasis( ...
                    solution=self, ...
                    evp=evp, ...
                    h=h, ...
                    modeNumber=modeNumber, ...
                    N2=@(z) self.N2(z), ...
                    rawVariableFunction=@(~,variable,z) IMExponentialGeostrophicAPVCatalog.rawVariable(modeData, self.N0, self.b, evp.g, variable, z), ...
                    rawUzFunction=@(~,z) IMExponentialGeostrophicAPVCatalog.rawFz(modeData, self.N0, self.b, evp.g, z), ...
                    normalization=options.normalization, ...
                    metadata=metadata, ...
                    modeSelectionDiagnostics=diagnostics);
                return;
            end
            f0 = evp.f0;
            gValue = evp.g;
            [h, roots, frequencies, modeNumber, modeKinds] = IMExponentialStratificationSolution.solveSpectrum(evp, self.N0, self.b, self.zDomain, options.nModes, f0, gValue);
            phaseSpeeds = sqrt(gValue*h);
            signFactors = IMExponentialStratificationSolution.surfaceSignFactors(self.N0, self.b, self.zDomain, frequencies, phaseSpeeds, gValue, modeKinds);

            modeData = struct();
            modeData.h = h;
            modeData.roots = roots;
            modeData.frequencies = frequencies;
            modeData.phaseSpeeds = phaseSpeeds;
            modeData.signFactors = signFactors;
            modeData.modeKinds = modeKinds;

            metadata = options.metadata;
            metadata.analyticalSolution = "exponentialStratification";
            metadata.roots = roots;
            metadata.frequencies = frequencies;
            metadata.phaseSpeeds = phaseSpeeds;
            metadata.signFactors = signFactors;
            metadata.modeKinds = modeKinds;

            basisSet = IMAnalyticalInternalModesBasis(solution=self, evp=evp, h=h, modeNumber=modeNumber, N2=@(z) self.N2(z), rawVariableFunction=@(basisSet,variable,z) self.rawVariable(modeData, evp, variable, z), rawUzFunction=@(basisSet,z) self.rawUz(modeData, evp, z), normalization=options.normalization, metadata=metadata);
        end

        function basisSet = geostrophicZeroAPVModesAtWavenumber(self,k,options)
            % Create exact canonical geostrophic zero-APV modes.
            %
            % The exact modified-Bessel columns solve the zero-APV equation
            % and are rescaled to unit endpoint responses under either the
            % free-surface or rigid-lid convention.
            %
            % - Topic: Compute geostrophic zero-APV modes
            % - Declaration: exactModes = geostrophicZeroAPVModesAtWavenumber(solution,k,options)
            % - Parameter k: positive horizontal wavenumbers
            % - Parameter options.endpoints: requested surface and bottom coordinates
            % - Parameter options.surfaceBoundary: `"freeSurface"` or `"rigidLid"`
            % - Parameter options.metadata: additional metadata
            % - Returns exactModes: exact boundary-normalized basis
            arguments
                self IMExponentialStratificationSolution
                k double {mustBeReal, mustBeFinite, mustBePositive}
                options.endpoints {mustBeText} = ["surface", "bottom"]
                options.surfaceBoundary {mustBeTextScalar} = "freeSurface"
                options.metadata struct = struct()
            end

            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=@(z) self.N2(z),zDomain=self.zDomain,f0=self.f0,g=self.g,k=k,endpoints=options.endpoints,surfaceBoundary=options.surfaceBoundary,metadata=options.metadata);
            modeData = self.zeroAPVModeData(problem);
            metadata = options.metadata;
            metadata.analyticalSolution = "exponentialStratification";
            metadata.analyticalFamily = "geostrophicZeroAPV";
            metadata.endpoints = problem.endpoints;
            metadata.surfaceBoundary = problem.surfaceBoundary;
            metadata.responseReciprocalCondition = modeData.responseReciprocalCondition;
            metadata.formulaFamily = "scaledModifiedBessel";
            basisSet = IMAnalyticalGeostrophicZeroAPVModesBasis(solution=self,problem=problem,FFunction=@(z) self.zeroAPVVariable(modeData,problem,"F",z),GFunction=@(z) self.zeroAPVVariable(modeData,problem,"G",z),metadata=metadata);
        end

        function summarize(self)
            % Print a readable solution-family summary.
            %
            % - Topic: Inspect analytical solutions
            % - Declaration: summarize(solution)
            arguments
                self IMExponentialStratificationSolution
            end

            summarize@IMAnalyticalSolution(self);
            fprintf("  N0: %g\n", self.N0);
            fprintf("  b: %g\n", self.b);
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
                error("IMExponentialStratificationSolution:InvalidEVP", "The EVP must be an IMInternalModes instance.");
            end
            IMExponentialStratificationSolution.validateEVPDomain(evp, self.zDomain);
        end

        function values = rawVariable(self, modeData, evp, variable, z)
            variable = string(variable);
            if variable ~= "G" && variable ~= "F"
                error("IMExponentialStratificationSolution:UnsupportedVariable", "Variable must be ""F"" or ""G"".");
            end

            z = z(:);
            values = zeros(length(z), length(modeData.h));
            for iMode = 1:length(modeData.h)
                if modeData.modeKinds(iMode) == "null"
                    G = zeros(size(z));
                    F = ones(size(z));
                else
                    [G, F] = IMExponentialStratificationSolution.rawModeValues(z, self.N0, self.b, self.zDomain, modeData.frequencies(iMode), modeData.phaseSpeeds(iMode), evp.g);
                end
                G = modeData.signFactors(iMode)*G;
                F = modeData.signFactors(iMode)*F;
                if variable == "G"
                    values(:,iMode) = G;
                else
                    values(:,iMode) = F;
                end
            end
        end

        function values = rawUz(self, modeData, evp, z)
            z = z(:);
            values = zeros(length(z), length(modeData.h));
            for iMode = 1:length(modeData.h)
                if modeData.modeKinds(iMode) == "null"
                    values(:,iMode) = zeros(size(z));
                    continue;
                end
                [G, F] = IMExponentialStratificationSolution.rawModeValues(z, self.N0, self.b, self.zDomain, modeData.frequencies(iMode), modeData.phaseSpeeds(iMode), evp.g);
                G = modeData.signFactors(iMode)*G;
                F = modeData.signFactors(iMode)*F;

                if evp.formulation == "G"
                    values(:,iMode) = F ./ modeData.h(iMode);
                else
                    values(:,iMode) = -(self.N2(z)/evp.g).*G;
                end
            end
        end

        function modeData = zeroAPVModeData(self,problem)
            nEndpoints = numel(problem.endpoints);
            nK = numel(problem.k);
            coefficients = zeros(2,nEndpoints,nK);
            eta = self.b*self.N0*problem.k/abs(self.f0);
            responseReciprocalCondition = zeros(1,nK);
            for iK = 1:nK
                [FSurface,FzSurface] = self.zeroAPVFundamental(eta(iK),problem.zDomain(2));
                [~,FzBottom] = self.zeroAPVFundamental(eta(iK),problem.zDomain(1));
                GSurface = -(self.g/self.N2(problem.zDomain(2)))*FzSurface;
                GBottom = -(self.g/self.N2(problem.zDomain(1)))*FzBottom;
                surfaceResponse = GSurface;
                if problem.surfaceBoundary == "freeSurface"
                    surfaceResponse = surfaceResponse-FSurface;
                end
                responseMatrix = [surfaceResponse; GBottom];
                responseReciprocalCondition(iK) = rcond(responseMatrix);
                if responseReciprocalCondition(iK) <= 100*eps
                    error("IMExponentialStratificationSolution:DegenerateEndpointResponse", "The exact endpoint-response matrix is numerically singular at k=%g.",problem.k(iK));
                end
                canonicalCoefficients = responseMatrix\eye(2);
                columnIndices = 1 + (problem.endpoints == "bottom");
                coefficients(:,:,iK) = canonicalCoefficients(:,columnIndices);
            end
            modeData = struct(eta=eta,coefficients=coefficients,responseReciprocalCondition=responseReciprocalCondition);
        end

        function values = zeroAPVVariable(self,modeData,problem,variable,z)
            z = z(:);
            values = zeros(numel(z),numel(problem.endpoints),numel(problem.k));
            N2Values = self.N2(z);
            for iK = 1:numel(problem.k)
                [fundamental,FzFundamental] = self.zeroAPVFundamental(modeData.eta(iK),z);
                if string(variable) == "F"
                    raw = fundamental;
                else
                    raw = -(self.g./N2Values).*FzFundamental;
                end
                values(:,:,iK) = raw*modeData.coefficients(:,:,iK);
            end
        end

        function [F,Fz] = zeroAPVFundamental(self,eta,z)
            z = z(:);
            exponential = exp(z/self.b);
            x = eta*exponential;
            xSurface = eta*exp(self.zDomain(2)/self.b);
            xBottom = eta*exp(self.zDomain(1)/self.b);
            IScale = exp(x-xSurface);
            KScale = exp(xBottom-x);
            I1 = IScale.*besseli(1,x,1);
            K1 = KScale.*besselk(1,x,1);
            I0 = IScale.*besseli(0,x,1);
            K0 = KScale.*besselk(0,x,1);
            F = [exponential.*I1 exponential.*K1];
            derivativeScale = (x/self.b).*exponential;
            Fz = [derivativeScale.*I0 -derivativeScale.*K0];
        end
    end

    methods (Static, Access = private)
        function validateEVPDomain(evp, zDomain)
            tolerance = 100*eps(max([1 abs(evp.zDomain) abs(zDomain)]));
            if max(abs(evp.zDomain - zDomain)) > tolerance
                error("IMExponentialStratificationSolution:DomainMismatch", "The analytical solution zDomain must match evp.zDomain.");
            end
        end

        function validateDomain(zDomain)
            tolerance = 100*eps(max(1,max(abs(zDomain))));
            if abs(zDomain(2)) > tolerance
                error("IMExponentialStratificationSolution:UnsupportedDomain", "Exponential-stratification analytical modes currently require the surface at z=0.");
            end
        end

        function [h, roots, frequencies, modeNumber, modeKinds] = solveSpectrum(evp, N0, b, zDomain, nModes, f0, gValue)
            surfaceBoundary = IMExponentialStratificationSolution.validateEVP(evp);
            evpName = string(evp.name);
            switch evpName
                case "waveModesAtWavenumber"
                    if ~isfield(evp.parameters, "k")
                        error("IMExponentialStratificationSolution:UnsupportedEVP", "A fixed-wavenumber EVP must include parameters.k.");
                    end
                    k = evp.parameters.k;
                    nInteriorModes = IMExponentialStratificationSolution.nInteriorModes(nModes, surfaceBoundary);
                    roots = IMExponentialStratificationSolution.rootsAtWavenumber(k, N0, b, zDomain, nInteriorModes, f0, gValue, surfaceBoundary);
                    h = (b*N0./roots).^2/gValue;
                    if surfaceBoundary == "free"
                        root0 = IMExponentialStratificationSolution.surfaceRootAtWavenumber(k, N0, b, zDomain, f0, gValue);
                        roots = [root0 roots];
                        h = [(b*N0/root0)^2/gValue h];
                    end
                    frequencies = sqrt(gValue*h*k*k + f0*f0);
                case "waveModesAtFrequency"
                    if ~isfield(evp.parameters, "omega")
                        error("IMExponentialStratificationSolution:UnsupportedEVP", "A fixed-frequency EVP must include parameters.omega.");
                    end
                    omega = evp.parameters.omega;
                    if omega >= N0 && surfaceBoundary == "free"
                        root0 = IMExponentialStratificationSolution.surfaceRootAtFrequency(omega, N0, b, zDomain, gValue);
                        roots = root0;
                        h = (b*N0/root0)^2/gValue;
                    else
                        nInteriorModes = IMExponentialStratificationSolution.nInteriorModes(nModes, surfaceBoundary);
                        [roots, eta] = IMExponentialStratificationSolution.rootsAtFrequency(omega, N0, b, zDomain, nInteriorModes, gValue, surfaceBoundary);
                        h = (b*eta./roots).^2/gValue;
                        if surfaceBoundary == "free"
                            root0 = IMExponentialStratificationSolution.surfaceRootAtFrequency(omega, N0, b, zDomain, gValue);
                            roots = [root0 roots];
                            h = [(b*N0/root0)^2/gValue h];
                        end
                    end
                    frequencies = omega*ones(size(h));
                case "hydrostaticGModes"
                    omega = 0;
                    nInteriorModes = IMExponentialStratificationSolution.nInteriorModes(nModes, surfaceBoundary);
                    [roots, eta] = IMExponentialStratificationSolution.rootsAtFrequency(omega, N0, b, zDomain, nInteriorModes, gValue, surfaceBoundary);
                    h = (b*eta./roots).^2/gValue;
                    if surfaceBoundary == "free"
                        root0 = IMExponentialStratificationSolution.surfaceRootAtFrequency(omega, N0, b, zDomain, gValue);
                        roots = [root0 roots];
                        h = [(b*N0/root0)^2/gValue h];
                    end
                    frequencies = zeros(size(h));
                case "hydrostaticFModes"
                    omega = 0;
                    nInteriorModes = max(nModes - 1,0);
                    [rootsInterior, eta] = IMExponentialStratificationSolution.rootsAtFrequency(omega, N0, b, zDomain, nInteriorModes, gValue, "rigid");
                    hInterior = (b*eta./rootsInterior).^2/gValue;
                    roots = [NaN rootsInterior];
                    h = [Inf hInterior];
                    frequencies = zeros(size(h));
                otherwise
                    error("IMExponentialStratificationSolution:UnsupportedEVP", "Exponential stratification supports fixed-wavenumber, fixed-frequency, hydrostatic G, and hydrostatic F EVPs.");
            end
            if evpName == "hydrostaticFModes"
                modeKinds = ["null" repmat("bessel",1,length(h)-1)];
            else
                modeKinds = repmat("bessel",1,length(h));
            end
            [roots, h, frequencies, modeKinds] = IMExponentialStratificationSolution.truncateModes(roots, h, frequencies, modeKinds, nModes);
            modeNumber = IMExponentialStratificationSolution.modeNumberFor(surfaceBoundary, length(h));
            if evpName == "hydrostaticFModes"
                modeNumber = 0:(length(h)-1);
            end
        end

        function surfaceBoundary = validateEVP(evp)
            surfaceBoundary = IMExponentialStratificationSolution.boundaryKind(evp.surfaceBoundary, evp.formulation);
            bottomBoundary = IMExponentialStratificationSolution.boundaryKind(evp.bottomBoundary, evp.formulation);
            if evp.formulation == "F"
                if evp.name ~= "hydrostaticFModes"
                    error("IMExponentialStratificationSolution:UnsupportedEVP", "Exponential stratification supports F-formulation EVPs only for hydrostatic F modes.");
                end
                if surfaceBoundary ~= "rigid" || bottomBoundary ~= "rigid"
                    error("IMExponentialStratificationSolution:UnsupportedBoundary", "Hydrostatic F exponential-stratification modes currently require rigid surface and rigid bottom boundaries.");
                end
                return;
            end
            if evp.formulation ~= "G"
                error("IMExponentialStratificationSolution:UnsupportedEVP", "Exponential stratification supports only F or G formulations.");
            end
            if bottomBoundary ~= "rigid"
                error("IMExponentialStratificationSolution:UnsupportedBoundary", "Exponential stratification currently supports only rigid bottom boundaries.");
            end
            if ~ismember(surfaceBoundary, ["rigid", "free"])
                error("IMExponentialStratificationSolution:UnsupportedBoundary", "Exponential stratification currently supports rigid or free surface boundaries.");
            end
        end

        function kind = boundaryKind(boundary, formulation)
            formulation = string(formulation);
            if formulation == "G"
                if IMExponentialStratificationSolution.isCondition(boundary, [1 0 0 0])
                    kind = "rigid";
                elseif IMExponentialStratificationSolution.isCondition(boundary, [0 1 1 0])
                    kind = "free";
                elseif IMExponentialStratificationSolution.isCondition(boundary, [0 1 0 0])
                    kind = "noSlip";
                else
                    kind = "custom";
                end
            else
                if IMExponentialStratificationSolution.isCondition(boundary, [0 1 0 0])
                    kind = "rigid";
                elseif IMExponentialStratificationSolution.isCondition(boundary, [1 0 0 0])
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

        function roots = rootsAtWavenumber(k, N0, b, zDomain, nModes, f0, gValue, surfaceBoundary)
            if nModes == 0
                roots = zeros(1,0);
                return;
            end
            depth = diff(zDomain);
            epsilon = f0/N0;
            lambda = k*b;
            x_lf = @(j,lambdaValue) (j - 1/4)*pi + lambdaValue*pi/2;
            x_hf = @(j,lambdaValue) lambdaValue.*(1 + 0.5*(3*pi*(4*j - 1)./(lambdaValue*8*sqrt(2))).^(2/3));
            if lambda < 2*(1 - 1/4)*1e-1
                lowerBound = x_lf(1, lambda);
            else
                lowerBound = x_hf(1, lambda);
            end
            if lambda < (IMExponentialStratificationSolution.nInitialSearchModes() - 1/4)
                upperBound = x_lf(1.1*IMExponentialStratificationSolution.nInitialSearchModes(), lambda);
            else
                upperBound = x_hf(5*IMExponentialStratificationSolution.nInitialSearchModes(), lambda);
            end

            nu = @(x) sqrt(epsilon*epsilon*x.*x + lambda*lambda);
            s = @(x) x;
            roots = IMExponentialStratificationSolution.findEnoughRoots(nu, s, [lowerBound upperBound], exp(-depth/b), nModes, surfaceBoundary, N0, b, gValue);
        end

        function [roots, eta] = rootsAtFrequency(omega, N0, b, zDomain, nModes, gValue, surfaceBoundary)
            depth = diff(zDomain);
            expMinusDOverB = exp(-depth/b);
            if omega >= N0
                error("IMExponentialStratificationSolution:UnsupportedFrequency", "Fixed-frequency exponential modes require omega < N0.");
            end
            if omega > N0*expMinusDOverB
                eta = (sqrt(N0*N0 - omega*omega) - omega*acos(omega/N0))/pi;
            else
                eta = (sqrt(N0*N0 - omega*omega) - sqrt(N0*N0*exp(-2*depth/b) - omega*omega) - omega*acos(omega/N0) + omega*acos(omega*exp(depth/b)/N0))/pi;
            end

            nu = @(x) omega*x/eta;
            s = @(x) N0*x/eta;
            if nModes == 0
                roots = zeros(1,0);
                return;
            end
            roots = IMExponentialStratificationSolution.findEnoughRoots(nu, s, [0.5 nModes+1], expMinusDOverB, nModes, surfaceBoundary, N0, b, gValue);
        end

        function root = surfaceRootAtWavenumber(k, N0, b, zDomain, f0, gValue)
            depth = diff(zDomain);
            kSafe = max(k,1e-15);
            rootEstimate = b*N0/sqrt(gValue*tanh(kSafe*depth)/kSafe);
            epsilon = f0/N0;
            lambda = k*b;
            nu = @(x) sqrt(epsilon*epsilon*x.*x + lambda*lambda);
            s = @(x) x;
            root = IMExponentialStratificationSolution.findFirstRoot(nu, s, [0.95 1.05]*rootEstimate, exp(-depth/b), "free", N0, b, gValue);
        end

        function root = surfaceRootAtFrequency(omega, N0, b, zDomain, gValue)
            depth = diff(zDomain);
            rootEstimate = max(b*N0*omega/gValue, b*N0/sqrt(gValue*depth));
            nu = @(x) omega*x/N0;
            s = @(x) x;
            root = IMExponentialStratificationSolution.findFirstRoot(nu, s, [0.95 1.5]*rootEstimate, exp(-depth/b), "free", N0, b, gValue);
        end

        function root = findFirstRoot(nu, s, bounds, expMinusDOverB, surfaceBoundary, N0, b, gValue)
            roots = IMExponentialStratificationSolution.findRootsInRange(nu, s, bounds, expMinusDOverB, surfaceBoundary, N0, b, gValue);
            if isempty(roots)
                error("IMExponentialStratificationSolution:RootSearchFailed", "Could not find the exponential free-surface root.");
            end
            root = roots(1);
        end

        function roots = findEnoughRoots(nu, s, bounds, expMinusDOverB, nModes, surfaceBoundary, N0, b, gValue)
            roots = IMExponentialStratificationSolution.findRootsInRange(nu, s, bounds, expMinusDOverB, surfaceBoundary, N0, b, gValue);
            iteration = 0;
            while length(roots) < nModes && iteration < 20
                iteration = iteration + 1;
                if length(roots) >= 2
                    dx = roots(end) - roots(end-1);
                else
                    dx = diff(bounds)/IMExponentialStratificationSolution.nInitialSearchModes();
                end
                dx = max(dx, eps(max(1,bounds(2))));
                oldUpperBound = bounds(2);
                nSearchModes = min(IMExponentialStratificationSolution.nInitialSearchModes(), nModes - length(roots) + 1);
                bounds = [oldUpperBound oldUpperBound + dx*nSearchModes];
                newRoots = IMExponentialStratificationSolution.findRootsInRange(nu, s, bounds, expMinusDOverB, surfaceBoundary, N0, b, gValue);
                roots = IMExponentialStratificationSolution.deduplicateRoots([roots; newRoots]);
            end
            if length(roots) < nModes
                error("IMExponentialStratificationSolution:RootSearchFailed", "Could not find %d exponential-stratification roots.", nModes);
            end
            roots = reshape(roots(1:nModes),1,[]);
        end

        function roots = findRootsInRange(nu, s, bounds, expMinusDOverB, surfaceBoundary, N0, b, gValue)
            x = linspace(bounds(1), bounds(2), IMExponentialStratificationSolution.nInitialSearchModes()).';
            intervals = zeros(0,2);

            bigNuIndices = find(nu(x) >= s(x)*expMinusDOverB);
            if ~isempty(bigNuIndices) && length(bigNuIndices) > 1
                intervals(end+1,:) = [x(min(bigNuIndices)) x(max(bigNuIndices))];
            end

            smallNuIndices = find(nu(x) < s(x)*expMinusDOverB);
            if ~isempty(smallNuIndices) && length(smallNuIndices) > 1
                if ~isempty(bigNuIndices)
                    intervals(end+1,:) = [x(max(bigNuIndices)) x(max(smallNuIndices))];
                else
                    intervals(end+1,:) = [x(min(smallNuIndices)) x(max(smallNuIndices))];
                end
            end

            intervalRoots = cell(size(intervals,1),1);
            for iInterval = 1:size(intervals,1)
                interval = intervals(iInterval,:);
                if interval(2) <= interval(1)
                    intervalRoots{iInterval} = zeros(0,1);
                    continue;
                end
                if iInterval == 1 && ~isempty(bigNuIndices)
                    residual = @(x) IMExponentialStratificationSolution.bigNuResidual(nu, s, expMinusDOverB, x, surfaceBoundary, N0, b, gValue);
                else
                    residual = @(x) IMExponentialStratificationSolution.smallNuResidual(nu, s, expMinusDOverB, x, surfaceBoundary, N0, b, gValue);
                end
                intervalRoots{iInterval} = IMExponentialStratificationSolution.scanRoots(residual, interval);
            end
            roots = vertcat(intervalRoots{:});
            roots = IMExponentialStratificationSolution.deduplicateRoots(roots);
        end

        function values = smallNuResidual(nu, s, expMinusDOverB, x, surfaceBoundary, N0, b, gValue)
            [A, B] = IMExponentialStratificationSolution.surfaceCoefficients(nu, s, x, surfaceBoundary, N0, b, gValue);
            values = A.*besselj(nu(x), expMinusDOverB*s(x)) + B.*bessely(nu(x), expMinusDOverB*s(x));
        end

        function values = bigNuResidual(nu, s, expMinusDOverB, x, surfaceBoundary, N0, b, gValue)
            [A, B] = IMExponentialStratificationSolution.surfaceCoefficients(nu, s, x, surfaceBoundary, N0, b, gValue);
            denominator = bessely(nu(x), expMinusDOverB*s(x));
            values = (A./denominator).*besselj(nu(x), expMinusDOverB*s(x)) + B;
        end

        function [A, B] = surfaceCoefficients(nu, s, x, surfaceBoundary, N0, b, gValue)
            order = nu(x);
            argument = s(x);
            switch surfaceBoundary
                case "rigid"
                    A = bessely(order, argument);
                    B = -besselj(order, argument);
                case "free"
                    surfaceCoefficient = b*N0*N0/(2*gValue);
                    A = bessely(order, argument) - (surfaceCoefficient./argument).*(bessely(order-1, argument) - bessely(order+1, argument));
                    B = -besselj(order, argument) + (surfaceCoefficient./argument).*(besselj(order-1, argument) - besselj(order+1, argument));
                otherwise
                    error("IMExponentialStratificationSolution:UnsupportedBoundary", "Unsupported exponential surface boundary ""%s"".", surfaceBoundary);
            end
        end

        function roots = scanRoots(residual, interval)
            nSamples = max(2048, ceil(64*diff(interval)));
            x = linspace(interval(1), interval(2), nSamples).';
            y = IMExponentialStratificationSolution.safeEvaluate(residual, x);
            roots = zeros(length(x)-1,1);
            nRoots = 0;
            finite = isfinite(y);
            branchScale = median(abs(y(finite)));
            if isempty(branchScale) || ~isfinite(branchScale) || branchScale == 0
                branchScale = 1;
            end
            residualTolerance = 1e-6*max(1,branchScale);

            for iSample = 1:(length(x)-1)
                if ~finite(iSample) || ~finite(iSample+1)
                    continue;
                end
                if y(iSample) == 0
                    nRoots = nRoots + 1;
                    roots(nRoots) = x(iSample);
                elseif sign(y(iSample)) == sign(y(iSample+1))
                    continue;
                else
                    try
                        root = fzero(residual, [x(iSample) x(iSample+1)]);
                        residualAtRoot = IMExponentialStratificationSolution.safeEvaluate(residual, root);
                        if isfinite(root) && isfinite(residualAtRoot) && abs(residualAtRoot) <= residualTolerance
                            nRoots = nRoots + 1;
                            roots(nRoots) = root;
                        end
                    catch
                    end
                end
            end
            roots = roots(1:nRoots);
            roots = IMExponentialStratificationSolution.deduplicateRoots(roots);
        end

        function values = safeEvaluate(residual, x)
            try
                rawValues = residual(x);
            catch
                rawValues = arrayfun(residual, x);
            end
            values = real(rawValues);
            values(abs(imag(rawValues)) > 1e-10*max(1,abs(values))) = NaN;
        end

        function roots = deduplicateRoots(roots)
            roots = sort(real(roots(:)));
            roots = roots(isfinite(roots));
            if isempty(roots)
                return;
            end
            tolerance = 1e-8*max(1,max(abs(roots)));
            keep = [true; abs(diff(roots)) > tolerance];
            roots = roots(keep);
        end

        function [G, F] = rawModeValues(z, N0, b, zDomain, omega, c, gValue)
            depth = diff(zDomain);
            nu = b*omega/c;
            argument = (b*N0/c)*exp(z(:)/b);
            bottomArgument = (b*N0/c)*exp(-depth/b);
            alpha = besselj(nu, bottomArgument)/bessely(nu, bottomArgument);
            if abs(alpha) < 1e-15
                G = besselj(nu, argument);
                F = (N0*exp(z(:)/b)*c/(2*gValue)).*(besselj(nu-1, argument) - besselj(nu+1, argument));
            else
                G = besselj(nu, argument) - alpha*bessely(nu, argument);
                F = (N0*exp(z(:)/b)*c/(2*gValue)).*((besselj(nu-1, argument) - besselj(nu+1, argument)) - alpha*(bessely(nu-1, argument) - bessely(nu+1, argument)));
            end
        end

        function signFactors = surfaceSignFactors(N0, b, zDomain, frequencies, phaseSpeeds, gValue, modeKinds)
            signFactors = ones(1,length(frequencies));
            for iMode = 1:length(frequencies)
                if modeKinds(iMode) == "null"
                    continue;
                end
                [~, FSurface] = IMExponentialStratificationSolution.rawModeValues(0, N0, b, zDomain, frequencies(iMode), phaseSpeeds(iMode), gValue);
                if FSurface < 0
                    signFactors(iMode) = -1;
                end
            end
        end

        function value = nInitialSearchModes()
            value = 128;
        end

        function nModes = nInteriorModes(totalModes, surfaceBoundary)
            if surfaceBoundary == "free"
                nModes = max(0,totalModes - 1);
            else
                nModes = totalModes;
            end
        end

        function [roots, h, frequencies, modeKinds] = truncateModes(roots, h, frequencies, modeKinds, nModes)
            keep = 1:min(nModes,length(h));
            roots = roots(keep);
            h = h(keep);
            frequencies = frequencies(keep);
            modeKinds = modeKinds(keep);
        end

        function modeNumber = modeNumberFor(surfaceBoundary, nModes)
            if surfaceBoundary == "free"
                modeNumber = [-1 1:(nModes-1)];
            else
                modeNumber = 1:nModes;
            end
        end
    end
end
