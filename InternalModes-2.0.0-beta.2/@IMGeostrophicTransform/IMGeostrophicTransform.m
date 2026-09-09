classdef IMGeostrophicTransform
    % Compose APV and zero-APV transforms at positive wavenumber.
    %
    % For APV eigendepth $$h_j$$ and horizontal wavenumber $$\kappa$$,
    %
    % $$
    % \mu_\kappa^j=\kappa^2+\frac{f_0^2}{g h_j}.
    % $$
    %
    % The APV endpoint-response columns are
    %
    % $$
    % \mathbf r_q^{\kappa j}=-\frac{f_0}{g\mu_\kappa^j}
    % \begin{bmatrix}B_{\mathrm s}^j\\G_j(z_b)\end{bmatrix},
    % $$
    %
    % restricted to active endpoints, where
    % $$B_{\mathrm s}=G(0)-F(0)$$ for a free surface and
    % $$B_{\mathrm s}=G(0)$$ for a rigid lid.
    %
    % Sampled APV and source inputs have leading dimensions
    % `nZ x nK x ...`. APV coefficient outputs have leading dimensions
    % `nAPVModes x nK x ...`, while endpoint and zero-APV arrays have
    % leading dimensions `nEndpoints x nK x ...`. Trailing field
    % dimensions and complex values are preserved.
    %
    % ```matlab
    % transform = IMGeostrophicTransform(apvTransform=apvTransform,zeroAPVModes=zeroAPVModes,g0=g0,gd=gd);
    % [Aq,A0] = transform.transformStateForward(APV=q,endpointAnomalies=b);
    % ```
    %
    % - Topic: Create geostrophic transforms
    % - Topic: Transform admissible states
    % - Topic: Project generic sources
    % - Topic: Inspect geostrophic transforms
    % - Declaration: classdef IMGeostrophicTransform

    properties (SetAccess = private)
        % Horizontal-wavenumber pages $$\kappa$$.
        % - Topic: Inspect geostrophic transforms
        k

        % Active endpoint coordinates in canonical order.
        % - Topic: Inspect geostrophic transforms
        activeEndpoints

        % APV endpoint-response pages $$\mathsf R_q^\kappa$$.
        %
        % Dimensions are nEndpoints x nAPVModes x nK.
        % - Topic: Inspect geostrophic transforms
        apvEndpointResponse

        % Compatibility and singularity diagnostics.
        % - Topic: Inspect geostrophic transforms
        compatibilityDiagnostics

        % Surface generalized-energy acceleration $$g_0$$.
        % - Topic: Inspect geostrophic transforms
        g0

        % Bottom generalized-energy acceleration $$g_d$$.
        % - Topic: Inspect geostrophic transforms
        gd

        % Coriolis parameter $$f_0$$.
        % - Topic: Inspect geostrophic transforms
        f0
    end

    properties (Access = private)
        apvTransform
        zeroAPVModes
        g
        depth
        endpointIndices
        requiredSourceEndpointSamples
        zeroFPairingMatrices
        zeroGPairingMatrices
        zeroSourceSolveMatrices
    end

    methods
        function self = IMGeostrophicTransform(options)
            % Create an APV and zero-APV composition transform.
            %
            % Finite endpoint accelerations, including zero, activate the
            % corresponding endpoint. Positive infinity makes it inactive.
            % Construction rejects modes satisfying
            %
            % $$
            % \frac{|\mu_\kappa^j|}
            % {\kappa^2+|f_0^2/(g h_j)|}\leq\texttt{muTolerance}.
            % $$
            %
            % - Topic: Create geostrophic transforms
            % - Declaration: transform = IMGeostrophicTransform(options)
            % - Parameter options.apvTransform: generalized-energy APV transform
            % - Parameter options.zeroAPVModes: canonical zero-APV basis
            % - Parameter options.g0: surface acceleration
            % - Parameter options.gd: bottom acceleration
            % - Parameter options.muTolerance: relative singularity tolerance
            % - Returns self: geostrophic composition transform
            arguments
                options.apvTransform (1,1) IMInternalModesDiscreteTransform
                options.zeroAPVModes (1,1) IMGeostrophicZeroAPVModesBasis
                options.g0 (1,1) double {mustBeReal}
                options.gd (1,1) double {mustBeReal}
                options.muTolerance (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = sqrt(eps)
            end

            self.validateEndpointAcceleration(options.g0,"surface");
            self.validateEndpointAcceleration(options.gd,"bottom");
            if isinf(options.g0) && isinf(options.gd)
                error("IMGeostrophicTransform:NoActiveEndpoint", "Both endpoint accelerations are Inf. Use the APV transform directly when no endpoint is active.");
            end
            if options.zeroAPVModes.rotationName ~= "boundaryNormalized"
                error("IMGeostrophicTransform:CanonicalZeroAPVModesRequired", "zeroAPVModes must use canonical boundary-normalized coordinates. Supply rotated coordinates to individual transform methods.");
            end

            apvTransform = options.apvTransform;
            zeroAPVModes = options.zeroAPVModes;
            metadata = apvTransform.problemMetadata;
            if ~isstruct(metadata) || ~all(isfield(metadata,["g0","gd","surfaceBoundary"]))
                error("IMGeostrophicTransform:InvalidAPVTransform", "apvTransform must preserve g0, gd, and surfaceBoundary metadata from geostrophicAPVModes.");
            end
            g0Matches = isequal(metadata.g0,options.g0);
            gdMatches = isequal(metadata.gd,options.gd);
            surfaceConventionMatches = string(metadata.surfaceBoundary) == zeroAPVModes.surfaceBoundary;
            if ~g0Matches || ~gdMatches
                error("IMGeostrophicTransform:EndpointParameterMismatch", "The supplied g0 and gd must exactly match apvTransform.problemMetadata.");
            end
            if ~surfaceConventionMatches
                error("IMGeostrophicTransform:SurfaceConventionMismatch", "The APV and zero-APV surfaceBoundary conventions must match.");
            end

            canonicalEndpoints = ["surface","bottom"];
            activeEndpoints = canonicalEndpoints([isfinite(options.g0),isfinite(options.gd)]);
            activeEndpointsMatch = isequal(zeroAPVModes.endpoints,activeEndpoints);
            if ~activeEndpointsMatch
                error("IMGeostrophicTransform:ActiveEndpointMismatch", "zeroAPVModes.endpoints must equal the finite-parameter subset of surface and bottom in canonical order.");
            end

            compatibilityTolerance = sqrt(eps);
            domainMismatch = self.relativeMismatch(apvTransform.zDomain,zeroAPVModes.zDomain);
            gravityMismatch = self.relativeMismatch(apvTransform.g,zeroAPVModes.g);
            if domainMismatch > compatibilityTolerance
                error("IMGeostrophicTransform:DomainMismatch", "The APV and zero-APV vertical domains are incompatible.");
            end
            if gravityMismatch > compatibilityTolerance
                error("IMGeostrophicTransform:GravityMismatch", "The APV and zero-APV gravitational accelerations are incompatible.");
            end
            zeroN2Values = zeroAPVModes.N2(apvTransform.z);
            zeroN2Values = zeroN2Values(:);
            if numel(zeroN2Values) ~= numel(apvTransform.z) || any(~isfinite(zeroN2Values)) || any(zeroN2Values <= 0)
                error("IMGeostrophicTransform:InvalidStratification", "zeroAPVModes.N2 must return one finite positive value for every APV sample point.");
            end
            stratificationRelativeMismatch = self.relativeMismatch(apvTransform.N2Values,zeroN2Values);
            if stratificationRelativeMismatch > compatibilityTolerance
                error("IMGeostrophicTransform:StratificationMismatch", "The APV and zero-APV sampled stratifications differ by more than sqrt(eps) relatively.");
            end

            f0Value = zeroAPVModes.f0;
            gValue = apvTransform.g;
            kValues = reshape(zeroAPVModes.k,1,[]);
            inverseDepths = 1./reshape(apvTransform.h,[],1);
            frequencyTerm = (f0Value^2/gValue)*inverseDepths;
            mu = frequencyTerm + kValues.^2;
            muDenominator = abs(frequencyTerm) + kValues.^2;
            relativeMuSeparation = abs(mu)./muDenominator;
            if any(~isfinite(mu),"all") || any(~isfinite(relativeMuSeparation),"all")
                error("IMGeostrophicTransform:InvalidAPVEigendepth", "APV eigendepths must produce finite mu and relative-mu separation values.");
            end
            if any(relativeMuSeparation <= options.muTolerance,"all")
                error("IMGeostrophicTransform:NearSingularMu", "A retained APV mode has relativeMuSeparation <= muTolerance. Reduce the retained band or change the parameters.");
            end

            FEndpoints = apvTransform.endpointValues(variable="F");
            GEndpoints = apvTransform.endpointValues(variable="G");
            surfaceResponse = GEndpoints(1,:);
            if zeroAPVModes.surfaceBoundary == "freeSurface"
                surfaceResponse = surfaceResponse-FEndpoints(1,:);
            end
            fullEndpointResponse = [surfaceResponse;GEndpoints(2,:)];
            endpointRows = zeros(1,numel(activeEndpoints));
            endpointRows(activeEndpoints == "surface") = 1;
            endpointRows(activeEndpoints == "bottom") = 2;
            selectedEndpointResponse = fullEndpointResponse(endpointRows,:);
            nEndpoints = numel(activeEndpoints);
            nAPVModes = numel(apvTransform.modeNumber);
            nK = numel(kValues);
            apvEndpointResponse = zeros(nEndpoints,nAPVModes,nK);
            for iK = 1:nK
                responseScale = reshape((-f0Value/gValue)./mu(:,iK),1,[]);
                apvEndpointResponse(:,:,iK) = selectedEndpointResponse.*responseScale;
            end

            endpointZ = [apvTransform.zDomain(2);apvTransform.zDomain(1)];
            endpointTolerance = sqrt(eps)*max(1,max(abs(apvTransform.zDomain)));
            endpointIndices = zeros(1,2);
            for iEndpoint = 1:2
                index = find(abs(apvTransform.z-endpointZ(iEndpoint)) <= endpointTolerance,1);
                if ~isempty(index)
                    endpointIndices(iEndpoint) = index;
                end
            end
            requiredSourceEndpointSamples = [isfinite(options.g0) && options.g0 ~= 0,isfinite(options.gd) && options.gd ~= 0];

            zeroF = zeroAPVModes.F(apvTransform.z);
            zeroG = zeroAPVModes.G(apvTransform.z);
            zeroFEndpoints = zeroAPVModes.F(endpointZ);
            zeroGEndpoints = zeroAPVModes.G(endpointZ);
            zeroFPairingMatrices = zeros(nEndpoints,numel(apvTransform.z),nK);
            zeroGPairingMatrices = zeros(nEndpoints,numel(apvTransform.z),nK);
            FWeights = reshape(apvTransform.weights/apvTransform.depth,1,[]);
            GWeights = reshape(apvTransform.weights.*apvTransform.N2Values/gValue,1,[]);
            for iK = 1:nK
                zeroFPairingMatrices(:,:,iK) = zeroF(:,:,iK).'.*FWeights;
                zeroGPairingMatrices(:,:,iK) = zeroG(:,:,iK).'.*GWeights;
                if isfinite(options.g0) && options.g0 ~= 0 && endpointIndices(1) > 0
                    zeroSurfaceResponse = zeroGEndpoints(1,:,iK);
                    if zeroAPVModes.surfaceBoundary == "freeSurface"
                        zeroSurfaceResponse = zeroSurfaceResponse-zeroFEndpoints(1,:,iK);
                    end
                    zeroGPairingMatrices(:,endpointIndices(1),iK) = zeroGPairingMatrices(:,endpointIndices(1),iK) + (options.g0/gValue)*zeroSurfaceResponse.';
                end
                if isfinite(options.gd) && options.gd ~= 0 && endpointIndices(2) > 0
                    zeroGPairingMatrices(:,endpointIndices(2),iK) = zeroGPairingMatrices(:,endpointIndices(2),iK) + (options.gd/gValue)*zeroGEndpoints(2,:,iK).';
                end
            end

            matrixG0 = options.g0;
            matrixGd = options.gd;
            if isinf(matrixG0)
                matrixG0 = 0;
            end
            if isinf(matrixGd)
                matrixGd = 0;
            end
            generalizedEnergy = zeroAPVModes.generalizedEnergyMatrix(g0=matrixG0,gd=matrixGd);
            zeroSourceSolveMatrices = zeros(nEndpoints,nEndpoints,nK);
            zeroAPVGramReciprocalCondition = zeros(1,nK);
            zeroAPVGramRelativeSeparation = zeros(1,nK);
            zeroGramTolerance = sqrt(eps);
            for iK = 1:nK
                normalizedGram = (2*kValues(iK)^2/apvTransform.depth)*generalizedEnergy(:,:,iK);
                formScale = norm(zeroAPVModes.energyMatrix(:,:,iK),2) ...
                    + abs(matrixG0)*norm(zeroAPVModes.surfaceBuoyancyMatrix(:,:,iK),2) ...
                    + abs(matrixGd)*norm(zeroAPVModes.bottomBuoyancyMatrix(:,:,iK),2);
                normalizedScale = (2*kValues(iK)^2/apvTransform.depth)*formScale;
                singularValues = svd(normalizedGram);
                zeroAPVGramReciprocalCondition(iK) = rcond(normalizedGram);
                zeroAPVGramRelativeSeparation(iK) = min(singularValues)/max(normalizedScale,realmin);
                if ~isfinite(zeroAPVGramReciprocalCondition(iK)) || zeroAPVGramReciprocalCondition(iK) <= zeroGramTolerance ...
                        || ~isfinite(zeroAPVGramRelativeSeparation(iK)) || zeroAPVGramRelativeSeparation(iK) <= zeroGramTolerance
                    error("IMGeostrophicTransform:SingularZeroAPVGram", "The normalized zero-APV generalized-energy Gram matrix is singular or ill-conditioned on page %d.",iK);
                end
                zeroSourceSolveMatrices(:,:,iK) = normalizedGram\eye(nEndpoints);
            end

            self.k = kValues;
            self.activeEndpoints = activeEndpoints;
            self.apvEndpointResponse = apvEndpointResponse;
            self.g0 = options.g0;
            self.gd = options.gd;
            self.f0 = f0Value;
            self.apvTransform = apvTransform;
            self.zeroAPVModes = zeroAPVModes;
            self.g = gValue;
            self.depth = apvTransform.depth;
            self.endpointIndices = endpointIndices;
            self.requiredSourceEndpointSamples = requiredSourceEndpointSamples;
            self.zeroFPairingMatrices = zeroFPairingMatrices;
            self.zeroGPairingMatrices = zeroGPairingMatrices;
            self.zeroSourceSolveMatrices = zeroSourceSolveMatrices;
            self.compatibilityDiagnostics = struct( ...
                stratificationRelativeMismatch=stratificationRelativeMismatch, ...
                domainRelativeMismatch=domainMismatch, ...
                gravityRelativeMismatch=gravityMismatch, ...
                endpointParameterMatches=struct(g0=g0Matches,gd=gdMatches), ...
                surfaceConventionMatches=surfaceConventionMatches, ...
                activeEndpointsMatch=activeEndpointsMatch, ...
                mu=mu, ...
                relativeMuSeparation=relativeMuSeparation, ...
                minimumRelativeMuSeparation=min(relativeMuSeparation,[],"all"), ...
                muTolerance=options.muTolerance, ...
                zeroAPVGramReciprocalCondition=zeroAPVGramReciprocalCondition, ...
                zeroAPVGramRelativeSeparation=zeroAPVGramRelativeSeparation, ...
                zeroAPVGramTolerance=zeroGramTolerance, ...
                sourceEndpointSamplesPresent=endpointIndices > 0);
        end

        function [APVCoefficients,zeroAPVCoefficients] = transformStateForward(self, options)
            % Transform an admissible APV and endpoint-anomaly state.
            %
            % $$
            % \mathbf A_0^\kappa=-\frac{g\kappa^2}{f_0}
            % (\mathbf b^\kappa-\mathsf R_q^\kappa\mathbf A_q^\kappa).
            % $$
            %
            % - Topic: Transform admissible states
            % - Declaration: [APVCoefficients,zeroAPVCoefficients] = transformStateForward(transform,options)
            % - Parameter options.APV: sampled APV pages
            % - Parameter options.endpointAnomalies: endpoint-anomaly pages
            % - Parameter options.zeroAPVCoordinates: output zero-APV coordinates
            % - Returns APVCoefficients: APV coefficient pages
            % - Returns zeroAPVCoefficients: zero-APV coefficient pages
            arguments
                self IMGeostrophicTransform
                options.APV double
                options.endpointAnomalies double
                options.zeroAPVCoordinates (1,1) IMGeostrophicZeroAPVModesBasis = self.zeroAPVModes
            end

            self.requireAPVChannel("F","transformStateForward");
            [APVPages,trailingShape] = self.reshapePages(options.APV,numel(self.apvTransform.z),"APV");
            [endpointPages,endpointTrailingShape] = self.reshapePages(options.endpointAnomalies,numel(self.activeEndpoints),"endpointAnomalies");
            self.requireMatchingTrailingShape(trailingShape,endpointTrailingShape,"APV","endpointAnomalies");
            coordinates = self.requireCoordinateBasis(options.zeroAPVCoordinates);
            nFields = size(APVPages,3);
            nAPVModes = numel(self.apvTransform.modeNumber);
            nK = numel(self.k);

            APVFlat = reshape(APVPages,numel(self.apvTransform.z),[]);
            APVCoefficients = self.apvTransform.transformForward(APVFlat,variable="F");
            APVCoefficientPages = reshape(APVCoefficients,nAPVModes,nK,nFields);
            canonicalZeroPages = zeros(numel(self.activeEndpoints),nK,nFields,"like",APVCoefficientPages+endpointPages(1,1,1));
            for iK = 1:nK
                Aq = reshape(APVCoefficientPages(:,iK,:),nAPVModes,nFields);
                b = reshape(endpointPages(:,iK,:),numel(self.activeEndpoints),nFields);
                canonicalZeroPages(:,iK,:) = reshape((-self.g*self.k(iK)^2/self.f0)*(b-self.apvEndpointResponse(:,:,iK)*Aq),numel(self.activeEndpoints),1,nFields);
            end
            zeroPages = self.coordinatesFromCanonical(canonicalZeroPages,coordinates);
            APVCoefficients = self.restorePages(APVCoefficientPages,nAPVModes,trailingShape);
            zeroAPVCoefficients = self.restorePages(zeroPages,numel(self.activeEndpoints),trailingShape);
        end

        function [APV,endpointAnomalies] = transformStateBack(self, options)
            % Reconstruct sampled APV and active endpoint anomalies.
            %
            % $$
            % q=A_{\mathrm i}^{F}\mathbf A_q,\qquad
            % \mathbf b^\kappa=\mathsf R_q^\kappa\mathbf A_q^\kappa
            % -\frac{f_0}{g\kappa^2}\mathbf A_0^\kappa.
            % $$
            %
            % - Topic: Transform admissible states
            % - Declaration: [APV,endpointAnomalies] = transformStateBack(transform,options)
            % - Parameter options.APVCoefficients: APV coefficient pages
            % - Parameter options.zeroAPVCoefficients: zero-APV coefficient pages
            % - Parameter options.zeroAPVCoordinates: input zero-APV coordinates
            % - Returns APV: sampled APV pages
            % - Returns endpointAnomalies: endpoint-anomaly pages
            arguments
                self IMGeostrophicTransform
                options.APVCoefficients double
                options.zeroAPVCoefficients double
                options.zeroAPVCoordinates (1,1) IMGeostrophicZeroAPVModesBasis = self.zeroAPVModes
            end

            nAPVModes = numel(self.apvTransform.modeNumber);
            [APVCoefficientPages,trailingShape] = self.reshapePages(options.APVCoefficients,nAPVModes,"APVCoefficients");
            [zeroPages,zeroTrailingShape] = self.reshapePages(options.zeroAPVCoefficients,numel(self.activeEndpoints),"zeroAPVCoefficients");
            self.requireMatchingTrailingShape(trailingShape,zeroTrailingShape,"APVCoefficients","zeroAPVCoefficients");
            coordinates = self.requireCoordinateBasis(options.zeroAPVCoordinates);
            canonicalZeroPages = self.coordinatesToCanonical(zeroPages,coordinates);
            nFields = size(APVCoefficientPages,3);
            nK = numel(self.k);
            endpointPages = zeros(numel(self.activeEndpoints),nK,nFields,"like",APVCoefficientPages+canonicalZeroPages(1,1,1));
            for iK = 1:nK
                Aq = reshape(APVCoefficientPages(:,iK,:),nAPVModes,nFields);
                A0 = reshape(canonicalZeroPages(:,iK,:),numel(self.activeEndpoints),nFields);
                endpointPages(:,iK,:) = reshape(self.apvEndpointResponse(:,:,iK)*Aq-(self.f0/(self.g*self.k(iK)^2))*A0,numel(self.activeEndpoints),1,nFields);
            end

            APVCoefficientsFlat = reshape(APVCoefficientPages,nAPVModes,[]);
            APV = self.apvTransform.transformBack(APVCoefficientsFlat,variable="F");
            APV = self.restorePages(reshape(APV,numel(self.apvTransform.z),nK,nFields),numel(self.apvTransform.z),trailingShape);
            endpointAnomalies = self.restorePages(endpointPages,numel(self.activeEndpoints),trailingShape);
        end

        function [APVSourceCoefficients,zeroAPVSourceCoefficients] = transformSourceForward(self, options)
            % Project generic vorticity and displacement sources.
            %
            % This method uses the APV transform's
            % `modeProjectionFunctional` operation rather than
            % `transformForward`. The source equations require raw modal
            % pairings before the sampled Gram solve.
            %
            % $$
            % S_q^j=\frac{1}{D}\int F_jS_\omega\,dz
            % -\frac{f_0}{D}\mathcal G_q^j[S_\eta].
            % $$
            %
            % Zero-APV source coefficients satisfy
            %
            % $$
            % \widehat{\mathsf H}_g^\kappa\mathbf S_0^\kappa
            % =\mathbf p_0^\kappa,\qquad
            % \widehat{\mathsf H}_g^\kappa=\frac{2\kappa^2}{D}\mathsf H_g^\kappa.
            % $$
            %
            % - Topic: Project generic sources
            % - Declaration: [APVSourceCoefficients,zeroAPVSourceCoefficients] = transformSourceForward(transform,options)
            % - Parameter options.vorticitySource: sampled vorticity-source pages
            % - Parameter options.displacementSource: sampled displacement-source pages
            % - Parameter options.zeroAPVCoordinates: output zero-APV coordinates
            % - Returns APVSourceCoefficients: APV source-coefficient pages
            % - Returns zeroAPVSourceCoefficients: zero-APV source-coefficient pages
            arguments
                self IMGeostrophicTransform
                options.vorticitySource double
                options.displacementSource double
                options.zeroAPVCoordinates (1,1) IMGeostrophicZeroAPVModesBasis = self.zeroAPVModes
            end

            self.requireAPVChannel("F","transformSourceForward");
            self.requireAPVChannel("G","transformSourceForward");
            missingEndpoints = self.requiredSourceEndpointSamples & self.endpointIndices == 0;
            if any(missingEndpoints)
                endpoints = ["surface","bottom"];
                error("IMGeostrophicTransform:MissingSourceEndpointSample", "transformSourceForward requires these endpoint samples in apvTransform.z: %s.",join(endpoints(missingEndpoints),", "));
            end
            [vorticityPages,trailingShape] = self.reshapePages(options.vorticitySource,numel(self.apvTransform.z),"vorticitySource");
            [displacementPages,displacementTrailingShape] = self.reshapePages(options.displacementSource,numel(self.apvTransform.z),"displacementSource");
            self.requireMatchingTrailingShape(trailingShape,displacementTrailingShape,"vorticitySource","displacementSource");
            coordinates = self.requireCoordinateBasis(options.zeroAPVCoordinates);
            nFields = size(vorticityPages,3);
            nK = numel(self.k);
            nAPVModes = numel(self.apvTransform.modeNumber);

            vorticityFlat = reshape(vorticityPages,numel(self.apvTransform.z),[]);
            displacementFlat = reshape(displacementPages,numel(self.apvTransform.z),[]);
            FPairings = self.apvTransform.modeProjectionFunctional(vorticityFlat,variable="F")/self.depth;
            GPairings = self.apvTransform.modeProjectionFunctional(displacementFlat,variable="G");
            APVSourcePages = reshape(FPairings-(self.f0/self.depth)*GPairings,nAPVModes,nK,nFields);

            canonicalZeroPages = zeros(numel(self.activeEndpoints),nK,nFields,"like",vorticityPages+displacementPages(1,1,1));
            for iK = 1:nK
                vorticity = reshape(vorticityPages(:,iK,:),numel(self.apvTransform.z),nFields);
                displacement = reshape(displacementPages(:,iK,:),numel(self.apvTransform.z),nFields);
                zeroPairings = self.zeroFPairingMatrices(:,:,iK)*vorticity-(self.f0/self.depth)*(self.zeroGPairingMatrices(:,:,iK)*displacement);
                canonicalZeroPages(:,iK,:) = reshape(self.zeroSourceSolveMatrices(:,:,iK)*zeroPairings,numel(self.activeEndpoints),1,nFields);
            end
            zeroPages = self.coordinatesFromCanonical(canonicalZeroPages,coordinates);
            APVSourceCoefficients = self.restorePages(APVSourcePages,nAPVModes,trailingShape);
            zeroAPVSourceCoefficients = self.restorePages(zeroPages,numel(self.activeEndpoints),trailingShape);
        end
    end

    methods (Access = private)
        function requireAPVChannel(self, variable, operation)
            if ~self.apvTransform.hasForwardTransform(variable=variable)
                reason = self.apvTransform.forwardTransformReason(variable=variable);
                error("IMGeostrophicTransform:UnavailableAPVChannel", "%s requires the APV %s forward channel. %s",operation,variable,reason);
            end
        end

        function coordinates = requireCoordinateBasis(self, coordinates)
            if ~coordinates.sharesCanonicalBasisWith(self.zeroAPVModes)
                error("IMGeostrophicTransform:IncompatibleZeroAPVCoordinates", "zeroAPVCoordinates must derive from the canonical zeroAPVModes supplied at construction.");
            end
        end

        function coordinatePages = coordinatesFromCanonical(self, canonicalPages, coordinates)
            coordinatePages = zeros(size(canonicalPages),"like",canonicalPages);
            nFields = size(canonicalPages,3);
            for iK = 1:numel(self.k)
                canonical = reshape(canonicalPages(:,iK,:),numel(self.activeEndpoints),nFields);
                coordinatePages(:,iK,:) = reshape(coordinates.rotationMatrix(:,:,iK)\canonical,numel(self.activeEndpoints),1,nFields);
            end
        end

        function canonicalPages = coordinatesToCanonical(self, coordinatePages, coordinates)
            canonicalPages = zeros(size(coordinatePages),"like",coordinatePages);
            nFields = size(coordinatePages,3);
            for iK = 1:numel(self.k)
                coordinatesAtK = reshape(coordinatePages(:,iK,:),numel(self.activeEndpoints),nFields);
                canonicalPages(:,iK,:) = reshape(coordinates.rotationMatrix(:,:,iK)*coordinatesAtK,numel(self.activeEndpoints),1,nFields);
            end
        end

        function [pages,trailingShape] = reshapePages(self, values, nLeading, name)
            if isempty(values) || size(values,1) ~= nLeading || size(values,2) ~= numel(self.k)
                error("IMGeostrophicTransform:InvalidArrayShape", "%s must have leading dimensions %d x %d.",name,nLeading,numel(self.k));
            end
            valueSize = size(values);
            trailingShape = valueSize(3:end);
            pages = reshape(values,nLeading,numel(self.k),prod(trailingShape));
        end

        function values = restorePages(self, pages, nLeading, trailingShape)
            if isempty(trailingShape)
                values = reshape(pages,nLeading,numel(self.k));
            else
                values = reshape(pages,[nLeading,numel(self.k),trailingShape]);
            end
        end

        function requireMatchingTrailingShape(~, firstShape, secondShape, firstName, secondName)
            if ~isequal(firstShape,secondShape)
                error("IMGeostrophicTransform:TrailingDimensionMismatch", "%s and %s must have identical trailing dimensions.",firstName,secondName);
            end
        end
    end

    methods (Static, Access = private)
        function validateEndpointAcceleration(value, endpoint)
            if isnan(value) || value == -Inf
                error("IMGeostrophicTransform:InvalidEndpointAcceleration", "%s acceleration must be signed finite, zero, or positive Inf.",endpoint);
            end
        end

        function mismatch = relativeMismatch(first, second)
            if ~isequal(size(first),size(second))
                mismatch = Inf;
                return;
            end
            scale = max([abs(first(:));abs(second(:))]);
            if isempty(scale) || scale == 0
                scale = 1;
            end
            mismatch = max(abs(first(:)-second(:)))/scale;
        end
    end
end
