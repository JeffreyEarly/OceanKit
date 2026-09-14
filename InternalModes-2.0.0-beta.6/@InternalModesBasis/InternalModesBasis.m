classdef InternalModesBasis < CAAnnotatedClass
    % Store solved vertical modes together with their component roles.
    %
    % `InternalModesBasis` is the canonical container for a solved vertical
    % mode basis. It records the sampled inverse modes $$F_j(z)$$ and
    % $$G_j(z)$$, the equivalent depths $$h_j$$, the quadrature grid, and
    % whether each component has a Sturm-Liouville forward projection.
    %
    % Geostrophic bases at $$\omega=0$$ have canonical F and G projections,
    %
    % $$
    % \mathcal{F}_g^j[u] = \gamma_j^{-1}\int_{-D}^{0} u F_g^j\,dz,
    % \qquad
    % \mathcal{G}_g^j[\eta] = \frac{1}{g}\int_{-D}^{0} N^2 \eta G_g^j\,dz.
    % $$
    %
    % Nonzero-$$\kappa$$ IGW bases have a canonical G projection,
    %
    % $$
    % \mathcal{G}_\kappa^j[\eta] =
    % \frac{1}{g}\int_{-D}^{0} \left(N^2-f_0^2\right)\eta G_\kappa^j\,dz,
    % $$
    %
    % but no independent canonical wave-F projection.
    %
    % ```matlab
    % basis = InternalModesBasis.fromSolverAtFrequency(im,0,nModes=32);
    % transform = basis.nativeTransform(component="G");
    % ```
    %
    % - Topic: Create vertical bases
    % - Topic: Inspect basis properties
    % - Topic: Build vertical transforms
    % - Topic: Analyze vertical spectra
    % - Topic: Persist vertical bases
    % - Topic: Developer topics
    % - Declaration: classdef InternalModesBasis < CAAnnotatedClass

    properties (SetAccess = private)
        % Depth grid where the inverse modes are sampled.
        %
        % The grid is a column vector in meters. Projection operators built
        % from this basis act on scalar fields sampled at these depths unless
        % an observation projection supplies its own observation operator.
        %
        % - Topic: Inspect basis properties
        z

        % Row coordinate used by annotated NetCDF persistence.
        %
        % - Topic: Developer topics
        % - Developer: true
        zIndex

        % F-mode index coordinate used by annotated NetCDF persistence.
        %
        % - Topic: Developer topics
        % - Developer: true
        modeF

        % G-mode index coordinate used by annotated NetCDF persistence.
        %
        % - Topic: Developer topics
        % - Developer: true
        modeG

        % Sampled F inverse modes with rows matching `z`.
        %
        % The columns are the solver-returned $$F_j(z)$$ modes. For
        % geostrophic transforms, methods may prepend the barotropic mode
        % $$F_0=1$$ when building an F transform.
        %
        % - Topic: Inspect basis properties
        F

        % Sampled G inverse modes with rows matching `z`.
        %
        % The columns are the solver-returned $$G_j(z)$$ modes.
        %
        % - Topic: Inspect basis properties
        G

        % Equivalent depths associated with the stored mode columns.
        %
        % - Topic: Inspect basis properties
        h

        % Buoyancy frequency squared sampled at `z`.
        %
        % The canonical G projection uses either $$N^2/g$$ or
        % $$(N^2-f_0^2)/g$$ as its vertical weight, depending on component role.
        %
        % - Topic: Inspect basis properties
        N2

        % Coriolis parameter used by the mode solve.
        %
        % - Topic: Inspect basis properties
        f0

        % Gravitational acceleration used by the mode solve.
        %
        % - Topic: Inspect basis properties
        g

        % Water-column depth $$D=z_{\max}-z_{\min}$$.
        %
        % The barotropic geostrophic F spectral weight is $$\gamma_0=D$$.
        %
        % - Topic: Inspect basis properties
        D

        % Horizontal wavenumber associated with the solve.
        %
        % - Topic: Inspect basis properties
        kappa

        % Frequency associated with the solve.
        %
        % - Topic: Inspect basis properties
        omega

        % Text label identifying the vertical eigenvalue problem.
        %
        % Common values are `"geostrophic"`, `"igwWavenumber"`, and
        % `"fixedFrequency"`.
        %
        % - Topic: Inspect basis properties
        problemType

        % Text label identifying the source solver or factory path.
        %
        % - Topic: Inspect basis properties
        sourceDescription

        % Component role label for the F component.
        %
        % Values such as `"eigenfunction"` and `"diagnostic"` explain
        % whether a canonical F projection exists.
        %
        % - Topic: Inspect basis properties
        componentRoleF

        % Component role label for the G component.
        %
        % - Topic: Inspect basis properties
        componentRoleG

        % Boolean flag indicating whether canonical F projection exists.
        %
        % For nonzero-$$\kappa$$ IGW bases this is false, because wave F
        % modes are diagnostic/evaluation modes rather than a separate
        % Sturm-Liouville projection basis.
        %
        % - Topic: Inspect basis properties
        forwardProjectionAvailableF

        % Boolean flag indicating whether canonical G projection exists.
        %
        % - Topic: Inspect basis properties
        forwardProjectionAvailableG

        % Orthogonality weight label for F projections.
        %
        % - Topic: Inspect basis properties
        orthogonalityWeightF

        % Orthogonality weight label for G projections.
        %
        % - Topic: Inspect basis properties
        orthogonalityWeightG
    end

    methods
        function self = InternalModesBasis(options)
            % Create a persisted vertical basis from canonical state.
            %
            % This constructor is intentionally a cheap state constructor for
            % annotated persistence. Scientific setup should usually use
            % `fromSolverAtFrequency`, `fromSolverAtWavenumber`, or
            % `fromSolvedModes`.
            %
            % - Topic: Create vertical bases
            % - Declaration: basis = InternalModesBasis(options)
            % - Parameter options.z: depth grid where modes are sampled
            % - Parameter options.F: sampled F inverse modes
            % - Parameter options.G: sampled G inverse modes
            % - Parameter options.h: equivalent-depth vector
            % - Parameter options.N2: buoyancy frequency squared sampled at `z`
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.D: water-column depth
            % - Parameter options.kappa: horizontal wavenumber metadata
            % - Parameter options.omega: frequency metadata
            % - Parameter options.problemType: text label for the EVP
            % - Parameter options.sourceDescription: text label for the source solver
            % - Parameter options.componentRoleF: F component role label
            % - Parameter options.componentRoleG: G component role label
            % - Parameter options.forwardProjectionAvailableF: true if F has a canonical projection
            % - Parameter options.forwardProjectionAvailableG: true if G has a canonical projection
            % - Parameter options.orthogonalityWeightF: F inner-product weight label
            % - Parameter options.orthogonalityWeightG: G inner-product weight label
            % - Returns basis: initialized InternalModesBasis instance
            arguments
                options.z (:,1) double = zeros(0,1)
                options.F double = zeros(0,0)
                options.G double = zeros(0,0)
                options.h double = zeros(1,0)
                options.N2 (:,1) double = zeros(0,1)
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
                options.D (1,1) double {mustBeNonnegative} = 0
                options.kappa (1,1) double = NaN
                options.omega (1,1) double = NaN
                options.problemType {mustBeTextScalar} = "unknown"
                options.sourceDescription {mustBeTextScalar} = "manual"
                options.componentRoleF {mustBeTextScalar} = "diagnostic"
                options.componentRoleG {mustBeTextScalar} = "eigenfunction"
                options.forwardProjectionAvailableF (1,1) logical = false
                options.forwardProjectionAvailableG (1,1) logical = true
                options.orthogonalityWeightF {mustBeTextScalar} = "none"
                options.orthogonalityWeightG {mustBeTextScalar} = "unknown"
            end

            self@CAAnnotatedClass();
            self.z = options.z(:);
            self.zIndex = (1:length(self.z)).';
            self.F = options.F;
            self.G = options.G;
            self.h = reshape(options.h,1,[]);
            self.N2 = options.N2(:);
            self.f0 = options.f0;
            self.g = options.g;
            self.D = options.D;
            self.kappa = options.kappa;
            self.omega = options.omega;
            self.problemType = string(options.problemType);
            self.sourceDescription = string(options.sourceDescription);
            self.componentRoleF = string(options.componentRoleF);
            self.componentRoleG = string(options.componentRoleG);
            self.forwardProjectionAvailableF = options.forwardProjectionAvailableF;
            self.forwardProjectionAvailableG = options.forwardProjectionAvailableG;
            self.orthogonalityWeightF = string(options.orthogonalityWeightF);
            self.orthogonalityWeightG = string(options.orthogonalityWeightG);
            self.modeF = (1:size(self.F,2)).';
            self.modeG = (1:size(self.G,2)).';
            self.validateState();
        end

        function transform = nativeTransform(self, options)
            % Build a transform on the basis' native vertical grid.
            %
            % The native transform maps between modal coefficients and
            % fields sampled at `basis.z`. Components with Dirichlet rows,
            % such as rigid-lid G modes, are solved on their active interior
            % rows and expanded back to full-grid matrices. With
            % `projectionMethod="auto"`, a square, well-conditioned active
            % inverse matrix uses a direct inverse; otherwise the transform
            % uses a weighted pseudoinverse.
            %
            % - Topic: Build vertical transforms
            % - Declaration: transform = nativeTransform(self,options)
            % - Parameter self: InternalModesBasis instance
            % - Parameter options.component: `"F"`, `"G"`, or `"both"`
            % - Parameter options.nModes: maximum number of modes retained
            % - Parameter options.projectionMethod: `"auto"`, `"directInverse"`, `"weightedPseudoinverse"`, or `"canonical"`
            % - Parameter options.allowNoncanonical: true to allow numerical wave-F projections
            % - Parameter options.maxConditionNumber: condition-number limit for direct inverses
            % - Returns transform: InternalModesTransform on the native grid
            arguments
                self InternalModesBasis
                options.component {mustBeTextScalar} = "both"
                options.nModes double = []
                options.projectionMethod {mustBeTextScalar} = "auto"
                options.allowNoncanonical (1,1) logical = false
                options.maxConditionNumber (1,1) double {mustBePositive} = 1e6
            end

            transform = self.buildTransform(component=options.component,nModes=options.nModes, ...
                projectionMethod=options.projectionMethod,allowNoncanonical=options.allowNoncanonical, ...
                maxConditionNumber=options.maxConditionNumber);
        end

        function transform = modelTransform(self, options)
            % Build a prefix-retained transform for model vertical grids.
            %
            % Projection resolvability and nonlinear de-aliasing are treated
            % as separate limits. If `nonlinearAliasingPolicy` is
            % `"quadratic"`, the transform automatically applies the
            % two-thirds retained-mode cap and records that diagnostic
            % separately from projection leakage.
            %
            % $$
            % n_{\mathrm{retained}} =
            % \min(n_{\mathrm{requested}},n_{\mathrm{projection}},
            % n_{\mathrm{nonlinear}}).
            % $$
            %
            % - Topic: Build vertical transforms
            % - Declaration: transform = modelTransform(self,options)
            % - Parameter self: InternalModesBasis instance
            % - Parameter options.component: `"F"`, `"G"`, or `"both"`
            % - Parameter options.nModes: requested retained mode count
            % - Parameter options.projectionTolerance: leakage tolerance for prefix selection
            % - Parameter options.maxConditionNumber: condition-number limit for prefix selection
            % - Parameter options.nTailCheck: number of rejected tail modes checked for leakage
            % - Parameter options.nonlinearAliasingPolicy: `"none"` or `"quadratic"`
            % - Parameter options.allowNoncanonical: true to allow numerical wave-F projections
            % - Returns transform: InternalModesTransform with prefix-retained modes
            arguments
                self InternalModesBasis
                options.component {mustBeTextScalar} = "both"
                options.nModes double = []
                options.projectionTolerance (1,1) double {mustBePositive} = 1e-2
                options.maxConditionNumber (1,1) double {mustBePositive} = 1e6
                options.nTailCheck (1,1) double {mustBeInteger,mustBeNonnegative} = 8
                options.nonlinearAliasingPolicy {mustBeMember(options.nonlinearAliasingPolicy,["none","quadratic"])} = "none"
                options.allowNoncanonical (1,1) logical = true
            end

            requested = self.requestedModeCount(options.component, options.nModes);
            projectionResolved = self.resolvedModeCount(component=options.component, ...
                nModes=requested,projectionTolerance=options.projectionTolerance, ...
                maxConditionNumber=options.maxConditionNumber,nTailCheck=options.nTailCheck);
            nonlinearLimit = self.nonlinearAliasLimit(options.component, requested, options.nonlinearAliasingPolicy);
            retained = min([requested projectionResolved nonlinearLimit]);

            transform = self.buildTransform(component=options.component,nModes=retained, ...
                projectionMethod="weightedPseudoinverse",allowNoncanonical=options.allowNoncanonical, ...
                maxConditionNumber=options.maxConditionNumber);
            transform = transform.withDiagnostics(selectionReason=self.selectionReason(requested, projectionResolved, nonlinearLimit, retained), ...
                projectionResolvedModes=projectionResolved,nonlinearAliasLimit=nonlinearLimit);
        end

        function transform = fixedGridTransform(self, referenceTransform, options)
            % Build a transform on another transform's vertical grid.
            %
            % This method evaluates the current basis on the reference grid
            % and builds a weighted pseudoinverse there. It is intended for
            % cases where a model fixes the hydrostatic quadrature grid but
            % needs nonzero-$$\kappa$$ IGW modes evaluated on that grid.
            %
            % - Topic: Build vertical transforms
            % - Declaration: transform = fixedGridTransform(self,referenceTransform,options)
            % - Parameter self: InternalModesBasis instance
            % - Parameter referenceTransform: InternalModesTransform whose grid and weights define the fixed grid
            % - Parameter options.component: `"F"`, `"G"`, or `"both"`
            % - Parameter options.projectionTolerance: leakage tolerance for retained prefix modes
            % - Parameter options.maxConditionNumber: condition-number limit
            % - Parameter options.nTailCheck: number of rejected tail modes checked for leakage
            % - Parameter options.preserveSize: true to keep rejected rows as zeros
            % - Returns transform: InternalModesTransform evaluated on the fixed grid
            arguments
                self InternalModesBasis
                referenceTransform InternalModesTransform
                options.component {mustBeTextScalar} = "both"
                options.projectionTolerance (1,1) double {mustBePositive} = 1e-2
                options.maxConditionNumber (1,1) double {mustBePositive} = 1e6
                options.nTailCheck (1,1) double {mustBeInteger,mustBeNonnegative} = 8
                options.preserveSize (1,1) logical = true
            end

            fixedBasis = self.interpolateToGrid(referenceTransform.z);
            requested = fixedBasis.requestedModeCount(options.component, []);
            projectionResolved = fixedBasis.resolvedModeCount(component=options.component, ...
                nModes=requested,projectionTolerance=options.projectionTolerance, ...
                maxConditionNumber=options.maxConditionNumber,nTailCheck=options.nTailCheck);
            retained = projectionResolved;
            transform = fixedBasis.buildTransform(component=options.component,nModes=retained, ...
                projectionMethod="weightedPseudoinverse",allowNoncanonical=true, ...
                maxConditionNumber=options.maxConditionNumber,preserveSize=options.preserveSize);
            transform = transform.withDiagnostics(selectionReason="projectionQualityLimit", ...
                projectionResolvedModes=projectionResolved,nonlinearAliasLimit=requested);
        end

        function projection = observationProjection(self, observation, options)
            % Build a vertical projection for an arbitrary observation grid.
            %
            % The observation problem uses the sampled matrix
            % $$B=H\Phi$$, where $$H$$ is either a point-sampling grid or an
            % explicit observation operator and $$\Phi$$ is the inverse mode
            % matrix. Rank-revealing QR selects a resolvable, possibly
            % non-contiguous set of modes.
            %
            % - Topic: Build vertical transforms
            % - Declaration: projection = observationProjection(self,observation,options)
            % - Parameter self: InternalModesBasis instance
            % - Parameter observation: observation depths or explicit observation matrix
            % - Parameter options.component: `"F"` or `"G"`
            % - Parameter options.nModes: number of candidate modes
            % - Parameter options.weights: observation weights
            % - Parameter options.rankTolerance: relative QR pivot tolerance
            % - Parameter options.maxConditionNumber: retained Gram condition-number limit
            % - Returns projection: InternalModesProjection for the observation operator
            arguments
                self InternalModesBasis
                observation double
                options.component {mustBeMember(options.component,["F","G"])} = "G"
                options.nModes double = []
                options.weights double = []
                options.rankTolerance (1,1) double {mustBePositive} = 5e-2
                options.maxConditionNumber (1,1) double {mustBePositive} = 10
            end

            component = string(options.component);
            [Phi,weights,modeNumbers,modeHeights,spectralWeights,componentRole,status] = self.componentMatrix(component, options.nModes);
            if isvector(observation)
                observationZ = observation(:);
                B = interp1(self.z, Phi, observationZ, "pchip", "extrap");
                observationMatrix = [];
            else
                observationZ = zeros(0,1);
                observationMatrix = observation;
                B = observationMatrix * Phi;
            end
            if isempty(options.weights)
                observationWeights = ones(size(B,1),1);
            else
                observationWeights = options.weights(:);
            end

            [retainedModes,selectionDiagnostics] = InternalModesProjection.selectResolvableModes(B, observationWeights, ...
                options.rankTolerance, options.maxConditionNumber);
            [projectionMatrix,resolutionMatrix,gramMatrix] = InternalModesProjection.weightedProjection(B, observationWeights, retainedModes);
            rejectedModes = setdiff(modeNumbers(:).', retainedModes);
            projection = InternalModesProjection(component=component,observationZ=observationZ, ...
                observationMatrix=observationMatrix,candidateModes=modeNumbers(:),retainedModes=retainedModes(:), ...
                rejectedModes=rejectedModes(:),projectionMatrix=projectionMatrix,reconstructionMatrix=Phi(:,retainedModes), ...
                resolutionMatrix=resolutionMatrix,aliasingMatrix=resolutionMatrix(:,rejectedModes), ...
                spectralWindow=abs(resolutionMatrix).^2,gramMatrix=gramMatrix,conditionNumber=cond(gramMatrix), ...
                spectralWeights=spectralWeights(:),transformStatus=status,componentRole=componentRole, ...
                isCanonicalComponent=self.isCanonicalComponent(component),diagnosticRankTolerance=options.rankTolerance, ...
                diagnosticFinalRelativePivot=selectionDiagnostics.finalRelativePivot);
        end

        function weights = spectralWeights(self, component, options)
            % Return Parseval weights for modal spectra.
            %
            % For canonical geostrophic F modes this returns
            % $$\gamma_0=D$$ and $$\gamma_j=h_g^j$$. For canonical G modes
            % the current normalization gives a modal spectrum weighted by
            % $$g$$. These Parseval spectrum weights are distinct from the
            % normalization used by canonical forward projections.
            %
            % - Topic: Analyze vertical spectra
            % - Declaration: weights = spectralWeights(self,component,options)
            % - Parameter self: InternalModesBasis instance
            % - Parameter component: `"F"` or `"G"`
            % - Parameter options.nModes: number of weights requested
            % - Returns weights: column vector of modal spectrum weights
            arguments
                self InternalModesBasis
                component {mustBeMember(component,["F","G"])}
                options.nModes double = []
            end

            [~,~,~,~,weights] = self.componentMatrix(string(component), options.nModes);
        end
    end

    methods (Access = private)
        function validateState(self)
            if ~isempty(self.F) && size(self.F,1) ~= length(self.z)
                error("InternalModesBasis:InvalidShape", "F must have one row for each z point.");
            end
            if ~isempty(self.G) && size(self.G,1) ~= length(self.z)
                error("InternalModesBasis:InvalidShape", "G must have one row for each z point.");
            end
            if ~isempty(self.N2) && length(self.N2) ~= length(self.z)
                error("InternalModesBasis:InvalidShape", "N2 must have one entry for each z point.");
            end
        end

        function transform = buildTransform(self, options)
            arguments
                self InternalModesBasis
                options.component {mustBeTextScalar} = "both"
                options.nModes double = []
                options.projectionMethod {mustBeTextScalar} = "auto"
                options.allowNoncanonical (1,1) logical = false
                options.maxConditionNumber (1,1) double {mustBePositive} = 1e6
                options.preserveSize (1,1) logical = false
            end

            component = string(options.component);
            includeF = component == "both" || component == "F";
            includeG = component == "both" || component == "G";

            forwardF = zeros(0,length(self.z));
            inverseF = zeros(length(self.z),0);
            retainedModesF = zeros(0,1);
            rejectedModesF = zeros(0,1);
            spectralWeightsF = zeros(0,1);
            transformStatusF = "diagnosticOnly";
            componentRoleF = self.componentRoleF;
            hF = zeros(0,1);
            weightsF = zeros(length(self.z),1);
            conditionNumberF = NaN;
            gramErrorF = NaN;

            if includeF
                if ~self.forwardProjectionAvailableF && ~options.allowNoncanonical
                    error("InternalModesBasis:NoncanonicalFProjection", ...
                        "This basis has no canonical F projection. Pass allowNoncanonical=true to build a numerical F transform.");
                end
                [PhiF,weightsF,modeNumbersF,hF,spectralWeightsF,componentRoleF,defaultStatusF,activeRowsF,projectionWeightsF] = self.componentMatrix("F", options.nModes);
                [forwardF,transformStatusF,conditionNumberF,gramErrorF] = self.forwardMatrix(PhiF, weightsF, projectionWeightsF, options.projectionMethod, defaultStatusF, options.maxConditionNumber, activeRowsF);
                inverseF = PhiF;
                inactiveRowsF = setdiff((1:length(self.z)).', activeRowsF);
                inverseF(inactiveRowsF,:) = 0;
                retainedModesF = modeNumbersF(:);
                rejectedModesF = setdiff((1:self.nAvailableModes("F")).', retainedModesF);
                if options.preserveSize
                    nAvailableF = self.nAvailableModes("F");
                    [~,~,~,~,spectralWeightsFAll] = self.componentMatrix("F", []);
                    forwardFPadded = zeros(nAvailableF,length(self.z));
                    inverseFPadded = zeros(length(self.z),nAvailableF);
                    forwardFPadded(retainedModesF,:) = forwardF;
                    inverseFPadded(:,retainedModesF) = inverseF;
                    forwardF = forwardFPadded;
                    inverseF = inverseFPadded;
                    spectralWeightsF = spectralWeightsFAll(:);
                end
            end

            forwardG = zeros(0,length(self.z));
            inverseG = zeros(length(self.z),0);
            retainedModesG = zeros(0,1);
            rejectedModesG = zeros(0,1);
            spectralWeightsG = zeros(0,1);
            transformStatusG = "diagnosticOnly";
            componentRoleG = self.componentRoleG;
            hG = zeros(0,1);
            weightsG = zeros(length(self.z),1);
            conditionNumberG = NaN;
            gramErrorG = NaN;

            if includeG
                [PhiG,weightsG,modeNumbersG,hG,spectralWeightsG,componentRoleG,defaultStatusG,activeRowsG,projectionWeightsG] = self.componentMatrix("G", options.nModes);
                [forwardG,transformStatusG,conditionNumberG,gramErrorG] = self.forwardMatrix(PhiG, weightsG, projectionWeightsG, options.projectionMethod, defaultStatusG, options.maxConditionNumber, activeRowsG);
                inverseG = PhiG;
                inactiveRowsG = setdiff((1:length(self.z)).', activeRowsG);
                inverseG(inactiveRowsG,:) = 0;
                retainedModesG = modeNumbersG(:);
                rejectedModesG = setdiff((1:self.nAvailableModes("G")).', retainedModesG);
                if options.preserveSize
                    nAvailableG = self.nAvailableModes("G");
                    [~,~,~,~,spectralWeightsGAll] = self.componentMatrix("G", []);
                    forwardGPadded = zeros(nAvailableG,length(self.z));
                    inverseGPadded = zeros(length(self.z),nAvailableG);
                    forwardGPadded(retainedModesG,:) = forwardG;
                    inverseGPadded(:,retainedModesG) = inverseG;
                    forwardG = forwardGPadded;
                    inverseG = inverseGPadded;
                    spectralWeightsG = spectralWeightsGAll(:);
                end
            end

            transform = InternalModesTransform(z=self.z,N2=self.N2,f0=self.f0,g=self.g,D=self.D,kappa=self.kappa,omega=self.omega, ...
                forwardF=forwardF,inverseF=inverseF,forwardG=forwardG,inverseG=inverseG, ...
                weightsF=weightsF,weightsG=weightsG,spectralWeightsF=spectralWeightsF,spectralWeightsG=spectralWeightsG, ...
                retainedModesF=retainedModesF,rejectedModesF=rejectedModesF,retainedModesG=retainedModesG,rejectedModesG=rejectedModesG, ...
                hF=hF,hG=hG,componentRoleF=componentRoleF,componentRoleG=componentRoleG, ...
                transformStatusF=transformStatusF,transformStatusG=transformStatusG, ...
                forwardProjectionAvailableF=self.forwardProjectionAvailableF,forwardProjectionAvailableG=self.forwardProjectionAvailableG, ...
                conditionNumberF=conditionNumberF,conditionNumberG=conditionNumberG,gramErrorF=gramErrorF,gramErrorG=gramErrorG, ...
                problemType=self.problemType,sourceDescription=self.sourceDescription);
        end

        function [Phi,weights,modeNumbers,modeHeights,spectrumWeights,componentRole,status,activeRows,projectionWeights] = componentMatrix(self, component, nModes)
            component = string(component);
            nAvailable = self.nAvailableModes(component);

            if component == "F"
                if self.isGeostrophic()
                    PhiAll = cat(2,ones(length(self.z),1),self.F);
                    hAll = [self.D reshape(self.h,1,[])];
                    spectrumAll = hAll;
                    projectionAll = hAll;
                else
                    PhiAll = self.F;
                    hAll = reshape(self.h,1,[]);
                    spectrumAll = hAll;
                    projectionAll = hAll;
                end
                componentRole = self.componentRoleF;
                status = self.defaultStatusForComponent("F");
                weights = self.quadratureWeights();
                activeRows = (1:size(PhiAll,1)).';
            else
                PhiAll = self.G;
                hAll = reshape(self.h,1,[]);
                spectrumAll = self.g * ones(1,size(PhiAll,2));
                projectionAll = ones(1,size(PhiAll,2));
                componentRole = self.componentRoleG;
                status = self.defaultStatusForComponent("G");
                weights = self.gQuadratureWeights();
                activeRows = self.activeProjectionRows(PhiAll);
            end

            nTransformable = min(nAvailable, length(activeRows));
            if isempty(nModes)
                nRetained = nTransformable;
            else
                nRetained = min(nTransformable, nModes);
            end

            Phi = PhiAll(:,1:nRetained);
            modeNumbers = (1:nRetained).';
            modeHeights = reshape(hAll(1:nRetained),[],1);
            spectrumWeights = reshape(spectrumAll(1:nRetained),[],1);
            projectionWeights = reshape(projectionAll(1:nRetained),[],1);
        end

        function n = nAvailableModes(self, component)
            component = string(component);
            if component == "F"
                n = size(self.F,2) + double(self.isGeostrophic());
            else
                n = size(self.G,2);
            end
        end

        function n = requestedModeCount(self, component, nModes)
            component = string(component);
            if component == "both"
                nAvailable = min(self.nAvailableModes("F"), self.nAvailableModes("G"));
            else
                nAvailable = self.nAvailableModes(component);
            end
            if isempty(nModes)
                n = nAvailable;
            else
                n = min(nAvailable, nModes);
            end
        end

        function nResolved = resolvedModeCount(self, options)
            arguments
                self InternalModesBasis
                options.component {mustBeTextScalar} = "G"
                options.nModes (1,1) double {mustBeInteger,mustBePositive}
                options.projectionTolerance (1,1) double {mustBePositive}
                options.maxConditionNumber (1,1) double {mustBePositive}
                options.nTailCheck (1,1) double {mustBeInteger,mustBeNonnegative}
            end

            component = string(options.component);
            if component == "both"
                if self.problemType == "igwWavenumber"
                    nResolved = self.resolvedModeCount(component="G",nModes=options.nModes,projectionTolerance=options.projectionTolerance, ...
                        maxConditionNumber=options.maxConditionNumber,nTailCheck=options.nTailCheck);
                    return;
                end
                nResolved = min([self.resolvedModeCount(component="F",nModes=options.nModes,projectionTolerance=options.projectionTolerance, ...
                    maxConditionNumber=options.maxConditionNumber,nTailCheck=options.nTailCheck), ...
                    self.resolvedModeCount(component="G",nModes=options.nModes,projectionTolerance=options.projectionTolerance, ...
                    maxConditionNumber=options.maxConditionNumber,nTailCheck=options.nTailCheck)]);
                return;
            end

            [Phi,weights,~,~,~,~,~,activeRows] = self.componentMatrix(component, []);
            PhiActive = Phi(activeRows,:);
            weightsActive = weights(activeRows);
            columnScale = max(abs(Phi),[],1);
            columnScale(columnScale == 0) = 1;
            Phi = Phi ./ columnScale;
            PhiActive = PhiActive ./ columnScale;
            nAvailable = min(options.nModes, size(Phi,2));
            iTailEnd = min(size(Phi,2), nAvailable + options.nTailCheck);
            nResolved = 0;
            for n = 1:nAvailable
                PhiN = PhiActive(:,1:n);
                gram = PhiN.' * (weightsActive .* PhiN);
                if cond(gram) > options.maxConditionNumber
                    break;
                end
                A = gram \ (weightsActive .* PhiN).';
                if n < iTailEnd
                    leakage = A * PhiActive(:,n+1:iTailEnd);
                    maxLeakage = max(vecnorm(leakage,2,1));
                else
                    maxLeakage = 0;
                end
                if maxLeakage > options.projectionTolerance
                    break;
                end
                nResolved = n;
            end
        end

        function rows = activeProjectionRows(~, Phi)
            if isempty(Phi)
                rows = zeros(0,1);
                return;
            end
            rowScale = max(abs(Phi),[],2);
            tolerance = sqrt(eps(max(1,max(rowScale))));
            rows = find(rowScale > tolerance);
            if isempty(rows)
                rows = (1:size(Phi,1)).';
            end
        end

        function limit = nonlinearAliasLimit(self, component, requested, policy)
            if string(policy) == "none"
                limit = requested;
                return;
            end
            component = string(component);
            if component == "F" || component == "both"
                limitF = floor(2*(self.nAvailableModes("F") - double(self.isGeostrophic()))/3);
                limitF = max(1, limitF);
            else
                limitF = requested;
            end
            if component == "G" || component == "both"
                limitG = max(1, floor(2*self.nAvailableModes("G")/3));
            else
                limitG = requested;
            end
            limit = min([requested limitF limitG]);
        end

        function reason = selectionReason(~, requested, projectionResolved, nonlinearLimit, retained)
            if retained == requested
                reason = "requestedModeCount";
            elseif retained == nonlinearLimit && nonlinearLimit <= projectionResolved
                reason = "nonlinearAliasingLimit";
            else
                reason = "projectionQualityLimit";
            end
        end

        function value = isGeostrophic(self)
            value = self.problemType == "geostrophic" || (isfinite(self.omega) && abs(self.omega) <= 10*eps(max(1,abs(self.f0))));
        end

        function value = isCanonicalComponent(self, component)
            component = string(component);
            if component == "F"
                value = self.forwardProjectionAvailableF;
            else
                value = self.forwardProjectionAvailableG;
            end
        end

        function status = defaultStatusForComponent(self, component)
            if self.isCanonicalComponent(component)
                status = "canonical";
            else
                status = "diagnosticOnly";
            end
        end

        function w = quadratureWeights(self)
            if isempty(self.z)
                w = zeros(0,1);
                return;
            end
            if self.isGeostrophic()
                FWithBarotropic = cat(2,ones(length(self.z),1),self.F);
            else
                FWithBarotropic = cat(2,ones(length(self.z),1),self.F);
            end
            b = zeros(size(FWithBarotropic,2),1);
            b(1) = self.D;
            if size(FWithBarotropic,1) == size(FWithBarotropic,2) && rank(FWithBarotropic.') == size(FWithBarotropic,2)
                w = (FWithBarotropic.') \ b;
            else
                w = pinv(FWithBarotropic.') * b;
            end
            if any(~isfinite(w)) || norm(w,1) == 0
                w = InternalModesBasis.trapezoidalWeights(self.z);
            end
        end

        function weights = gQuadratureWeights(self)
            w = self.quadratureWeights();
            if self.problemType == "igwWavenumber"
                weights = w .* (self.N2 - self.f0*self.f0) / self.g;
            else
                weights = w .* self.N2 / self.g;
            end
        end

        function [forward,status,conditionNumber,gramError] = forwardMatrix(~, Phi, weights, projectionWeights, projectionMethod, defaultStatus, maxConditionNumber, activeRows)
            projectionMethod = string(projectionMethod);
            status = string(defaultStatus);
            if nargin < 8 || isempty(activeRows)
                activeRows = (1:size(Phi,1)).';
            end
            PhiActive = Phi(activeRows,:);
            weightsActive = weights(activeRows);
            if isempty(PhiActive) || size(PhiActive,2) == 0
                forward = zeros(0,size(Phi,1));
                conditionNumber = 0;
                gramError = 0;
                return;
            end
            conditionNumber = cond(PhiActive);
            if projectionMethod == "auto"
                if size(PhiActive,1) == size(PhiActive,2) && conditionNumber <= maxConditionNumber
                    projectionMethod = "directInverse";
                else
                    projectionMethod = "weightedPseudoinverse";
                end
            end

            if projectionMethod == "directInverse"
                forwardActive = inv(PhiActive);
                status = "directInverse";
            elseif projectionMethod == "canonical"
                forwardActive = (PhiActive.' .* weightsActive.') ./ projectionWeights;
                status = "canonical";
            elseif projectionMethod == "weightedPseudoinverse"
                gram = PhiActive.' * (weightsActive .* PhiActive);
                forwardActive = gram \ (weightsActive .* PhiActive).';
                status = "weightedPseudoinverse";
                conditionNumber = cond(gram);
            else
                error("InternalModesBasis:InvalidProjectionMethod", ...
                    "projectionMethod must be auto, directInverse, weightedPseudoinverse, or canonical.");
            end
            forward = zeros(size(forwardActive,1),size(Phi,1));
            forward(:,activeRows) = forwardActive;
            gramError = norm(forward*Phi - eye(size(Phi,2)),'fro')/max(1,norm(eye(size(Phi,2)),'fro'));
        end

        function basis = interpolateToGrid(self, zTarget)
            FTarget = interp1(self.z,self.F,zTarget,"pchip","extrap");
            GTarget = interp1(self.z,self.G,zTarget,"pchip","extrap");
            N2Target = interp1(self.z,self.N2,zTarget,"pchip","extrap");
            basis = InternalModesBasis.fromSolvedModes(zTarget,FTarget,GTarget,self.h,N2=N2Target, ...
                f0=self.f0,g=self.g,kappa=self.kappa,omega=self.omega,problemType=self.problemType, ...
                sourceDescription=self.sourceDescription + " interpolated to fixed grid", ...
                componentRoleF=self.componentRoleF,componentRoleG=self.componentRoleG, ...
                forwardProjectionAvailableF=self.forwardProjectionAvailableF,forwardProjectionAvailableG=self.forwardProjectionAvailableG, ...
                orthogonalityWeightF=self.orthogonalityWeightF,orthogonalityWeightG=self.orthogonalityWeightG);
        end
    end

    methods (Static)
        function basis = fromSolverAtFrequency(solver, omega, options)
            % Solve modes at fixed frequency and return an annotated basis.
            %
            % For `omega=0`, the returned basis is marked geostrophic and
            % both F and G components are canonical vertical projection
            % bases.
            %
            % - Topic: Create vertical bases
            % - Declaration: basis = InternalModesBasis.fromSolverAtFrequency(solver,omega,options)
            % - Parameter solver: InternalModesBase solver instance
            % - Parameter omega: fixed frequency in radians per second
            % - Parameter options.nModes: number of mode columns retained
            % - Parameter options.useModeAdaptedGrid: true to request mode-adapted quadrature points when available
            % - Parameter options.nQuadraturePoints: number of quadrature points for mode-adapted solves
            % - Parameter options.g: gravitational acceleration used when solver state does not expose it publicly
            % - Returns basis: InternalModesBasis containing solved modes and component role
            arguments
                solver InternalModesBase
                omega (1,1) double
                options.nModes double = []
                options.useModeAdaptedGrid (1,1) logical = false
                options.nQuadraturePoints double = []
                options.g (1,1) double {mustBePositive} = 9.81
            end

            if options.useModeAdaptedGrid && isa(solver,"InternalModesSpectral")
                if isempty(options.nQuadraturePoints)
                    if isempty(options.nModes)
                        nQuadraturePoints = max(2, solver.nModes + 1);
                    else
                        nQuadraturePoints = options.nModes + 1;
                    end
                else
                    nQuadraturePoints = options.nQuadraturePoints;
                end
                [F,G,h,z] = solver.modesAtQuadraturePoints(nPoints=nQuadraturePoints,omega=omega);
            else
                [F,G,h] = solver.modesAtFrequency(omega);
                z = solver.z;
            end
            if ~isempty(options.nModes)
                nModes = min(options.nModes, size(F,2));
                F = F(:,1:nModes);
                G = G(:,1:nModes);
                h = h(1:nModes);
            end
            N2 = InternalModesBasis.sampleN2(solver,z);
            isGeostrophic = abs(omega) <= 10*eps(max(1,abs(solver.f0)));
            if isGeostrophic
                problemType = "geostrophic";
                componentRoleF = "eigenfunction";
                forwardF = true;
                weightF = "1";
                weightG = "N2/g";
            else
                problemType = "fixedFrequency";
                componentRoleF = "eigenfunction";
                forwardF = true;
                weightF = "1";
                weightG = "N2/g";
            end
            basis = InternalModesBasis.fromSolvedModes(z,F,G,h,N2=N2,f0=solver.f0,g=options.g,kappa=NaN,omega=omega, ...
                problemType=problemType,sourceDescription=class(solver),componentRoleF=componentRoleF, ...
                componentRoleG="eigenfunction",forwardProjectionAvailableF=forwardF,forwardProjectionAvailableG=true, ...
                orthogonalityWeightF=weightF,orthogonalityWeightG=weightG);
        end

        function basis = fromSolverAtWavenumber(solver, kappa, options)
            % Solve modes at fixed horizontal wavenumber.
            %
            % For nonzero $$\kappa$$, the returned basis marks G as the
            % canonical Sturm-Liouville projection component and F as
            % diagnostic/evaluation-only.
            %
            % - Topic: Create vertical bases
            % - Declaration: basis = InternalModesBasis.fromSolverAtWavenumber(solver,kappa,options)
            % - Parameter solver: InternalModesBase solver instance
            % - Parameter kappa: horizontal wavenumber in radians per meter
            % - Parameter options.nModes: number of mode columns retained
            % - Parameter options.useModeAdaptedGrid: true to request mode-adapted quadrature points when available
            % - Parameter options.nQuadraturePoints: number of quadrature points for mode-adapted solves
            % - Parameter options.g: gravitational acceleration used when solver state does not expose it publicly
            % - Returns basis: InternalModesBasis containing solved modes and component role
            arguments
                solver InternalModesBase
                kappa (1,1) double {mustBeNonnegative}
                options.nModes double = []
                options.useModeAdaptedGrid (1,1) logical = false
                options.nQuadraturePoints double = []
                options.g (1,1) double {mustBePositive} = 9.81
            end

            if options.useModeAdaptedGrid && isa(solver,"InternalModesSpectral")
                if isempty(options.nQuadraturePoints)
                    if isempty(options.nModes)
                        nQuadraturePoints = max(2, solver.nModes + 1);
                    else
                        nQuadraturePoints = options.nModes + 1;
                    end
                else
                    nQuadraturePoints = options.nQuadraturePoints;
                end
                [F,G,h,z] = solver.modesAtQuadraturePoints(nPoints=nQuadraturePoints,k=kappa);
                omega = sqrt(options.g*h*kappa*kappa + solver.f0*solver.f0);
            else
                [F,G,h,omega] = solver.modesAtWavenumber(kappa);
                z = solver.z;
            end
            if ~isempty(options.nModes)
                nModes = min(options.nModes, size(F,2));
                F = F(:,1:nModes);
                G = G(:,1:nModes);
                h = h(1:nModes);
                omega = omega(1:nModes);
            end
            N2 = InternalModesBasis.sampleN2(solver,z);
            isHydrostaticLimit = kappa == 0;
            if isHydrostaticLimit
                problemType = "geostrophic";
                componentRoleF = "eigenfunction";
                forwardF = true;
                weightF = "1";
                weightG = "N2/g";
            else
                problemType = "igwWavenumber";
                componentRoleF = "diagnostic";
                forwardF = false;
                weightF = "none";
                weightG = "(N2-f0^2)/g";
            end
            basis = InternalModesBasis.fromSolvedModes(z,F,G,h,N2=N2,f0=solver.f0,g=options.g,kappa=kappa,omega=mean(omega,"omitnan"), ...
                problemType=problemType,sourceDescription=class(solver),componentRoleF=componentRoleF, ...
                componentRoleG="eigenfunction",forwardProjectionAvailableF=forwardF,forwardProjectionAvailableG=true, ...
                orthogonalityWeightF=weightF,orthogonalityWeightG=weightG);
        end

        function basis = fromSolvedModes(z, F, G, h, options)
            % Create a basis from already-solved mode arrays.
            %
            % Use this factory when another workflow has already computed
            % sampled inverse modes and equivalent depths.
            %
            % - Topic: Create vertical bases
            % - Declaration: basis = InternalModesBasis.fromSolvedModes(z,F,G,h,options)
            % - Parameter z: depth grid
            % - Parameter F: sampled F inverse modes
            % - Parameter G: sampled G inverse modes
            % - Parameter h: equivalent-depth vector
            % - Parameter options.N2: buoyancy frequency squared sampled at `z`
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.kappa: horizontal wavenumber metadata
            % - Parameter options.omega: frequency metadata
            % - Parameter options.problemType: text label for the EVP
            % - Returns basis: InternalModesBasis containing the supplied modes
            arguments
                z (:,1) double
                F double
                G double
                h double
                options.N2 (:,1) double = zeros(size(z))
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
                options.kappa (1,1) double = NaN
                options.omega (1,1) double = NaN
                options.problemType {mustBeTextScalar} = "unknown"
                options.sourceDescription {mustBeTextScalar} = "manual"
                options.componentRoleF {mustBeTextScalar} = "diagnostic"
                options.componentRoleG {mustBeTextScalar} = "eigenfunction"
                options.forwardProjectionAvailableF (1,1) logical = false
                options.forwardProjectionAvailableG (1,1) logical = true
                options.orthogonalityWeightF {mustBeTextScalar} = "none"
                options.orthogonalityWeightG {mustBeTextScalar} = "unknown"
            end

            basis = InternalModesBasis(z=z,F=F,G=G,h=h,N2=options.N2,f0=options.f0,g=options.g,D=max(z)-min(z), ...
                kappa=options.kappa,omega=options.omega,problemType=options.problemType,sourceDescription=options.sourceDescription, ...
                componentRoleF=options.componentRoleF,componentRoleG=options.componentRoleG, ...
                forwardProjectionAvailableF=options.forwardProjectionAvailableF,forwardProjectionAvailableG=options.forwardProjectionAvailableG, ...
                orthogonalityWeightF=options.orthogonalityWeightF,orthogonalityWeightG=options.orthogonalityWeightG);
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CADimensionProperty('zIndex','','depth-grid row index');
            propertyAnnotations(end+1) = CADimensionProperty('modeF','','F-mode index');
            propertyAnnotations(end+1) = CADimensionProperty('modeG','','G-mode index');
            propertyAnnotations(end+1) = CANumericProperty('z',{'zIndex'},'m','depth grid');
            propertyAnnotations(end+1) = CANumericProperty('F',{'zIndex','modeF'},'','sampled F inverse modes');
            propertyAnnotations(end+1) = CANumericProperty('G',{'zIndex','modeG'},'','sampled G inverse modes');
            propertyAnnotations(end+1) = CANumericProperty('h',{'modeG'},'m','equivalent depths');
            propertyAnnotations(end+1) = CANumericProperty('N2',{'zIndex'},'s^-2','buoyancy frequency squared');
            propertyAnnotations(end+1) = CANumericProperty('f0',{},'s^-1','Coriolis parameter');
            propertyAnnotations(end+1) = CANumericProperty('g',{},'m s^-2','gravitational acceleration');
            propertyAnnotations(end+1) = CANumericProperty('D',{},'m','water-column depth');
            propertyAnnotations(end+1) = CANumericProperty('kappa',{},'rad m^-1','horizontal wavenumber');
            propertyAnnotations(end+1) = CANumericProperty('omega',{},'s^-1','frequency');
            propertyAnnotations(end+1) = CAPropertyAnnotation('problemType','vertical eigenvalue problem label');
            propertyAnnotations(end+1) = CAPropertyAnnotation('sourceDescription','source solver or factory description');
            propertyAnnotations(end+1) = CAPropertyAnnotation('componentRoleF','F component role');
            propertyAnnotations(end+1) = CAPropertyAnnotation('componentRoleG','G component role');
            propertyAnnotations(end+1) = CANumericProperty('forwardProjectionAvailableF',{},'','canonical F projection availability');
            propertyAnnotations(end+1) = CANumericProperty('forwardProjectionAvailableG',{},'','canonical G projection availability');
            propertyAnnotations(end+1) = CAPropertyAnnotation('orthogonalityWeightF','F orthogonality weight label');
            propertyAnnotations(end+1) = CAPropertyAnnotation('orthogonalityWeightG','G orthogonality weight label');
        end

        function names = classRequiredPropertyNames()
            names = {'z','F','G','h','N2','f0','g','D','kappa','omega','problemType','sourceDescription', ...
                'componentRoleF','componentRoleG','forwardProjectionAvailableF','forwardProjectionAvailableG', ...
                'orthogonalityWeightF','orthogonalityWeightG'};
        end
    end

    methods (Static, Access = private)
        function N2 = sampleN2(solver, z)
            if isprop(solver,"N2_function") && isa(solver.N2_function,"function_handle")
                N2 = solver.N2_function(z);
            elseif length(solver.z) == length(z) && max(abs(solver.z(:) - z(:))) < 1e-10
                N2 = solver.N2(:);
            else
                N2 = interp1(solver.z(:),solver.N2(:),z(:),"pchip","extrap");
            end
        end

        function w = trapezoidalWeights(z)
            z = z(:);
            w = zeros(size(z));
            if length(z) == 1
                return;
            end
            dz = diff(z);
            w(1) = dz(1)/2;
            w(end) = dz(end)/2;
            if length(z) > 2
                w(2:end-1) = (dz(1:end-1) + dz(2:end))/2;
            end
            w = abs(w);
        end
    end
end
