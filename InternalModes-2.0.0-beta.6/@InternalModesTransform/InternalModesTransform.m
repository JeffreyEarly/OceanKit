classdef InternalModesTransform < CAAnnotatedClass
    % Apply vertical modal transforms and compute vertical modal spectra.
    %
    % `InternalModesTransform` stores forward and inverse vertical
    % operators for one `InternalModesBasis`. It works only on vertical
    % scalar fields or modal coefficients; full wave-vortex state
    % projectors remain outside this package.
    %
    % Public matrices use the full transform grid `z`. Some components are
    % solved on component-specific active rows before being expanded back to
    % that grid. For example, rigid-lid G transforms omit the zero Dirichlet
    % endpoint rows during the solve, then store zero endpoint columns in the
    % forward matrix and zero endpoint rows in the inverse matrix.
    %
    % For modal coefficients $$p_j$$ and $$q_j$$, the vertical
    % cross-spectrum is
    %
    % $$
    % S[p,q]_j = m\,s_j\,\Re\{p_j q_j^*\},
    % $$
    %
    % where $$s_j$$ is the component spectral weight and $$m$$ is an
    % optional horizontal multiplicity supplied by the caller.
    %
    % ```matlab
    % coefficients = transform.projectVertical(profile,component="G");
    % spectrum = transform.spectrum(coefficients,component="G");
    % ```
    %
    % - Topic: Inspect transform properties
    % - Topic: Apply vertical transforms
    % - Topic: Analyze vertical spectra
    % - Topic: Persist vertical transforms
    % - Topic: Developer topics
    % - Declaration: classdef InternalModesTransform < CAAnnotatedClass

    properties (SetAccess = private)
        % Depth grid for vertical fields acted on by this transform.
        %
        % - Topic: Inspect transform properties
        z

        % Row coordinate used by annotated NetCDF persistence.
        %
        % - Topic: Developer topics
        % - Developer: true
        zIndex

        % F-mode coordinate used by annotated NetCDF persistence.
        %
        % - Topic: Developer topics
        % - Developer: true
        modeF

        % G-mode coordinate used by annotated NetCDF persistence.
        %
        % - Topic: Developer topics
        % - Developer: true
        modeG

        % Buoyancy frequency squared sampled at `z`.
        %
        % - Topic: Inspect transform properties
        N2

        % Coriolis parameter used by the source mode solve.
        %
        % - Topic: Inspect transform properties
        f0

        % Gravitational acceleration used by the source mode solve.
        %
        % - Topic: Inspect transform properties
        g

        % Water-column depth used by the barotropic F spectral weight.
        %
        % - Topic: Inspect transform properties
        D

        % Horizontal wavenumber metadata.
        %
        % - Topic: Inspect transform properties
        kappa

        % Frequency metadata.
        %
        % - Topic: Inspect transform properties
        omega

        % Forward F projection matrix from samples to modal coefficients.
        %
        % The rows are retained F modes and columns are vertical grid
        % samples.
        %
        % - Topic: Inspect transform properties
        forwardF

        % Inverse F reconstruction matrix from coefficients to samples.
        %
        % - Topic: Inspect transform properties
        inverseF

        % Forward G projection matrix from samples to modal coefficients.
        %
        % For rigid-lid G modes, zero endpoint columns may be present so the
        % matrix remains compatible with fields sampled on the full `z`
        % grid while the projection itself uses the active interior rows.
        %
        % - Topic: Inspect transform properties
        forwardG

        % Inverse G reconstruction matrix from coefficients to samples.
        %
        % For rigid-lid G modes, zero endpoint rows may be present so the
        % matrix reconstructs onto the full `z` grid.
        %
        % - Topic: Inspect transform properties
        inverseG

        % Quadrature weights used to build F numerical projections.
        %
        % - Topic: Inspect transform properties
        weightsF

        % Quadrature weights used to build G numerical projections.
        %
        % - Topic: Inspect transform properties
        weightsG

        % Parseval weights used for F modal spectra.
        %
        % For geostrophic F modes this stores $$\gamma_j$$.
        %
        % - Topic: Analyze vertical spectra
        spectralWeightsF

        % Parseval weights used for G modal spectra.
        %
        % For canonical G modes this stores $$g$$ for each retained mode.
        %
        % - Topic: Analyze vertical spectra
        spectralWeightsG

        % Retained F mode numbers.
        %
        % - Topic: Inspect transform properties
        retainedModesF

        % Rejected F mode numbers.
        %
        % - Topic: Inspect transform properties
        rejectedModesF

        % Retained G mode numbers.
        %
        % - Topic: Inspect transform properties
        retainedModesG

        % Rejected G mode numbers.
        %
        % - Topic: Inspect transform properties
        rejectedModesG

        % Equivalent depths for retained F modes.
        %
        % - Topic: Inspect transform properties
        hF

        % Equivalent depths for retained G modes.
        %
        % - Topic: Inspect transform properties
        hG

        % F component role label.
        %
        % - Topic: Inspect transform properties
        componentRoleF

        % G component role label.
        %
        % - Topic: Inspect transform properties
        componentRoleG

        % F transform construction status.
        %
        % Values include `"canonical"`, `"directInverse"`,
        % `"weightedPseudoinverse"`, and `"diagnosticOnly"`.
        %
        % - Topic: Inspect transform properties
        transformStatusF

        % G transform construction status.
        %
        % - Topic: Inspect transform properties
        transformStatusG

        % True when a canonical F projection exists for this basis.
        %
        % - Topic: Inspect transform properties
        forwardProjectionAvailableF

        % True when a canonical G projection exists for this basis.
        %
        % - Topic: Inspect transform properties
        forwardProjectionAvailableG

        % Condition-number diagnostic for the F forward operator.
        %
        % - Topic: Inspect transform properties
        conditionNumberF

        % Condition-number diagnostic for the G forward operator.
        %
        % - Topic: Inspect transform properties
        conditionNumberG

        % Relative round-trip error for the retained F modes.
        %
        % - Topic: Inspect transform properties
        gramErrorF

        % Relative round-trip error for the retained G modes.
        %
        % - Topic: Inspect transform properties
        gramErrorG

        % Prefix-selection reason for model-style transforms.
        %
        % - Topic: Inspect transform properties
        selectionReason

        % Number of modes accepted by the projection-quality diagnostic.
        %
        % - Topic: Inspect transform properties
        projectionResolvedModes

        % Number of modes allowed by the nonlinear aliasing policy.
        %
        % - Topic: Inspect transform properties
        nonlinearAliasLimit

        % Text label identifying the source vertical eigenvalue problem.
        %
        % - Topic: Inspect transform properties
        problemType

        % Text label identifying the source solver or factory path.
        %
        % - Topic: Inspect transform properties
        sourceDescription
    end

    methods
        function self = InternalModesTransform(options)
            % Create a vertical transform from canonical persisted state.
            %
            % This constructor stores already-built operators. Use
            % `InternalModesBasis` factories to solve modes and build new
            % transforms from stratification.
            %
            % - Topic: Persist vertical transforms
            % - Declaration: transform = InternalModesTransform(options)
            % - Parameter options.z: depth grid
            % - Parameter options.forwardF: F forward matrix
            % - Parameter options.inverseF: F inverse matrix
            % - Parameter options.forwardG: G forward matrix
            % - Parameter options.inverseG: G inverse matrix
            % - Parameter options.spectralWeightsF: F spectrum weights
            % - Parameter options.spectralWeightsG: G spectrum weights
            % - Returns transform: initialized InternalModesTransform instance
            arguments
                options.z (:,1) double = zeros(0,1)
                options.N2 (:,1) double = zeros(0,1)
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
                options.D (1,1) double {mustBeNonnegative} = 0
                options.kappa (1,1) double = NaN
                options.omega (1,1) double = NaN
                options.forwardF double = zeros(0,0)
                options.inverseF double = zeros(0,0)
                options.forwardG double = zeros(0,0)
                options.inverseG double = zeros(0,0)
                options.weightsF (:,1) double = zeros(0,1)
                options.weightsG (:,1) double = zeros(0,1)
                options.spectralWeightsF (:,1) double = zeros(0,1)
                options.spectralWeightsG (:,1) double = zeros(0,1)
                options.retainedModesF (:,1) double = zeros(0,1)
                options.rejectedModesF (:,1) double = zeros(0,1)
                options.retainedModesG (:,1) double = zeros(0,1)
                options.rejectedModesG (:,1) double = zeros(0,1)
                options.hF (:,1) double = zeros(0,1)
                options.hG (:,1) double = zeros(0,1)
                options.componentRoleF {mustBeTextScalar} = "diagnostic"
                options.componentRoleG {mustBeTextScalar} = "eigenfunction"
                options.transformStatusF {mustBeTextScalar} = "diagnosticOnly"
                options.transformStatusG {mustBeTextScalar} = "diagnosticOnly"
                options.forwardProjectionAvailableF (1,1) logical = false
                options.forwardProjectionAvailableG (1,1) logical = false
                options.conditionNumberF (1,1) double = NaN
                options.conditionNumberG (1,1) double = NaN
                options.gramErrorF (1,1) double = NaN
                options.gramErrorG (1,1) double = NaN
                options.selectionReason {mustBeTextScalar} = "notSelected"
                options.projectionResolvedModes (1,1) double = NaN
                options.nonlinearAliasLimit (1,1) double = NaN
                options.problemType {mustBeTextScalar} = "unknown"
                options.sourceDescription {mustBeTextScalar} = "manual"
            end

            self@CAAnnotatedClass();
            self.z = options.z(:);
            self.zIndex = (1:length(self.z)).';
            self.N2 = options.N2(:);
            self.f0 = options.f0;
            self.g = options.g;
            self.D = options.D;
            self.kappa = options.kappa;
            self.omega = options.omega;
            self.forwardF = options.forwardF;
            self.inverseF = options.inverseF;
            self.forwardG = options.forwardG;
            self.inverseG = options.inverseG;
            self.weightsF = options.weightsF(:);
            self.weightsG = options.weightsG(:);
            self.spectralWeightsF = options.spectralWeightsF(:);
            self.spectralWeightsG = options.spectralWeightsG(:);
            self.retainedModesF = options.retainedModesF(:);
            self.rejectedModesF = options.rejectedModesF(:);
            self.retainedModesG = options.retainedModesG(:);
            self.rejectedModesG = options.rejectedModesG(:);
            self.hF = options.hF(:);
            self.hG = options.hG(:);
            self.componentRoleF = string(options.componentRoleF);
            self.componentRoleG = string(options.componentRoleG);
            self.transformStatusF = string(options.transformStatusF);
            self.transformStatusG = string(options.transformStatusG);
            self.forwardProjectionAvailableF = options.forwardProjectionAvailableF;
            self.forwardProjectionAvailableG = options.forwardProjectionAvailableG;
            self.conditionNumberF = options.conditionNumberF;
            self.conditionNumberG = options.conditionNumberG;
            self.gramErrorF = options.gramErrorF;
            self.gramErrorG = options.gramErrorG;
            self.selectionReason = string(options.selectionReason);
            self.projectionResolvedModes = options.projectionResolvedModes;
            self.nonlinearAliasLimit = options.nonlinearAliasLimit;
            self.problemType = string(options.problemType);
            self.sourceDescription = string(options.sourceDescription);
            self.modeF = (1:size(self.inverseF,2)).';
            self.modeG = (1:size(self.inverseG,2)).';
        end

        function matrix = forward(self, options)
            % Return a forward vertical projection matrix.
            %
            % The returned matrix maps samples on `z` to vertical modal
            % coefficients. Noncanonical wave-F projections are returned
            % only when `allowNoncanonical=true`. Components with inactive
            % rows keep full-grid matrix shapes; those inactive columns are
            % zero in the returned forward matrix.
            %
            % - Topic: Apply vertical transforms
            % - Declaration: matrix = forward(self,options)
            % - Parameter self: InternalModesTransform instance
            % - Parameter options.component: `"F"` or `"G"`
            % - Parameter options.allowNoncanonical: true to allow numerical wave-F forward matrices
            % - Returns matrix: forward projection matrix
            arguments
                self InternalModesTransform
                options.component {mustBeMember(options.component,["F","G"])} = "G"
                options.allowNoncanonical (1,1) logical = false
            end

            component = string(options.component);
            if component == "F"
                if ~self.forwardProjectionAvailableF && ~options.allowNoncanonical
                    error("InternalModesTransform:NoncanonicalFProjection", ...
                        "This transform has no canonical F projection. Pass allowNoncanonical=true to retrieve the numerical F matrix.");
                end
                matrix = self.forwardF;
            else
                matrix = self.forwardG;
            end
        end

        function matrix = inverse(self, options)
            % Return an inverse vertical reconstruction matrix.
            %
            % The returned matrix maps vertical modal coefficients to field
            % samples on `z`.
            %
            % - Topic: Apply vertical transforms
            % - Declaration: matrix = inverse(self,options)
            % - Parameter self: InternalModesTransform instance
            % - Parameter options.component: `"F"` or `"G"`
            % - Returns matrix: inverse reconstruction matrix
            arguments
                self InternalModesTransform
                options.component {mustBeMember(options.component,["F","G"])} = "G"
            end

            if string(options.component) == "F"
                matrix = self.inverseF;
            else
                matrix = self.inverseG;
            end
        end

        function coefficients = projectVertical(self, field, options)
            % Project vertical samples onto modal coefficients.
            %
            % If `field` has multiple columns, each column is projected
            % independently with the same vertical operator.
            %
            % - Topic: Apply vertical transforms
            % - Declaration: coefficients = projectVertical(self,field,options)
            % - Parameter self: InternalModesTransform instance
            % - Parameter field: vertical samples with rows matching `z`
            % - Parameter options.component: `"F"` or `"G"`
            % - Parameter options.allowNoncanonical: true to allow numerical wave-F projection
            % - Returns coefficients: modal coefficients
            arguments
                self InternalModesTransform
                field double
                options.component {mustBeMember(options.component,["F","G"])} = "G"
                options.allowNoncanonical (1,1) logical = false
            end

            coefficients = self.forward(component=options.component,allowNoncanonical=options.allowNoncanonical) * field;
        end

        function field = reconstructVertical(self, coefficients, options)
            % Reconstruct vertical samples from modal coefficients.
            %
            % - Topic: Apply vertical transforms
            % - Declaration: field = reconstructVertical(self,coefficients,options)
            % - Parameter self: InternalModesTransform instance
            % - Parameter coefficients: modal coefficients
            % - Parameter options.component: `"F"` or `"G"`
            % - Returns field: reconstructed vertical samples on `z`
            arguments
                self InternalModesTransform
                coefficients double
                options.component {mustBeMember(options.component,["F","G"])} = "G"
            end

            field = self.inverse(component=options.component) * coefficients;
        end

        function matrix = crossTransformTo(self, targetTransform, options)
            % Map coefficients from a target transform into this transform.
            %
            % The matrix is
            %
            % $$
            % C = P_{\mathrm{self}}\Phi_{\mathrm{target}},
            % $$
            %
            % where `P_self` is this transform's forward matrix and
            % `Phi_target` is the target transform's inverse matrix for the
            % selected component.
            %
            % - Topic: Apply vertical transforms
            % - Declaration: matrix = crossTransformTo(self,targetTransform,options)
            % - Parameter self: InternalModesTransform instance
            % - Parameter targetTransform: transform whose coefficients should be evaluated in this transform
            % - Parameter options.component: `"F"` or `"G"`
            % - Parameter options.allowNoncanonical: true to allow numerical F matrices
            % - Returns matrix: coefficient-to-coefficient transform matrix
            arguments
                self InternalModesTransform
                targetTransform InternalModesTransform
                options.component {mustBeMember(options.component,["F","G"])} = "G"
                options.allowNoncanonical (1,1) logical = false
            end

            matrix = self.forward(component=options.component,allowNoncanonical=options.allowNoncanonical) * targetTransform.inverse(component=options.component);
        end

        function S = spectrum(self, coefficients, options)
            % Compute a vertical modal auto-spectrum.
            %
            % This is equivalent to `crossSpectrum(coefficients,coefficients)`.
            %
            % - Topic: Analyze vertical spectra
            % - Declaration: S = spectrum(self,coefficients,options)
            % - Parameter self: InternalModesTransform instance
            % - Parameter coefficients: vertical modal coefficients
            % - Parameter options.component: `"F"` or `"G"`
            % - Parameter options.horizontalMultiplicity: multiplicity factor, often 1 or 2
            % - Parameter options.requireCanonical: true to reject noncanonical component spectra
            % - Returns S: modal auto-spectrum
            arguments
                self InternalModesTransform
                coefficients double
                options.component {mustBeMember(options.component,["F","G"])} = "G"
                options.horizontalMultiplicity double = 1
                options.requireCanonical (1,1) logical = true
            end

            S = self.crossSpectrum(coefficients,coefficients,component=options.component, ...
                horizontalMultiplicity=options.horizontalMultiplicity,requireCanonical=options.requireCanonical);
        end

        function S = crossSpectrum(self, coefficientsA, coefficientsB, options)
            % Compute a vertical modal cross-spectrum.
            %
            % For component spectral weights $$s_j$$, coefficients $$a_j$$
            % and $$b_j$$, and multiplicity $$m$$, this method returns
            %
            % $$
            % S_j = m\,s_j\,\Re\{a_j b_j^*\}.
            % $$
            %
            % - Topic: Analyze vertical spectra
            % - Declaration: S = crossSpectrum(self,coefficientsA,coefficientsB,options)
            % - Parameter self: InternalModesTransform instance
            % - Parameter coefficientsA: first modal coefficient array
            % - Parameter coefficientsB: second modal coefficient array
            % - Parameter options.component: `"F"` or `"G"`
            % - Parameter options.horizontalMultiplicity: multiplicity factor, often 1 or 2
            % - Parameter options.requireCanonical: true to reject noncanonical component spectra
            % - Returns S: modal cross-spectrum
            arguments
                self InternalModesTransform
                coefficientsA double
                coefficientsB double
                options.component {mustBeMember(options.component,["F","G"])} = "G"
                options.horizontalMultiplicity double = 1
                options.requireCanonical (1,1) logical = true
            end

            component = string(options.component);
            if options.requireCanonical && ~self.isCanonicalComponent(component)
                error("InternalModesTransform:NoncanonicalSpectrum", ...
                    "A canonical %s spectrum is not available for this transform.", component);
            end
            weights = self.componentSpectralWeights(component);
            multiplicity = options.horizontalMultiplicity;
            if isscalar(multiplicity)
                multiplicity = multiplicity * ones(size(weights));
            else
                multiplicity = multiplicity(:);
            end
            S = (multiplicity .* weights) .* real(coefficientsA .* conj(coefficientsB));
        end

        function transform = withDiagnostics(self, options)
            % Return a copy with updated selection diagnostics.
            %
            % - Topic: Developer topics
            % - Declaration: transform = withDiagnostics(self,options)
            % - Parameter self: InternalModesTransform instance
            % - Parameter options.selectionReason: reason retained modes were selected
            % - Parameter options.projectionResolvedModes: projection-quality mode count
            % - Parameter options.nonlinearAliasLimit: nonlinear aliasing mode limit
            % - Returns transform: copied transform with diagnostics updated
            % - Developer: true
            arguments
                self InternalModesTransform
                options.selectionReason {mustBeTextScalar} = self.selectionReason
                options.projectionResolvedModes (1,1) double = self.projectionResolvedModes
                options.nonlinearAliasLimit (1,1) double = self.nonlinearAliasLimit
            end

            transform = InternalModesTransform(z=self.z,N2=self.N2,f0=self.f0,g=self.g,D=self.D,kappa=self.kappa,omega=self.omega, ...
                forwardF=self.forwardF,inverseF=self.inverseF,forwardG=self.forwardG,inverseG=self.inverseG, ...
                weightsF=self.weightsF,weightsG=self.weightsG,spectralWeightsF=self.spectralWeightsF,spectralWeightsG=self.spectralWeightsG, ...
                retainedModesF=self.retainedModesF,rejectedModesF=self.rejectedModesF,retainedModesG=self.retainedModesG,rejectedModesG=self.rejectedModesG, ...
                hF=self.hF,hG=self.hG,componentRoleF=self.componentRoleF,componentRoleG=self.componentRoleG, ...
                transformStatusF=self.transformStatusF,transformStatusG=self.transformStatusG, ...
                forwardProjectionAvailableF=self.forwardProjectionAvailableF,forwardProjectionAvailableG=self.forwardProjectionAvailableG, ...
                conditionNumberF=self.conditionNumberF,conditionNumberG=self.conditionNumberG,gramErrorF=self.gramErrorF,gramErrorG=self.gramErrorG, ...
                selectionReason=options.selectionReason,projectionResolvedModes=options.projectionResolvedModes, ...
                nonlinearAliasLimit=options.nonlinearAliasLimit,problemType=self.problemType,sourceDescription=self.sourceDescription);
        end

        function ncfile = writeToFile(self, path, properties, options)
            % Write this transform to an annotated NetCDF file.
            %
            % Empty component blocks are omitted so a G-only or F-only
            % transform can still be persisted without creating zero-length
            % NetCDF dimensions.
            %
            % - Topic: Persist vertical transforms
            % - Declaration: ncfile = writeToFile(self,path,properties,options)
            % - Parameter self: InternalModesTransform instance
            % - Parameter path: destination NetCDF path
            % - Parameter properties: optional additional annotated properties
            % - Parameter options.shouldOverwriteExisting: true to replace an existing file
            % - Parameter options.shouldAddRequiredProperties: true to include the persisted transform state
            % - Parameter options.attributes: additional NetCDF attributes
            % - Returns ncfile: NetCDFFile written by CAAnnotatedClass
            arguments (Input)
                self InternalModesTransform
                path char {mustBeNonempty}
            end
            arguments (Input,Repeating)
                properties char
            end
            arguments (Input)
                options.shouldOverwriteExisting logical = false
                options.shouldAddRequiredProperties logical = true
                options.attributes = configureDictionary("string","string")
            end
            arguments (Output)
                ncfile
            end

            if options.shouldAddRequiredProperties
                properties = union(properties, self.nonemptyRequiredPropertyNames());
            end
            ncfile = writeToFile@CAAnnotatedClass(self,path,properties{:}, ...
                shouldOverwriteExisting=options.shouldOverwriteExisting, ...
                shouldAddRequiredProperties=false,attributes=options.attributes);
        end
    end

    methods (Access = private)
        function names = nonemptyRequiredPropertyNames(self)
            names = InternalModesTransform.classRequiredPropertyNames();
            if isempty(self.modeF)
                names = setdiff(names, {'forwardF','inverseF','weightsF','spectralWeightsF','retainedModesF','hF', ...
                    'componentRoleF','transformStatusF','forwardProjectionAvailableF','conditionNumberF','gramErrorF'});
            end
            if isempty(self.modeG)
                names = setdiff(names, {'forwardG','inverseG','weightsG','spectralWeightsG','retainedModesG','hG', ...
                    'componentRoleG','transformStatusG','forwardProjectionAvailableG','conditionNumberG','gramErrorG'});
            end
        end

        function weights = componentSpectralWeights(self, component)
            if string(component) == "F"
                weights = self.spectralWeightsF(:);
            else
                weights = self.spectralWeightsG(:);
            end
        end

        function value = isCanonicalComponent(self, component)
            if string(component) == "F"
                value = self.forwardProjectionAvailableF;
            else
                value = self.forwardProjectionAvailableG;
            end
        end
    end

    methods (Static)
        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CADimensionProperty('zIndex','','depth-grid row index');
            propertyAnnotations(end+1) = CADimensionProperty('modeF','','F-mode index');
            propertyAnnotations(end+1) = CADimensionProperty('modeG','','G-mode index');
            propertyAnnotations(end+1) = CANumericProperty('z',{'zIndex'},'m','depth grid');
            propertyAnnotations(end+1) = CANumericProperty('N2',{'zIndex'},'s^-2','buoyancy frequency squared');
            propertyAnnotations(end+1) = CANumericProperty('f0',{},'s^-1','Coriolis parameter');
            propertyAnnotations(end+1) = CANumericProperty('g',{},'m s^-2','gravitational acceleration');
            propertyAnnotations(end+1) = CANumericProperty('D',{},'m','water-column depth');
            propertyAnnotations(end+1) = CANumericProperty('kappa',{},'rad m^-1','horizontal wavenumber');
            propertyAnnotations(end+1) = CANumericProperty('omega',{},'s^-1','frequency');
            propertyAnnotations(end+1) = CANumericProperty('forwardF',{'modeF','zIndex'},'','F forward projection matrix');
            propertyAnnotations(end+1) = CANumericProperty('inverseF',{'zIndex','modeF'},'','F inverse reconstruction matrix');
            propertyAnnotations(end+1) = CANumericProperty('forwardG',{'modeG','zIndex'},'','G forward projection matrix');
            propertyAnnotations(end+1) = CANumericProperty('inverseG',{'zIndex','modeG'},'','G inverse reconstruction matrix');
            propertyAnnotations(end+1) = CANumericProperty('weightsF',{'zIndex'},'','F quadrature weights');
            propertyAnnotations(end+1) = CANumericProperty('weightsG',{'zIndex'},'','G quadrature weights');
            propertyAnnotations(end+1) = CANumericProperty('spectralWeightsF',{'modeF'},'','F spectral weights');
            propertyAnnotations(end+1) = CANumericProperty('spectralWeightsG',{'modeG'},'','G spectral weights');
            propertyAnnotations(end+1) = CANumericProperty('retainedModesF',{'modeF'},'','retained F mode numbers');
            propertyAnnotations(end+1) = CANumericProperty('retainedModesG',{'modeG'},'','retained G mode numbers');
            propertyAnnotations(end+1) = CANumericProperty('hF',{'modeF'},'m','F equivalent depths');
            propertyAnnotations(end+1) = CANumericProperty('hG',{'modeG'},'m','G equivalent depths');
            propertyAnnotations(end+1) = CANumericProperty('conditionNumberF',{},'','F condition number');
            propertyAnnotations(end+1) = CANumericProperty('conditionNumberG',{},'','G condition number');
            propertyAnnotations(end+1) = CANumericProperty('gramErrorF',{},'','F round-trip error');
            propertyAnnotations(end+1) = CANumericProperty('gramErrorG',{},'','G round-trip error');
            propertyAnnotations(end+1) = CANumericProperty('projectionResolvedModes',{},'','projection-resolved mode count');
            propertyAnnotations(end+1) = CANumericProperty('nonlinearAliasLimit',{},'','nonlinear alias mode limit');
            propertyAnnotations(end+1) = CAPropertyAnnotation('componentRoleF','F component role');
            propertyAnnotations(end+1) = CAPropertyAnnotation('componentRoleG','G component role');
            propertyAnnotations(end+1) = CAPropertyAnnotation('transformStatusF','F transform construction status');
            propertyAnnotations(end+1) = CAPropertyAnnotation('transformStatusG','G transform construction status');
            propertyAnnotations(end+1) = CANumericProperty('forwardProjectionAvailableF',{},'','canonical F projection availability');
            propertyAnnotations(end+1) = CANumericProperty('forwardProjectionAvailableG',{},'','canonical G projection availability');
            propertyAnnotations(end+1) = CAPropertyAnnotation('selectionReason','mode-selection reason');
            propertyAnnotations(end+1) = CAPropertyAnnotation('problemType','vertical eigenvalue problem label');
            propertyAnnotations(end+1) = CAPropertyAnnotation('sourceDescription','source solver or factory description');
        end

        function names = classRequiredPropertyNames()
            names = {'z','N2','f0','g','D','kappa','omega','forwardF','inverseF','forwardG','inverseG', ...
                'weightsF','weightsG','spectralWeightsF','spectralWeightsG','retainedModesF','retainedModesG', ...
                'hF','hG','componentRoleF','componentRoleG','transformStatusF','transformStatusG', ...
                'forwardProjectionAvailableF','forwardProjectionAvailableG','conditionNumberF','conditionNumberG', ...
                'gramErrorF','gramErrorG','projectionResolvedModes','nonlinearAliasLimit','selectionReason','problemType','sourceDescription'};
        end
    end
end
