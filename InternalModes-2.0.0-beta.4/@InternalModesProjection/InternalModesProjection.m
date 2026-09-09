classdef InternalModesProjection < CAAnnotatedClass
    % Project arbitrary vertical observations onto resolvable modes.
    %
    % `InternalModesProjection` stores the weighted least-squares operator
    % for an observation matrix $$B=H\Phi$$. Unlike model transforms, an
    % observation grid may identify a non-contiguous subset of modes. The
    % resolution matrix
    %
    % $$
    % R = A B
    % $$
    %
    % describes how true candidate-mode coefficients appear in the
    % recovered retained coefficients.
    %
    % ```matlab
    % projection = basis.observationProjection(zObs,component="G");
    % etaHat = projection.project(etaObs);
    % ```
    %
    % - Topic: Inspect projection properties
    % - Topic: Apply observation projections
    % - Topic: Analyze observation spectra
    % - Topic: Persist observation projections
    % - Topic: Developer topics
    % - Declaration: classdef InternalModesProjection < CAAnnotatedClass

    properties (SetAccess = private)
        % Component label, either `"F"` or `"G"`.
        %
        % - Topic: Inspect projection properties
        component

        % Observation depths for point-sampled projections.
        %
        % This is empty when the projection was built from an explicit
        % observation matrix.
        %
        % - Topic: Inspect projection properties
        observationZ

        % Observation index coordinate used by annotated persistence.
        %
        % - Topic: Developer topics
        % - Developer: true
        observationIndex

        % Candidate-mode coordinate used by annotated persistence.
        %
        % - Topic: Developer topics
        % - Developer: true
        candidateModeIndex

        % Retained-mode coordinate used by annotated persistence.
        %
        % - Topic: Developer topics
        % - Developer: true
        retainedModeIndex

        % Explicit observation matrix, when supplied by the caller.
        %
        % - Topic: Inspect projection properties
        observationMatrix

        % Candidate mode numbers considered by the projection.
        %
        % - Topic: Inspect projection properties
        candidateModes

        % Retained mode numbers selected as resolvable.
        %
        % - Topic: Inspect projection properties
        retainedModes

        % Rejected candidate mode numbers.
        %
        % - Topic: Inspect projection properties
        rejectedModes

        % Forward projection matrix from observations to retained coefficients.
        %
        % - Topic: Inspect projection properties
        projectionMatrix

        % Reconstruction matrix from retained coefficients to native samples.
        %
        % - Topic: Inspect projection properties
        reconstructionMatrix

        % Resolution matrix from true candidate coefficients to recovered rows.
        %
        % - Topic: Inspect projection properties
        resolutionMatrix

        % Alias matrix formed from rejected-mode columns of the resolution matrix.
        %
        % - Topic: Inspect projection properties
        aliasingMatrix

        % Squared resolution matrix used as an expected spectral window.
        %
        % For uncorrelated modal coefficients, this matrix maps a true
        % modal spectrum to the expected recovered retained-mode spectrum.
        %
        % - Topic: Analyze observation spectra
        spectralWindow

        % Retained Gram matrix used in the weighted least-squares solve.
        %
        % - Topic: Inspect projection properties
        gramMatrix

        % Condition number of `gramMatrix`.
        %
        % - Topic: Inspect projection properties
        conditionNumber

        % Candidate-mode spectral weights.
        %
        % - Topic: Analyze observation spectra
        spectralWeights

        % Transform construction status.
        %
        % - Topic: Inspect projection properties
        transformStatus

        % Component role label.
        %
        % - Topic: Inspect projection properties
        componentRole

        % True when this projection has a canonical modal spectrum.
        %
        % - Topic: Inspect projection properties
        isCanonicalComponent

        % QR rank tolerance used for mode selection.
        %
        % - Topic: Inspect projection properties
        diagnosticRankTolerance

        % Final retained relative QR pivot.
        %
        % - Topic: Inspect projection properties
        diagnosticFinalRelativePivot
    end

    methods
        function self = InternalModesProjection(options)
            % Create an observation projection from canonical persisted state.
            %
            % - Topic: Persist observation projections
            % - Declaration: projection = InternalModesProjection(options)
            % - Parameter options.component: `"F"` or `"G"`
            % - Parameter options.projectionMatrix: matrix from observations to retained coefficients
            % - Parameter options.resolutionMatrix: matrix from candidate coefficients to recovered coefficients
            % - Parameter options.spectralWindow: spectrum transfer matrix
            % - Returns projection: initialized InternalModesProjection instance
            arguments
                options.component {mustBeTextScalar} = "G"
                options.observationZ (:,1) double = zeros(0,1)
                options.observationMatrix double = zeros(0,0)
                options.candidateModes (:,1) double = zeros(0,1)
                options.retainedModes (:,1) double = zeros(0,1)
                options.rejectedModes (:,1) double = zeros(0,1)
                options.projectionMatrix double = zeros(0,0)
                options.reconstructionMatrix double = zeros(0,0)
                options.resolutionMatrix double = zeros(0,0)
                options.aliasingMatrix double = zeros(0,0)
                options.spectralWindow double = zeros(0,0)
                options.gramMatrix double = zeros(0,0)
                options.conditionNumber (1,1) double = NaN
                options.spectralWeights (:,1) double = zeros(0,1)
                options.transformStatus {mustBeTextScalar} = "weightedPseudoinverse"
                options.componentRole {mustBeTextScalar} = "eigenfunction"
                options.isCanonicalComponent (1,1) logical = true
                options.diagnosticRankTolerance (1,1) double = NaN
                options.diagnosticFinalRelativePivot (1,1) double = NaN
            end

            self@CAAnnotatedClass();
            self.component = string(options.component);
            self.observationZ = options.observationZ(:);
            self.observationMatrix = options.observationMatrix;
            self.candidateModes = options.candidateModes(:);
            self.retainedModes = options.retainedModes(:);
            self.rejectedModes = options.rejectedModes(:);
            self.projectionMatrix = options.projectionMatrix;
            self.reconstructionMatrix = options.reconstructionMatrix;
            self.resolutionMatrix = options.resolutionMatrix;
            self.aliasingMatrix = options.aliasingMatrix;
            self.spectralWindow = options.spectralWindow;
            self.gramMatrix = options.gramMatrix;
            self.conditionNumber = options.conditionNumber;
            self.spectralWeights = options.spectralWeights(:);
            self.transformStatus = string(options.transformStatus);
            self.componentRole = string(options.componentRole);
            self.isCanonicalComponent = options.isCanonicalComponent;
            self.diagnosticRankTolerance = options.diagnosticRankTolerance;
            self.diagnosticFinalRelativePivot = options.diagnosticFinalRelativePivot;
            self.observationIndex = (1:size(self.projectionMatrix,2)).';
            self.candidateModeIndex = (1:length(self.candidateModes)).';
            self.retainedModeIndex = (1:length(self.retainedModes)).';
        end

        function coefficients = project(self, observations)
            % Project observations onto retained modal coefficients.
            %
            % - Topic: Apply observation projections
            % - Declaration: coefficients = project(self,observations)
            % - Parameter self: InternalModesProjection instance
            % - Parameter observations: observation values with rows matching the observation operator
            % - Returns coefficients: retained modal coefficients
            arguments
                self InternalModesProjection
                observations double
            end

            coefficients = self.projectionMatrix * observations;
        end

        function values = reconstruct(self, coefficients)
            % Reconstruct sampled values from retained coefficients.
            %
            % - Topic: Apply observation projections
            % - Declaration: values = reconstruct(self,coefficients)
            % - Parameter self: InternalModesProjection instance
            % - Parameter coefficients: retained modal coefficients
            % - Returns values: reconstructed values on the projection's reconstruction grid
            arguments
                self InternalModesProjection
                coefficients double
            end

            values = self.reconstructionMatrix * coefficients;
        end

        function S = spectrum(self, coefficients, options)
            % Compute a retained-mode observation spectrum.
            %
            % This method uses the retained rows of the candidate spectral
            % weights. For canonical G modes this is the potential-energy
            % spectrum implied by the basis normalization.
            %
            % - Topic: Analyze observation spectra
            % - Declaration: S = spectrum(self,coefficients,options)
            % - Parameter self: InternalModesProjection instance
            % - Parameter coefficients: retained modal coefficients
            % - Parameter options.horizontalMultiplicity: multiplicity factor, often 1 or 2
            % - Parameter options.requireCanonical: true to reject noncanonical spectra
            % - Returns S: retained-mode spectrum
            arguments
                self InternalModesProjection
                coefficients double
                options.horizontalMultiplicity double = 1
                options.requireCanonical (1,1) logical = true
            end

            S = self.crossSpectrum(coefficients,coefficients,horizontalMultiplicity=options.horizontalMultiplicity, ...
                requireCanonical=options.requireCanonical);
        end

        function S = crossSpectrum(self, coefficientsA, coefficientsB, options)
            % Compute a retained-mode observation cross-spectrum.
            %
            % - Topic: Analyze observation spectra
            % - Declaration: S = crossSpectrum(self,coefficientsA,coefficientsB,options)
            % - Parameter self: InternalModesProjection instance
            % - Parameter coefficientsA: first retained coefficient array
            % - Parameter coefficientsB: second retained coefficient array
            % - Parameter options.horizontalMultiplicity: multiplicity factor, often 1 or 2
            % - Parameter options.requireCanonical: true to reject noncanonical spectra
            % - Returns S: retained-mode cross-spectrum
            arguments
                self InternalModesProjection
                coefficientsA double
                coefficientsB double
                options.horizontalMultiplicity double = 1
                options.requireCanonical (1,1) logical = true
            end

            if options.requireCanonical && ~self.isCanonicalComponent
                error("InternalModesProjection:NoncanonicalSpectrum", ...
                    "A canonical spectrum is not available for this observation projection.");
            end
            retainedWeights = self.spectralWeights(self.retainedModes);
            multiplicity = options.horizontalMultiplicity;
            if isscalar(multiplicity)
                multiplicity = multiplicity * ones(size(retainedWeights));
            else
                multiplicity = multiplicity(:);
            end
            S = (multiplicity .* retainedWeights) .* real(coefficientsA .* conj(coefficientsB));
        end

        function expectedSpectrum = expectedRecoveredSpectrum(self, trueSpectrum)
            % Apply the spectral window to a true candidate-mode spectrum.
            %
            % For uncorrelated modal coefficients, `spectralWindow` maps
            % the candidate spectrum into the expected retained recovered
            % spectrum.
            %
            % - Topic: Analyze observation spectra
            % - Declaration: expectedSpectrum = expectedRecoveredSpectrum(self,trueSpectrum)
            % - Parameter self: InternalModesProjection instance
            % - Parameter trueSpectrum: candidate-mode spectrum
            % - Returns expectedSpectrum: expected retained recovered spectrum
            arguments
                self InternalModesProjection
                trueSpectrum double
            end

            expectedSpectrum = self.spectralWindow * trueSpectrum;
        end
    end

    methods (Static)
        function [retainedModes,diagnostics] = selectResolvableModes(B, weights, rankTolerance, maxConditionNumber)
            % Select observation-resolvable columns with pivoted QR.
            %
            % - Topic: Developer topics
            % - Declaration: [retainedModes,diagnostics] = InternalModesProjection.selectResolvableModes(B,weights,rankTolerance,maxConditionNumber)
            % - Parameter B: sampled observation matrix
            % - Parameter weights: observation weights
            % - Parameter rankTolerance: relative QR pivot tolerance
            % - Parameter maxConditionNumber: maximum retained Gram condition number
            % - Returns retainedModes: selected column numbers
            % - Returns diagnostics: structure with QR diagnostics
            % - Developer: true
            weightedB = sqrt(weights(:)) .* B;
            columnNorms = vecnorm(weightedB,2,1);
            columnNorms(columnNorms == 0) = 1;
            normalizedB = weightedB ./ columnNorms;
            [~,R,pivotOrder] = qr(normalizedB,0);
            if isempty(R)
                retainedModes = zeros(1,0);
                diagnostics.relativePivots = zeros(0,1);
                diagnostics.finalRelativePivot = NaN;
                diagnostics.conditionNumber = NaN;
                return;
            end
            relativePivots = abs(diag(R))/abs(R(1,1));
            nRank = find(relativePivots >= rankTolerance,1,'last');
            if isempty(nRank)
                nRank = 0;
            end

            selectedInPivotOrder = [];
            for iPivot = 1:nRank
                candidateModes = sort([selectedInPivotOrder pivotOrder(iPivot)]);
                candidateConditionNumber = cond(normalizedB(:,candidateModes).'*normalizedB(:,candidateModes));
                if candidateConditionNumber <= maxConditionNumber
                    selectedInPivotOrder = [selectedInPivotOrder pivotOrder(iPivot)]; %#ok<AGROW>
                end
            end

            retainedModes = sort(selectedInPivotOrder);
            diagnostics.relativePivots = relativePivots;
            if isempty(selectedInPivotOrder)
                diagnostics.finalRelativePivot = NaN;
                diagnostics.conditionNumber = NaN;
            else
                diagnostics.finalRelativePivot = relativePivots(length(selectedInPivotOrder));
                diagnostics.conditionNumber = cond(normalizedB(:,retainedModes).'*normalizedB(:,retainedModes));
            end
        end

        function [A,resolutionMatrix,gramMatrix] = weightedProjection(B, weights, retainedModes)
            % Build a weighted projection for selected observation columns.
            %
            % - Topic: Developer topics
            % - Declaration: [A,resolutionMatrix,gramMatrix] = InternalModesProjection.weightedProjection(B,weights,retainedModes)
            % - Parameter B: sampled observation matrix
            % - Parameter weights: observation weights
            % - Parameter retainedModes: retained column numbers
            % - Returns A: weighted projection matrix
            % - Returns resolutionMatrix: coefficient resolution matrix
            % - Returns gramMatrix: retained weighted Gram matrix
            % - Developer: true
            BS = B(:,retainedModes);
            weightedBS = weights(:) .* BS;
            gramMatrix = BS.' * weightedBS;
            A = gramMatrix \ weightedBS.';
            resolutionMatrix = A * B;
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CADimensionProperty('observationIndex','','observation row index');
            propertyAnnotations(end+1) = CADimensionProperty('candidateModeIndex','','candidate mode index');
            propertyAnnotations(end+1) = CADimensionProperty('retainedModeIndex','','retained mode index');
            propertyAnnotations(end+1) = CAPropertyAnnotation('component','component label');
            propertyAnnotations(end+1) = CANumericProperty('observationZ',{'observationIndex'},'m','observation depths');
            propertyAnnotations(end+1) = CANumericProperty('observationMatrix',{'observationIndex','candidateModeIndex'},'','explicit observation matrix');
            propertyAnnotations(end+1) = CANumericProperty('candidateModes',{'candidateModeIndex'},'','candidate mode numbers');
            propertyAnnotations(end+1) = CANumericProperty('retainedModes',{'retainedModeIndex'},'','retained mode numbers');
            propertyAnnotations(end+1) = CANumericProperty('projectionMatrix',{'retainedModeIndex','observationIndex'},'','forward projection matrix');
            propertyAnnotations(end+1) = CANumericProperty('resolutionMatrix',{'retainedModeIndex','candidateModeIndex'},'','resolution matrix');
            propertyAnnotations(end+1) = CANumericProperty('spectralWindow',{'retainedModeIndex','candidateModeIndex'},'','spectral window');
            propertyAnnotations(end+1) = CANumericProperty('gramMatrix',{'retainedModeIndex','retainedModeIndex'},'','retained Gram matrix');
            propertyAnnotations(end+1) = CANumericProperty('conditionNumber',{},'','Gram condition number');
            propertyAnnotations(end+1) = CANumericProperty('spectralWeights',{'candidateModeIndex'},'','candidate spectral weights');
            propertyAnnotations(end+1) = CAPropertyAnnotation('transformStatus','transform construction status');
            propertyAnnotations(end+1) = CAPropertyAnnotation('componentRole','component role');
            propertyAnnotations(end+1) = CANumericProperty('isCanonicalComponent',{},'','canonical spectrum availability');
            propertyAnnotations(end+1) = CANumericProperty('diagnosticRankTolerance',{},'','QR rank tolerance');
            propertyAnnotations(end+1) = CANumericProperty('diagnosticFinalRelativePivot',{},'','final retained relative QR pivot');
        end

        function names = classRequiredPropertyNames()
            names = {'component','observationZ','observationMatrix','candidateModes','retainedModes', ...
                'projectionMatrix','resolutionMatrix','spectralWindow','gramMatrix','conditionNumber','spectralWeights','transformStatus', ...
                'componentRole','isCanonicalComponent','diagnosticRankTolerance','diagnosticFinalRelativePivot'};
        end
    end
end
