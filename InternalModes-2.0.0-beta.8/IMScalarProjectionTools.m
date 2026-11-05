classdef (Hidden, Sealed) IMScalarProjectionTools
    % Share fixed-rule physical scalar adapters across capable family bases.

    methods (Static)
        function physical = projectionOnGrid(basis,z,weights,options)
            arguments (Input)
                basis
                z (:,1) double {mustBeReal,mustBeFinite}
                weights (:,1) double {mustBeReal,mustBeFinite}
                options.variable (1,1) string = ""
                options.columns (1,:) double {mustBeInteger,mustBePositive} = 1:basis.nColumns
                options.prepareRecovery (1,1) logical = false
                options.failureMode (1,1) string {mustBeMember(options.failureMode,["diagnostic","error"])} = "diagnostic"
                options.rankTolerance (1,1) double {mustBeReal,mustBeFinite,mustBePositive} = 1e-12
            end
            arguments (Output)
                physical (1,1) IMProjection
            end
            IMScalarProjectionTools.requirePhysicalVariable(basis,options.variable);
            if isprop(basis,"evp") && isa(basis.evp,"IMInternalModes")
                spec = basis.evp.innerProduct(options.variable);
                for term = reshape([spec.surfaceWeights;spec.bottomWeights],1,[])
                    if term.coefficient ~= 0 && term.d ~= 0
                        error("IMProjection:UnsupportedDerivativeTrace","The physical pairing requires independent endpoint derivative observations; scalar values do not supply them.");
                    end
                end
                for term = reshape(spec.endpointInnerProductTerms,1,[])
                    if term.coefficient ~= 0 && string(term.variable) ~= options.variable
                        error("IMProjection:UnsupportedCompanionTrace","The physical pairing requires independent companion endpoint observations; scalar values do not supply them.");
                    end
                end
            end
            selected = basis.selectColumns(options.columns);
            if isa(selected,"IMNeutralLabelModesBasis")
                physical = IMScalarProjectionTools.neutralProjection(selected,z,weights);
            else
                recipe = selected.projectionRecipe(variable=options.variable);
                physical = IMProjection.fromRecipe(recipe,z,weights,prepareRecovery=options.prepareRecovery,failureMode=options.failureMode,rankTolerance=options.rankTolerance);
            end
        end

        function transform = scalarTransform(basis,z,weights,options)
            arguments (Input)
                basis
                z (:,1) double {mustBeReal,mustBeFinite}
                weights (:,1) double {mustBeReal,mustBeFinite}
                options.variable (1,1) string = ""
                options.recovery = []
                options.policy (1,1) string {mustBeMember(options.policy,["correctedPhysicalDual","sampledLeastSquares"])} = "correctedPhysicalDual"
                options.coefficientScale (:,:) double {mustBeFinite} = []
                options.observationFactor (:,:) double {mustBeFinite} = []
                options.fitFactor (:,:) double {mustBeFinite} = []
                options.rankTolerance (1,1) double {mustBeReal,mustBeFinite,mustBePositive} = 1e-12
                options.failureMode (1,1) string {mustBeMember(options.failureMode,["diagnostic","error"])} = "diagnostic"
            end
            arguments (Output)
                transform (1,1) IMScalarDiscreteTransform
            end
            physical = IMScalarProjectionTools.projectionOnGrid(basis,z,weights,variable=options.variable,prepareRecovery=false);
            synthesis = basis.evaluate(z,variable=options.variable);
            if isempty(options.recovery)
                if options.policy == "sampledLeastSquares" && isempty(options.fitFactor) && isempty(options.observationFactor)
                    error("IMScalarProjection:MissingFitFactor","sampledLeastSquares requires an explicit fitFactor or observationFactor; physical quadrature does not declare a fitting or noise norm.");
                end
                coefficientScale = options.coefficientScale;
                if isempty(coefficientScale)
                    active = physical.activeColumnMask;
                    majorant = physical.majorantGramMatrix(active,active);
                    diagonalScale = 1./sqrt(diag(majorant));
                    if isempty(diagonalScale)
                        coefficientScale = zeros(0);
                    else
                        [factor,notPositive] = chol(diagonalScale.*majorant.*diagonalScale.');
                        if notPositive ~= 0
                            error("IMScalarProjection:InvalidMajorant","The active continuous scalar majorant must be positive definite to define coefficient units.");
                        end
                        coefficientScale = diag(diagonalScale)/factor;
                    end
                end
                recovery = IMDiscreteTransform.fromObservationMatrices(synthesis,policy=options.policy,columnLabels=basis.columnLabels,activeColumnMask=physical.activeColumnMask,coefficientScale=coefficientScale,observationFactor=options.observationFactor,fitFactor=options.fitFactor,physicalDual=physical.physicalDualMatrix,rankTolerance=options.rankTolerance,failureMode=options.failureMode,metadata=struct(construction="scalarTransform",variable=options.variable,z=z,normalization=string(basis.normalization),physicalProjectionProvenance=physical.provenance));
            else
                if ~isa(options.recovery,"IMDiscreteTransform") || ~isscalar(options.recovery)
                    error("IMScalarProjection:InvalidPreparedRecovery","recovery must be one engine-prepared IMDiscreteTransform.");
                end
                if ~isempty(options.coefficientScale) || ~isempty(options.observationFactor) || ~isempty(options.fitFactor) || options.policy ~= "correctedPhysicalDual" || options.rankTolerance ~= 1e-12 || options.failureMode ~= "diagnostic"
                    error("IMScalarProjection:ConflictingRecoveryOptions","Supply either a prepared recovery or recovery preparation options.");
                end
                recovery = options.recovery;
            end
            if ~isequal(recovery.synthesisMatrix,synthesis)
                error("IMScalarProjection:IncompatiblePreparedRecovery","The prepared recovery must use the exact requested scalar synthesis and observation order.");
            end
            if isfield(recovery.metadata,"variable") && string(recovery.metadata.variable) ~= options.variable
                error("IMScalarProjection:IncompatiblePreparedRecovery","The prepared recovery variable must match the physical scalar channel.");
            end
            legacyModes = [];
            if isprop(basis,"modeNumber"), legacyModes = basis.modeNumber; end
            transform = IMScalarDiscreteTransform(physical,recovery,z=z,weights=weights,normalization=string(basis.normalization),modeNumber=legacyModes);
        end
    end

    methods (Static, Access = private)
        function requirePhysicalVariable(basis,variable)
            if variable == ""
                error("IMScalarProjection:MissingVariable","Specify the independent scalar observations with variable=.");
            end
            if isa(basis,"IMNeutralLabelModesBasis")
                supported = variable == "APV";
            elseif isa(basis,"IMInternalModesBasis") || isa(basis,"IMAnalyticalInternalModesBasis")
                name = string(basis.evp.name);
                supported = (basis.evp.modeFamily == "hydrostatic" || ismember(name,["hydrostaticFModes","hydrostaticGModes","geostrophicAPVModes"])) && ismember(variable,["F","G"]);
                supported = supported || name == "meanDensityAnomalyModes" && variable == "G";
            else
                supported = isa(basis,"IMBasisSet") && variable == "u";
            end
            if ~supported
                error("IMScalarProjection:UnavailableProjection","This family does not define an independent physical scalar %s projector. Evaluable profiles may still support explicitly prepared observation recovery.",variable);
            end
        end

        function projection = neutralProjection(basis,z,weights)
            if basis.f0 == 0
                error("IMNeutralLabelModesBasis:ZeroRotationProjection","APV is identically zero at f0=0 and cannot identify neutral-label amplitudes.");
            end
            if numel(z) < 2 || any(diff(z) <= 0)
                error("IMNeutralLabelModesBasis:InvalidGrid","z must contain at least two strictly increasing points.");
            end
            if numel(weights) ~= numel(z) || ~any(weights ~= 0)
                error("IMNeutralLabelModesBasis:InvalidWeights","Supply one weight per point and at least one nonzero weight.");
            end
            derivatives = basis.Gz(z,normalization="rmsG");
            factors = basis.normalizationFactors();
            pairing = (-derivatives.'./(basis.f0*diff(basis.zDomain)*basis.m_n(:).^2)).*factors(:).*weights.';
            if any(~isfinite(pairing),"all")
                error("IMNeutralLabelModesBasis:UnrepresentableProjection","The APV dual exceeds floating-point range for this rotation and geometry.");
            end
            recipe = basis.projectionRecipe(variable="APV");
            majorant = recipe.majorantGramMatrix;
            if any(diag(majorant) <= 0) || any(~isfinite(majorant),"all")
                error("IMNeutralLabelModesBasis:UnrepresentableProjection","The physical APV norm is outside floating-point range for this rotation and geometry.");
            end
            provenance = recipe.provenance;
            provenance.continuousAPVGramMatrix = majorant;
            projection = IMProjection.fromPrescribedDual(pairing,eye(basis.nColumns),eye(basis.nColumns),majorantGramMatrix=majorant,columnLabels=basis.columnLabels,provenance=provenance);
        end
    end
end
