function [transform, assessment] = certifiedDiscreteTransform(self, options)
% Select and independently fit a certified scalar modal band.
%
% Each candidate count receives a fresh quadrature-weight fit to exactly
% that band. The search therefore never judges a short prefix using weights
% fitted to a larger, poorly resolved family. The largest Gram-certified
% count is selected first. Enabled leakage or quadratic policies may then
% reduce the count; every reduction is refitted and reassessed until the
% fitted and retained counts agree.
%
% `nPoints` asks this basis to design a mode-root grid. Explicit `z` may be
% accompanied by the `gridDesign` returned by another family's
% `modeRootGrid`, making shared-grid provenance inspectable.
%
% - Topic: Build discrete transforms
% - Declaration: [transform,assessment] = certifiedDiscreteTransform(basisSet,options)
% - Parameter options.nPoints: exact requested mode-root point count
% - Parameter options.z: explicit increasing physical sample points
% - Parameter options.gridDesign: optional provenance returned by `modeRootGrid`
% - Parameter options.gramTolerance: normalized-Gram operator-error tolerance
% - Parameter options.leakageTolerance: optional rejected-mode leakage tolerance
% - Parameter options.quadraticAliasingTolerance: optional quadratic-aliasing tolerance
% - Parameter options.nCheckModes: optional rejected-mode check count
% - Returns transform: independently fitted certified transform
% - Returns assessment: final diagnostics, grid provenance, and count-search table
arguments
    self IMBasisSet
    options.nPoints double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.z (:,1) double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.gridDesign struct = struct.empty
    options.gramTolerance (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative} = 1e-2
    options.leakageTolerance double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.quadraticAliasingTolerance double {mustBeReal, mustBeFinite} = zeros(0,1)
    options.nCheckModes double {mustBeReal, mustBeFinite} = zeros(0,1)
end

IMDiscreteTransformTools.validateOptionalPositiveInteger(options.nPoints,"nPoints");
IMDiscreteTransformTools.validateOptionalPositiveInteger(options.nCheckModes,"nCheckModes");
IMDiscreteTransformTools.validateOptionalPositiveTolerance(options.leakageTolerance,"leakageTolerance");
IMDiscreteTransformTools.validateOptionalPositiveTolerance(options.quadraticAliasingTolerance,"quadraticAliasingTolerance");
hasPointCount = ~isempty(options.nPoints);
hasExplicitPoints = ~isempty(options.z);
if hasPointCount == hasExplicitPoints
    error("IMBasisSet:InvalidDiscretePointSpecification", "Specify exactly one of nPoints or z.");
end
if hasPointCount && ~isempty(options.gridDesign)
    error("IMBasisSet:InvalidDiscretePointSpecification", "gridDesign accompanies explicit z; nPoints creates its own mode-root provenance.");
end

nAvailableModes = size(self.nativeModes,2);
if hasPointCount
    [z,maximumModeCount,gridDesign] = IMDiscreteTransformTools.pointsForExactCount(self,options.nPoints,nAvailableModes);
else
    z = options.z(:);
    maximumModeCount = min(nAvailableModes,length(z));
    gridDesign = IMDiscreteTransformTools.validatedGridDesign(self,z,options.gridDesign);
end

[searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage] = emptySearchColumns();
transform = [];
assessment = [];
for modeCount = maximumModeCount:-1:1
    try
        [candidateTransform,candidateAssessment] = self.discreteTransform(z=z,nModes=modeCount, ...
            gramTolerance=options.gramTolerance,allowRetainedPrefix=true);
        candidateGramError = candidateAssessment.gramPolicy.error(end);
        accepted = candidateAssessment.gramPolicy.accepted(end);
        [searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage] = ...
            appendSearch(searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage, ...
            "gramSearch",modeCount,true,accepted,candidateAssessment.retainedModeCount,candidateGramError,"","");
        if accepted
            transform = candidateTransform;
            assessment = candidateAssessment;
            break
        end
    catch exception
        if ~isRecoverableSearchFailure(exception)
            rethrow(exception)
        end
        [searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage] = ...
            appendSearch(searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage, ...
            "gramSearch",modeCount,false,false,0,NaN,string(exception.identifier),string(exception.message));
    end
end
if isempty(assessment)
    error("IMBasisSet:NoCertifiedDiscreteTransform", "No independently fitted mode count passed the Gram tolerance on the supplied points.");
end

if ~isempty(options.leakageTolerance) || ~isempty(options.quadraticAliasingTolerance)
    modeCount = assessment.retainedModeCount;
    while true
        try
            [transform,assessment] = self.discreteTransform(z=z,nModes=modeCount, ...
                gramTolerance=options.gramTolerance,leakageTolerance=options.leakageTolerance, ...
                quadraticAliasingTolerance=options.quadraticAliasingTolerance,nCheckModes=options.nCheckModes, ...
                allowRetainedPrefix=true);
            accepted = assessment.retainedModeCount == modeCount;
            [searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage] = ...
                appendSearch(searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage, ...
                "policyRefit",modeCount,true,accepted,assessment.retainedModeCount,assessment.gramPolicy.error(end),"","");
            if accepted
                break
            end
            modeCount = assessment.retainedModeCount;
        catch exception
            if ~isRecoverableSearchFailure(exception) || modeCount == 1
                rethrow(exception)
            end
            [searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage] = ...
                appendSearch(searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage, ...
                "policyRefit",modeCount,false,false,0,NaN,string(exception.identifier),string(exception.message));
            modeCount = modeCount-1;
        end
    end
end

certificationSearch = table(searchStage,searchedModeCount,fitSucceeded,candidateAccepted,searchedRetainedModeCount,searchedGramError,failureIdentifier,failureMessage, ...
    VariableNames=["stage","modeCount","fitSucceeded","accepted","retainedModeCount","gramError","failureIdentifier","failureMessage"]);
assessment = assessment.withCertificationMetadata(gridDesign,certificationSearch);
if assessment.weightFitModeCount ~= assessment.retainedModeCount
    error("IMBasisSet:InconsistentCertifiedTransform", "Certified weights must be fitted to the final retained mode count.");
end
end

function [stage,modeCount,fitSucceeded,accepted,retainedModeCount,gramError,failureIdentifier,failureMessage] = emptySearchColumns()
stage = strings(0,1);
modeCount = zeros(0,1);
fitSucceeded = false(0,1);
accepted = false(0,1);
retainedModeCount = zeros(0,1);
gramError = zeros(0,1);
failureIdentifier = strings(0,1);
failureMessage = strings(0,1);
end

function [stage,modeCount,fitSucceeded,accepted,retainedModeCount,gramError,failureIdentifier,failureMessage] = appendSearch( ...
        stage,modeCount,fitSucceeded,accepted,retainedModeCount,gramError,failureIdentifier,failureMessage, ...
        newStage,newModeCount,newFitSucceeded,newAccepted,newRetainedModeCount,newGramError,newFailureIdentifier,newFailureMessage)
stage(end+1,1) = newStage;
modeCount(end+1,1) = newModeCount;
fitSucceeded(end+1,1) = newFitSucceeded;
accepted(end+1,1) = newAccepted;
retainedModeCount(end+1,1) = newRetainedModeCount;
gramError(end+1,1) = newGramError;
failureIdentifier(end+1,1) = newFailureIdentifier;
failureMessage(end+1,1) = newFailureMessage;
end

function tf = isRecoverableSearchFailure(exception)
tf = ismember(string(exception.identifier),["IMBasisSet:DiscreteTransformPolicyFailed", ...
    "IMBasisSet:NoAcceptableDiscreteTransformPrefix","IMBasisSet:QuadratureWeightFitFailed"]);
end
