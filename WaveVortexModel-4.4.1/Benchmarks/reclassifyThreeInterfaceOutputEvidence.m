function [comparison,recovery] = reclassifyThreeInterfaceOutputEvidence(comparison,relativeTolerance,absoluteTolerance)
% Reclassify stored numerical payload evidence at matched method tolerances.
arguments
    comparison (1,1) struct
    relativeTolerance (1,1) double {mustBePositive}
    absoluteTolerance (1,1) double {mustBeNonnegative}
end
if ~isfield(comparison,"outputGraph") || ~isfield(comparison.outputGraph,"categories") || ~isfield(comparison.outputGraph,"differences")
    error("WaveVortexBenchmark:IncompleteOutputEvidence","Stored output-graph evidence is incomplete.")
end
graph = comparison.outputGraph;
categories = graph.categories;
reclassifiedCategories = strings(0,1);
for iCategory = 1:numel(categories)
    category = categories(iCategory);
    finiteEvidence = isfinite(category.maximumAbsoluteError) && isfinite(category.maximumRelativeError);
    safelyWithinTolerance = finiteEvidence && (category.maximumAbsoluteError<=absoluteTolerance || category.maximumRelativeError<=relativeTolerance);
    if ~logical(category.passed) && safelyWithinTolerance
        categories(iCategory).passed = true;
        reclassifiedCategories(end+1,1) = string(category.name); %#ok<AGROW>
    end
end
if ischar(graph.differences)
    differences = string(graph.differences);
else
    differences = string(graph.differences(:));
end
payloadDifference = endsWith(differences," payload differs");
if all([categories.passed])
    differences = differences(~payloadDifference);
end
graph.categories = categories;
graph.differences = differences;
graph.passed = all([categories.passed]) && isempty(differences);
comparison.outputGraph = graph;
comparison.outputAgreementPassed = graph.passed;
requiredFlags = [logical(comparison.integratorAgreementPassed) logical(comparison.adaptiveWorkAgreementPassed) logical(comparison.absoluteToleranceFingerprintAgreementPassed) logical(comparison.memoryAgreementPassed)];
if isfield(comparison,"endpointTrajectoryAgreementPassed")
    requiredFlags(end+1) = logical(comparison.endpointTrajectoryAgreementPassed);
end
comparison.matchedContractPassed = all(requiredFlags) && graph.passed;
recovery = struct("relativeTolerance",relativeTolerance,"absoluteTolerance",absoluteTolerance,"reclassifiedCategories",reclassifiedCategories,"remainingDifferences",differences,"passed",graph.passed);
end
