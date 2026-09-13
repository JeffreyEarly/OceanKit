function decision = compiledPreviewBenchmarkDecision(comparison)
% Apply the public compiled-preview availability gate.
arguments
    comparison (:,1) struct
end
required = ["status" "compiledSpeedup" "maximumRelativeError" "libraryIdentityPassed" "nativeExecutionPassed" "noFallback" "lifecyclePassed" "planCount" "persistentFullHermitianBytes"];
for name = required
    if ~all(arrayfun(@(item)isfield(item,name),comparison))
        error("WaveVortexModel:CompiledPreviewBenchmarkDecision","Comparison records are missing %s.",name);
    end
end
complete = string({comparison.status})' == "complete";
correct = [comparison.maximumRelativeError]' <= 1e-12;
native = [comparison.libraryIdentityPassed]' & [comparison.nativeExecutionPassed]' & [comparison.noFallback]';
lifecycle = [comparison.lifecyclePassed]' & [comparison.planCount]' == 17 & [comparison.persistentFullHermitianBytes]' == 0;
speed = [comparison.compiledSpeedup]' >= 1.25;
available = ~isempty(comparison) && all(complete & correct & native & lifecycle & speed);
decision = struct( ...
    "status",conditional(available,"PREVIEW-AVAILABLE","PREVIEW-NOT-AVAILABLE"), ...
    "available",available, ...
    "correctnessPassed",all(correct), ...
    "nativeExecutionPassed",all(native), ...
    "lifecyclePassed",all(lifecycle), ...
    "speedPassed",all(speed), ...
    "minimumSpeedup",min([comparison.compiledSpeedup]), ...
    "maximumRelativeError",max([comparison.maximumRelativeError]), ...
    "maximumExactRetainedRatio",max([comparison.exactRetainedRatio]), ...
    "maximumOperationPeakRSSRatio",max([comparison.operationPeakRSSRatio]));
end

function value = conditional(condition,trueValue,falseValue)
if condition
    value = trueValue;
else
    value = falseValue;
end
end
