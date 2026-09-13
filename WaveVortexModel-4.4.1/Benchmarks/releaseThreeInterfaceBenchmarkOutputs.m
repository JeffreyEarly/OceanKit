function [runs,releasedBytes] = releaseThreeInterfaceBenchmarkOutputs(runs,comparison,workFolder)
% Require all repeat gates before deleting any completed worker output.
if isempty(comparison) || ~all([comparison.matchedContractPassed]) || ~all([comparison.outputAgreementPassed]) || ~all([comparison.endpointTrajectoryAgreementPassed])
    error("WaveVortexBenchmark:ThreeInterfaceRepeatComparison","Repeat comparison failed; all current-repeat outputs are retained for investigation.");
end
releasedBytes = 0;
for iRun = 1:numel(runs)
    if string(runs(iRun).status) ~= "complete"
        continue
    end
    pathname = string(runs(iRun).output.path);
    runs(iRun).output.retention = "released-after-repeat-correctness";
    runs(iRun).output.releasedBytes = 0;
    if pathname == "" || ~isfile(pathname)
        continue
    end
    if ~startsWith(pathname,workFolder+filesep)
        error("WaveVortexBenchmark:UnsafeOutputRelease","Refusing to release a validated output outside the benchmark temporary directory: %s",pathname);
    end
    information = dir(pathname);
    runs(iRun).output.releasedBytes = information.bytes;
    releasedBytes = releasedBytes+information.bytes;
    delete(pathname);
end
end

