function scores = waveVortexBenchmarkScores(referenceSeconds,measuredSeconds,familyIds)
% Calculate case, family, and suite benchmark scores.
arguments
    referenceSeconds (:,1) double {mustBePositive}
    measuredSeconds (:,1) double {mustBePositive}
    familyIds (:,1) string
end
if numel(referenceSeconds) ~= numel(measuredSeconds) || numel(referenceSeconds) ~= numel(familyIds)
    error("WaveVortexBenchmark:ScoreSizeMismatch","Reference times, measured times, and family IDs must have equal lengths.");
end
caseScores = 100*referenceSeconds./measuredSeconds;
families = unique(familyIds,"stable");
familyScores = NaN(size(families));
for iFamily = 1:numel(families)
    familyScores(iFamily) = exp(mean(log(caseScores(familyIds == families(iFamily)))));
end
scores = struct("caseScores",caseScores,"familyIds",families,"familyScores",familyScores,"suiteScore",exp(mean(log(familyScores))));
end
