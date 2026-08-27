function [results,archive] = revalidateThreeInterfaceBenchmarkArtifact(rawArtifactPath,options)
% Revalidate a completed v2 worker matrix after a correctness-policy fix.
arguments
    rawArtifactPath (1,1) string {mustBeFile}
    options.outputDirectory (1,1) string = ""
    options.archiveDirectory (1,1) string = ""
end
repositoryRoot = string(fileparts(fileparts(mfilename("fullpath"))));
raw = jsondecode(fileread(rawArtifactPath));
if string(raw.schemaVersion)~="three-interface-benchmark-v2" || logical(raw.source.isDirty) || string(raw.status)~="failed"
    error("WaveVortexBenchmark:ArtifactNotRevalidatable","Revalidation requires failed, clean three-interface-benchmark-v2 evidence.")
end
allowedFailures = ["WaveVortexBenchmark:NumericalMismatch" "WaveVortexBenchmark:OutputGraphMismatch" "WaveVortexBenchmark:MatchedContractFailed"];
if ~isfield(raw,"failure") || string(raw.failure.stage)~="correctness" || ~ismember(string(raw.failure.identifier),allowedFailures)
    error("WaveVortexBenchmark:ArtifactNotRevalidatable","Only a final correctness-policy failure can be revalidated without rerunning workers.")
end
expectedRunCount = double(raw.configuration.processRunCount)*3*numel(raw.cases);
if numel(raw.runs)~=expectedRunCount || any(string({raw.runs.status})~="complete")
    error("WaveVortexBenchmark:ArtifactNotRevalidatable","The artifact does not contain the complete fresh-process worker matrix.")
end
[validatorCommit,validatorTree,validatorDirty] = gitIdentity(repositoryRoot);
if validatorDirty
    error("WaveVortexBenchmark:DirtyRevalidationSource","Commit the validation-policy fix before revalidating canonical evidence.")
end
originalArtifactSHA256 = sha256File(rawArtifactPath);
definitions = raw.cases;
recoveries = repmat(struct("caseId","","relativeTolerance",NaN,"absoluteTolerance",NaN,"reclassifiedCategories",strings(0,1),"remainingDifferences",strings(0,1),"passed",false),numel(definitions),1);
for iCase = 1:numel(definitions)
    definitions(iCase).outputRelativeTolerance = double(definitions(iCase).relativeTolerance);
    definitions(iCase).outputAbsoluteTolerance = double(definitions(iCase).absoluteTolerance);
    comparisonIndex = find(string({raw.comparison.id})==string(definitions(iCase).id),1);
    if isempty(comparisonIndex)
        error("WaveVortexBenchmark:ArtifactNotRevalidatable","Case %s has no aggregate comparison.",string(definitions(iCase).id))
    end
    [raw.comparison(comparisonIndex),recovery] = reclassifyThreeInterfaceOutputEvidence(raw.comparison(comparisonIndex),definitions(iCase).outputRelativeTolerance,definitions(iCase).outputAbsoluteTolerance);
    recoveries(iCase) = mergeStruct(struct("caseId",string(definitions(iCase).id)),recovery);
    for iRepeat = 1:numel(raw.repeatComparisonEvidence)
        repeatIndex = find(string({raw.repeatComparisonEvidence(iRepeat).comparison.id})==string(definitions(iCase).id),1);
        if isempty(repeatIndex)
            error("WaveVortexBenchmark:ArtifactNotRevalidatable","Repeat %d has no comparison for %s.",iRepeat,string(definitions(iCase).id))
        end
        raw.repeatComparisonEvidence(iRepeat).comparison(repeatIndex) = reclassifyThreeInterfaceOutputEvidence(raw.repeatComparisonEvidence(iRepeat).comparison(repeatIndex),definitions(iCase).outputRelativeTolerance,definitions(iCase).outputAbsoluteTolerance);
    end
end
if any(~[recoveries.passed])
    error("WaveVortexBenchmark:OutputGraphMismatch","Stored evidence still contains an output mismatch outside matched method tolerances or a structural mismatch.")
end
raw.cases = definitions;
raw.configuration.outputAgreementTolerancePolicy = "matched method RelTol plus base component AbsTol";
raw.revalidation = struct("policy","conservative stored-evidence reclassification; every previously failing category maximum absolute or relative error must fit the matched method tolerance, and structural differences remain failures","originalStatus",string(raw.status),"originalFailure",raw.failure,"originalArtifactSHA256",originalArtifactSHA256,"validatedAtUTC",utcTimestamp,"validatorCommit",validatorCommit,"validatorTree",validatorTree,"validatorDirty",false,"cases",recoveries);
raw.status = "complete";
raw.failure = emptyFailure;
validateThreeInterfaceBenchmarkContract(raw);
if options.outputDirectory==""
    options.outputDirectory = fullfile(fileparts(rawArtifactPath),"revalidated");
end
if isfolder(options.outputDirectory)
    error("WaveVortexBenchmark:ThreeInterfaceOutputExists","Output already exists: %s",options.outputDirectory)
end
mkdir(options.outputDirectory);
outputPath = fullfile(options.outputDirectory,"three-interface-benchmark.json");
writeText(outputPath,jsonencode(raw,PrettyPrint=true));
writeText(fullfile(options.outputDirectory,"summary.md"),summaryMarkdown(raw));
if options.archiveDirectory==""
    options.archiveDirectory = fullfile(fileparts(repositoryRoot),"wave-vortex-model-benchmark-artifacts","three-interface");
end
if ~isfolder(options.archiveDirectory), mkdir(options.archiveDirectory); end
archiveName = string(raw.runId)+"-three-interface-benchmark.json.gz";
generated = gzip(outputPath,options.archiveDirectory);
generatedPath = string(generated{1});
archivePath = fullfile(options.archiveDirectory,archiveName);
if generatedPath~=archivePath, movefile(generatedPath,archivePath,"f"); end
information = dir(archivePath);
archive = struct("fileName",archiveName,"sha256",sha256File(archivePath),"compressedBytes",information.bytes,"location","external author archive; not distributed with source");
results = raw;
end

function [commit,tree,isDirty] = gitIdentity(repositoryRoot)
[status,commit] = system("git -C "+shellQuote(repositoryRoot)+" rev-parse HEAD");
if status~=0, error("WaveVortexBenchmark:GitIdentity","Unable to read the validator commit."); end
[status,tree] = system("git -C "+shellQuote(repositoryRoot)+" rev-parse HEAD^{tree}");
if status~=0, error("WaveVortexBenchmark:GitIdentity","Unable to read the validator tree."); end
[status,changes] = system("git -C "+shellQuote(repositoryRoot)+" status --porcelain --untracked-files=no");
if status~=0, error("WaveVortexBenchmark:GitIdentity","Unable to inspect the validator worktree."); end
commit = strtrim(string(commit));
tree = strtrim(string(tree));
isDirty = strlength(strtrim(string(changes)))>0;
end

function markdown = summaryMarkdown(results)
lines = ["# Matched three-interface benchmark"; ""; "Status: `"+results.status+"`."; ""; "The complete worker matrix was revalidated from stored numerical evidence after applying matched method tolerances. No performance worker was rerun."; ""; "| Case | Interface | Integration (s) | Peak RSS (GiB) |"; "|---|---|---:|---:|"];
for benchmarkCase = reshape(results.comparison,1,[])
    for item = reshape(benchmarkCase.interfaces,1,[])
        lines(end+1) = sprintf('| %s | %s | %.6f | %.3f |',benchmarkCase.id,item.id,item.integrationSeconds,item.totalPeakRSSBytes/2^30); %#ok<AGROW>
    end
end
markdown = strjoin(lines,newline)+newline;
end

function value = mergeStruct(left,right)
value = left;
for name = string(fieldnames(right))'
    value.(name) = right.(name);
end
end

function value = sha256File(pathname)
[status,output] = system("/usr/bin/shasum -a 256 "+shellQuote(pathname));
if status~=0, error("WaveVortexBenchmark:HashFailed","%s",output); end
value = string(extractBefore(strtrim(output),65));
end

function value = shellQuote(value)
singleQuoteEscape = char([39 34 39 34 39]);
value = "'"+replace(string(value),"'",singleQuoteEscape)+"'";
end

function value = utcTimestamp
value = string(datetime("now","TimeZone","UTC","Format","yyyy-MM-dd'T'HH:mm:ss.SSS'Z'"));
end

function value = emptyFailure
value = struct("stage","","identifier","","message","","report","");
end

function writeText(pathname,value)
fileId = fopen(pathname,"w");
if fileId<0, error("WaveVortexBenchmark:WriteFailed","Unable to write %s.",pathname); end
cleanup = onCleanup(@()fclose(fileId));
fprintf(fileId,"%s",value);
clear cleanup
end
