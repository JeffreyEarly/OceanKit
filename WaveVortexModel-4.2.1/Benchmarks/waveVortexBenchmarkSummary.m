function summary = waveVortexBenchmarkSummary(results)
% Render a WaveVortex benchmark result as a Markdown summary.
arguments
    results (1,1) struct
end
lines = ["# WaveVortex benchmark";"";"- Status: `" + results.status + "`";"- Run: `" + results.runId + "`";"- MATLAB: `" + results.environment.matlabRelease + "`";"- Architecture: `" + results.environment.architecture + "`";"";"## Suite scores";"";"| Suite | Backend | Score |";"|---|---|---:|"];
for iSuite = 1:numel(results.suites)
    if ~isModelOperationSuite(results.suites(iSuite))
        continue
    end
    for iScore = 1:numel(results.suites(iSuite).suiteScores)
        score = results.suites(iSuite).suiteScores(iScore);
        lines(end+1) = sprintf("| %s | %s | %.3f |",results.suites(iSuite).id,score.backendId,score.score); %#ok<AGROW>
    end
end
lines = [lines;"";"## Family scores";"";"| Suite | Family | Backend | Score |";"|---|---|---|---:|"];
for iSuite = 1:numel(results.suites)
    if ~isModelOperationSuite(results.suites(iSuite))
        continue
    end
    for iScore = 1:numel(results.suites(iSuite).familyScores)
        score = results.suites(iSuite).familyScores(iScore);
        lines(end+1) = sprintf("| %s | %s | %s | %.3f |",results.suites(iSuite).id,score.id,score.backendId,score.score); %#ok<AGROW>
    end
end
lines = [lines;"";"## Timing and scores";"";"| Suite | Case | Transform | Backend | Median (ms) | Reference score | Same-host speedup | Error |";"|---|---|---|---|---:|---:|---:|---:|"];
for iSuite = 1:numel(results.suites)
    suite = results.suites(iSuite);
    if ~isModelOperationSuite(suite)
        continue
    end
    for iCase = 1:numel(suite.cases)
        benchmarkCase = suite.cases(iCase);
        for iBackend = 1:numel(benchmarkCase.backends)
            backend = benchmarkCase.backends(iBackend);
            lines(end+1) = sprintf("| %s | %s | %s | %s | %.3f | %.3f | %.3f | %.3g |",suite.id,benchmarkCase.id,benchmarkCase.transformId,backend.id,1e3*backend.medianSeconds,backend.caseScore,backend.sameHostSpeedup,backend.relativeError); %#ok<AGROW>
        end
    end
end
lines = [lines;"";"## Construction and cache diagnostics";"";"| Suite | Case | Backend | Construction (s) | First call (ms) | Same-state cache hit (ms) |";"|---|---|---|---:|---:|---:|"];
for iSuite = 1:numel(results.suites)
    suite = results.suites(iSuite);
    if ~isModelOperationSuite(suite)
        continue
    end
    for iCase = 1:numel(suite.cases)
        benchmarkCase = suite.cases(iCase);
        for iBackend = 1:numel(benchmarkCase.backends)
            backend = benchmarkCase.backends(iBackend);
            lines(end+1) = sprintf("| %s | %s | %s | %.3f | %.3f | %.3f |",suite.id,benchmarkCase.id,backend.id,backend.constructionSeconds,1e3*backend.firstCallSeconds,1e3*backend.sameStateCacheHitSeconds); %#ok<AGROW>
        end
    end
end
lines = [lines;"";"## Memory";"";"| Suite | Case | Backend | Persistent increment (MiB) | Peak increment (MiB) | Provider |";"|---|---|---|---:|---:|---|"];
for iSuite = 1:numel(results.suites)
    suite = results.suites(iSuite);
    if ~isModelOperationSuite(suite)
        continue
    end
    for iCase = 1:numel(suite.cases)
        benchmarkCase = suite.cases(iCase);
        for iBackend = 1:numel(benchmarkCase.backends)
            backend = benchmarkCase.backends(iBackend);
            lines(end+1) = sprintf("| %s | %s | %s | %.3f | %.3f | %s |",suite.id,benchmarkCase.id,backend.id,backend.memory.persistentIncrementBytes/2^20,backend.memory.peakIncrementBytes/2^20,backend.memory.provider); %#ok<AGROW>
        end
    end
end
for iSuite = 1:numel(results.suites)
    suite = results.suites(iSuite);
    if isfield(suite,"kind") && string(suite.kind) == "transform-layout"
        lines = [lines;transformLayoutLines(suite)]; %#ok<AGROW>
    end
end
summary = strjoin(lines,newline) + newline;
end

function tf = isModelOperationSuite(suite)
tf = ~isfield(suite,"kind") || string(suite.kind) == "model-operation";
end

function lines = transformLayoutLines(suite)
lines = ["";"## Transform-layout diagnostic";"";"The strict winner is the smallest median. The current `wv-sorted-linear` path remains preferred whenever it is within 3% of that median.";"";"### Extraction and complete-forward winners";"";"| Case | Antialias | Operation | Current (ms) | Strict fastest | Strict (ms) | Preferred | Current / fastest |";"|---|---:|---|---:|---|---:|---|---:|"];
lines = [lines;selectionRows(suite,["extract" "forward-complete"])];
lines = [lines;"";"### Insertion and complete-inverse winners";"";"| Case | Antialias | Operation | Current (ms) | Strict fastest | Strict (ms) | Preferred | Current / fastest |";"|---|---:|---|---:|---|---:|---|---:|"];
lines = [lines;selectionRows(suite,["insert-primary" "insert-conjugate" "insert-complete" "inverse-complete"])];
lines = [lines;"";"### Mapping-array and working-array storage";"";"| Case | Strategy | Mapping arrays (MiB) | Persistent full buffer (MiB) | Real input (MiB) | Full spectrum (MiB) | WV source (MiB) | Timed WV result (MiB) | Timed real result (MiB) |";"|---|---|---:|---:|---:|---:|---:|---:|---:|"];
for iCase = 1:numel(suite.cases)
    benchmarkCase = suite.cases(iCase);
    if benchmarkCase.status ~= "complete"
        continue
    end
    realBytes = ledgerBytes(benchmarkCase.storage,"realInput");
    fullBytes = ledgerBytes(benchmarkCase.storage,"fullSpectrum");
    wvBytes = ledgerBytes(benchmarkCase.storage,"wvArray");
    for iStrategy = 1:numel(benchmarkCase.strategies)
        strategy = benchmarkCase.strategies(iStrategy);
        wvResultBytes = operationBytes(strategy,"forward-complete");
        realResultBytes = operationBytes(strategy,"inverse-complete");
        lines(end+1) = sprintf("| %s | %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |",benchmarkCase.id,strategy.id,strategy.mappingBytes/2^20,strategy.buffer.bytes/2^20,realBytes/2^20,fullBytes/2^20,wvBytes/2^20,wvResultBytes/2^20,realResultBytes/2^20); %#ok<AGROW>
    end
end
lines = [lines;"";"### Correctness and observable copy semantics";"";"| Case | Strategy | Maximum error | Source arrays unchanged | Buffer reused | Timed clearing | Copy status |";"|---|---|---:|---|---|---|---|"];
for iCase = 1:numel(suite.cases)
    benchmarkCase = suite.cases(iCase);
    if benchmarkCase.status ~= "complete"
        lines(end+1) = "| " + benchmarkCase.id + " | failed | NaN | — | — | — | " + benchmarkCase.failure.identifier + " |"; %#ok<AGROW>
        continue
    end
    for iStrategy = 1:numel(benchmarkCase.strategies)
        strategy = benchmarkCase.strategies(iStrategy);
        maximumError = max([strategy.operations.relativeError]);
        lines(end+1) = sprintf("| %s | %s | %.3g | %s | %s | %s | %s |",benchmarkCase.id,strategy.id,maximumError,yesNo(strategy.sourceArraysUnchanged),yesNo(strategy.persistentBufferReused),yesNo(strategy.timedBufferClearing),strategy.copyObservation); %#ok<AGROW>
    end
end
end

function rows = selectionRows(suite,operationIds)
rows = strings(0,1);
for iCase = 1:numel(suite.cases)
    benchmarkCase = suite.cases(iCase);
    if benchmarkCase.status ~= "complete"
        continue
    end
    for operationId = operationIds
        iSelection = find(string({benchmarkCase.selections.operationId}) == operationId,1);
        selection = benchmarkCase.selections(iSelection);
        rows(end+1,1) = sprintf("| %s | %d | %s | %.3f | %s | %.3f | %s | %.3f |",benchmarkCase.id,benchmarkCase.shouldAntialias,operationId,1e3*selection.currentSeconds,selection.strictFastestStrategy,1e3*selection.strictFastestSeconds,selection.preferredStrategy,selection.currentRelativeToFastest); %#ok<AGROW>
    end
end
end

function bytes = ledgerBytes(ledger,name)
iEntry = find(string({ledger.name}) == name,1);
bytes = ledger(iEntry).bytes;
end

function bytes = operationBytes(strategy,operationId)
iOperation = find(string({strategy.operations.id}) == operationId,1);
bytes = strategy.operations(iOperation).resultStorage.bytes;
end

function value = yesNo(tf)
if tf
    value = "yes";
else
    value = "no";
end
end
