function analysis = profileCodeHotspots(runFcn, options)
% Profile MATLAB code and rank actionable hotspot targets.
%
% Use `profileCodeHotspots` to profile a function-handle workload and turn
% the profiler output into ranked raw tables plus project-focused hotspot
% views that are easier to act on first.
%
% ```matlab
% analysis = profileCodeHotspots(@() myWorkload(), projectRoots=pwd);
% ```

arguments
    runFcn (1,1) function_handle
    options.projectRoots = string(pwd)
    options.maxFunctions (1,1) double {mustBeInteger, mustBePositive} = 25
    options.maxLines (1,1) double {mustBeInteger, mustBePositive} = 25
    options.label = ""
    options.shouldPrintReport = []
    options.shouldIncludeToolboxInPriorityTargets (1,1) logical = false
end

shouldPrintReport = options.shouldPrintReport;
if isempty(shouldPrintReport)
    shouldPrintReport = (nargout == 0);
end

profile clear
profile on
cleanupProfile = onCleanup(@() profile("off"));

startTime = tic;
runFcn();
elapsedTime = toc(startTime);
profileInfo = profile("info");

analysis = analyzeProfileInfo(profileInfo, ...
    projectRoots=options.projectRoots, ...
    maxFunctions=options.maxFunctions, ...
    maxLines=options.maxLines, ...
    label=options.label, ...
    elapsedTime=elapsedTime, ...
    shouldIncludeToolboxInPriorityTargets=options.shouldIncludeToolboxInPriorityTargets);

if shouldPrintReport
    printAnalysisSummary(analysis);
end
end

function printAnalysisSummary(analysis)
disp("Profile hotspot summary:")
for iSummary = 1:numel(analysis.summary)
    disp("  - " + analysis.summary(iSummary))
end

disp("Project hotspots by self time:")
disp(selectFunctionDisplayColumns(analysis.topProjectBySelfTime))

disp("Project hotspots by total time:")
disp(selectFunctionDisplayColumns(analysis.topProjectByTotalTime))

disp("Project hotspots by number of calls:")
disp(selectFunctionDisplayColumns(analysis.topProjectByNumCalls))

if ~isempty(analysis.topActionableLines)
    disp("Actionable compute lines:")
    disp(selectLineDisplayColumns(analysis.topActionableLines))
end

if ~isempty(analysis.topCallsiteLines)
    disp("Callsite lines:")
    disp(selectLineDisplayColumns(analysis.topCallsiteLines))
end

if ~isempty(analysis.topExternalBySelfTime)
    disp("External or dependency hotspots by self time:")
    disp(selectFunctionDisplayColumns(analysis.topExternalBySelfTime))
end

if ~isempty(analysis.priorityTargets)
    disp("Priority targets:")
    disp(analysis.priorityTargets)
end

disp("Raw functions by self time:")
disp(selectFunctionDisplayColumns(analysis.topBySelfTime))

if ~isempty(analysis.topLines)
    disp("Raw executed lines:")
    disp(selectLineDisplayColumns(analysis.topLines))
end
end

function displayTable = selectFunctionDisplayColumns(functionTable)
if isempty(functionTable)
    displayTable = functionTable;
    return
end

availableNames = string(functionTable.Properties.VariableNames);
preferredNames = ["selfTime", "totalTime", "numCalls", "selfTimePerCall", ...
    "totalTimePerCall", "completeName"];
displayNames = preferredNames(ismember(preferredNames, availableNames));
displayTable = functionTable(:, displayNames);
end

function displayTable = selectLineDisplayColumns(lineTable)
if isempty(lineTable)
    displayTable = lineTable;
    return
end

availableNames = string(lineTable.Properties.VariableNames);
preferredNames = ["lineTime", "lineFractionOfFunction", "lineKind", ...
    "calledTargetName", "lineNumber", "functionName", "fileName"];
displayNames = preferredNames(ismember(preferredNames, availableNames));
displayTable = lineTable(:, displayNames);
end
