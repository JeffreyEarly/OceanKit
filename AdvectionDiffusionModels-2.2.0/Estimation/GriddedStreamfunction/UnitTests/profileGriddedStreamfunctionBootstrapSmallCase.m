function results = profileGriddedStreamfunctionBootstrapSmallCase(options)
% Profile the small gridded-streamfunction bootstrap benchmark case.
%
% ```matlab
% results = profileGriddedStreamfunctionBootstrapSmallCase();
% ```

arguments
    options.nBootstraps (1,1) double {mustBeInteger, mustBePositive} = 4
    options.randomSeed (1,1) double {mustBeInteger, mustBeNonnegative} = 5
    options.maxFunctions (1,1) double {mustBeInteger, mustBePositive} = 25
    options.shouldPrintScores (1,1) logical = true
end

repoRoot = fileparts(fileparts(fileparts(fileparts(mfilename("fullpath")))));
workspaceRoot = fileparts(repoRoot);
dependencyRoot = fullfile(workspaceRoot, "spline-core");
distributionRoot = fullfile(workspaceRoot, "distributions");

addPackageToPath(repoRoot);

if ~isfolder(dependencyRoot)
    error("profileGriddedStreamfunctionBootstrapSmallCase:MissingDependencyRepo", ...
        "Expected sibling dependency repo at %s.", dependencyRoot);
end
if ~isfolder(distributionRoot)
    error("profileGriddedStreamfunctionBootstrapSmallCase:MissingDependencyRepo", ...
        "Expected sibling dependency repo at %s.", distributionRoot);
end

addPackageToPath(distributionRoot);
addPackageToPath(dependencyRoot);

[~, ~, ~, ~, trajectories] = GriddedStreamfunctionBootstrapUnitTests.synchronousLinearFieldData();

profile clear
profile on
cleanupProfile = onCleanup(@() profile("off")); %#ok<NASGU>

bootstrap = GriddedStreamfunctionBootstrap.fromTrajectories( ...
    trajectories, ...
    nBootstraps=options.nBootstraps, ...
    randomSeed=options.randomSeed);
profileInfo = profile("info");
topFunctions = topFunctionTable(profileInfo.FunctionTable, options.maxFunctions);
selfLikeFunctions = selfLikeFunctionTable(profileInfo.FunctionTable, options.maxFunctions);
leafFunctions = leafFunctionTable(profileInfo.FunctionTable, options.maxFunctions);

if options.shouldPrintScores
    scoreTable = table( ...
        bootstrap.scores.uv(:), ...
        bootstrap.scores.strain(:), ...
        bootstrap.scores.zeta(:), ...
        bootstrap.scores.joint(:), ...
        'VariableNames', ["uv", "strain", "zeta", "joint"]);
    disp(scoreTable)
end

results = struct( ...
    "bootstrap", bootstrap, ...
    "profileInfo", profileInfo, ...
    "topFunctions", topFunctions, ...
    "selfLikeFunctions", selfLikeFunctions, ...
    "leafFunctions", leafFunctions);

if nargout == 0
    disp("Top functions by total time:")
    disp(topFunctions)
    disp("Likely self-time hotspots:")
    disp(selfLikeFunctions)
    disp("Likely leaf hotspots by total time:")
    disp(leafFunctions)
end
end

function topFunctions = topFunctionTable(functionTable, maxFunctions)
if isempty(functionTable)
    topFunctions = table('Size', [0 3], ...
        'VariableTypes', ["double", "double", "string"], ...
        'VariableNames', ["totalTime", "numCalls", "functionName"]);
    return
end

totalTime = transpose([functionTable.TotalTime]);
numCalls = transpose([functionTable.NumCalls]);
functionName = transpose(string({functionTable.FunctionName}));

[~, order] = sort(totalTime, "descend");
order = order(1:min(maxFunctions, numel(order)));
topFunctions = table(totalTime(order), numCalls(order), functionName(order), ...
    'VariableNames', ["totalTime", "numCalls", "functionName"]);
end

function leafFunctions = leafFunctionTable(functionTable, maxFunctions)
if isempty(functionTable)
    leafFunctions = table('Size', [0 3], ...
        'VariableTypes', ["double", "double", "string"], ...
        'VariableNames', ["totalTime", "numCalls", "functionName"]);
    return
end

isLeaf = arrayfun(@(entry) isempty(entry.Children), functionTable);
leafEntries = functionTable(isLeaf);

if isempty(leafEntries)
    leafFunctions = table('Size', [0 3], ...
        'VariableTypes', ["double", "double", "string"], ...
        'VariableNames', ["totalTime", "numCalls", "functionName"]);
    return
end

totalTime = transpose([leafEntries.TotalTime]);
numCalls = transpose([leafEntries.NumCalls]);
functionName = transpose(string({leafEntries.FunctionName}));

[~, order] = sort(totalTime, "descend");
order = order(1:min(maxFunctions, numel(order)));
leafFunctions = table(totalTime(order), numCalls(order), functionName(order), ...
    'VariableNames', ["totalTime", "numCalls", "functionName"]);
end

function selfLikeFunctions = selfLikeFunctionTable(functionTable, maxFunctions)
if isempty(functionTable)
    selfLikeFunctions = table('Size', [0 3], ...
        'VariableTypes', ["double", "double", "string"], ...
        'VariableNames', ["selfLikeTime", "numCalls", "functionName"]);
    return
end

selfLikeTime = nan(numel(functionTable), 1);
for iFunction = 1:numel(functionTable)
    childTime = childTotalTime(functionTable(iFunction).Children);
    if ~isnan(childTime)
        selfLikeTime(iFunction) = max(functionTable(iFunction).TotalTime - childTime, 0);
    end
end

valid = isfinite(selfLikeTime);
if ~any(valid)
    selfLikeFunctions = table('Size', [0 3], ...
        'VariableTypes', ["double", "double", "string"], ...
        'VariableNames', ["selfLikeTime", "numCalls", "functionName"]);
    return
end

numCalls = transpose([functionTable(valid).NumCalls]);
functionName = transpose(string({functionTable(valid).FunctionName}));
selfLikeTime = selfLikeTime(valid);

[~, order] = sort(selfLikeTime, "descend");
order = order(1:min(maxFunctions, numel(order)));
selfLikeFunctions = table(selfLikeTime(order), numCalls(order), functionName(order), ...
    'VariableNames', ["selfLikeTime", "numCalls", "functionName"]);
end

function totalTime = childTotalTime(children)
if isempty(children)
    totalTime = 0;
    return
end

if isstruct(children)
    if isfield(children, "TotalTime")
        totalTime = sum([children.TotalTime]);
        return
    end
    if isfield(children, "Time")
        totalTime = sum([children.Time]);
        return
    end
end

totalTime = NaN;
end

function addPackageToPath(repoRoot)
metadataPath = fullfile(repoRoot, "resources", "mpackage.json");
pathEntries = string(strsplit(path, pathsep));

if isfile(metadataPath)
    metadata = jsondecode(fileread(metadataPath));
    for iFolder = 1:length(metadata.folders)
        folderPath = fullfile(repoRoot, metadata.folders(iFolder).path);
        if ~any(pathEntries == folderPath)
            addpath(folderPath);
            pathEntries(end + 1) = folderPath; %#ok<AGROW>
        end
    end
end

if ~any(pathEntries == repoRoot)
    addpath(repoRoot);
end
end
