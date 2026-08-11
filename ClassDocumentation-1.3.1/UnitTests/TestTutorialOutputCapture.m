repoRoot = fileparts(fileparts(mfilename('fullpath')));
cleanupRepoPath = addPathIfNeeded(repoRoot); %#ok<NASGU>

helperRoot = tempname;
mkdir(helperRoot);
cleanupHelperRoot = onCleanup(@() removeDirectoryIfPresent(helperRoot)); %#ok<NASGU>
cleanupHelperPath = addPathIfNeeded(helperRoot); %#ok<NASGU>

tutorialPath = fullfile(helperRoot, 'OutputDemo.m');
buildFolder = fullfile(helperRoot, 'build');

writeTutorialSource(tutorialPath);

tutorialDocumentation = TutorialDocumentation.documentationFromSourceFiles( ...
    tutorialPath, ...
    buildFolder=buildFolder, ...
    sourceRoot=helperRoot, ...
    executionPaths=string(helperRoot));
tutorialDocumentation.writeToFile();

markdownPath = fullfile(buildFolder, 'tutorials', 'output-demo.md');
markdownText = fileread(markdownPath);

assert(contains(markdownText, "```text"), ...
    'tutorialOutputCapture should render a fenced text block in the generated markdown.');
assert(contains(markdownText, "count") && contains(markdownText, "label") && contains(markdownText, "alpha"), ...
    'The generated markdown should contain the captured table output.');
assert(~contains(markdownText, "<strong>"), ...
    'Captured output should strip MATLAB rich-display markup tags.');
assert(contains(markdownText, "*Table output from `disp`.*"), ...
    'Captured output captions should be written beneath the fenced block.');

firstCodeIndex = strfind(markdownText, "```matlab" + newline + "values = [1 2 3];");
outputIndex = strfind(markdownText, "```text" + newline);
secondCodeIndex = strfind(markdownText, "```matlab" + newline + "x = linspace(0, 1, 5);");
figureIndex = strfind(markdownText, "![A small figure for output-capture testing.");
assert(~isempty(firstCodeIndex) && ~isempty(outputIndex) && ~isempty(secondCodeIndex) && ~isempty(figureIndex), ...
    'The markdown should contain the surrounding code, output, and figure blocks.');
assert(firstCodeIndex(1) < outputIndex(1) && outputIndex(1) < secondCodeIndex(1) && secondCodeIndex(1) < figureIndex(1), ...
    'Captured output should preserve ordering between surrounding code and figure blocks.');

assert(isfile(fullfile(buildFolder, 'tutorials', 'output-demo', 'demo-figure.png')), ...
    'The tutorial figure should still be generated alongside the captured output block.');

function writeTutorialSource(tutorialPath)
lines = { ...
    '%% Tutorial Metadata', ...
    '% Title: Output demo', ...
    '% Slug: output-demo', ...
    '% Description: Exercise explicit tutorial output capture.', ...
    '', ...
    '%% Output capture', ...
    '% This tutorial exists only to exercise explicit output capture.', ...
    'values = [1 2 3];', ...
    'summaryTable = table(values.'', ["alpha"; "beta"; "gamma"], VariableNames=["count", "label"]);', ...
    'if exist("tutorialOutputCapture", "var") && isa(tutorialOutputCapture, "function_handle"), tutorialOutputCapture(@() disp(summaryTable), Caption="Table output from `disp`."); end', ...
    'x = linspace(0, 1, 5);', ...
    'figure', ...
    'plot(x, x.^2, LineWidth=1.5)', ...
    'xlabel("x")', ...
    'ylabel("x^2")', ...
    'if exist("tutorialFigureCapture", "var") && isa(tutorialFigureCapture, "function_handle"), tutorialFigureCapture("demo-figure", Caption="A small figure for output-capture testing."); end'};
writeTextFile(tutorialPath, lines);
end

function cleanup = addPathIfNeeded(pathToAdd)
if contains(path, [pathsep pathToAdd pathsep]) || startsWith(path, [pathToAdd pathsep]) ...
        || endsWith(path, [pathsep pathToAdd]) || strcmp(path, pathToAdd)
    cleanup = onCleanup(@() []);
else
    addpath(pathToAdd, '-begin');
    cleanup = onCleanup(@() rmpath(pathToAdd));
end
end

function writeTextFile(path, text)
if ischar(text)
    lines = {text};
elseif isstring(text)
    lines = cellstr(text);
else
    lines = text;
end

fileID = fopen(path, 'w');
assert(fileID >= 0, 'Unable to open %s for writing.', path);
cleanup = onCleanup(@() fclose(fileID)); %#ok<NASGU>
for iLine = 1:numel(lines)
    fprintf(fileID, '%s', lines{iLine});
    if iLine < numel(lines)
        fprintf(fileID, '\n');
    end
end
end

function removeDirectoryIfPresent(pathToRemove)
if isfolder(pathToRemove)
    rmdir(pathToRemove, 's');
end
end
