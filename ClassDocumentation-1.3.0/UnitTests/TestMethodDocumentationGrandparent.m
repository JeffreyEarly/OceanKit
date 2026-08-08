repoRoot = fileparts(fileparts(mfilename('fullpath')));
cleanupRepoPath = addPathIfNeeded(repoRoot);

helperRoot = tempname;
mkdir(helperRoot);
cleanupHelperRoot = onCleanup(@() removeDirectoryIfPresent(helperRoot));
cleanupHelperPath = addPathIfNeeded(helperRoot);

className = 'DocumentationGrandparentTestClass';
writeTextFile(fullfile(helperRoot, strcat(className, '.m')), { ...
    sprintf('classdef %s', className), ...
    '    % Class page summary.', ...
    '    %', ...
    '    % - Topic: Inspect values', ...
    sprintf('    %% - Declaration: classdef %s', className), ...
    '    methods', ...
    '        function value = sampleMethod(~)', ...
    '            % Return a sample value.', ...
    '            %', ...
    '            % - Topic: Inspect values', ...
    '            value = 1;', ...
    '        end', ...
    '    end', ...
    'end'});

rehash
clear(className)

customDocumentation = ClassDocumentation(className, buildFolder=helperRoot, websiteFolder="custom", parent="Transforms", grandparent="Class documentation", methodGrandparent="Transforms");
customDocumentation.writeToFile();
customClassText = fileread(fullfile(customDocumentation.pathOfClassFolderOnHardDrive, "index.md"));
customMethodText = fileread(fullfile(customDocumentation.pathOfClassFolderOnHardDrive, "samplemethod.md"));
assert(contains(customClassText, "parent: Transforms") && contains(customClassText, "grand_parent: Class documentation"), ...
    "Class pages should retain their configured parent and grandparent.");
assert(contains(customMethodText, sprintf("parent: %s", className)) && contains(customMethodText, "grand_parent: Transforms"), ...
    "Method pages should use the configured method grandparent.");

defaultDocumentation = ClassDocumentation(className, buildFolder=helperRoot, websiteFolder="default");
defaultDocumentation.writeToFile();
defaultMethodText = fileread(fullfile(defaultDocumentation.pathOfClassFolderOnHardDrive, "samplemethod.md"));
assert(contains(defaultMethodText, "grand_parent: Classes"), ...
    "Method pages should preserve the existing default grandparent.");

batchDocumentation = ClassDocumentation.classDocumentationFromClassNames({className}, buildFolder=helperRoot, websiteFolder="batch", parent="Transforms", grandparent="Class documentation", methodGrandparent="Transforms");
batchDocumentation(1).writeToFile();
batchMethodText = fileread(fullfile(batchDocumentation(1).pathOfClassFolderOnHardDrive, "samplemethod.md"));
assert(contains(batchMethodText, "grand_parent: Transforms"), ...
    "Batch construction should forward the configured method grandparent.");

function cleanup = addPathIfNeeded(pathToAdd)
if contains(path, [pathsep pathToAdd pathsep]) || startsWith(path, [pathToAdd pathsep]) || ...
        endsWith(path, [pathsep pathToAdd]) || strcmp(path, pathToAdd)
    cleanup = onCleanup(@() []);
else
    addpath(pathToAdd, '-begin');
    cleanup = onCleanup(@() rmpath(pathToAdd));
end
end

function writeTextFile(path, text)
fileID = fopen(path, 'w');
assert(fileID >= 0, 'Unable to open %s for writing.', path);
cleanup = onCleanup(@() fclose(fileID));
for iLine = 1:numel(text)
    fprintf(fileID, '%s', text{iLine});
    if iLine < numel(text)
        fprintf(fileID, '\n');
    end
end
end

function removeDirectoryIfPresent(pathToRemove)
if isfolder(pathToRemove)
    rmdir(pathToRemove, 's');
end
end
