repoRoot = fileparts(fileparts(mfilename('fullpath')));
cleanupRepoPath = addPathIfNeeded(repoRoot);

helperRoot = tempname;
mkdir(helperRoot);
cleanupHelperRoot = onCleanup(@() removeDirectoryIfPresent(helperRoot));
cleanupHelperPath = addPathIfNeeded(helperRoot);

className = 'DocumentationIndexCollisionTestClass';
writeTextFile(fullfile(helperRoot, strcat(className, '.m')), { ...
    sprintf('classdef %s', className), ...
    '    % Class page summary.', ...
    '    %', ...
    '    % - Topic: Inspect values', ...
    sprintf('    %% - Declaration: classdef %s', className), ...
    '    properties', ...
    '        % Observed index values.', ...
    '        %', ...
    '        % - Topic: Inspect values', ...
    '        index', ...
    '    end', ...
    'end'});

rehash
clear(className)

classDocumentation = ClassDocumentation(className, buildFolder=helperRoot);
classDocumentation.writeToFile();

classIndexPath = fullfile(classDocumentation.pathOfClassFolderOnHardDrive, 'index.md');
propertyIndexPath = fullfile(classDocumentation.pathOfClassFolderOnHardDrive, 'index_.md');

assert(isfile(classIndexPath), 'Class index page should be written to index.md.');
assert(isfile(propertyIndexPath), 'API item named index should be written to index_.md.');

classIndexText = fileread(classIndexPath);
propertyIndexText = fileread(propertyIndexPath);
assert(contains(classIndexText, sprintf('title: %s', className)), ...
    'Class index page should not be overwritten by the index property page.');
assert(contains(classIndexText, 'index_.html'), ...
    'Class topic links should point to the escaped index property page.');
assert(contains(propertyIndexText, 'title: index'), ...
    'Escaped index property page should still document the index item.');

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
