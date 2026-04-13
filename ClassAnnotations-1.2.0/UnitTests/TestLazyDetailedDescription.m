repoRoot = fileparts(fileparts(mfilename('fullpath')));
cleanupRepoPath = addPathIfNeeded(repoRoot);

helperRoot = tempname;
mkdir(helperRoot);
cleanupHelperRoot = onCleanup(@() removeDirectoryIfPresent(helperRoot)); %#ok<NASGU>
cleanupHelperPath = addPathIfNeeded(helperRoot); %#ok<NASGU>
clear buildPropertyAnnotation buildMethodAnnotation

writeTextFile(fullfile(helperRoot, 'buildPropertyAnnotation.m'), { ...
    'function annotation = buildPropertyAnnotation(options)', ...
    'arguments', ...
    '    options.detailedDescription char = ''''', ...
    'end', ...
    'annotation = CAPropertyAnnotation(''sidecarProperty'', ''sidecar-backed property'', detailedDescription=options.detailedDescription);', ...
    'end'});
writeTextFile(fullfile(helperRoot, 'buildMethodAnnotation.m'), { ...
    'function annotation = buildMethodAnnotation(options)', ...
    'arguments', ...
    '    options.detailedDescription char = ''''', ...
    'end', ...
    'annotation = CAMethodAnnotation(''sidecarMethod'', ''sidecar-backed method'', detailedDescription=options.detailedDescription);', ...
    'end'});
rehash

propertySidecarPath = fullfile(helperRoot, 'sidecarProperty.md');
methodSidecarPath = fullfile(helperRoot, 'sidecarMethod.md');
writeTextFile(propertySidecarPath, 'Property detailed description from sidecar.');
writeTextFile(methodSidecarPath, 'Method detailed description from sidecar.');

propertyAnnotation = buildPropertyAnnotation();
assert(isempty(propertyAnnotation.detailedDescription), ...
    'Property sidecars should not load during ordinary annotation construction.');

methodAnnotation = buildMethodAnnotation();
assert(isempty(methodAnnotation.detailedDescription), ...
    'Method sidecars should not load during ordinary annotation construction.');

outerCleanup = CAEnableDetailedDescriptionSidecars(); %#ok<NASGU>
scopedPropertyAnnotation = buildPropertyAnnotation();
assert(strcmp(scopedPropertyAnnotation.detailedDescription, 'Property detailed description from sidecar.'), ...
    'Property sidecars should load when documentation builders enable sidecar loading.');

scopedMethodAnnotation = buildMethodAnnotation();
assert(strcmp(scopedMethodAnnotation.detailedDescription, 'Method detailed description from sidecar.'), ...
    'Method sidecars should load when documentation builders enable sidecar loading.');

explicitPropertyAnnotation = buildPropertyAnnotation(detailedDescription='Explicit property detailed description.');
assert(strcmp(explicitPropertyAnnotation.detailedDescription, 'Explicit property detailed description.'), ...
    'Explicit property detailedDescription text should override any sidecar text.');

explicitMethodAnnotation = buildMethodAnnotation(detailedDescription='Explicit method detailed description.');
assert(strcmp(explicitMethodAnnotation.detailedDescription, 'Explicit method detailed description.'), ...
    'Explicit method detailedDescription text should override any sidecar text.');

innerCleanup = CAEnableDetailedDescriptionSidecars(); %#ok<NASGU>
nestedPropertyAnnotation = buildPropertyAnnotation();
assert(strcmp(nestedPropertyAnnotation.detailedDescription, 'Property detailed description from sidecar.'), ...
    'Nested sidecar-loading scopes should preserve the enabled state.');
clear innerCleanup

afterInnerCleanupPropertyAnnotation = buildPropertyAnnotation();
assert(strcmp(afterInnerCleanupPropertyAnnotation.detailedDescription, 'Property detailed description from sidecar.'), ...
    'Leaving an inner scope should restore the previous enabled state.');
clear outerCleanup

restoredPropertyAnnotation = buildPropertyAnnotation();
assert(isempty(restoredPropertyAnnotation.detailedDescription), ...
    'Sidecar-loading state should restore to disabled after the outer scope ends.');

function cleanup = addPathIfNeeded(pathToAdd)
if contains(path, [pathsep pathToAdd pathsep]) || startsWith(path, [pathToAdd pathsep]) || endsWith(path, [pathsep pathToAdd]) || strcmp(path, pathToAdd)
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
