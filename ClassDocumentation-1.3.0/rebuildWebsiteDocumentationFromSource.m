function rebuildWebsiteDocumentationFromSource(sourceFolder, buildFolder, preservedRelativePaths, options)
arguments
    sourceFolder
    buildFolder
    preservedRelativePaths string = string.empty(0,1)
    options.preservedRelativeDirectories string = string.empty(0,1)
end

sourceFolder = string(sourceFolder);
buildFolder = string(buildFolder);
preservedRelativePaths = normalizeRelativeEntries(preservedRelativePaths, "file");
preservedRelativeDirectories = normalizeRelativeEntries(options.preservedRelativeDirectories, "directory");

if ~isscalar(sourceFolder) || ~isscalar(buildFolder)
    error("rebuildWebsiteDocumentationFromSource:InvalidFolderInput", ...
        "sourceFolder and buildFolder must each be a single path.");
end

if ~isfolder(sourceFolder)
    error("rebuildWebsiteDocumentationFromSource:MissingSourceFolder", ...
        "Website documentation source folder does not exist: %s", sourceFolder);
end

temporaryPreserveFolder = "";
cleanupTemporaryPreserveFolder = []; %#ok<NASGU>
if isfolder(buildFolder) && (~isempty(preservedRelativePaths) || ~isempty(preservedRelativeDirectories))
    temporaryPreserveFolder = string(tempname);
    mkdir(temporaryPreserveFolder);
    cleanupTemporaryPreserveFolder = onCleanup(@() removeTemporaryFolder(temporaryPreserveFolder)); %#ok<NASGU>

    for iPath = 1:numel(preservedRelativePaths)
        relativePath = preservedRelativePaths(iPath);
        sourcePath = fullfile(buildFolder, relativePath);
        if ~isfile(sourcePath)
            error("rebuildWebsiteDocumentationFromSource:MissingPreservedAsset", ...
                "Preserved website asset is missing from the existing docs build: %s", sourcePath);
        end

        temporaryPath = fullfile(temporaryPreserveFolder, relativePath);
        ensureParentFolderExists(temporaryPath);
        copyfile(sourcePath, temporaryPath, "f");
    end

    for iDirectory = 1:numel(preservedRelativeDirectories)
        relativeDirectory = preservedRelativeDirectories(iDirectory);
        sourceDirectory = fullfile(buildFolder, relativeDirectory);
        if ~isfolder(sourceDirectory)
            continue;
        end

        temporaryDirectory = fullfile(temporaryPreserveFolder, relativeDirectory);
        ensureParentFolderExists(temporaryDirectory);
        copyfile(sourceDirectory, temporaryDirectory);
    end
end

if isfolder(buildFolder)
    rmdir(buildFolder, "s");
end

copyfile(sourceFolder, buildFolder);

if temporaryPreserveFolder == ""
    return;
end

for iPath = 1:numel(preservedRelativePaths)
    relativePath = preservedRelativePaths(iPath);
    temporaryPath = fullfile(temporaryPreserveFolder, relativePath);
    targetPath = fullfile(buildFolder, relativePath);
    ensureParentFolderExists(targetPath);
    copyfile(temporaryPath, targetPath, "f");
end

for iDirectory = 1:numel(preservedRelativeDirectories)
    relativeDirectory = preservedRelativeDirectories(iDirectory);
    temporaryDirectory = fullfile(temporaryPreserveFolder, relativeDirectory);
    if ~isfolder(temporaryDirectory)
        continue;
    end

    targetDirectory = fullfile(buildFolder, relativeDirectory);
    if isfolder(targetDirectory)
        rmdir(targetDirectory, "s");
    end
    ensureParentFolderExists(targetDirectory);
    copyfile(temporaryDirectory, targetDirectory);
end
end

function entries = normalizeRelativeEntries(entries, entryKind)
entries = reshape(string(entries), [], 1);
entries = strtrim(entries);
entries = entries(strlength(entries) > 0);

for iEntry = 1:numel(entries)
    entry = replace(entries(iEntry), "\", "/");
    pathSegments = split(entry, "/");
    if startsWith(entry, "/") || ~isempty(regexp(char(entry), "^[A-Za-z]:", "once"))
        error("rebuildWebsiteDocumentationFromSource:AbsolutePreservedPath", ...
            "Preserved %s '%s' must be relative to the docs build folder.", entryKind, entries(iEntry));
    end
    if any(pathSegments == "." | pathSegments == "..")
        error("rebuildWebsiteDocumentationFromSource:InvalidPreservedPath", ...
            "Preserved %s '%s' must not contain '.' or '..' segments.", entryKind, entries(iEntry));
    end

    entries(iEntry) = join(pathSegments, "/");
end

entries = unique(entries, "stable");
end

function ensureParentFolderExists(filePath)
parentFolder = fileparts(filePath);
if strlength(parentFolder) == 0 || isfolder(parentFolder)
    return;
end

mkdir(parentFolder);
end

function removeTemporaryFolder(folderPath)
if folderPath ~= "" && isfolder(folderPath)
    rmdir(folderPath, "s");
end
end
