function metadata = writeMetadata(metadataPath,options)
arguments
    metadataPath (1,1) string
    options.mode (1,1) string {mustBeMember(options.mode,["legacy","pilot"])}
    options.workflowRef (1,1) string
    options.workflowSha (1,1) string
    options.sourceStartSha (1,1) string
    options.repositoryRoot (1,1) string {mustBeFolder}
    options.oldVersion (1,1) string
    options.bumpType (1,1) string
    options.releaseDate (1,1) string
    options.documentationChecked (1,1) logical
    options.documentationBuilt (1,1) logical
    options.unreleasedPromoted (1,1) logical
    options.shouldPackageForDistribution (1,1) logical
    options.exportRoot (1,1) string = ""
    options.releaseBodyPath (1,1) string = ""
end

manifest = jsondecode(fileread(fullfile(options.repositoryRoot,"resources","mpackage.json")));
packageName = string(manifest.name);
version = string(manifest.version);
snapshotFolder = packageName + "-" + version;
distributionCreated = false;
exportPath = string(missing);
snapshotValue = string(missing);
if options.shouldPackageForDistribution
    candidate = fullfile(options.exportRoot,snapshotFolder);
    if ~isfolder(candidate)
        error("OceanKitRelease:DistributionNotCreated", "Requested distribution was not created at %s.", candidate);
    end
    exportManifestPath = fullfile(candidate,"resources","mpackage.json");
    if ~isfile(exportManifestPath)
        error("OceanKitRelease:DistributionManifestNotFound", "The distribution has no package manifest at %s.", exportManifestPath);
    end
    exportManifest = jsondecode(fileread(exportManifestPath));
    if string(exportManifest.name) ~= packageName || string(exportManifest.version) ~= version
        error("OceanKitRelease:DistributionMetadataMismatch", ...
            "The exported manifest does not match %s %s.", packageName, version);
    end
    distributionCreated = true;
    exportPath = canonicalPath(string(candidate));
    snapshotValue = snapshotFolder;
end

bodyPath = string(missing);
if options.releaseBodyPath ~= ""
    if ~isfile(options.releaseBodyPath)
        error("OceanKitRelease:ReleaseBodyNotFound", "Release body was not created at %s.", options.releaseBodyPath);
    end
    bodyPath = canonicalPath(options.releaseBodyPath);
end

authoringPaths = oceankitrelease.changedAuthoringPaths(options.repositoryRoot,options.mode);
metadata = struct( ...
    "schemaVersion",1, ...
    "mode",options.mode, ...
    "workflowRef",options.workflowRef, ...
    "workflowSha",options.workflowSha, ...
    "sourceStartSha",options.sourceStartSha, ...
    "packageName",packageName, ...
    "oldVersion",options.oldVersion, ...
    "version",version, ...
    "bumpType",options.bumpType, ...
    "releaseDate",options.releaseDate, ...
    "versionChanged",version ~= options.oldVersion, ...
    "documentationChecked",options.documentationChecked, ...
    "documentationBuilt",options.documentationBuilt, ...
    "unreleasedPromoted",options.unreleasedPromoted, ...
    "distributionCreated",distributionCreated, ...
    "snapshotFolder",snapshotValue, ...
    "exportPath",exportPath, ...
    "releaseBodyPath",bodyPath);
metadata.authoringPaths = cellstr(authoringPaths);

metadataFolder = string(fileparts(metadataPath));
if metadataFolder ~= "" && ~isfolder(metadataFolder)
    mkdir(metadataFolder);
end
temporaryPath = string(tempname(metadataFolder));
temporaryCleanup = onCleanup(@()deleteIfPresent(temporaryPath));
fileID = fopen(temporaryPath,"w");
if fileID < 0
    error("OceanKitRelease:MetadataWriteFailed", "Unable to write %s.", temporaryPath);
end
fileCleanup = onCleanup(@()fclose(fileID));
fwrite(fileID,jsonencode(metadata,PrettyPrint=true));
clear fileCleanup
movefile(temporaryPath,metadataPath,"f");
clear temporaryCleanup
fprintf("Release metadata: %s\n",metadataPath);
end

function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end

function path = canonicalPath(path)
path = string(java.io.File(char(path)).getCanonicalPath());
end
