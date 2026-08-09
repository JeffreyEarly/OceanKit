function body = promoteUnreleased(changelogPath,version,releaseDate,releaseBodyPath)
arguments
    changelogPath (1,1) string {mustBeFile}
    version (1,1) string
    releaseDate (1,1) string
    releaseBodyPath (1,1) string
end

if isempty(regexp(version,'^\d+\.\d+\.\d+$','once'))
    error("OceanKitRelease:InvalidVersion", "Release version must use x.y.z semantic versioning.");
end
if isempty(regexp(releaseDate,'^\d{4}-\d{2}-\d{2}$','once'))
    error("OceanKitRelease:InvalidReleaseDate", "Release date must use YYYY-MM-DD.");
end

body = oceankitrelease.unreleasedBody(changelogPath);
source = fileread(changelogPath);
[headingStart,headingEnd] = regexp(source,'(?m)^## \[Unreleased\][ \t]*\r?\n','start','end','once');
followingText = source(headingEnd+1:end);
nextHeading = regexp(followingText,'(?m)^## \[','once','start');
if isempty(nextHeading)
    suffix = "";
else
    suffix = string(followingText(nextHeading:end));
end
prefix = string(source(1:headingStart-1));

promoted = prefix + "## [Unreleased]" + newline + newline + ...
    "## [" + version + "] - " + releaseDate + newline + newline + ...
    body + newline + newline + suffix;

writeTextAtomically(releaseBodyPath, body + newline);
writeTextAtomically(changelogPath, promoted);
end

function writeTextAtomically(destination,text)
destinationFolder = string(fileparts(destination));
if destinationFolder ~= "" && ~isfolder(destinationFolder)
    mkdir(destinationFolder);
end
temporaryPath = string(tempname(destinationFolder));
temporaryCleanup = onCleanup(@()deleteIfPresent(temporaryPath));
fileID = fopen(temporaryPath,"w");
if fileID < 0
    error("OceanKitRelease:WriteFailed", "Unable to write %s.", temporaryPath);
end
fileCleanup = onCleanup(@()fclose(fileID));
fwrite(fileID,text);
clear fileCleanup
movefile(temporaryPath,destination,"f");
clear temporaryCleanup
end

function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end
