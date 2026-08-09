function paths = changedAuthoringPaths(repositoryRoot,mode)
arguments
    repositoryRoot (1,1) string {mustBeFolder}
    mode (1,1) string {mustBeMember(mode,["legacy","pilot"])}
end

[exitCode,output] = system(sprintf('git -C "%s" status --porcelain=v1 --untracked-files=all',repositoryRoot));
if exitCode ~= 0
    error("OceanKitRelease:GitStatusFailed", "Unable to inspect %s:\n%s", repositoryRoot, output);
end
lines = splitlines(string(output));
lines(lines == "") = [];
paths = strings(numel(lines),1);
for iLine = 1:numel(lines)
    line = lines(iLine);
    if strlength(line) < 4
        error("OceanKitRelease:UnexpectedGitStatus", "Unexpected git status entry: %s", line);
    end
    path = extractAfter(line,3);
    if contains(path," -> ")
        path = extractAfter(path," -> ");
    end
    paths(iLine) = path;
end
paths = unique(paths,'sorted');

if mode == "legacy"
    paths(startsWith(paths,"dist/")) = [];
end

if mode == "pilot"
    allowed = ["resources/mpackage.json";"CHANGELOG.md";"docs/version-history.md"];
else
    allowed = ["resources/mpackage.json";"CHANGELOG.md"];
end
for path = paths'
    isAllowed = any(path == allowed) || (mode == "legacy" && startsWith(path,"docs/"));
    if ~isAllowed
        error("OceanKitRelease:UnexpectedAuthoringDrift", ...
            "Release preparation changed an unexpected path: %s", path);
    end
end
end
