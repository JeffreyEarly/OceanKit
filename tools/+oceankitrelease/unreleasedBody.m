function body = unreleasedBody(changelogPath)
arguments
    changelogPath (1,1) string {mustBeFile}
end

source = fileread(changelogPath);
[headingStarts,headingEnds] = regexp(source,'(?m)^## \[Unreleased\][ \t]*\r?\n','start','end');
if numel(headingStarts) ~= 1
    error("OceanKitRelease:InvalidUnreleasedSection", ...
        "CHANGELOG.md must contain exactly one '## [Unreleased]' section.");
end

followingText = source(headingEnds+1:end);
nextHeading = regexp(followingText,'(?m)^## \[','once','start');
if isempty(nextHeading)
    rawBody = followingText;
else
    rawBody = followingText(1:nextHeading-1);
end

body = string(regexprep(rawBody,'^\s*|\s*$',''));
if body == ""
    error("OceanKitRelease:EmptyUnreleasedSection", ...
        "CHANGELOG.md must contain a nonempty '## [Unreleased]' section.");
end
end
