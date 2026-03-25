function build_website_documentation(options)
arguments
    options.rootDir = ".."
end

rootDir = char(java.io.File(char(options.rootDir)).getCanonicalPath());
buildFolder = fullfile(rootDir, "docs");
sourceFolder = fullfile(rootDir, "Documentation", "WebsiteDocumentation");

if ~isfolder(sourceFolder)
    error("build_website_documentation:SourceMissing", ...
        "Could not find source documentation at %s", sourceFolder);
end

if isfolder(buildFolder)
    rmdir(buildFolder, "s");
end

copyfile(sourceFolder, buildFolder);

changelogPath = fullfile(rootDir, "CHANGELOG.md");
if isfile(changelogPath)
    header = "---" + newline + ...
             "layout: default" + newline + ...
             "title: Version History" + newline + ...
             "nav_order: 100" + newline + ...
             "---" + newline + newline;
    versionHistoryText = header + fileread(changelogPath);
    versionHistoryFilePath = fullfile(buildFolder, "version-history.md");
    fid = fopen(versionHistoryFilePath, "w");
    assert(fid ~= -1, "Could not open version-history.md for writing");
    fwrite(fid, versionHistoryText);
    fclose(fid);
end
end
