function ci_release(options)
arguments
    options.rootDir = ".."
    options.bumpType = "none"
    options.notes string = ""
    options.notesFile (1,1) string = ""
    options.shouldBuildWebsiteDocumentation (1,1) logical = false
    options.shouldPackageForDistribution (1,1) logical = false
    options.shouldPromoteUnreleased (1,1) logical = false
    options.shouldRequireDocumentationBuilder (1,1) logical = false
    options.releaseDate (1,1) string = ""
    options.releaseBodyPath (1,1) string = ""
    options.outputRoot (1,1) string = ""
    options.dist_folder (1,1) string = "dist"
    options.excluded_dist_folders string = [".git", ".github", "docs", "tools", "Documentation", "OceanKit"]
end
%CI_RELEASE CI entry point for MPM release.
%   CI_RELEASE(options) where
%       options.bumpType is "patch", "minor", or "major". If it is left
%
%   Steps:
%     1) Bump version in resources/mpackage.json using matlab.mpm.Package
%     2) Run custom documentation build
%     3) Export package root to dist/<name>-<version> for MPM repo
%
%   This script assumes that options.rootDir points at the *package root*,
%   i.e. the folder that contains the resources/mpackage.json file.

% 1) Read package metadata and bump version (semantic x.y.z) if requested

notes = options.notes;
if options.notesFile ~= ""
    if ~isfile(options.notesFile)
        error("ci_release:notesFileNotFound", "Could not find release notes at %s", options.notesFile);
    end
    notes = string(fileread(options.notesFile));
end

mpkgPath = fullfile(options.rootDir, "resources", "mpackage.json");
if ~isfile(mpkgPath)
    error("ci_release:mpackageNotFound", "Could not find mpackage.json at %s", mpkgPath);
end

% Create a matlab.mpm.Package object for this root folder. Modifying
% properties on this object updates mpackage.json for us.
pkg = matlab.mpm.Package(options.rootDir);

pkgName       = string(pkg.Name);
currentVerObj = pkg.Version;
oldVer        = string(currentVerObj);

bumpType = string(options.bumpType);
bumpType = lower(bumpType);

switch bumpType
    case "major"
        newVerObj = matlab.mpm.Version(currentVerObj.Major + 1, 0, 0);
    case "minor"
        newVerObj = matlab.mpm.Version(currentVerObj.Major, currentVerObj.Minor + 1, 0);
    case "patch"
        newVerObj = matlab.mpm.Version(currentVerObj.Major, currentVerObj.Minor, currentVerObj.Patch + 1);
    otherwise
        % "none" or anything else: keep existing version
        newVerObj = currentVerObj;
end

newVer = string(newVerObj);

isVersionBump = any(bumpType == ["major" "minor" "patch"]);
if options.shouldPromoteUnreleased
    if ~isVersionBump
        error("ci_release:promotionRequiresVersionBump", "Promoting Unreleased requires a major, minor, or patch bump.");
    end
    if options.releaseDate == "" || options.releaseBodyPath == ""
        error("ci_release:promotionOutputRequired", "Promoting Unreleased requires a release date and release-body path.");
    end
    changelogPath = fullfile(options.rootDir, "CHANGELOG.md");
    oceankitrelease.unreleasedBody(changelogPath);
end

% If we are actually bumping, assign back to the package so that
% mpackage.json is rewritten by MATLAB Package Manager.
if isVersionBump
    pkg.Version = newVerObj;
    fprintf('Bumping version: %s -> %s (%s)\n', oldVer, string(newVerObj), bumpType);
else
    fprintf('Not bumping version (current version %s)\n', oldVer);
end



% If we bumped the version and have release notes, update the changelog.
if isVersionBump && options.shouldPromoteUnreleased
    changelogPath = fullfile(options.rootDir, "CHANGELOG.md");
    oceankitrelease.promoteUnreleased(changelogPath, newVer, options.releaseDate, options.releaseBodyPath);
elseif isVersionBump && strlength(strtrim(notes)) > 0
    changelogPath = fullfile(options.rootDir, "CHANGELOG.md");
    update_changelog(changelogPath, notes, newVer);
end

%% 2) Run your custom documentation build
if options.shouldBuildWebsiteDocumentation
    % Replace this with your actual doc build entry point
    % e.g. waveVortexDiagnostics_build_docs, or build_docs
    if exist("build_website_documentation","file")
        fprintf('Running documentation builder\n');
        build_website_documentation(rootDir=options.rootDir);
    elseif options.shouldRequireDocumentationBuilder
        error("ci_release:documentationBuilderNotFound", "Requested documentation builder build_website_documentation was not found.");
    end
end

%% 3) Export package root to dist/<name>-<version> for MPM repo

if options.shouldPackageForDistribution == true
    if options.outputRoot == ""
        distDir = fullfile(options.rootDir, options.dist_folder);
    else
        distDir = options.outputRoot;
    end
    if ~isfolder(distDir)
        mkdir(distDir);
    end

    pkgFolderName = pkgName + "-" + newVer;
    targetRoot    = fullfile(distDir, pkgFolderName);

    % Clean any stale output
    if isfolder(targetRoot)
        rmdir(targetRoot, "s");
    end

    fprintf('Exporting package root to %s\n', targetRoot);
    copyfile(options.rootDir, targetRoot);

    % Strip CI-only junk from the exported package
    % (best-effort: ignore errors if these don't exist)
    for iFolder = 1:numel(options.excluded_dist_folders)
        try
            rmdir(fullfile(targetRoot, options.excluded_dist_folders(iFolder)), "s");
        catch
        end
    end

    try
        rmdir(fullfile(targetRoot, "dist"), "s");
    catch
    end

    % Write a small metadata file for the GitHub Action
    metaPath = fullfile(distDir, "mpm_release_metadata.txt");
    fid = fopen(metaPath, "w");
    assert(fid ~= -1, "Could not open metadata file for writing");
    fprintf(fid, "NAME=%s\nVERSION=%s\nFOLDER=%s\n", pkgName, newVer, pkgFolderName);
    fclose(fid);

    fprintf('ci_release complete: %s %s\n', pkgName, newVer);
end
end
