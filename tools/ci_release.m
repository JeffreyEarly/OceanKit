function ci_release(options)
arguments
    options.rootDir = ".."
    options.bumpType = "none"
    options.notes string = ""
    options.shouldBuildWebsiteDocumentation (1,1) logical = false
    options.shouldPackageForDistribution (1,1) logical = false
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

% If we are actually bumping, assign back to the package
% so that mpackage.json is rewritten by MATLAB Package Manager.
if bumpType == "major" || bumpType == "minor" || bumpType == "patch"
    try
        pkg.Version = newVerObj;
        fprintf('Bumping version: %s -> %s (%s)\n', oldVer, string(newVerObj), bumpType);
    catch ME
        % Dump extended info to the Action log
        fprintf(2, '*** Error setting pkg.Version in ci_release.m ***\n');
        fprintf(2, 'Identifier: %s\n', ME.identifier);
        fprintf(2, '%s\n', getReport(ME, "extended", "hyperlinks", "off"));

        % Optionally also dump some environment info:
        fprintf(2, '\n--- License in use ---\n');
        try
            li = license("inuse");
            disp(li);
        catch
        end

        fprintf(2, '\n--- MATLAB version ---\n');
        try
            ver
        catch
        end

        rethrow(ME);  % keep failing the build, but with a lot more detail
    end
else
    fprintf('Not bumping version (current version %s)\n', oldVer);
end

newVer = string(newVerObj);

% If we bumped the version and have release notes, update the changelog.
if (bumpType == "major" || bumpType == "minor" || bumpType == "patch") && strlength(strtrim(options.notes)) > 0
    changelogPath = fullfile(options.rootDir, "CHANGELOG.md");
    update_changelog(changelogPath, options.notes, newVer);
end

%% 2) Run your custom documentation build
if options.shouldBuildWebsiteDocumentation
    % Replace this with your actual doc build entry point
    % e.g. waveVortexDiagnostics_build_docs, or build_docs
    if exist("build_website_documentation","file")
        fprintf('Running documentation builder\n');
        build_website_documentation(rootDir=options.rootDir);
    end
end

%% 3) Export package root to dist/<name>-<version> for MPM repo

if options.shouldPackageForDistribution == true
    distDir = fullfile(options.rootDir, options.dist_folder);
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