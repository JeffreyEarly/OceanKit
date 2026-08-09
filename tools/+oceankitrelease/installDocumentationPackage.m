function dependency = installDocumentationPackage(specifier,options)
arguments
    specifier (1,1) string
    options.repositoryRoot (1,1) string = ""
end

tokens = regexp(specifier,'^([^@]+)@(\d+\.\d+\.\d+)$','tokens','once');
if isempty(tokens)
    error("OceanKitRelease:InvalidDocumentationPackageSpecifier", ...
        "Documentation package specifier must use Name@x.y.z.");
end
expectedName = string(tokens{1});
expectedVersion = string(tokens{2});

packages = mpmlist;
matches = find(string([packages.Name]) == expectedName);
if isempty(matches)
    mpminstall(specifier,Prompt=false);
    packages = mpmlist;
    matches = find(string([packages.Name]) == expectedName);
end
if numel(matches) ~= 1
    error("OceanKitRelease:DocumentationPackageNotResolved", ...
        "Expected exactly one installed %s package.", expectedName);
end
resolved = packages(matches);
actualVersion = string(resolved.Version);
if actualVersion ~= expectedVersion
    error("OceanKitRelease:DocumentationPackageVersionMismatch", ...
        "Requested %s, but version %s resolved from %s.", specifier, actualVersion, string(resolved.PackageRoot));
end

actualRoot = canonicalPath(string(resolved.PackageRoot));
if options.repositoryRoot ~= ""
    expectedRoot = canonicalPath(fullfile(options.repositoryRoot,expectedName + "-" + expectedVersion));
    if actualRoot ~= expectedRoot
        error("OceanKitRelease:DocumentationPackagePathMismatch", ...
            "Requested %s from %s, but the installed package resolves from %s.", ...
            specifier, expectedRoot, actualRoot);
    end
end

dependency = struct("Name",expectedName,"Version",actualVersion,"Root",actualRoot);
fprintf("Documentation dependency: %s %s (%s)\n",dependency.Name,dependency.Version,dependency.Root);
end

function path = canonicalPath(path)
path = string(java.io.File(char(path)).getCanonicalPath());
end
