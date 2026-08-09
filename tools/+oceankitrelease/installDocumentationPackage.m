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

if options.repositoryRoot ~= ""
    sourceManifestPath = fullfile(options.repositoryRoot,expectedName + "-" + expectedVersion,"resources","mpackage.json");
    if ~isfile(sourceManifestPath)
        error("OceanKitRelease:DocumentationPackageSourceNotFound", ...
            "Requested package source was not found at %s.", sourceManifestPath);
    end
    sourceManifest = jsondecode(fileread(sourceManifestPath));
    if string(sourceManifest.name) ~= expectedName || string(sourceManifest.version) ~= expectedVersion
        error("OceanKitRelease:DocumentationPackageSourceMismatch", ...
            "The source manifest at %s does not match %s.", sourceManifestPath, specifier);
    end
end

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
if ~isfolder(actualRoot)
    error("OceanKitRelease:DocumentationPackageRootNotFound", ...
        "The installed package root does not exist: %s.", actualRoot);
end

dependency = struct("Name",expectedName,"Version",actualVersion,"Root",actualRoot);
fprintf("Documentation dependency: %s %s (%s)\n",dependency.Name,dependency.Version,dependency.Root);
end

function path = canonicalPath(path)
path = string(java.io.File(char(path)).getCanonicalPath());
end
