function runDocumentationCheck(repositoryRoot,taskName)
arguments
    repositoryRoot (1,1) string {mustBeFolder}
    taskName (1,1) string
end

buildFile = fullfile(repositoryRoot,"buildfile.m");
if ~isfile(buildFile)
    error("OceanKitRelease:BuildFileNotFound", "Documentation checking requires %s.", buildFile);
end
before = repositoryStatus(repositoryRoot);
if before ~= ""
    error("OceanKitRelease:DirtyCheckoutBeforeDocumentationCheck", ...
        "Documentation checking requires a clean checkout:\n%s", before);
end

buildtool(taskName,"-buildFile",buildFile);

after = repositoryStatus(repositoryRoot);
if after ~= ""
    error("OceanKitRelease:DocumentationCheckModifiedCheckout", ...
        "Documentation checking modified the checkout:\n%s", after);
end
end

function status = repositoryStatus(repositoryRoot)
command = sprintf('git -C "%s" status --porcelain=v1 --untracked-files=all',repositoryRoot);
[exitCode,output] = system(command);
if exitCode ~= 0
    error("OceanKitRelease:GitStatusFailed", "Unable to inspect %s:\n%s", repositoryRoot, output);
end
status = strtrim(string(output));
end
