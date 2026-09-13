function wasLoaded = wvCompiledBackendResolveModule(moduleName,modulePath)
% Put one exact MEX module on the path without calling foreign same-name code.
arguments
    moduleName (1,1) string
    modulePath (1,1) string
end
expectedPath = canonicalPath(modulePath);
[~,mexFiles] = inmem("-completenames");
loadedPaths = string(mexFiles);
loadedPaths = loadedPaths(endsWith(loadedPaths,filesep+moduleName+"."+mexext));
wasLoaded = ~isempty(loadedPaths);
if wasLoaded && (numel(loadedPaths) ~= 1 || canonicalPath(loadedPaths(1)) ~= expectedPath)
    error("WaveVortexModel:CompiledBackendModuleConflict","A different %s MEX module is already loaded. Delete its compiled owners, clear the MEX module, and move any legacy package-local copy out of the MATLAB path.",moduleName);
end
moduleDirectory = string(fileparts(expectedPath));
if ~wasLoaded
    addpath(moduleDirectory,"-begin");
end
resolvedPath = string(which(moduleName));
if resolvedPath == "" || canonicalPath(resolvedPath) ~= expectedPath
    error("WaveVortexModel:CompiledBackendModuleConflict","MATLAB resolved %s instead of the qualified module %s. Change out of a shadowing working directory and move any legacy package-local MEX copy out of the MATLAB path.",resolvedPath,expectedPath);
end
end

function pathname = canonicalPath(pathname)
[status,output] = system("/bin/realpath "+shellQuote(pathname));
if status ~= 0
    error("WaveVortexModel:CompiledBackendPath","Unable to resolve %s.",pathname);
end
pathname = string(strtrim(output));
end

function value = shellQuote(value)
value = "'"+replace(string(value),"'","'""'""'")+"'";
end
