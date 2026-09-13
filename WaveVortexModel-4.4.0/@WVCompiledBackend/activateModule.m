function moduleName = activateModule(capabilities)
% Resolve the exact qualified MEX module for a compiled transform owner.
arguments
    capabilities (1,1) struct
end
if ~isfield(capabilities,"isAvailable") || ~capabilities.isAvailable || ...
        ~isfield(capabilities,"module") || ~isfield(capabilities.module,"identityValidated") || ...
        ~capabilities.module.identityValidated
    error("WaveVortexModel:CompiledBackendUnavailable","A qualified compiled backend is required before activating its MEX module.");
end
module = capabilities.module;
if ~isfield(module,"name") || ~isfield(module,"path") || ~isfield(module,"sha256") || ...
        string(module.name) == "" || string(module.path) == "" || string(module.sha256) == ""
    error("WaveVortexModel:CompiledBackendCapabilitySchema","The compiled-backend capability record has no complete module identity.");
end
moduleName = string(module.name);
modulePath = string(module.path);
if ~isfile(modulePath)
    error("WaveVortexModel:CompiledBackendNotBuilt","The qualified MEX module is missing: %s",modulePath);
end
wvCompiledBackendResolveModule(moduleName,modulePath);
if sha256File(modulePath) ~= string(module.sha256)
    error("WaveVortexModel:CompiledBackendIdentityMismatch","The resolved MEX module does not match the qualified binary identity.");
end
end

function hash = sha256File(pathname)
[status,output] = system("/usr/bin/shasum -a 256 "+shellQuote(pathname));
if status ~= 0
    error("WaveVortexModel:CompiledBackendHash","Unable to hash %s.",pathname);
end
hash = extractBefore(string(strtrim(output))," ");
end

function value = shellQuote(value)
value = "'"+replace(string(value),"'","'""'""'")+"'";
end
