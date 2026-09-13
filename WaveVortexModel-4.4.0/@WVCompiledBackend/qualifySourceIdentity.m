function qualification = qualifySourceIdentity()
% Verify the installed binaries and current source tree against the build record.
capabilities = WVCompiledBackend.capabilities();
if ~capabilities.isAvailable
    error("WaveVortexModel:CompiledBackendSourceIdentity", ...
        "The compiled backend binary identity is not valid: %s",capabilities.failure.message);
end
recorded = capabilities.sourceIdentity;
if ~isfield(recorded,"manifest") || ~isfield(recorded.manifest,"aggregateSHA256")
    error("WaveVortexModel:CompiledBackendSourceIdentity", ...
        "The validated build record does not contain a source manifest.");
end
classFile = string(which("WVCompiledBackend"));
packageRoot = string(fileparts(fileparts(classFile)));
current = wvCompiledBackendNewSourceManifest(packageRoot);
if string(current.aggregateSHA256) ~= string(recorded.manifest.aggregateSHA256)
    error("WaveVortexModel:CompiledBackendSourceIdentityMismatch", ...
        "The current compiled-backend sources do not match the validated installed module.");
end
qualification = recorded;
qualification.currentManifest = current;
qualification.matches = true;
end
