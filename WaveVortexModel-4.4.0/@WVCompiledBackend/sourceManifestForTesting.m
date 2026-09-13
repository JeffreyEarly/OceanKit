function manifest = sourceManifestForTesting(overrides)
% Exercise deterministic source-manifest construction without building.
arguments
    overrides (1,1) struct
end
if ~isfield(overrides,"PackageRoot")
    error("WaveVortexModel:CompiledBackendTestOverride","PackageRoot is required.");
end
manifest = wvCompiledBackendNewSourceManifest(string(overrides.PackageRoot));
end
