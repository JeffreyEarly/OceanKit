function capabilities = capabilities()
% Inspect native compiled-backend support without downloading or building.
%
% Expected unavailability is returned as structured data. This method does
% not download source, invoke a compiler, create cache output, warn, or throw
% because the provider is absent or unsupported.
%
% - Topic: Inspect compiled support
% - Declaration: capabilities = WVCompiledBackend.capabilities()
% - Returns capabilities: schema `1.0.0` support, identity, validation, storage, build-attempt, and failure information

capabilities = wvCompiledBackendCapabilities(struct());
end
