function capabilities = capabilitiesForTesting(overrides)
% Inspect capabilities with injected platform and filesystem state.
%
% - Topic: Test compiled support
% - Declaration: capabilities = WVCompiledBackend.capabilitiesForTesting(overrides)
% - Parameter overrides: scalar struct of test-only environment overrides
% - Returns capabilities: schema `1.0.0` capability record
% - Developer: true

arguments
    overrides (1,1) struct
end
capabilities = wvCompiledBackendCapabilities(overrides);
end
