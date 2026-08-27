function capabilities = buildForTesting(overrides)
% Exercise the build transaction with injected test state.
%
% - Topic: Test compiled support
% - Declaration: capabilities = WVCompiledBackend.buildForTesting(overrides)
% - Parameter overrides: scalar struct of test-only environment overrides
% - Returns capabilities: schema `1.0.0` capability record
% - Developer: true

arguments
    overrides (1,1) struct
end
capabilities = wvCompiledBackendBuild(overrides);
end
