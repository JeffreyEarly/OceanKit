classdef WVCompiledBackend
    % Inspect and build the source-only compiled constant-stratification backend.
    %
    % `WVCompiledBackend` is a developer-facing capability and build surface.
    % It does not select a computational backend for model objects. Detection
    % is side-effect free with respect to downloads and compilation; building
    % is always an explicit action.
    %
    % ```matlab
    % capabilities = WVCompiledBackend.capabilities();
    % capabilities = WVCompiledBackend.build();
    % ```
    %
    % - Topic: Inspect compiled support
    % - Topic: Build compiled support
    % - Topic: Test compiled support
    % - Declaration: classdef WVCompiledBackend

    methods (Static)
        capabilities = capabilities()
        capabilities = build()
    end

    methods (Static, Hidden)
        capabilities = capabilitiesForTesting(overrides)
        capabilities = buildForTesting(overrides)
    end
end
