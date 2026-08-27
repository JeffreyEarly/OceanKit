function capabilities = build()
% Compile, validate, and transactionally install the pinned native provider.
%
% The build downloads FFTW 3.3.11 from the official archive, verifies its
% pinned SHA-256 digest, builds only the NEON/pthreads shared provider, stages
% and validates the MEX module, and installs it beside the package root.
% Failures are returned in the capability schema and retained as the latest
% build attempt.
%
% - Topic: Build compiled support
% - Declaration: capabilities = WVCompiledBackend.build()
% - Returns capabilities: schema `1.0.0` support, identity, validation, storage, build-attempt, and failure information

capabilities = wvCompiledBackendBuild(struct());
end
