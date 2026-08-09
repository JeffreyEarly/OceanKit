function backends = waveVortexBenchmarkBackends(backendIds)
% Return registered WaveVortex benchmark backend descriptors.
arguments
    backendIds (1,:) string = "builtin"
end
registry = struct("id","builtin","description","MATLAB builtin WaveVortex transform implementation");
unknownIds = setdiff(backendIds,string({registry.id}));
if ~isempty(unknownIds)
    error("WaveVortexBenchmark:UnknownBackend","Unknown benchmark backend: %s.",strjoin(unknownIds,", "));
end
backends = registry(ismember(string({registry.id}),backendIds));
end
