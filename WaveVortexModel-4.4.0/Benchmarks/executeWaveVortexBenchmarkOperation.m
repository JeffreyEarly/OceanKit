function outputs = executeWaveVortexBenchmarkOperation(wvt,operationId)
% Execute one registered benchmark operation and return comparable outputs.
arguments
    wvt WVTransform
    operationId (1,1) string
end
switch operationId
    case "nonlinearAdvection"
        if isa(wvt,"WVTransformBarotropicQG") || isa(wvt,"WVTransformStratifiedQG")
            outputs = {wvt.nonlinearFlux()};
        else
            [Fp,Fm,F0] = wvt.nonlinearFlux();
            outputs = {Fp,Fm,F0};
        end
    otherwise
        error("WaveVortexBenchmark:UnknownOperation","Unknown benchmark operation %s.",operationId);
end
end
