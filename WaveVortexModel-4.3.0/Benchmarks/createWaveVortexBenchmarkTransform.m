function wvt = createWaveVortexBenchmarkTransform(benchmarkCase,backendId)
% Construct one registered transform/backend benchmark candidate.
arguments
    benchmarkCase (1,1) struct
    backendId (1,1) string = "builtin"
end
if backendId == "compiled" && ~startsWith(string(benchmarkCase.transformId),"constant-")
    error("WaveVortexBenchmark:UnsupportedBackend","The compiled preview supports only constant-stratification benchmark cases.");
end
computationalBackend = conditional(backendId=="compiled","compiled","matlab");

switch benchmarkCase.transformId
    case "constant-nonhydrostatic"
        wvt = WVTransformConstantStratification(benchmarkCase.Lxyz,benchmarkCase.Nxyz,isHydrostatic=false,shouldAntialias=benchmarkCase.shouldAntialias,computationalBackend=computationalBackend);
    case "constant-hydrostatic"
        wvt = WVTransformConstantStratification(benchmarkCase.Lxyz,benchmarkCase.Nxyz,isHydrostatic=true,shouldAntialias=benchmarkCase.shouldAntialias,computationalBackend=computationalBackend);
    case "hydrostatic"
        wvt = WVTransformHydrostatic(benchmarkCase.Lxyz,benchmarkCase.Nxyz,N2=@benchmarkN2,shouldAntialias=benchmarkCase.shouldAntialias);
    case "boussinesq"
        wvt = WVTransformBoussinesq(benchmarkCase.Lxyz,benchmarkCase.Nxyz,N2=@benchmarkN2,shouldAntialias=benchmarkCase.shouldAntialias);
    case "stratified-qg"
        wvt = WVTransformStratifiedQG(benchmarkCase.Lxyz,benchmarkCase.Nxyz,N2=@benchmarkN2,shouldAntialias=benchmarkCase.shouldAntialias);
    case "barotropic-qg"
        wvt = WVTransformBarotropicQG(benchmarkCase.Lxyz,benchmarkCase.Nxyz,shouldAntialias=benchmarkCase.shouldAntialias);
    otherwise
        error("WaveVortexBenchmark:UnknownTransform","Unknown transform ID %s.",benchmarkCase.transformId);
end

function value = conditional(condition,trueValue,falseValue)
if condition
    value = trueValue;
else
    value = falseValue;
end
end
end

function N2 = benchmarkN2(z)
N0 = 5.2e-3;
N2 = N0*N0*ones(size(z));
end
