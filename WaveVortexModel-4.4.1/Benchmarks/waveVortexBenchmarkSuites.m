function suites = waveVortexBenchmarkSuites(suiteIds)
% Return validated, versioned WaveVortex benchmark suite definitions.
arguments
    suiteIds (1,:) string = ["core-v1"]
end

registry = suiteRegistry();
unknownIds = setdiff(suiteIds,string({registry.id}));
if ~isempty(unknownIds)
    error("WaveVortexBenchmark:UnknownSuite","Unknown benchmark suite: %s.",strjoin(unknownIds,", "));
end
suites = registry(ismember(string({registry.id}),suiteIds));
for iSuite = 1:numel(suites)
    validateSuite(suites(iSuite));
end
end

function suites = suiteRegistry()
suites = [smokeSuite(),coreSuite(),standardScalingSuite(),largeScalingSuite(),transformLayoutSuite()];
end

function suite = smokeSuite()
suite = baseSuite("smoke-v1","Small correctness and harness-validation cases",false);
suite.cases = [ ...
    make3DCase("smoke-constant-nonhydrostatic","constant-nonhydrostatic",[16 16 9],false,2,1,1101), ...
    make3DCase("smoke-constant-hydrostatic","constant-hydrostatic",[16 16 9],true,2,1,1102), ...
    make3DCase("smoke-hydrostatic","hydrostatic",[16 16 9],true,2,1,1103), ...
    make3DCase("smoke-boussinesq","boussinesq",[16 16 9],false,2,1,1104), ...
    make3DCase("smoke-stratified-qg","stratified-qg",[16 16 9],true,2,1,1105), ...
    make2DCase("smoke-barotropic-qg","barotropic-qg",[32 32],2,1,1106)];
end

function suite = coreSuite()
suite = baseSuite("core-v1","Canonical constant-stratification nonlinear-advection score",true);
suite.cases = [ ...
    make3DCase("constant-nonhydrostatic-256x256x65","constant-nonhydrostatic",[256 256 65],false,7,2,2101), ...
    make3DCase("constant-hydrostatic-256x256x65","constant-hydrostatic",[256 256 65],true,7,2,2102), ...
    make3DCase("constant-nonhydrostatic-512x512x129","constant-nonhydrostatic",[512 512 129],false,3,2,2103), ...
    make3DCase("constant-hydrostatic-512x512x129","constant-hydrostatic",[512 512 129],true,3,2,2104)];
end

function suite = standardScalingSuite()
suite = baseSuite("scaling-standard-v1","Standard horizontal and vertical scaling across transform families",true);
families = ["constant-nonhydrostatic" "constant-hydrostatic" "hydrostatic" "boussinesq" "stratified-qg"];
hydrostatic = [false true true false true];
horizontalSizes = [64 64 65; 128 128 65; 256 256 65];
verticalSizes = [128 128 33; 128 128 129; 128 128 257];
cases = emptyCases();
seed = 3000;
for iFamily = 1:numel(families)
    sizes = [horizontalSizes;verticalSizes];
    for iSize = 1:size(sizes,1)
        seed = seed + 1;
        caseId = families(iFamily) + "-" + join(string(sizes(iSize,:)),"x");
        cases(end+1) = make3DCase(caseId,families(iFamily),sizes(iSize,:),hydrostatic(iFamily),7,2,seed); %#ok<AGROW>
    end
end
for nxy = [128 256 512 1024]
    seed = seed + 1;
    cases(end+1) = make2DCase("barotropic-qg-" + nxy + "x" + nxy,"barotropic-qg",[nxy nxy],7,2,seed); %#ok<AGROW>
end
suite.cases = cases;
end

function suite = largeScalingSuite()
suite = baseSuite("scaling-large-v1","Large-memory scaling across transform families",true);
families = ["constant-nonhydrostatic" "constant-hydrostatic" "hydrostatic" "boussinesq" "stratified-qg"];
hydrostatic = [false true true false true];
commonSizes = [256 256 129; 512 512 129; 512 512 257];
cases = emptyCases();
seed = 4000;
for iFamily = 1:numel(families)
    sizes = commonSizes;
    if startsWith(families(iFamily),"constant-")
        sizes(end+1,:) = [1024 1024 129]; %#ok<AGROW>
    end
    for iSize = 1:size(sizes,1)
        seed = seed + 1;
        caseId = families(iFamily) + "-" + join(string(sizes(iSize,:)),"x");
        cases(end+1) = make3DCase(caseId,families(iFamily),sizes(iSize,:),hydrostatic(iFamily),3,2,seed); %#ok<AGROW>
    end
end
for nxy = [2048 4096]
    seed = seed + 1;
    cases(end+1) = make2DCase("barotropic-qg-" + nxy + "x" + nxy,"barotropic-qg",[nxy nxy],3,2,seed); %#ok<AGROW>
end
suite.cases = cases;
end

function suite = transformLayoutSuite()
suite = baseSuite("transform-layout-v1","Full-complex horizontal mapping and Fourier diagnostics",false);
suite.kind = "transform-layout";
suite.operation = "transform-layout";
sizes = [64 48 17;65 63 17;256 256 65;512 512 129];
cases = emptyCases();
seed = 6900;
for iSize = 1:size(sizes,1)
    for shouldAntialias = [false true]
        seed = seed+1;
        antialiasId = "antialias-" + string(double(shouldAntialias));
        caseId = "full-layout-" + join(string(sizes(iSize,:)),"x") + "-" + antialiasId;
        sampleCount = 7;
        if sizes(iSize,1) == 512
            sampleCount = 3;
        end
        benchmarkCase = baseCase(caseId,"layout-only",sizes(iSize,:),false,sampleCount,2,seed);
        benchmarkCase.Lxyz = [15e3 15e3];
        benchmarkCase.operation = "transform-layout";
        benchmarkCase.shouldAntialias = shouldAntialias;
        cases(end+1) = benchmarkCase; %#ok<AGROW>
    end
end
suite.cases = cases;
end

function suite = baseSuite(id,description,isScored)
suite = struct("id",id,"version",1,"kind","model-operation","description",description,"operation","nonlinearAdvection","isScored",isScored,"selectionIsComplete",true,"cases",emptyCases());
end

function benchmarkCase = make3DCase(id,transformId,Nxyz,isHydrostatic,sampleCount,warmupCount,seed)
benchmarkCase = baseCase(id,transformId,Nxyz,isHydrostatic,sampleCount,warmupCount,seed);
benchmarkCase.Lxyz = [15e3 15e3 1300];
end

function benchmarkCase = make2DCase(id,transformId,Nxy,sampleCount,warmupCount,seed)
benchmarkCase = baseCase(id,transformId,Nxy,true,sampleCount,warmupCount,seed);
benchmarkCase.Lxyz = [15e3 15e3];
end

function benchmarkCase = baseCase(id,transformId,Nxyz,isHydrostatic,sampleCount,warmupCount,seed)
benchmarkCase = struct("id",string(id),"transformId",string(transformId),"scoreFamily",string(transformId),"operation","nonlinearAdvection","Lxyz",[],"Nxyz",double(Nxyz),"isHydrostatic",logical(isHydrostatic),"shouldAntialias",true,"seed",double(seed),"warmupCount",double(warmupCount),"sampleCount",double(sampleCount));
end

function cases = emptyCases()
cases = struct("id",{},"transformId",{},"scoreFamily",{},"operation",{},"Lxyz",{},"Nxyz",{},"isHydrostatic",{},"shouldAntialias",{},"seed",{},"warmupCount",{},"sampleCount",{});
end

function validateSuite(suite)
if strlength(suite.id) == 0 || isempty(suite.cases)
    error("WaveVortexBenchmark:InvalidSuite","Suite definitions require an ID and at least one case.");
end
caseIds = string({suite.cases.id});
if numel(unique(caseIds)) ~= numel(caseIds)
    error("WaveVortexBenchmark:DuplicateCase","Suite %s contains duplicate case IDs.",suite.id);
end
knownTransforms = ["constant-nonhydrostatic" "constant-hydrostatic" "hydrostatic" "boussinesq" "stratified-qg" "barotropic-qg" "layout-only"];
if ~all(ismember(string({suite.cases.transformId}),knownTransforms))
    error("WaveVortexBenchmark:InvalidTransform","Suite %s contains an unknown transform ID.",suite.id);
end
end
