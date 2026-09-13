function definitions = matchedModelBenchmarkCaseDefinitions(options)
% Construct the two-case matched-model RK78 contract without running it.
arguments
    options (1,1) struct
end
required = ["Nxyz" "Lxyz" "deltaT" "relativeTolerance" "absoluteTolerance" "adaptiveInitialStep" "modelConfigurations" "pilotFinalTime"];
if ~all(isfield(options,required)) || numel(options.modelConfigurations)~=1
    error("WaveVortexBenchmark:ModelConfiguration","Matched-model case construction requires one model configuration and complete integration controls.")
end
modelConfiguration = string(options.modelConfigurations);
if ~ismember(modelConfiguration,["constant-nonhydrostatic" "hydrostatic-exponential" "boussinesq-exponential"])
    error("WaveVortexBenchmark:ModelConfiguration","Unsupported model configuration %s.",modelConfiguration)
end
finalTime = 7168;
if ~isnan(options.pilotFinalTime), finalTime = options.pilotFinalTime; end
physicalConfiguration = "nonhydrostatic";
if modelConfiguration=="hydrostatic-exponential", physicalConfiguration = "hydrostatic"; end
common = struct("Nxyz",options.Nxyz,"Lxyz",options.Lxyz,"deltaT",options.deltaT,"finalTime",finalTime,"relativeTolerance",options.relativeTolerance,"absoluteTolerance",options.absoluteTolerance,"outputRelativeTolerance",options.relativeTolerance,"outputAbsoluteTolerance",options.absoluteTolerance,"initialStep",options.adaptiveInitialStep,"maximumStep",[],"maximumStepPolicy","matlab-default","forcing","default WVNonlinearAdvection","shouldAntialias",true,"seed",4001,"operation","model-continuation","integrationStepCount",NaN,"denseOutputPointsPerStep",3,"denseOutputStartTime",[],"denseOutputEndTime",[],"denseOutputRecordCount",0,"denseOutputIntegrationRecordCount",0,"modelConfiguration",modelConfiguration);
definitions = repmat(merge(common,struct("id","","requestedIntegrator","","physicalConfiguration","","isHydrostatic",false,"workload","","outputInterval",NaN,"observerGraph","","outputScheduleSeconds",[])),2,1);
workloads = ["coefficient-endpoint" "composite-dense-output"];
for iCase = 1:2
    workload = workloads(iCase);
    isDense = workload=="composite-dense-output";
    observerGraph = "coefficient-only state with endpoint-only output";
    if isDense, observerGraph = "one persisted restart record plus four first-step records containing u, two 3-D particles with tracked u, one 3-D tracer, and one source-linked mooring; coefficients at endpoints"; end
    schedule = struct("denseOutputStartTime",conditional(isDense,0,[]),"denseOutputEndTime",conditional(isDense,128,[]),"denseOutputRecordCount",conditional(isDense,5,0),"denseOutputIntegrationRecordCount",conditional(isDense,4,0),"outputScheduleSeconds",conditional(isDense,[0 32 64 96 128],[0 finalTime]));
    identity = struct("id",modelConfiguration+"--adaptive-rk78--"+workload,"requestedIntegrator","adaptive-rk78","physicalConfiguration",physicalConfiguration,"isHydrostatic",physicalConfiguration=="hydrostatic","workload",workload,"outputInterval",conditional(isDense,32,finalTime),"observerGraph",observerGraph);
    definitions(iCase) = merge(common,merge(identity,schedule));
end
end

function value = merge(first,second)
value = first;
names = fieldnames(second);
for iName = 1:numel(names), value.(names{iName}) = second.(names{iName}); end
end

function value = conditional(condition,trueValue,falseValue)
if condition, value = trueValue; else, value = falseValue; end
end
