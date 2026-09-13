function waveVortexObserverCostBenchmarkWorker(configPath,outputPath)
% Fresh-process entry point used by runWaveVortexObserverCostBenchmark.
config = jsondecode(fileread(configPath));
try
    value = config.configuration;
    result = runWaveVortexObserverCostBenchmark(caseIds=string(config.caseId), ...
        Nxyz=value.Nxyz,Lxyz=value.Lxyz,deltaT=value.deltaT, ...
        integrationStepCount=value.integrationStepCount,denseOutputPointsPerStep=value.denseOutputPointsPerStep, ...
        seed=value.seed,shouldUseFreshProcess=false,phasePath=string(config.phasePath),plateauSeconds=config.plateauSeconds);
    payload = struct("status","complete","case",result.cases(1),"failure","");
catch exception
    payload = struct("status","failed","case",struct,"failure",string(getReport(exception,"extended","hyperlinks","off")));
end
fileId = fopen(outputPath,"w");
cleanup = onCleanup(@()fclose(fileId));
fprintf(fileId,"%s",jsonencode(payload));
clear cleanup
end
