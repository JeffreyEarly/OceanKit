function modulePath = buildConstantAdoptionMex(sourceRoot,providerRoot,outputFolder,moduleName,options)
% Build an isolated author-only MEX with the production compiler/provider flags.
arguments
    sourceRoot (1,1) string
    providerRoot (1,1) string
    outputFolder (1,1) string
    moduleName (1,1) string
    options.compact (1,1) logical = false
    options.horizontalWorkers (1,1) double {mustBeInteger,mustBePositive} = 12
    options.pointwiseWorkers (1,1) double {mustBeInteger,mustBePositive} = 12
end
assert(~isfolder(outputFolder),"Refusing to overwrite an isolated module build.");
mkdir(outputFolder);
adapter = fullfile(sourceRoot,"CompiledKernel","adapters","native-fftw");
core = fullfile(sourceRoot,"CompiledKernel","src",["WVKernelTypes.cpp","WVPreparedVerticalOperator.cpp","WVTransformStratifiedQGKernel.cpp","WVTransformHydrostaticKernel.cpp","WVTransformBoussinesqKernel.cpp","WVRetainedHorizontalOperator.cpp","WVTransformBarotropicQGKernel.cpp","WVTransformConstantStratificationKernel.cpp"]);
runtime = fullfile(sourceRoot,"PortableRuntime","src",["WVExtensionCatalog.cpp","WVPortableImplementationContract.cpp","WVForcingContracts.cpp","WVForcingEngine.cpp","WVBarotropicQGForcingEngine.cpp","WVBarotropicQGIntegrationSystem.cpp","WVIntegrationState.cpp","WVRungeKutta.cpp","WVObserverAdapter.cpp","WVObserverContracts.cpp","WVObservation.cpp","WVPortableTypedRecord.cpp","WVOutputSchedule.cpp","WVFieldEvaluationService.cpp","WVBarotropicQGFieldEvaluationAdapter.cpp","WVConstantStratificationIntegrationSystem.cpp","WVModelTransformAdapters.cpp","WVModel.cpp"]);
runtime = [runtime,fullfile(sourceRoot,"PortableRuntime","src",["WVStratifiedQGForcingEngine.cpp","WVHydrostaticForcingEngine.cpp","WVBoussinesqForcingEngine.cpp","WVStratifiedQGIntegrationSystem.cpp","WVHydrostaticIntegrationSystem.cpp","WVBoussinesqIntegrationSystem.cpp","WVStratifiedFieldEvaluationAdapter.cpp","WVDiagnosticFieldPlan.cpp","WVDensityEventEvaluation.cpp","WVNoMotionProfile.cpp","WVNoMotionProfileRecovery.cpp","WVNoMotionDensityMoments.cpp"])];
sources = [fullfile(adapter,"wv_compiled_backend_mex.cpp"),core,runtime,fullfile(adapter,"WVNativeFFTWEngine.cpp")];
compiler = "CXXFLAGS=$CXXFLAGS -std=c++17 -pthread -O3 -mcpu=native -mmacosx-version-min=13.3 -DWV_KERNEL_NATIVE_OPTIMIZATION=1 -DWV_KERNEL_COEFFICIENT_WORKERS=2 -DWV_MODEL_ENABLE_OUTPUT=0";
compiler = compiler+" -DWV_KERNEL_COMPACT_CONSTANT_CANDIDATE="+double(options.compact)+" -DWV_KERNEL_COMPACT_HORIZONTAL_WORKERS="+options.horizontalWorkers+" -DWV_KERNEL_COMPACT_POINTWISE_WORKERS="+options.pointwiseWorkers;
linker = "LDFLAGS=$LDFLAGS -pthread -mmacosx-version-min=13.3 -Wl,-rpath,"+fullfile(providerRoot,"lib");
sourceArguments = cellstr(sources);
mex("-R2018a",compiler,sourceArguments{:},"-I"+fullfile(sourceRoot,"CompiledKernel","include"),"-I"+fullfile(sourceRoot,"PortableRuntime","include"),"-I"+fullfile(sourceRoot,"PortableRuntime","src"),"-I"+adapter,"-I"+fullfile(providerRoot,"include"),linker,fullfile(providerRoot,"lib","libfftw3_threads.dylib"),fullfile(providerRoot,"lib","libfftw3.dylib"),"-outdir",outputFolder,"-output",moduleName);
modulePath = fullfile(outputFolder,moduleName+"."+mexext);
record = struct(modulePath=modulePath,moduleSHA256=portableCompatibilitySHA256(modulePath),sourceRoot=sourceRoot,providerRoot=providerRoot,compilerFlags=compiler,linkerFlags=linker,sources=sources,options=options,matlabVersion=string(version));
file = fopen(fullfile(outputFolder,"build.json"),"w"); assert(file>=0); cleanup = onCleanup(@()fclose(file));
fwrite(file,jsonencode(record,PrettyPrint=true));
end
