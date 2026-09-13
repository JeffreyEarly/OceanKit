classdef TestWVCompiledTransformBackend < matlab.unittest.TestCase
    methods (Test,TestTags="smoke")
        function compiledSourceIdentityIsDistinct(testCase)
            profile = @(z) 1e-4*exp(z/700);
            first = WVTransformHydrostatic([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=false);
            second = WVTransformHydrostatic([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=false);
            cleanup = onCleanup(@()deleteTransforms(first,second));
            testCase.verifyFalse(first.compiledSourceIdentity == second.compiledSourceIdentity);
            clear cleanup
        end

        function compiledSourceIdentityIsFreshAfterMatReload(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            profile = @(z) 1e-4*exp(z/700);
            original = WVTransformHydrostatic([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=false);
            cleanup = onCleanup(@()deleteTransforms(original));
            original.initWithRandomFlow(uvMax=0.01);
            original.t0 = 2;
            original.t = 7;
            originalIdentity = original.compiledSourceIdentity;
            path = fullfile(fixture.Folder,"transform.mat");
            saved = original;
            save(path,"saved","-v7");
            clear saved
            loadedData = load(path,"saved");
            loaded = loadedData.saved;
            loadedData = load(path,"saved");
            loadedAgain = loadedData.saved;
            loadedCleanup = onCleanup(@()deleteTransforms(loaded,loadedAgain));
            testCase.verifyFalse(loaded.compiledSourceIdentity == originalIdentity);
            testCase.verifyFalse(loadedAgain.compiledSourceIdentity == originalIdentity);
            testCase.verifyFalse(loaded.compiledSourceIdentity == loadedAgain.compiledSourceIdentity);
            for property = ["Lx" "Ly" "Lz" "Nx" "Ny" "Nz" "t0" "t"]
                testCase.verifyEqual(loaded.(property),original.(property));
                testCase.verifyEqual(loadedAgain.(property),original.(property));
            end
            for property = ["Ap" "Am" "A0" "kMode_wv" "lMode_wv"]
                testCase.verifyEqual(loaded.(property),original.(property));
                testCase.verifyEqual(loadedAgain.(property),original.(property));
            end
            clear loadedCleanup cleanup
        end

        function compiledSourceIdentitySurvivesClearFunctions(testCase)
            profile = @(z) 1e-4*exp(z/700);
            first = WVTransformHydrostatic([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=false);
            cleanup = onCleanup(@()deleteTransforms(first));
            firstIdentity = first.compiledSourceIdentity;
            clear functions
            second = WVTransformHydrostatic([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=false);
            secondCleanup = onCleanup(@()deleteTransforms(second));
            testCase.verifyTrue(firstIdentity == first.compiledSourceIdentity);
            testCase.verifyFalse(firstIdentity == second.compiledSourceIdentity);
            clear secondCleanup cleanup
        end

        function invalidInputsAreRejected(testCase)
            testCase.verifyError(@()WVCompiledTransformBackend.create(1),"MATLAB:validation:UnableToConvert");
        end
    end

    methods (Test,TestTags="optional")
        function allSupportedFamiliesCreateNativeSessions(testCase)
            capabilities = WVCompiledBackend.capabilities();
            testCase.assumeTrue(capabilities.isAvailable,capabilities.failure.message);
            WVCompiledBackend.activateModule(capabilities);
            definitions = supportedConfigurations();
            baseline = wv_compiled_backend_mex('moduleMetrics');
            for index = 1:numel(definitions)
                transform = definitions(index).create();
                transformCleanup = onCleanup(@()deleteTransforms(transform));
                backend = WVCompiledTransformBackend.create(transform);
                backendCleanup = onCleanup(@()deleteBackends(backend));
                metadata = backend.metadata();
                testCase.verifyEqual(string(metadata.transformClass),definitions(index).transformClass,definitions(index).name);
                testCase.verifyGreaterThan(metadata.engineBytes,0,definitions(index).name);
                testCase.verifyEqual(metadata.duplicateExecutions,0,definitions(index).name);
                expectedFFTThreads = 1;
                if ismember(definitions(index).transformClass,["WVTransformConstantStratification" "WVTransformBarotropicQG"])
                    expectedFFTThreads = capabilities.contract.threadCount;
                end
                testCase.verifyEqual(metadata.effectiveFFTThreads,expectedFFTThreads,definitions(index).name);
                clear backendCleanup transformCleanup
            end
            final = wv_compiled_backend_mex('moduleMetrics');
            testCase.verifyEqual(final.kernelCount,baseline.kernelCount);
            testCase.verifyEqual(final.matlabTransformCount,baseline.matlabTransformCount);
            testCase.verifyEqual(final.activePlans,baseline.activePlans);
        end

        function compiledAndLegacyHandlesShareLifecycle(testCase)
            capabilities = WVCompiledBackend.capabilities();
            testCase.assumeTrue(capabilities.isAvailable,capabilities.failure.message);
            WVCompiledBackend.activateModule(capabilities);
            profile = @(z) 1e-4*exp(z/700);
            hydrostatic = WVTransformHydrostatic([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=false);
            constant = WVTransformConstantStratification([4000 3000 1000],[6 6 5],isHydrostatic=true,shouldAntialias=false);
            transformCleanup = onCleanup(@()deleteTransforms(hydrostatic,constant));
            baseline = wv_compiled_backend_mex('moduleMetrics');
            legacy = WVCompiledConstantStratificationBackend.create(constant);
            modern = WVCompiledTransformBackend.create(hydrostatic);
            ownersCleanup = onCleanup(@()deleteBackends(legacy,modern));
            live = wv_compiled_backend_mex('moduleMetrics');
            testCase.verifyEqual(live.kernelCount,baseline.kernelCount+2);
            testCase.verifyEqual(live.matlabTransformCount,baseline.matlabTransformCount+1);
            clear ownersCleanup
            afterFirst = wv_compiled_backend_mex('moduleMetrics');
            testCase.verifyEqual(afterFirst.kernelCount,baseline.kernelCount);
            testCase.verifyEqual(afterFirst.matlabTransformCount,baseline.matlabTransformCount);
            modern = WVCompiledTransformBackend.create(hydrostatic);
            legacy = WVCompiledConstantStratificationBackend.create(constant);
            ownersCleanup = onCleanup(@()deleteBackends(legacy,modern));
            live = wv_compiled_backend_mex('moduleMetrics');
            testCase.verifyEqual(live.kernelCount,baseline.kernelCount+2);
            testCase.verifyEqual(live.matlabTransformCount,baseline.matlabTransformCount+1);
            clear ownersCleanup transformCleanup
            final = wv_compiled_backend_mex('moduleMetrics');
            testCase.verifyEqual(final.kernelCount,baseline.kernelCount);
            testCase.verifyEqual(final.matlabTransformCount,baseline.matlabTransformCount);
            testCase.verifyEqual(final.activePlans,baseline.activePlans);
        end

        function hydrostaticCallScopeReusesFieldsAndRecovers(testCase)
            capabilities = WVCompiledBackend.capabilities();
            testCase.assumeTrue(capabilities.isAvailable,capabilities.failure.message);
            WVCompiledBackend.activateModule(capabilities);
            profile = @(z) 1e-4*exp(z/700);
            wvt = WVTransformHydrostatic([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=false);
            other = WVTransformHydrostatic([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=false);
            cleanup = onCleanup(@()deleteTransforms(wvt,other));
            rng(503,"twister");
            wvt.initWithRandomFlow(uvMax=0.01);
            wvt.t0 = 2; wvt.t = 7;
            moduleBaseline = wv_compiled_backend_mex('moduleMetrics');
            backend = WVCompiledTransformBackend.create(wvt);
            backendCleanup = onCleanup(@()delete(backend));
            moduleLive = wv_compiled_backend_mex('moduleMetrics');
            testCase.verifyEqual(moduleLive.matlabTransformCount,moduleBaseline.matlabTransformCount+1);
            testCase.verifyEqual(moduleLive.kernelCount,moduleBaseline.kernelCount+1);
            testCase.verifyTrue(logical(moduleLive.moduleLocked));
            names = ["u","v","w","eta"];
            backend.prepare(names);
            result = backend.evaluate(wvt,names,shouldEvaluateFlux=true);
            nameCells = cellstr(names);
            expected = cell(1,numel(names)); [expected{:}] = wvt.variableWithName(nameCells{:});
            for index = 1:numel(names)
                testCase.verifyLessThanOrEqual(relativeError(result.values{index},expected{index}),1e-12);
            end
            [fp,fm,f0] = wvt.nonlinearFlux();
            testCase.verifyLessThanOrEqual(relativeError(result.flux{1},fp),1e-12);
            testCase.verifyLessThanOrEqual(relativeError(result.flux{2},fm),1e-12);
            testCase.verifyLessThanOrEqual(relativeError(result.flux{3},f0),1e-12);
            firstMetrics = result.metrics;
            repeatedMixed = backend.evaluate(wvt,names,shouldEvaluateFlux=true);
            testCase.verifyEqual(repeatedMixed.values,result.values,"Repeated field batch changed.");
            testCase.verifyEqual(repeatedMixed.flux,result.flux,"Repeated flux batch changed.");
            testCase.verifyEqual(repeatedMixed.metrics.stateValidations-firstMetrics.stateValidations,1);
            testCase.verifyEqual(repeatedMixed.metrics.phasePreparations-firstMetrics.phasePreparations,1);
            testCase.verifyEqual(repeatedMixed.metrics.duplicateExecutions,0);
            testCase.verifyEqual(repeatedMixed.metrics.liveEvaluationBytes,0);
            wvt.t = 11;
            changed = backend.evaluate(wvt,names);
            expectedChanged = cell(1,numel(names)); [expectedChanged{:}] = wvt.variableWithName(nameCells{:});
            for index = 1:numel(names)
                testCase.verifyLessThanOrEqual(relativeError(changed.values{index},expectedChanged{index}),1e-12);
            end
            testCase.verifyEqual(changed.metrics.stateValidations-repeatedMixed.metrics.stateValidations,1);
            testCase.verifyEqual(changed.metrics.phasePreparations-repeatedMixed.metrics.phasePreparations,1);
            testCase.verifyEqual(changed.metrics.duplicateExecutions,0);
            testCase.verifyEqual(changed.metrics.liveEvaluationBytes,0);
            testCase.verifyError(@()backend.evaluate(other,names),"WaveVortexModel:CompiledTransformMismatch");
            testCase.verifyError(@()backend.prepare(42),"WaveVortexModel:CompiledTransformVariables");
            testCase.verifyError(@()backend.prepare("unsupported-variable"),"WaveVortexModel:CompiledTransformExecution");
            backend.prepare(names);
            recovered = backend.evaluate(wvt,names);
            for index = 1:numel(names)
                testCase.verifyLessThanOrEqual(relativeError(recovered.values{index},expectedChanged{index}),1e-12);
            end
            clear backendCleanup
            moduleAfterDelete = wv_compiled_backend_mex('moduleMetrics');
            testCase.verifyEqual(moduleAfterDelete.matlabTransformCount,moduleBaseline.matlabTransformCount);
            testCase.verifyEqual(moduleAfterDelete.kernelCount,moduleBaseline.kernelCount);
            testCase.verifyEqual(moduleAfterDelete.activePlans,moduleBaseline.activePlans);
            testCase.verifyEqual(moduleAfterDelete.outstandingPlanningBytes,0);
            clear cleanup
        end
    end
end

function definitions = supportedConfigurations
profile = @(z)1e-4*exp(z/700);
definitions = [ ...
    struct("name","constant-hydrostatic","transformClass","WVTransformConstantStratification","create",@()WVTransformConstantStratification([4000 3000 1000],[6 6 5],isHydrostatic=true,shouldAntialias=false)) ...
    struct("name","constant-nonhydrostatic","transformClass","WVTransformConstantStratification","create",@()WVTransformConstantStratification([4000 3000 1000],[6 6 5],isHydrostatic=false,shouldAntialias=false)) ...
    struct("name","hydrostatic","transformClass","WVTransformHydrostatic","create",@()WVTransformHydrostatic([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=false)) ...
    struct("name","boussinesq","transformClass","WVTransformBoussinesq","create",@()WVTransformBoussinesq([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=false)) ...
    struct("name","stratified-qg","transformClass","WVTransformStratifiedQG","create",@()WVTransformStratifiedQG([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=false)) ...
    struct("name","barotropic-qg","transformClass","WVTransformBarotropicQG","create",@()WVTransformBarotropicQG([4000 3000],[6 6],h=1000,j=1,shouldAntialias=false))];
end

function deleteTransforms(varargin)
for transform = varargin
    if isa(transform{1},"handle") && isvalid(transform{1})
        delete(transform{1});
    end
end
end

function deleteBackends(varargin)
for backend = varargin
    if ~isempty(backend{1}) && isvalid(backend{1})
        delete(backend{1});
    end
end
end

function value = relativeError(actual,expected)
value = max(abs(actual(:)-expected(:)),[],'omitmissing')/max(max(abs(expected(:)),[],'omitmissing'),realmin);
end
