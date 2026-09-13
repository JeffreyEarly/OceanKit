classdef TestWVCompiledEvaluationScope < matlab.unittest.TestCase
    % Event-scope lifecycle and cache contracts for compiled transforms.
    methods (Test,TestTags="optional")
        function allFamiliesShareOneScopedStateAndCache(testCase)
            for definition = configurations()
                wvt = definition.create(); transformCleanup = onCleanup(@()delete(wvt));
                rng(506,"twister"); wvt.initWithRandomFlow(uvMax=0.01); wvt.t0 = 5; wvt.t = 13;
                backend = WVCompiledTransformBackend.create(wvt); backendCleanup = onCleanup(@()delete(backend));
                backend.prepare({'u','v'});

                baseline = backend.metadata();
                outer = backend.scopedEvaluation(wvt);
                started = backend.metadata();
                testCase.verifyEqual(started.scopedEvaluations-baseline.scopedEvaluations,1,definition.name);
                testCase.verifyEqual(started.stateInputCopyBytes-baseline.stateInputCopyBytes,started.stateStorageBytes,definition.name);
                testCase.verifyEqual(started.evaluationActive,1,definition.name);

                firstU = backend.evaluateVariables(wvt,{'u'}); afterU = backend.metadata();
                firstV = backend.evaluateVariables(wvt,{'v'}); afterV = backend.metadata();
                repeatedU = backend.evaluateVariables(wvt,{'u'}); afterRepeatedU = backend.metadata();
                verifyNear(testCase,firstU.values{1},wvt.u,definition.name+" u");
                verifyNear(testCase,firstV.values{1},wvt.v,definition.name+" v");
                testCase.verifyEqual(repeatedU.values{1},firstU.values{1},definition.name);
                testCase.verifyGreaterThan(afterU.producerExecutions,started.producerExecutions,definition.name);
                testCase.verifyGreaterThanOrEqual(afterV.producerExecutions,afterU.producerExecutions,definition.name);
                testCase.verifyEqual(afterRepeatedU.producerExecutions,afterV.producerExecutions,definition.name);
                testCase.verifyGreaterThan(afterRepeatedU.cacheHits,afterV.cacheHits,definition.name);
                testCase.verifyEqual(afterRepeatedU.duplicateExecutions-baseline.duplicateExecutions,0,definition.name);

                forcing = backend.evaluateVariables(wvt,forcingNames(wvt)); afterForcing = backend.metadata();
                testCase.verifyEqual(numel(forcing.values),numel(forcingNames(wvt)),definition.name);
                for value = forcing.values'
                    testCase.verifyTrue(all(isfinite(value{1}(:))),definition.name);
                end
                testCase.verifyGreaterThan(afterForcing.producerExecutions,afterRepeatedU.producerExecutions,definition.name);
                clear outer
                closed = backend.metadata();
                testCase.verifyEqual(closed.evaluationActive,0,definition.name);
                testCase.verifyEqual(closed.liveEvaluationBytes,0,definition.name);
                testCase.verifyEqual(closed.stateValidations-baseline.stateValidations,1,definition.name);
                if isQG(wvt)
                    testCase.verifyEqual(closed.phasePreparations-baseline.phasePreparations,0,definition.name);
                else
                    testCase.verifyEqual(closed.phasePreparations-baseline.phasePreparations,1,definition.name);
                end

                nestedBaseline = backend.metadata();
                outer = backend.scopedEvaluation(wvt); nested = backend.scopedEvaluation(wvt);
                testCase.verifyEqual(backend.metadata().scopedEvaluations-nestedBaseline.scopedEvaluations,1,definition.name);
                clear nested
                testCase.verifyEqual(backend.metadata().evaluationActive,1,definition.name);
                backend.evaluateVariables(wvt,{'u'});
                clear outer
                testCase.verifyEqual(backend.metadata().evaluationActive,0,definition.name);

                mutationScope = backend.scopedEvaluation(wvt);
                mutateCoefficient(wvt);
                testCase.verifyError(@()backend.evaluateVariables(wvt,{'u'}),"WaveVortexModel:CompiledTransformStateChanged",definition.name);
                testCase.verifyEqual(backend.metadata().evaluationActive,0,definition.name);
                clear mutationScope
                fresh = backend.scopedEvaluation(wvt); backend.evaluateVariables(wvt,{'u'}); clear fresh

                timeScope = backend.scopedEvaluation(wvt);
                wvt.t0 = wvt.t0+1;
                testCase.verifyError(@()backend.evaluateVariables(wvt,{'u'}),"WaveVortexModel:CompiledTransformStateChanged",definition.name);
                testCase.verifyEqual(backend.metadata().evaluationActive,0,definition.name);
                clear timeScope
                fresh = backend.scopedEvaluation(wvt); backend.evaluateVariables(wvt,{'u'}); clear fresh

                failureScope = backend.scopedEvaluation(wvt);
                cachedU = backend.evaluateVariables(wvt,{'u'}); beforeFailure = backend.metadata();
                testCase.verifyError(@()backend.evaluateVariables(wvt,{'not_a_portable_variable'}),"WaveVortexModel:CompiledTransformExecution",definition.name);
                afterFailure = backend.metadata();
                testCase.verifyEqual(afterFailure.producerExecutions,beforeFailure.producerExecutions,definition.name);
                testCase.verifyEqual(afterFailure.evaluationActive,1,definition.name);
                backend.evaluateVariables(wvt,{'v'}); afterSuccessfulV = backend.metadata();
                cachedAgain = backend.evaluateVariables(wvt,{'u'}); afterCachedAgain = backend.metadata();
                testCase.verifyEqual(cachedAgain.values{1},cachedU.values{1},definition.name);
                testCase.verifyEqual(afterCachedAgain.producerExecutions,afterSuccessfulV.producerExecutions,definition.name);
                testCase.verifyGreaterThan(afterCachedAgain.cacheHits,afterSuccessfulV.cacheHits,definition.name);
                clear failureScope
                testCase.verifyEqual(backend.metadata().liveEvaluationBytes,0,definition.name);

                pendingCleanup = backend.scopedEvaluation(wvt);
                delete(backend);
                testCase.verifyWarningFree(@()delete(pendingCleanup),definition.name);
                clear pendingCleanup backendCleanup transformCleanup
            end
        end
    end
end

function definitions = configurations()
profile = @(z)1e-4*exp(z/700);
definitions = [ ...
    struct("name","constant-hydrostatic","create",@()WVTransformConstantStratification([4000 3000 1000],[6 6 5],N0=5.2e-3,isHydrostatic=true,shouldAntialias=false)) ...
    struct("name","constant-nonhydrostatic","create",@()WVTransformConstantStratification([4000 3000 1000],[6 6 5],N0=5.2e-3,isHydrostatic=false,shouldAntialias=false)) ...
    struct("name","hydrostatic","create",@()WVTransformHydrostatic([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=false)) ...
    struct("name","boussinesq","create",@()WVTransformBoussinesq([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=false)) ...
    struct("name","stratified-qg","create",@()WVTransformStratifiedQG([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=false)) ...
    struct("name","barotropic-qg","create",@()WVTransformBarotropicQG([4000 3000],[6 6],h=1000,j=1,shouldAntialias=false))];
end

function names = forcingNames(wvt)
if isQG(wvt)
    names = {'Fqgpv_nonlinear_advection'};
else
    names = {'Fu_nonlinear_advection','Fv_nonlinear_advection','Feta_nonlinear_advection'};
end
end

function value = isQG(wvt)
value = isa(wvt,"WVTransformStratifiedQG") || isa(wvt,"WVTransformBarotropicQG");
end

function mutateCoefficient(wvt)
column = find(wvt.kMode_wv == 1 & wvt.lMode_wv == 1,1);
row = min(2,wvt.Nj);
index = sub2ind(size(wvt.A0),row,column);
if isQG(wvt)
    value = wvt.A0; value(index) = value(index)+complex(1e-7,-2e-7); wvt.A0 = value;
else
    value = wvt.Ap; value(index) = value(index)+complex(1e-7,-2e-7); wvt.Ap = value;
end
end

function verifyNear(testCase,actual,expected,diagnostic)
testCase.verifySize(actual,size(expected),diagnostic);
testCase.assertTrue(all(isfinite(actual(:))) && all(isfinite(expected(:))),diagnostic);
error = norm(actual(:)-expected(:))/max(norm(expected(:)),realmin);
testCase.verifyLessThanOrEqual(error,1e-12,diagnostic);
end
