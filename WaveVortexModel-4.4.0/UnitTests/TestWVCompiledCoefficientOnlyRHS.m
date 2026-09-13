classdef TestWVCompiledCoefficientOnlyRHS < matlab.unittest.TestCase
    % Qualify the sealed native coefficient-only WVModel right-hand side.
    methods (Test,TestTags="optional")
        function allSixFamiliesUseCoefficientOnlyLeaf(testCase)
            for definition = configurations()
                matlabWVT = definition.create("matlab");
                compiledWVT = definition.create("compiled");
                cleanup = onCleanup(@()deleteTransforms(matlabWVT,compiledWVT));
                seedState(matlabWVT,compiledWVT);
                matlabModel = WVModel(matlabWVT);
                compiledModel = WVModel(compiledWVT);
                testCase.verifyTrue(compiledWVT.canUseCompiledCoefficientOnlyRightHandSide(),definition.name);
                matlabState = matlabModel.initialConditionsCellArray();
                compiledState = compiledModel.initialConditionsCellArray();
                before = compiledWVT.computationalBackendMetadata.runtimeMetrics;
                expected = matlabModel.fluxAtTimeCellArray(13,matlabState);
                actual = compiledModel.fluxAtTimeCellArray(13,compiledState);
                verifyCellNear(testCase,actual,expected,definition.name+" coefficient-only RHS");
                after = compiledWVT.computationalBackendMetadata.runtimeMetrics;
                testCase.verifyEqual(after.coefficientOnlyEvaluations-before.coefficientOnlyEvaluations,1,definition.name);
                testCase.verifyEqual(after.stateInputCopyBytes-before.stateInputCopyBytes,0,definition.name);
                testCase.verifyEqual(after.planPreparations-before.planPreparations,0,definition.name);
                testCase.verifyEqual(after.plannedEvaluationBytes,before.plannedEvaluationBytes,definition.name);
                testCase.verifyEqual(after.scopedEvaluations-before.scopedEvaluations,0,definition.name);
                testCase.verifyEqual(after.cacheHits-before.cacheHits,0,definition.name);
                testCase.verifyEqual(after.duplicateExecutions-before.duplicateExecutions,0,definition.name);
                if ~isQG(compiledWVT)
                    testCase.verifyEqual(after.stateValidations-before.stateValidations,1,definition.name);
                    testCase.verifyEqual(after.phasePreparations-before.phasePreparations,1,definition.name);
                end
                ordinary = compiledWVT.compiledVariables("u");
                testCase.verifyTrue(all(isfinite(ordinary{1}(:))),definition.name+" later generic query");
                clear cleanup
            end
        end

        function nonsealedWorkloadsKeepEstablishedModelPath(testCase)
            WVCompiledBackend.activateModule(WVCompiledBackend.capabilities());
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fileparts(mfilename('fullpath'))));
            scenarios = ["damping" "tracer" "particle" "extra-observer" "coefficient-subclass" "density-correction"];
            definition = configurations();
            definition = definition(3);
            for scenario = scenarios
                matlabWVT = definition.create("matlab");
                compiledWVT = definition.create("compiled");
                cleanup = onCleanup(@()deleteTransforms(matlabWVT,compiledWVT));
                seedState(matlabWVT,compiledWVT);
                if scenario == "damping"
                    matlabWVT.setForcing([WVNonlinearAdvection(matlabWVT),WVAdaptiveDamping(matlabWVT)]);
                    compiledWVT.setForcing([WVNonlinearAdvection(compiledWVT),WVAdaptiveDamping(compiledWVT)]);
                elseif scenario == "density-correction"
                    matlabForcing = matlabWVT.forcingWithName("nonlinear advection");
                    compiledForcing = compiledWVT.forcingWithName("nonlinear advection");
                    matlabForcing.dLnN2 = matlabForcing.dLnN2 + 1e-3;
                    compiledForcing.dLnN2 = compiledForcing.dLnN2 + 1e-3;
                end
                matlabModel = WVModel(matlabWVT);
                compiledModel = WVModel(compiledWVT);
                switch scenario
                    case "tracer"
                        phi = reshape(1:prod(matlabWVT.spatialMatrixSize),matlabWVT.spatialMatrixSize);
                        matlabModel.addTracer(phi,"dye");
                        compiledModel.addTracer(phi,"dye");
                    case "particle"
                        matlabModel.setFloatPositions(500,400,-250,'u');
                        compiledModel.setFloatPositions(500,400,-250,'u');
                    case "extra-observer"
                        matlabModel.addFluxedObservingSystem(WVTestPortablePointDiagnostic(matlabModel,x=500,y=400,z=-250,fieldName="u"));
                        compiledModel.addFluxedObservingSystem(WVTestPortablePointDiagnostic(compiledModel,x=500,y=400,z=-250,fieldName="u"));
                    case "coefficient-subclass"
                        matlabModel.removeFluxedObservingSystem(matlabModel.wvCoefficientFluxedObservingSystem());
                        compiledModel.removeFluxedObservingSystem(compiledModel.wvCoefficientFluxedObservingSystem());
                        matlabModel.addFluxedCoefficients(WVTestOverriddenCoefficients(matlabModel));
                        compiledModel.addFluxedCoefficients(WVTestOverriddenCoefficients(compiledModel));
                end
                before = compiledWVT.computationalBackendMetadata.runtimeMetrics;
                if scenario == "damping" || scenario == "density-correction"
                    testCase.verifyFalse(compiledWVT.canUseCompiledCoefficientOnlyRightHandSide(),scenario);
                end
                expected = matlabModel.fluxAtTimeCellArray(13,matlabModel.initialConditionsCellArray());
                actual = compiledModel.fluxAtTimeCellArray(13,compiledModel.initialConditionsCellArray());
                verifyCellNear(testCase,actual,expected,scenario+" established RHS");
                after = compiledWVT.computationalBackendMetadata.runtimeMetrics;
                testCase.verifyEqual(after.coefficientOnlyEvaluations,before.coefficientOnlyEvaluations,scenario);
                clear cleanup
            end
        end

        function explicitEvaluationLeasesRejectFreshLeaf(testCase)
            definition = configurations();
            wvt = definition(3).create("compiled");
            cleanup = onCleanup(@()delete(wvt));
            rng(511,"twister");
            wvt.initWithRandomFlow(uvMax=.01);
            backend = WVCompiledTransformBackend.create(wvt);
            backendCleanup = onCleanup(@()delete(backend));
            scope = backend.scopedEvaluation(wvt);
            testCase.verifyFalse(backend.canBeginCoefficientOnlyEvaluation());
            [used,result] = backend.tryCoefficientOnlyRightHandSide(wvt);
            testCase.verifyFalse(used);
            testCase.verifyEmpty(result);
            testCase.verifyError(@()backend.coefficientOnlyRightHandSide(wvt),"WaveVortexModel:CompiledCoefficientOnlyWorkload");
            clear scope
            scope = wvt.scopedEvaluation();
            [used,result] = backend.tryCoefficientOnlyRightHandSide(wvt);
            testCase.verifyFalse(used);
            testCase.verifyEmpty(result);
            clear scope
            scope = backend.scopedEvaluation(wvt);
            wvt.t = wvt.t + 1;
            testCase.verifyError(@()backend.tryCoefficientOnlyRightHandSide(wvt),"WaveVortexModel:CompiledTransformStateChanged");
            testCase.verifyError(@()backend.coefficientOnlyRightHandSide(wvt),"WaveVortexModel:CompiledTransformStateChanged");
            testCase.verifyError(@()backend.coefficientOnlyRightHandSide(wvt),"WaveVortexModel:CompiledTransformStateChanged");
            clear scope backendCleanup cleanup
        end
    end
end

function definitions = configurations()
profile = @(z)1e-4*exp(z/700);
definitions = [ ...
    struct("name","constant-hydrostatic","create",@(backend)WVTransformConstantStratification([4000 3000 1000],[6 6 5],N0=5.2e-3,isHydrostatic=true,shouldAntialias=true,computationalBackend=backend)), ...
    struct("name","constant-nonhydrostatic","create",@(backend)WVTransformConstantStratification([4000 3000 1000],[6 6 5],N0=5.2e-3,isHydrostatic=false,shouldAntialias=true,computationalBackend=backend)), ...
    struct("name","hydrostatic","create",@(backend)WVTransformHydrostatic([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=true,computationalBackend=backend)), ...
    struct("name","boussinesq","create",@(backend)WVTransformBoussinesq([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=true,computationalBackend=backend)), ...
    struct("name","stratified-qg","create",@(backend)WVTransformStratifiedQG([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=true,computationalBackend=backend)), ...
    struct("name","barotropic-qg","create",@(backend)WVTransformBarotropicQG([4000 3000],[6 6],h=1000,j=1,shouldAntialias=true,computationalBackend=backend))];
end

function seedState(matlabWVT,compiledWVT)
rng(511,"twister");
matlabWVT.initWithRandomFlow(uvMax=.01);
compiledWVT.A0 = matlabWVT.A0;
if ~isQG(matlabWVT)
    compiledWVT.Ap = matlabWVT.Ap;
    compiledWVT.Am = matlabWVT.Am;
end
matlabWVT.t0 = 5; matlabWVT.t = 13;
compiledWVT.t0 = 5; compiledWVT.t = 13;
end

function value = isQG(wvt)
value = isa(wvt,"WVTransformStratifiedQG") || isa(wvt,"WVTransformBarotropicQG");
end

function verifyCellNear(testCase,actual,expected,diagnostic)
testCase.verifyEqual(numel(actual),numel(expected),diagnostic);
for index = 1:numel(expected)
    testCase.verifySize(actual{index},size(expected{index}),diagnostic);
    testCase.verifyTrue(all(isfinite(actual{index}(:))) && all(isfinite(expected{index}(:))),diagnostic);
    relativeError = norm(actual{index}(:)-expected{index}(:))/max(norm(expected{index}(:)),realmin);
    testCase.verifyLessThanOrEqual(relativeError,1e-12,diagnostic);
end
end

function deleteTransforms(varargin)
for transform = varargin
    if isvalid(transform{1}), delete(transform{1}); end
end
end
