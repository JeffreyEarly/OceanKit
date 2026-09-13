classdef TestWVCompiledConsumers < matlab.unittest.TestCase
    % Focused consumer tests for compiled state scopes and MATLAB observers.
    methods (Test,TestTags="optional")
        function nonlinearAndSpeedSharePhysicalProducersInEitherOrder(testCase)
            definitions = configurations();
            for definition = definitions(1:4)
                expected = definition.create("matlab");
                actual = definition.create("compiled");
                cleanup = onCleanup(@()deleteTransforms(expected,actual));
                seedState(expected,actual);
                expectedFlux = cell(1,3);
                [expectedFlux{:}] = expected.nonlinearFlux();
                for speedFirst = [false true]
                    before = actual.computationalBackendMetadata.runtimeMetrics;
                    scope = actual.scopedEvaluation();
                    if speedFirst, actual.compiledVariables("uvMax"); end
                    flux = cell(1,3); [flux{:}] = actual.nonlinearFlux();
                    speed = actual.compiledVariables("uvMax");
                    repeated = actual.compiledVariables("uvMax");
                    verifyNear(testCase,speed{1},expected.uvMax,definition.name+" speed reduction");
                    testCase.verifyEqual(repeated,speed,definition.name);
                    for channel = 1:3
                        verifyNear(testCase,flux{channel},expectedFlux{channel},definition.name+" speed/flux ordering");
                    end
                    clear scope
                    after = actual.computationalBackendMetadata.runtimeMetrics;
                    testCase.verifyEqual(after.physicalReconstructionExecutions-before.physicalReconstructionExecutions,4,definition.name);
                    testCase.verifyEqual(after.nonlinearProducerExecutions-before.nonlinearProducerExecutions,1,definition.name);
                    testCase.verifyEqual(after.duplicateExecutions-before.duplicateExecutions,0,definition.name);
                end
                clear cleanup
            end
        end

        function modifiedNonlinearDensityCorrectionPreservesMatlabBehavior(testCase)
            definitions = configurations();
            for definition = definitions(1:4)
                expected = definition.create("matlab");
                actual = definition.create("compiled");
                cleanup = onCleanup(@()deleteTransforms(expected,actual));
                seedState(expected,actual);
                for wvt = [expected actual]
                    forcing = wvt.forcingWithName("nonlinear advection");
                    forcing.dLnN2 = forcing.dLnN2 + 1e-3;
                end
                expectedFlux = cell(1,3); actualFlux = cell(1,3);
                [expectedFlux{:}] = expected.nonlinearFlux();
                before = actual.computationalBackendMetadata.runtimeMetrics;
                [actualFlux{:}] = actual.nonlinearFlux();
                for channel = 1:3
                    verifyNear(testCase,actualFlux{channel},expectedFlux{channel},definition.name+" modified density correction");
                end
                after = actual.computationalBackendMetadata.runtimeMetrics;
                testCase.verifyEqual(after.nonlinearProducerExecutions,before.nonlinearProducerExecutions,definition.name);
                clear cleanup
            end
        end

        function nonlinearCoefficientBoundarySharesRawAndPhysicalFields(testCase)
            definitions = configurations();
            for definition = definitions(1:4)
                expected = definition.create("matlab");
                actual = definition.create("compiled");
                cleanup = onCleanup(@()deleteTransforms(expected,actual));
                seedState(expected,actual);
                expectedFlux = cell(1,3);
                [expectedFlux{:}] = expected.nonlinearFlux();
                names = {'Fu_nonlinear_advection','Fv_nonlinear_advection','Feta_nonlinear_advection'};
                expectedRaw = cell(1,4);
                forcing = expected.forcingWithName("nonlinear advection");
                if isa(actual,'WVTransformBoussinesq') || (isa(actual,'WVTransformConstantStratification') && ~actual.isHydrostatic)
                    names{end+1} = 'Fw_nonlinear_advection';
                    [expectedRaw{1},expectedRaw{2},expectedRaw{4},expectedRaw{3}] = forcing.addNonhydrostaticSpatialForcing(expected,0,0,0,0);
                else
                    [expectedRaw{1:3}] = forcing.addHydrostaticSpatialForcing(expected,0,0,0);
                end
                for rawFirst = [false true]
                    before = actual.computationalBackendMetadata.runtimeMetrics;
                    scope = actual.scopedEvaluation();
                    if rawFirst, actual.compiledVariables(names); end
                    first = cell(1,3); second = cell(1,3);
                    [first{:}] = actual.nonlinearFlux();
                    [second{:}] = actual.nonlinearFlux();
                    raw = actual.compiledVariables(names);
                    for channel = 1:3
                        verifyNear(testCase,first{channel},expectedFlux{channel},definition.name+" coefficient boundary");
                        testCase.verifyEqual(second{channel},first{channel},definition.name);
                    end
                    for channel = 1:numel(names)
                        verifyNear(testCase,raw{channel},expectedRaw{channel},definition.name+" raw tendency");
                    end
                    beforeFields = actual.computationalBackendMetadata.runtimeMetrics;
                    fieldNames = ["u","v","w","eta"];
                    fields = actual.compiledVariables(fieldNames);
                    for fieldIndex = 1:numel(fieldNames)
                        name = fieldNames(fieldIndex);
                        verifyNear(testCase,fields{fieldIndex},expected.(name),definition.name+" shared "+name);
                    end
                    afterFields = actual.computationalBackendMetadata.runtimeMetrics;
                    testCase.verifyEqual(afterFields.physicalReconstructionExecutions,beforeFields.physicalReconstructionExecutions,definition.name);
                    clear scope
                    after = actual.computationalBackendMetadata.runtimeMetrics;
                    testCase.verifyEqual(after.nonlinearProducerExecutions-before.nonlinearProducerExecutions,1,definition.name);
                    testCase.verifyEqual(after.stateValidations-before.stateValidations,1,definition.name);
                    testCase.verifyEqual(after.phasePreparations-before.phasePreparations,1,definition.name);
                    testCase.verifyEqual(after.duplicateExecutions-before.duplicateExecutions,0,definition.name);
                end
                clear cleanup
            end
        end

        function nonlinearAdaptiveFluxMatchesAcrossAllFamilies(testCase)
            for definition = configurations()
                matlabWVT = definition.create("matlab");
                compiledWVT = definition.create("compiled");
                cleanup = onCleanup(@()deleteTransforms(matlabWVT,compiledWVT));
                seedState(matlabWVT,compiledWVT);
                for wvt = [matlabWVT compiledWVT]
                    wvt.setForcing([WVNonlinearAdvection(wvt),WVAdaptiveDamping(wvt)]);
                end
                matlabModel = WVModel(matlabWVT);
                compiledModel = WVModel(compiledWVT);
                matlabModel.setupIntegrator(integratorType="adaptive",integrator=@ode45,absTolerance=1e-10,relTolerance=1e-9,shouldShowIntegrationStats=0);
                compiledModel.setupIntegrator(integratorType="adaptive",integrator=@ode45,absTolerance=1e-10,relTolerance=1e-9,shouldShowIntegrationStats=0);
                expectedDiagnostics = matlabWVT.fluxForForcing();
                actualDiagnostics = compiledWVT.fluxForForcing();
                testCase.verifyEqual(sort(string(expectedDiagnostics.keys)),sort(string(actualDiagnostics.keys)),definition.name);
                for name = string(expectedDiagnostics.keys)
                    expectedValue = expectedDiagnostics{name}; actualValue = actualDiagnostics{name};
                    if isstruct(expectedValue)
                        for channel = string(fieldnames(expectedValue))'
                            verifyNear(testCase,actualValue.(channel),expectedValue.(channel),definition.name+" forcing diagnostic "+name+" "+channel);
                        end
                    else
                        verifyNear(testCase,actualValue,expectedValue,definition.name+" forcing diagnostic "+name);
                    end
                end
                matlabY = matlabModel.initialConditionsCellArray();
                compiledY = compiledModel.initialConditionsCellArray();
                expected = matlabModel.fluxAtTimeCellArray(13,matlabY);
                before = compiledWVT.computationalBackendMetadata.runtimeMetrics;
                actual = compiledModel.fluxAtTimeCellArray(13,compiledY);
                testCase.verifyEqual(numel(actual),numel(expected),definition.name);
                for index = 1:numel(expected)
                    verifyNear(testCase,actual{index},expected{index},definition.name+" forcing flux");
                end
                metrics = compiledWVT.computationalBackendMetadata.runtimeMetrics;
                testCase.verifyEqual(metrics.duplicateExecutions,0,definition.name);
                testCase.verifyEqual(metrics.stateValidations-before.stateValidations,1,definition.name);
                testCase.verifyEqual(metrics.scopedEvaluations-before.scopedEvaluations,1,definition.name);
                testCase.verifyGreaterThan(metrics.producerExecutions,0,definition.name);
                testCase.verifyTrue(any(cellfun(@(value)any(abs(value(:))>0),actual)),definition.name);
                clear cleanup
            end
        end

        function sampledObserversShareCompiledStateScope(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            definitions = configurations();
            references = cell(1,2);
            for iBackend = 1:2
                backends = ["matlab","compiled"];
                wvt = definitions(3).create(backends(iBackend));
                cleanup = onCleanup(@()delete(wvt));
                rng(511,"twister"); wvt.initWithRandomFlow(uvMax=.01);
                wvt.setForcing([WVNonlinearAdvection(wvt),WVAdaptiveDamping(wvt)]);
                model = WVModel(wvt);
                path = fullfile(fixture.Folder,backends(iBackend)+"-consumers.nc");
                output = model.createNetCDFFileForModelOutput(path,outputInterval=.5,shouldOverwriteExisting=true);
                fileCleanup = onCleanup(@()model.closeNetCDFFile());
                model.eulerianObservingSystem.addNetCDFOutputVariables("u","v");
                model.setFloatPositions([500 1500],[400 1200],[-250 -750],"u",trackedVarInterpolation="linear",advectionInterpolation="spline");
                model.addTracer(reshape(1:prod(wvt.spatialMatrixSize),wvt.spatialMatrixSize),"dye");
                shared = output.addNewEvenlySpacedOutputGroup("shared",outputInterval=.5);
                shared.addObservingSystem([model.fluxedObservingSystemWithName("float"),model.fluxedObservingSystemWithName("dye")]);
                mooring = WVMooring(model,name="mooring",x=[0 1000],y=[0 900],trackedFieldNames={'u'});
                output.outputGroupWithName(model.defaultOutputGroupName()).addObservingSystem(mooring);
                model.setupIntegrator(integratorType="fixed",deltaT=.25);
                state = model.initialConditionsCellArray();
                if iBackend==2, before = wvt.computationalBackendMetadata.runtimeMetrics; end
                rhs = model.fluxAtTimeCellArray(0,state);
                if iBackend==2
                    after = wvt.computationalBackendMetadata.runtimeMetrics;
                    testCase.verifyEqual(after.scopedEvaluations-before.scopedEvaluations,1);
                    testCase.verifyEqual(after.stateValidations-before.stateValidations,1);
                    testCase.verifyEqual(after.phasePreparations-before.phasePreparations,1);
                    testCase.verifyEqual(after.stateInputCopyBytes-before.stateInputCopyBytes,3*wvt.Nj*wvt.Nkl*16);
                    testCase.verifyEqual(after.duplicateExecutions,0);
                    for i = 1:numel(rhs), verifyNear(testCase,rhs{i},references{1}.rhs{i},"observer RHS"); end
                end
                model.integrateToTime(.5,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
                [x,y,z,tracked] = model.floatPositions();
                references{iBackend} = struct("rhs",{rhs},"u",wvt.u,"A0",wvt.A0,"Ap",wvt.Ap,"Am",wvt.Am,"tracer",model.tracer("dye"),"x",x,"y",y,"z",z,"trackedU",tracked.u);
                if iBackend==2
                    for field = ["u","A0","Ap","Am","tracer","x","y","z","trackedU"]
                        verifyNear(testCase,references{2}.(field),references{1}.(field),"integrated "+field);
                    end
                    testCase.verifyEqual(wvt.computationalBackendMetadata.runtimeMetrics.duplicateExecutions,0);
                end
                model.closeNetCDFFile(); clear fileCleanup
                file = NetCDFFile(path,shouldReadOnly=true); readerCleanup = onCleanup(@()closeIfOpen(file));
                storedU = file.readVariables('wave-vortex/u');
                verifyNear(testCase,storedU(:,:,:,end),references{iBackend}.u,"stored Eulerian velocity");
                file.close(); clear readerCleanup cleanup
            end
        end

        function coefficientSubclassesPreserveCallbacks(testCase)
            WVCompiledBackend.activateModule(WVCompiledBackend.capabilities());
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fileparts(mfilename('fullpath'))));
            definitions = configurations();
            for className = ["WVTestInheritedCoefficients","WVTestOverriddenCoefficients"]
                references = cell(1,2);
                for iBackend = 1:2
                    backends = ["matlab","compiled"];
                    wvt = definitions(3).create(backends(iBackend)); cleanup = onCleanup(@()delete(wvt));
                    rng(514,"twister"); wvt.initWithRandomFlow(uvMax=.01);
                    model = WVModel(wvt);
                    model.removeFluxedObservingSystem(model.wvCoefficientFluxedObservingSystem());
                    observer = feval(className,model); model.addFluxedCoefficients(observer);
                    model.addTracer(reshape(1:prod(wvt.spatialMatrixSize),wvt.spatialMatrixSize),"dye");
                    before = wvt.computationalBackendMetadata.runtimeMetrics;
                    references{iBackend} = model.fluxAtTimeCellArray(3,model.initialConditionsCellArray());
                    if className=="WVTestOverriddenCoefficients", testCase.verifyEqual(observer.fluxCalls,1); end
                    if iBackend==2
                        for i=1:numel(references{2}), verifyNear(testCase,references{2}{i},references{1}{i},className+" RHS"); end
                        after = wvt.computationalBackendMetadata.runtimeMetrics;
                        testCase.verifyEqual(after.duplicateExecutions,0);
                        if className=="WVTestInheritedCoefficients"
                            testCase.verifyEqual(after.scopedEvaluations-before.scopedEvaluations,1);
                        end
                    end
                    clear cleanup
                end
            end
        end

        function registryMutationInvalidatesOuterScope(testCase)
            definitions = configurations();
            wvt = definitions(3).create("compiled");
            cleanup = onCleanup(@()delete(wvt));
            wvt.t0 = 0; wvt.t = 0;
            rng(511,"twister");
            wvt.initWithRandomFlow(uvMax=.01);
            scope = wvt.scopedEvaluation();
            wvt.addOperation(WVOperation("consumer_scope_sentinel",WVVariableAnnotation("consumer_scope_sentinel",{}, "1","scope sentinel"),@(~)37),shouldSuppressWarning=true);
            testCase.verifyError(@()wvt.variableWithName("u"),"WaveVortexModel:CompiledTransformStateChanged");
            clear scope
            testCase.verifyEqual(wvt.consumer_scope_sentinel,37);
            testCase.verifyWarningFree(@()wvt.variableWithName("u"));
            scope = wvt.scopedEvaluation();
            wvt.addForcing(WVAdaptiveDamping(wvt));
            testCase.verifyError(@()wvt.variableWithName("v"),"WaveVortexModel:CompiledTransformStateChanged");
            clear scope
            testCase.verifyWarningFree(@()wvt.variableWithName("v"));
            clear cleanup
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

function verifyNear(testCase,actual,expected,diagnostic)
testCase.verifySize(actual,size(expected),diagnostic);
testCase.verifyTrue(all(isfinite(actual(:))) && all(isfinite(expected(:))),diagnostic);
testCase.verifyLessThanOrEqual(norm(actual(:)-expected(:))/max(norm(expected(:)),realmin),1e-12,diagnostic);
end

function deleteTransforms(varargin)
for transform = varargin
    if isvalid(transform{1}), delete(transform{1}); end
end
end

function closeIfOpen(file)
if isvalid(file) && ~isempty(file.id), file.close(); end
end
