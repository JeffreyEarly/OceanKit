classdef TestWVCompiledPublicAdapters < matlab.unittest.TestCase
    % Public MATLAB behavior backed by the shared compiled transform session.
    methods (Test,TestTags="optional")
        function verticalCalculusPreservesMatlabLayouts(testCase)
            definitions = configurations();
            for definition = definitions(1:5)
                expected = definition.create("matlab"); actual = definition.create("compiled");
                cleanup = onCleanup(@()deleteTransforms(expected,actual));
                rng(510,"twister");
                volume = rand(expected.spatialMatrixSize);
                for family = ["F","G"]
                    derivative = "diffZ"+family; integral = "intZ"+family;
                    for order = 1:4
                        verifyNear(testCase,feval(derivative,actual,volume,'n',order),feval(derivative,expected,volume,'n',order),definition.name+" "+derivative+" order "+order);
                    end
                    verifyNear(testCase,feval(integral,actual,volume),feval(integral,expected,volume),definition.name+" "+integral+" volume");
                    for columns = [0 1 5]
                        matrix = complex(rand(expected.Nz,columns),rand(expected.Nz,columns));
                        verifyNear(testCase,feval(integral,actual,matrix),feval(integral,expected,matrix),definition.name+" "+integral+" matrix");
                    end
                    complexVolume = complex(volume,.5*flip(volume,3));
                    verifyNear(testCase,feval(derivative,actual,complexVolume),feval(derivative,expected,complexVolume),definition.name+" "+derivative+" complex");
                end
                clear cleanup
            end
        end

        function densityReferencesShareOneImmutableState(testCase)
            definitions = configurations();
            for definition = definitions(1:4)
                expected = definition.create("matlab"); actual = definition.create("compiled");
                cleanup = onCleanup(@()deleteTransforms(expected,actual));
                seedMatchingState(expected,actual);
                mda = expected.J>0 & expected.Kh==0;
                expected.A0(mda) = .003*sin(find(mda)); actual.A0 = expected.A0;
                names = {'rho_nm','eta_true','ape','apv'};
                scope = actual.scopedEvaluation();
                before = actual.computationalBackendMetadata.runtimeMetrics;
                actualEta = [];
                for useActual = [true false true]
                    expected.shouldUseTrueNoMotionProfile = useActual;
                    actual.shouldUseTrueNoMotionProfile = useActual;
                    reference = cell(1,4); native = cell(1,4);
                    [reference{:}] = expected.variableWithName(names{:});
                    [native{:}] = actual.variableWithName(names{:});
                    span = max(expected.rho_total,[],"all")-min(expected.rho_total,[],"all");
                    % These are the existing qualified inverse-density bounds.
                    tolerances = [1e-6*span,1e-6*expected.Lz,1e-5*max(abs(reference{3}),[],"all")+1e-12,1e-6*max(abs(reference{4}),[],"all")+1e-12];
                    for i=1:4
                        testCase.verifySize(native{i},size(reference{i}),definition.name+" "+names{i});
                        testCase.assertTrue(all(isfinite(native{i}),"all") && all(isfinite(reference{i}),"all"));
                        testCase.verifyLessThanOrEqual(max(abs(native{i}-reference{i}),[],"all"),tolerances(i),definition.name+" "+names{i});
                    end
                    if useActual
                        actualEta = native{2};
                    else
                        testCase.verifyGreaterThan(max(abs(native{2}-actualEta),[],"all"),1e-8,definition.name+" density references are distinct");
                    end
                end
                after = actual.computationalBackendMetadata.runtimeMetrics;
                testCase.verifyEqual(after.stateValidations-before.stateValidations,0,definition.name);
                testCase.verifyEqual(after.duplicateExecutions,0,definition.name);
                clear scope
                report = actual.compiledDensityRecoveryReport();
                testCase.verifyEqual(string(report.solver),"dampedLeastSquares");
                testCase.verifyTrue(report.qualified);
                clear cleanup
            end
        end

        function allPublicFamiliesMatchMatlab(testCase)
            for definition = configurations()
                matlabWVT = definition.create("matlab");
                compiledWVT = definition.create("compiled");
                cleanup = onCleanup(@()deleteTransforms(matlabWVT,compiledWVT));
                seedMatchingState(matlabWVT,compiledWVT);
                diagnostic = definition.name;

                metadata = compiledWVT.computationalBackendMetadata;
                testCase.verifyEqual(compiledWVT.computationalBackend,"compiled",diagnostic);
                testCase.verifyEqual(metadata.requestedBackend,"compiled",diagnostic);
                testCase.verifyEqual(metadata.activeBackend,"compiled",diagnostic);
                testCase.verifyGreaterThan(metadata.runtimeMetrics.engineBytes,0,diagnostic);
                testCase.verifyEqual(metadata.runtimeMetrics.duplicateExecutions,0,diagnostic);

                verifyRawPublicMethods(testCase,matlabWVT,compiledWVT,diagnostic);
                verifyForwardTransform(testCase,matlabWVT,compiledWVT,diagnostic);
                verifyStateVariables(testCase,matlabWVT,compiledWVT,diagnostic);

                expectedFlux = flux(matlabWVT); actualFlux = flux(compiledWVT);
                testCase.verifyEqual(numel(actualFlux),numel(expectedFlux),diagnostic);
                for index = 1:numel(expectedFlux)
                    verifyNear(testCase,actualFlux{index},expectedFlux{index},diagnostic+" nonlinear flux");
                end

                custom = WVVariableAnnotation('compiled_public_callback',{},'1','compiled public callback sentinel');
                compiledWVT.addOperation(WVOperation('compiled_public_callback',custom,@(~)37),shouldSuppressWarning=true);
                testCase.verifyEqual(compiledWVT.compiled_public_callback,37,diagnostic);

                matlabReplacementAnnotation = WVVariableAnnotation('u',matlabWVT.spatialDimensionNames,'m s-1','replacement velocity callback');
                replacement = WVOperation('u',matlabReplacementAnnotation,@(wvt)0.125*ones(wvt.spatialMatrixSize));
                matlabWVT.addOperation(replacement,shouldOverwriteExisting=true,shouldSuppressWarning=true);
                compiledReplacementAnnotation = WVVariableAnnotation('u',compiledWVT.spatialDimensionNames,'m s-1','replacement velocity callback');
                compiledReplacement = WVOperation('u',compiledReplacementAnnotation,@(wvt)0.125*ones(wvt.spatialMatrixSize));
                compiledWVT.addOperation(compiledReplacement,shouldOverwriteExisting=true,shouldSuppressWarning=true);
                testCase.verifyEqual(compiledWVT.u,0.125*ones(compiledWVT.spatialMatrixSize),diagnostic);
                expectedFlux = flux(matlabWVT); actualFlux = flux(compiledWVT);
                for index = 1:numel(expectedFlux)
                    verifyNear(testCase,actualFlux{index},expectedFlux{index},diagnostic+" custom-u flux");
                end

                scope = compiledWVT.scopedEvaluation();
                mutateState(compiledWVT);
                testCase.verifyError(@()compiledWVT.variableWithName('v'),"WaveVortexModel:CompiledTransformStateChanged",diagnostic);
                clear scope
                testCase.verifyWarningFree(@()compiledWVT.variableWithName('v'),diagnostic);
                clear cleanup
            end
        end
    end
end

function definitions = configurations()
profile = @(z)1e-4*exp(z/700);
definitions = [ ...
    struct("name","constant-hydrostatic","create",@(backend)WVTransformConstantStratification([4000 3000 1000],[6 6 5],N0=5.2e-3,isHydrostatic=true,shouldAntialias=true,computationalBackend=backend)) ...
    struct("name","constant-nonhydrostatic","create",@(backend)WVTransformConstantStratification([4000 3000 1000],[6 6 5],N0=5.2e-3,isHydrostatic=false,shouldAntialias=true,computationalBackend=backend)) ...
    struct("name","hydrostatic","create",@(backend)WVTransformHydrostatic([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=true,computationalBackend=backend)) ...
    struct("name","boussinesq","create",@(backend)WVTransformBoussinesq([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=true,computationalBackend=backend)) ...
    struct("name","stratified-qg","create",@(backend)WVTransformStratifiedQG([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=true,computationalBackend=backend)) ...
    struct("name","barotropic-qg","create",@(backend)WVTransformBarotropicQG([4000 3000],[6 6],h=1000,j=1,shouldAntialias=true,computationalBackend=backend))];
end

function seedMatchingState(matlabWVT,compiledWVT)
rng(507,"twister"); matlabWVT.initWithRandomFlow(uvMax=0.01);
if ~isQG(matlabWVT)
    compiledWVT.Ap = matlabWVT.Ap; compiledWVT.Am = matlabWVT.Am;
end
compiledWVT.A0 = matlabWVT.A0;
matlabWVT.t0 = 5; matlabWVT.t = 13; compiledWVT.t0 = 5; compiledWVT.t = 13;
end

function verifyRawPublicMethods(testCase,matlabWVT,compiledWVT,diagnostic)
rng(508,"twister"); volume = rand(matlabWVT.spatialMatrixSize);
expected = matlabWVT.transformFromSpatialDomainWithFourier(volume);
actual = compiledWVT.transformFromSpatialDomainWithFourier(volume);
verifyNear(testCase,actual,expected,diagnostic+" Fourier forward");
verifyNear(testCase,compiledWVT.transformToSpatialDomainWithFourier(complex(expected)),matlabWVT.transformToSpatialDomainWithFourier(expected),diagnostic+" Fourier inverse");
for order = [1 2 4]
    verifyNear(testCase,compiledWVT.diffX(volume,n=order),matlabWVT.diffX(volume,n=order),diagnostic+" diffX");
    verifyNear(testCase,compiledWVT.diffY(volume,n=order),matlabWVT.diffY(volume,n=order),diagnostic+" diffY");
end
if isa(matlabWVT,"WVTransformBarotropicQG"), return, end
grid = complex(rand(matlabWVT.Nz,matlabWVT.Nkl),rand(matlabWVT.Nz,matlabWVT.Nkl));
verifyNear(testCase,compiledWVT.transformFromSpatialDomainWithFg(grid),matlabWVT.transformFromSpatialDomainWithFg(grid),diagnostic+" F projection");
verifyNear(testCase,compiledWVT.transformFromSpatialDomainWithGg(grid),matlabWVT.transformFromSpatialDomainWithGg(grid),diagnostic+" G projection");
column = complex(rand(matlabWVT.Nz,1),rand(matlabWVT.Nz,1));
verifyNear(testCase,compiledWVT.transformFromSpatialDomainWithFio(column),matlabWVT.transformFromSpatialDomainWithFio(column),diagnostic+" Fio projection");
end

function verifyForwardTransform(testCase,matlabWVT,compiledWVT,diagnostic)
if isQG(matlabWVT), return, end
rng(509,"twister");
u = rand(matlabWVT.spatialMatrixSize); v = rand(matlabWVT.spatialMatrixSize); eta = rand(matlabWVT.spatialMatrixSize);
expected = cell(1,3); actual = cell(1,3);
if isa(matlabWVT,"WVTransformBoussinesq") || (isa(matlabWVT,"WVTransformConstantStratification") && ~matlabWVT.isHydrostatic)
    w = rand(matlabWVT.spatialMatrixSize);
    [expected{:}] = matlabWVT.transformUVWEtaToWaveVortex(u,v,w,eta);
    [actual{:}] = compiledWVT.transformUVWEtaToWaveVortex(u,v,w,eta);
else
    [expected{:}] = matlabWVT.transformUVEtaToWaveVortex(u,v,eta);
    [actual{:}] = compiledWVT.transformUVEtaToWaveVortex(u,v,eta);
end
for index = 1:3
    verifyNear(testCase,actual{index},expected{index},diagnostic+" forward transform");
end
end

function verifyStateVariables(testCase,matlabWVT,compiledWVT,diagnostic)
if isQG(matlabWVT)
    stateNames = {'A0t'}; fieldNames = {'u','v'};
else
    stateNames = {'phase','conjPhase','Apt','Amt','A0t'}; fieldNames = {'u','v','w','eta'};
end
verifyVariables(testCase,matlabWVT,compiledWVT,stateNames,diagnostic+" current state");
componentCandidates = {'u_g','v_g','w_g','eta_g','u_w','v_w','w_w','eta_w','u_io','v_io','w_io','eta_io','u_mda','v_mda','w_mda','eta_mda'};
componentNames = componentCandidates(cellfun(@(name)matlabWVT.hasVariableWithName(name) && compiledWVT.hasVariableWithName(name),componentCandidates));
scope = compiledWVT.scopedEvaluation();
before = compiledWVT.computationalBackendMetadata.runtimeMetrics;
verifyVariables(testCase,matlabWVT,compiledWVT,fieldNames,diagnostic+" fields");
verifyVariables(testCase,matlabWVT,compiledWVT,componentNames,diagnostic+" components");
middle = compiledWVT.computationalBackendMetadata.runtimeMetrics;
verifyVariables(testCase,matlabWVT,compiledWVT,fieldNames,diagnostic+" repeated fields");
after = compiledWVT.computationalBackendMetadata.runtimeMetrics;
testCase.verifyGreaterThan(middle.producerExecutions,before.producerExecutions,diagnostic);
testCase.verifyEqual(after.producerExecutions,middle.producerExecutions,diagnostic);
testCase.verifyEqual(after.duplicateExecutions,0,diagnostic);
clear scope
testCase.verifyEqual(compiledWVT.computationalBackendMetadata.runtimeMetrics.liveEvaluationBytes,0,diagnostic);
end

function verifyVariables(testCase,matlabWVT,compiledWVT,names,diagnostic)
if isempty(names), return, end
expected = cell(size(names)); actual = cell(size(names));
[expected{:}] = matlabWVT.variableWithName(names{:});
[actual{:}] = compiledWVT.variableWithName(names{:});
for index = 1:numel(names)
    verifyNear(testCase,actual{index},expected{index},diagnostic+" "+names{index});
end
end

function values = flux(wvt)
if isQG(wvt)
    values = {wvt.nonlinearFlux()};
else
    values = cell(1,3); [values{:}] = wvt.nonlinearFlux();
end
end

function value = isQG(wvt)
value = isa(wvt,"WVTransformStratifiedQG") || isa(wvt,"WVTransformBarotropicQG");
end

function mutateState(wvt)
if isQG(wvt)
    value = wvt.A0; value(2) = value(2)+1e-8; wvt.A0 = value;
else
    value = wvt.Ap; value(2) = value(2)+1e-8; wvt.Ap = value;
end
end

function verifyNear(testCase,actual,expected,diagnostic)
testCase.verifySize(actual,size(expected),diagnostic);
testCase.assertTrue(all(isfinite(actual(:))) && all(isfinite(expected(:))),diagnostic);
error = norm(actual(:)-expected(:))/max(norm(expected(:)),realmin);
testCase.verifyLessThanOrEqual(error,1e-12,diagnostic);
end

function deleteTransforms(varargin)
for transform = varargin
    if isvalid(transform{1}), delete(transform{1}), end
end
end
