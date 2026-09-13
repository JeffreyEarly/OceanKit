classdef TestWVCompiledAdapterRestore < matlab.unittest.TestCase
    % Factory and reconstruction ownership for compiled public transforms.

    properties (SetAccess=private)
        temporaryFolder (1,1) string
    end

    methods (TestMethodSetup)
        function createTemporaryFolder(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.temporaryFolder = string(fixture.Folder);
        end
    end

    methods (Test,TestTags="optional")
        function allFamiliesRestoreWithFreshRuntimeOwnership(testCase)
            capabilities = WVCompiledBackend.capabilities();
            testCase.assumeTrue(capabilities.isAvailable,capabilities.failure.message);
            definitions = configurations();
            for index = 1:numel(definitions)
                definition = definitions(index);
                source = definition.create("compiled");
                cleanup = onCleanup(@()deleteTransforms(source));
                seedState(source,5100+index);
                path = fullfile(testCase.temporaryFolder,definition.name+".nc");
                file = source.writeToFile(char(path),shouldOverwriteExisting=true);
                file.close();

                matlabRestored = WVTransform.waveVortexTransformFromFile(char(path));
                cleanupMatlab = onCleanup(@()deleteTransforms(matlabRestored));
                testCase.verifyEqual(matlabRestored.computationalBackend,"matlab",definition.name);
                verifyWritableOpen(testCase,path);

                [compiledRestored,reader] = WVTransform.waveVortexTransformFromFile(char(path),computationalBackend="compiled");
                cleanupCompiled = onCleanup(@()deleteTransforms(compiledRestored));
                cleanupReader = onCleanup(@()closeIfOpen(reader));
                testCase.verifyNotEmpty(reader.id,definition.name);
                testCase.verifyClass(compiledRestored,class(source),definition.name);
                verifyCompiled(testCase,compiledRestored,definition.name);
                verifyFreshIdentity(testCase,source,matlabRestored,compiledRestored,definition.name);
                verifyState(testCase,source,matlabRestored,compiledRestored,definition.name);
                verifyFieldsAndFlux(testCase,matlabRestored,compiledRestored,definition.name+" restore");

                groupFactory = char(string(class(source))+".transformFromGroup");
                groupRestored = feval(groupFactory,reader,'computationalBackend','compiled');
                cleanupGroup = onCleanup(@()deleteTransforms(groupRestored));
                verifyCompiled(testCase,groupRestored,definition.name+" group");
                testCase.verifyFalse(groupRestored.compiledSourceIdentity == source.compiledSourceIdentity,definition.name);
                clear cleanupGroup

                reader.close();
                clear cleanupReader
                verifyWritableOpen(testCase,path);

                matlabResolution = matlabRestored.waveVortexTransformWithResolution(definition.resolution);
                compiledResolution = compiledRestored.waveVortexTransformWithResolution(definition.resolution);
                cleanupResolution = onCleanup(@()deleteTransforms(matlabResolution,compiledResolution));
                verifyCompiled(testCase,compiledResolution,definition.name+" resolution");
                testCase.verifyFalse(compiledResolution.compiledSourceIdentity == compiledRestored.compiledSourceIdentity,definition.name);
                verifyForcing(testCase,compiledRestored,compiledResolution,definition.name+" resolution");
                verifyFieldsAndFlux(testCase,matlabResolution,compiledResolution,definition.name+" resolution");
                clear cleanupResolution

                if definition.explicitAntialias
                    matlabExplicit = matlabRestored.waveVortexTransformWithExplicitAntialiasing();
                    compiledExplicit = compiledRestored.waveVortexTransformWithExplicitAntialiasing();
                    cleanupExplicit = onCleanup(@()deleteTransforms(matlabExplicit,compiledExplicit));
                    verifyCompiled(testCase,compiledExplicit,definition.name+" explicit antialias");
                    testCase.verifyFalse(compiledExplicit.compiledSourceIdentity == compiledRestored.compiledSourceIdentity,definition.name);
                    verifyForcing(testCase,matlabExplicit,compiledExplicit,definition.name+" explicit antialias");
                    testCase.verifyTrue(compiledExplicit.hasForcingWithName("antialias filter"),definition.name);
                    verifyFieldsAndFlux(testCase,matlabExplicit,compiledExplicit,definition.name+" explicit antialias");
                    clear cleanupExplicit
                end

                clear cleanupCompiled cleanupMatlab cleanup
            end
        end
    end
end

function definitions = configurations()
profile = @(z)1e-4*exp(z/700);
definitions = [ ...
    struct("name","constant-hydrostatic","create",@(backend)WVTransformConstantStratification([4000 3000 1000],[6 6 5],N0=5.2e-3,isHydrostatic=true,shouldAntialias=true,computationalBackend=backend),"resolution",[8 6 6],"explicitAntialias",true) ...
    struct("name","constant-nonhydrostatic","create",@(backend)WVTransformConstantStratification([4000 3000 1000],[6 6 5],N0=5.2e-3,isHydrostatic=false,shouldAntialias=true,computationalBackend=backend),"resolution",[8 6 6],"explicitAntialias",true) ...
    struct("name","hydrostatic","create",@(backend)WVTransformHydrostatic([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=true,computationalBackend=backend),"resolution",[8 6 6],"explicitAntialias",true) ...
    struct("name","boussinesq","create",@(backend)WVTransformBoussinesq([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=true,computationalBackend=backend),"resolution",[8 6 6],"explicitAntialias",true) ...
    struct("name","stratified-qg","create",@(backend)WVTransformStratifiedQG([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=true,computationalBackend=backend),"resolution",[8 6 6],"explicitAntialias",false) ...
    struct("name","barotropic-qg","create",@(backend)WVTransformBarotropicQG([4000 3000],[6 6],h=1000,j=1,shouldAntialias=true,computationalBackend=backend),"resolution",[8 6],"explicitAntialias",false)];
end

function seedState(wvt,seed)
rng(seed,"twister");
wvt.initWithRandomFlow(uvMax=0.01);
wvt.t0 = 5;
wvt.t = 13;
end

function verifyCompiled(testCase,wvt,diagnostic)
metadata = wvt.computationalBackendMetadata;
testCase.verifyEqual(wvt.computationalBackend,"compiled",diagnostic);
testCase.verifyEqual(metadata.requestedBackend,"compiled",diagnostic);
testCase.verifyEqual(metadata.activeBackend,"compiled",diagnostic);
testCase.verifyTrue(metadata.module.identityValidated,diagnostic);
testCase.verifyGreaterThan(metadata.runtimeMetrics.engineBytes,0,diagnostic);
end

function verifyFreshIdentity(testCase,source,matlabRestored,compiledRestored,diagnostic)
testCase.verifyFalse(source.compiledSourceIdentity == matlabRestored.compiledSourceIdentity,diagnostic);
testCase.verifyFalse(source.compiledSourceIdentity == compiledRestored.compiledSourceIdentity,diagnostic);
testCase.verifyFalse(matlabRestored.compiledSourceIdentity == compiledRestored.compiledSourceIdentity,diagnostic);
end

function verifyState(testCase,source,matlabRestored,compiledRestored,diagnostic)
testCase.verifyEqual(matlabRestored.t,source.t,diagnostic);
testCase.verifyEqual(matlabRestored.t0,source.t0,diagnostic);
testCase.verifyEqual(compiledRestored.t,source.t,diagnostic);
testCase.verifyEqual(compiledRestored.t0,source.t0,diagnostic);
testCase.verifyEqual(matlabRestored.A0,source.A0,diagnostic);
testCase.verifyEqual(compiledRestored.A0,source.A0,diagnostic);
if ~isQG(source)
    testCase.verifyEqual(matlabRestored.Ap,source.Ap,diagnostic);
    testCase.verifyEqual(matlabRestored.Am,source.Am,diagnostic);
    testCase.verifyEqual(compiledRestored.Ap,source.Ap,diagnostic);
    testCase.verifyEqual(compiledRestored.Am,source.Am,diagnostic);
end
verifyForcing(testCase,source,matlabRestored,diagnostic);
verifyForcing(testCase,source,compiledRestored,diagnostic);
end

function verifyForcing(testCase,expected,actual,diagnostic)
expectedForcing = expected.forcing;
actualForcing = actual.forcing;
testCase.verifyEqual(string({actualForcing.name}),string({expectedForcing.name}),diagnostic);
testCase.verifyEqual(string(arrayfun(@class,actualForcing,UniformOutput=false)),string(arrayfun(@class,expectedForcing,UniformOutput=false)),diagnostic);
end

function verifyFieldsAndFlux(testCase,expected,actual,diagnostic)
names = {'u','v'};
if ~isQG(expected), names = {'u','v','w','eta'}; end
expectedFields = cell(size(names)); actualFields = cell(size(names));
[expectedFields{:}] = expected.variableWithName(names{:});
[actualFields{:}] = actual.variableWithName(names{:});
for index = 1:numel(names)
    verifyNear(testCase,actualFields{index},expectedFields{index},diagnostic+" "+names{index});
end
expectedFlux = flux(expected); actualFlux = flux(actual);
for index = 1:numel(expectedFlux)
    verifyNear(testCase,actualFlux{index},expectedFlux{index},diagnostic+" flux");
end
end

function values = flux(wvt)
if isQG(wvt)
    values = {wvt.nonlinearFlux()};
else
    values = cell(1,3);
    [values{:}] = wvt.nonlinearFlux();
end
end

function value = isQG(wvt)
value = isa(wvt,"WVTransformStratifiedQG") || isa(wvt,"WVTransformBarotropicQG");
end

function verifyNear(testCase,actual,expected,diagnostic)
testCase.verifySize(actual,size(expected),diagnostic);
testCase.assertTrue(all(isfinite(actual(:))) && all(isfinite(expected(:))),diagnostic);
error = norm(actual(:)-expected(:))/max(norm(expected(:)),realmin);
testCase.verifyLessThanOrEqual(error,2e-11,diagnostic);
end

function verifyWritableOpen(testCase,path)
file = NetCDFFile(path,shouldReadOnly=false);
cleanup = onCleanup(@()closeIfOpen(file));
testCase.verifyNotEmpty(file.id);
file.close();
clear cleanup
end

function closeIfOpen(file)
if ~isempty(file) && isvalid(file) && ~isempty(file.id), file.close(); end
end

function deleteTransforms(varargin)
for value = varargin
    if ~isempty(value{1}) && isvalid(value{1}), delete(value{1}); end
end
end
