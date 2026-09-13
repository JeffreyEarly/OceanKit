classdef TestWVCompiledConfiguration < matlab.unittest.TestCase
    % Immutability of scientific configuration owned by a compiled transform.

    methods (Test,TestTags="optional")
        function compiledConfigurationRejectsPublicMutation(testCase)
            capabilities = WVCompiledBackend.capabilities();
            testCase.assumeTrue(capabilities.isAvailable,capabilities.failure.message);

            transforms = configurations("compiled");
            cleanup = onCleanup(@()deleteTransforms(transforms{:})); %#ok<NASGU>
            verifyLockedProperties(testCase,transforms{1},boussinesqProperties());
            verifyLockedProperties(testCase,transforms{2},constantProperties());
            verifyLockedProperties(testCase,transforms{3},{'h'});

            verifyIndexedMutationRejected(testCase,transforms{1},"PFpm");
            verifyIndexedMutationRejected(testCase,transforms{1},"UAp");
            verifyIndexedMutationRejected(testCase,transforms{1},"UA0");
            verifyIndexedMutationRejected(testCase,transforms{1},"Apm_TE_factor");
            verifyIndexedMutationRejected(testCase,transforms{2},"DCT");
            verifyIndexedMutationRejected(testCase,transforms{2},"cos_alpha");
        end

        function rejectedMutationPreservesActiveEvaluation(testCase)
            capabilities = WVCompiledBackend.capabilities();
            testCase.assumeTrue(capabilities.isAvailable,capabilities.failure.message);

            transforms = configurations("compiled");
            cleanup = onCleanup(@()deleteTransforms(transforms{:})); %#ok<NASGU>
            wvt = transforms{1};
            scope = wvt.scopedEvaluation(); %#ok<NASGU>
            before = wvt.computationalBackendMetadata.runtimeMetrics;
            original = wvt.PF0;
            testCase.verifyError(@()assignProperty(wvt,"PF0",differentValue(original)),configurationLockedIdentifier());
            testCase.verifyEqual(wvt.PF0,original);
            testCase.verifyWarningFree(@()wvt.variableWithName('u'));
            after = wvt.computationalBackendMetadata.runtimeMetrics;
            testCase.verifyEqual(after.scopedEvaluations,before.scopedEvaluations);
        end

        function privateBackendsLeaseMatlabConfiguration(testCase)
            capabilities = WVCompiledBackend.capabilities();
            testCase.assumeTrue(capabilities.isAvailable,capabilities.failure.message);

            wvt = boussinesqConfiguration("matlab");
            cleanup = onCleanup(@()deleteTransforms(wvt)); %#ok<NASGU>
            original = wvt.PF0; replacement = differentValue(original);
            first = WVCompiledTransformBackend.create(wvt);
            firstCleanup = onCleanup(@()delete(first)); %#ok<NASGU>
            testCase.verifyError(@()assignProperty(wvt,"PF0",replacement),configurationLockedIdentifier());
            second = WVCompiledTransformBackend.create(wvt);
            secondCleanup = onCleanup(@()delete(second)); %#ok<NASGU>
            delete(first);
            testCase.verifyError(@()assignProperty(wvt,"PF0",replacement),configurationLockedIdentifier());
            delete(second);
            assignProperty(wvt,"PF0",replacement);
            testCase.verifyEqual(wvt.PF0,replacement);
        end

        function failedPrivateBackendCreationDoesNotLease(testCase)
            capabilities = WVCompiledBackend.capabilities();
            testCase.assumeTrue(capabilities.isAvailable,capabilities.failure.message);

            wvt = boussinesqConfiguration("matlab");
            cleanup = onCleanup(@()deleteTransforms(wvt)); %#ok<NASGU>
            original = wvt.PF0;
            wvt.PF0 = 1;
            didFail = false;
            try
                unexpected = WVCompiledTransformBackend.create(wvt);
                delete(unexpected);
            catch
                didFail = true;
            end
            testCase.verifyTrue(didFail);
            assignProperty(wvt,"PF0",original);
            testCase.verifyEqual(wvt.PF0,original);
        end

        function matlabConfigurationRemainsWritable(testCase)
            transforms = configurations("matlab");
            cleanup = onCleanup(@()deleteTransforms(transforms{:})); %#ok<NASGU>

            boussinesq = transforms{1};
            replacement = differentValue(boussinesq.PF0);
            assignProperty(boussinesq,"PF0",replacement);
            testCase.verifyEqual(boussinesq.PF0,replacement);
            replacement = differentValue(boussinesq.UAp);
            assignFirstElement(boussinesq,"UAp",replacement(1));
            testCase.verifyEqual(boussinesq.UAp(1),replacement(1));

            constant = transforms{2};
            replacement = differentValue(constant.cos_alpha);
            assignFirstElement(constant,"cos_alpha",replacement(1));
            testCase.verifyEqual(constant.cos_alpha(1),replacement(1));

            barotropic = transforms{3};
            replacement = differentValue(barotropic.h);
            assignProperty(barotropic,"h",replacement);
            testCase.verifyEqual(barotropic.h,replacement);
        end
    end
end

function transforms = configurations(backend)
transforms = { ...
    boussinesqConfiguration(backend), ...
    WVTransformConstantStratification([4000 3000 1000],[6 6 5],N0=5.2e-3,isHydrostatic=false,shouldAntialias=true,computationalBackend=backend), ...
    WVTransformBarotropicQG([4000 3000],[6 6],h=1000,j=1,shouldAntialias=true,computationalBackend=backend)};
end

function wvt = boussinesqConfiguration(backend)
profile = @(z)1e-4*exp(z/700);
wvt = WVTransformBoussinesq([4000 3000 1000],[6 6 5],Nj=3,N2Function=profile,shouldAntialias=true,computationalBackend=backend);
end

function names = boussinesqProperties()
names = { ...
    'planetaryRadius','rotationRate','latitude','g', ...
    'dLnN2','PF0inv','QG0inv','PF0','QG0','h_0','h_pm','P0','Q0', ...
    'K2unique','iK2unique','K2uniqueK2Map','PFpmInv','QGpmInv','PFpm','QGpm','QGwg','Ppm','Qpm', ...
    'UAp','VAp','WAp','NAp','UAm','VAm','WAm','NAm','ApmD','ApmN','iOmega', ...
    'UA0','VA0','NA0','PA0','A0Z','A0N', ...
    'Apm_TE_factor','A0_TE_factor','A0_TZ_factor','A0_QGPV_factor','A0_Psi_factor','A0_KE_factor','A0_PE_factor'};
end

function names = constantProperties()
names = { ...
    'planetaryRadius','rotationRate','latitude','g', ...
    'N0','isHydrostatic','h_0','h_pm','F_g','G_g','F_wg','G_wg','DCT','iDCT','DST','iDST', ...
    'UAp','VAp','WAp','NAp','UAm','VAm','WAm','NAm','ApmD','ApmN','iOmega', ...
    'UA0','VA0','NA0','PA0','A0Z','A0N', ...
    'Apm_TE_factor','A0_TE_factor','A0_TZ_factor','A0_QGPV_factor','A0_Psi_factor','A0_KE_factor','A0_PE_factor', ...
    'cos_alpha','sin_alpha','ApmD_scaled','ApmW_scaled'};
end

function verifyLockedProperties(testCase,wvt,names)
for index = 1:numel(names)
    name = string(names{index});
    original = wvt.(name);
    testCase.verifyError(@()assignProperty(wvt,name,differentValue(original)),configurationLockedIdentifier(),name);
    testCase.verifyEqual(wvt.(name),original,name);
end
end

function verifyIndexedMutationRejected(testCase,wvt,name)
original = wvt.(name);
replacement = differentValue(original);
testCase.verifyError(@()assignFirstElement(wvt,name,replacement(1)),configurationLockedIdentifier(),name);
testCase.verifyEqual(wvt.(name),original,name);
end

function assignProperty(wvt,name,value)
wvt.(name) = value;
end

function assignFirstElement(wvt,name,value)
wvt.(name)(1) = value;
end

function value = differentValue(value)
if iscell(value)
    if isempty(value)
        value = {1};
    else
        value{1} = differentValue(value{1});
    end
elseif isempty(value)
    value = 1;
elseif islogical(value)
    value(1) = ~value(1);
else
    value(1) = value(1) + 1;
end
end

function identifier = configurationLockedIdentifier()
identifier = "WaveVortexModel:CompiledTransformConfigurationLocked";
end

function deleteTransforms(varargin)
for index = 1:numel(varargin)
    transform = varargin{index};
    if ~isempty(transform) && isvalid(transform), delete(transform); end
end
end
