classdef TestPublicInterpolationContract < matlab.unittest.TestCase

    properties
        wvt
    end

    properties (TestParameter)
        invalidInterpolationMethod = {'exact','finufft','nearest'}
        supportedInterpolationMethod = {'linear','spline'}
    end

    methods (TestClassSetup)
        function classSetup(testCase)
            Lxyz = [1000 800 500];
            Nxyz = [8 8 7];
            N2 = @(z) (5.2e-3)^2*(1 + 0.2*z/Lxyz(3));
            testCase.wvt = WVTransformHydrostatic(Lxyz,Nxyz,N2=N2,shouldAntialias=false);
        end
    end

    methods (Test, TestTags = "full")
        function testTransformRejectsUnsupportedMethods(testCase,invalidInterpolationMethod)
            testCase.verifyError(@() testCase.wvt.variableAtPositionWithName(0,0,[],'u',interpolationMethod=invalidInterpolationMethod),'MATLAB:validators:mustBeMember');
        end

        function testLagrangianParticlesRejectUnsupportedAdvection(testCase,invalidInterpolationMethod)
            model = testCase.model();
            testCase.verifyError(@() WVLagrangianParticles(model,name='invalid',x=0,y=0,z=-100,trackedFieldNames={},advectionInterpolation=invalidInterpolationMethod),'MATLAB:validators:mustBeMember');
        end

        function testLagrangianParticlesRejectUnsupportedTrackedFields(testCase,invalidInterpolationMethod)
            model = testCase.model();
            testCase.verifyError(@() WVLagrangianParticles(model,name='invalid',x=0,y=0,z=-100,trackedFieldNames={},trackedVarInterpolation=invalidInterpolationMethod),'MATLAB:validators:mustBeMember');
        end

        function testAddParticlesRejectsUnsupportedMethods(testCase,invalidInterpolationMethod)
            model = testCase.model();
            testCase.verifyError(@() model.addParticles('invalid-advection',false,0,0,-100,advectionInterpolation=invalidInterpolationMethod),'MATLAB:validators:mustBeMember');
            testCase.verifyError(@() model.addParticles('invalid-tracked',false,0,0,-100,trackedVarInterpolation=invalidInterpolationMethod),'MATLAB:validators:mustBeMember');
        end

        function testSetFloatPositionsRejectsUnsupportedMethods(testCase,invalidInterpolationMethod)
            model = testCase.model();
            testCase.verifyError(@() model.setFloatPositions(0,0,-100,advectionInterpolation=invalidInterpolationMethod),'MATLAB:validators:mustBeMember');
            testCase.verifyError(@() model.setFloatPositions(0,0,-100,trackedVarInterpolation=invalidInterpolationMethod),'MATLAB:validators:mustBeMember');
        end

        function testSetDrifterPositionsRejectsUnsupportedMethods(testCase,invalidInterpolationMethod)
            model = testCase.model();
            testCase.verifyError(@() model.setDrifterPositions(0,0,[],advectionInterpolation=invalidInterpolationMethod),'MATLAB:validators:mustBeMember');
            testCase.verifyError(@() model.setDrifterPositions(0,0,[],trackedVarInterpolation=invalidInterpolationMethod),'MATLAB:validators:mustBeMember');
        end

        function testSupportedMethodsPropagateThroughParticleAPIs(testCase,supportedInterpolationMethod)
            directModel = testCase.model();
            direct = WVLagrangianParticles(directModel,name='direct',x=0,y=0,z=-100,trackedFieldNames={},advectionInterpolation=supportedInterpolationMethod,trackedVarInterpolation=supportedInterpolationMethod);
            testCase.verifyEqual(direct.advectionInterpolation,supportedInterpolationMethod);
            testCase.verifyEqual(direct.trackedVarInterpolation,supportedInterpolationMethod);

            particleModel = testCase.model();
            particleModel.addParticles('particles',false,0,0,-100,advectionInterpolation=supportedInterpolationMethod,trackedVarInterpolation=supportedInterpolationMethod);
            particles = particleModel.fluxedObservingSystemWithName('particles');
            testCase.verifyEqual(particles.advectionInterpolation,supportedInterpolationMethod);
            testCase.verifyEqual(particles.trackedVarInterpolation,supportedInterpolationMethod);

            floatModel = testCase.model();
            floatModel.setFloatPositions(0,0,-100,advectionInterpolation=supportedInterpolationMethod,trackedVarInterpolation=supportedInterpolationMethod);
            floats = floatModel.fluxedObservingSystemWithName('float');
            testCase.verifyEqual(floats.advectionInterpolation,supportedInterpolationMethod);
            testCase.verifyEqual(floats.trackedVarInterpolation,supportedInterpolationMethod);

            drifterModel = testCase.model();
            drifterModel.setDrifterPositions(0,0,[],advectionInterpolation=supportedInterpolationMethod,trackedVarInterpolation=supportedInterpolationMethod);
            drifters = drifterModel.fluxedObservingSystemWithName('drifter');
            testCase.verifyEqual(drifters.advectionInterpolation,supportedInterpolationMethod);
            testCase.verifyEqual(drifters.trackedVarInterpolation,supportedInterpolationMethod);
        end

        function testParticleInterpolationDefaults(testCase)
            directModel = testCase.model();
            direct = WVLagrangianParticles(directModel,name='direct',x=0,y=0,z=-100,trackedFieldNames={});
            testCase.verifyEqual(direct.advectionInterpolation,'linear');
            testCase.verifyEqual(direct.trackedVarInterpolation,'linear');

            particleModel = testCase.model();
            particleModel.addParticles('particles',false,0,0,-100);
            particles = particleModel.fluxedObservingSystemWithName('particles');
            testCase.verifyEqual(particles.advectionInterpolation,'linear');
            testCase.verifyEqual(particles.trackedVarInterpolation,'spline');

            floatModel = testCase.model();
            floatModel.setFloatPositions(0,0,-100);
            floats = floatModel.fluxedObservingSystemWithName('float');
            testCase.verifyEqual(floats.advectionInterpolation,'linear');
            testCase.verifyEqual(floats.trackedVarInterpolation,'linear');

            drifterModel = testCase.model();
            drifterModel.setDrifterPositions(0,0,[]);
            drifters = drifterModel.fluxedObservingSystemWithName('drifter');
            testCase.verifyEqual(drifters.advectionInterpolation,'linear');
            testCase.verifyEqual(drifters.trackedVarInterpolation,'linear');
        end
    end

    methods (Access=private)
        function model = model(testCase)
            model = WVModel(testCase.wvt,shouldUseLinearDynamics=true);
        end
    end

end
