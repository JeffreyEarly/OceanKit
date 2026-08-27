classdef TestObservingSystems < matlab.unittest.TestCase
    properties (TestParameter)
        integratorType = struct(fixed="fixed",adaptive="adaptive")
    end

    methods (Test, TestTags = "full")
        function modelRegistryUsesIdentityAndAtomicMutations(testCase)
            model = TestObservingSystems.linearModel3D();
            particles = TestObservingSystems.particles(model,"particles",0);
            model.addFluxedObservingSystem(particles);

            model.addFluxedObservingSystem(particles);
            testCase.verifyEqual(numel(model.fluxedObservingSystems),1)
            testCase.verifyTrue(model.fluxedObservingSystems(1) == particles)

            duplicate = TestObservingSystems.particles(model,"particles",100);
            testCase.verifyError(@()model.addFluxedObservingSystem(duplicate),'')
            TestObservingSystems.verifyRegistry(testCase,model,particles,{1:3})

            otherModel = TestObservingSystems.linearModel3D();
            foreign = TestObservingSystems.particles(otherModel,"foreign",0);
            testCase.verifyError(@()model.addFluxedObservingSystem(foreign),'')
            testCase.verifyError(@()model.removeFluxedObservingSystem(foreign),'')
            TestObservingSystems.verifyRegistry(testCase,model,particles,{1:3})

            testCase.verifyError(@()model.removeFluxedObservingSystem(duplicate),'')
            testCase.verifyError(@()model.fluxedObservingSystemWithName("missing"),'')
            TestObservingSystems.verifyRegistry(testCase,model,particles,{1:3})

            model.removeFluxedObservingSystem(particles);
            testCase.verifyEmpty(model.fluxedObservingSystems)
            testCase.verifyEqual(model.nFluxComponents,0)
            testCase.verifyEmpty(model.indicesForFluxedSystem)
        end

        function coefficientsRemainFirstAndIndicesRemainContiguous(testCase)
            model = TestObservingSystems.linearModel3D();
            particles = TestObservingSystems.particles(model,"particles",0);
            tracer = WVTracer(model,name="tracer",phi=zeros(model.wvt.spatialMatrixSize));
            coefficients = WVCoefficients(model);

            model.addFluxedObservingSystem(particles);
            model.addFluxedObservingSystem(tracer);
            model.addFluxedObservingSystem(coefficients);
            TestObservingSystems.verifyRegistry(testCase,model,[coefficients particles tracer],{1:3,4:6,7})

            duplicateCoefficients = WVCoefficients(model);
            testCase.verifyError(@()model.addFluxedCoefficients(duplicateCoefficients),'')
            TestObservingSystems.verifyRegistry(testCase,model,[coefficients particles tracer],{1:3,4:6,7})

            model.removeFluxedObservingSystem(particles);
            TestObservingSystems.verifyRegistry(testCase,model,[coefficients tracer],{1:3,4})

            replacementParticles = TestObservingSystems.particles(model,"replacement",100);
            model.addFluxedObservingSystem(replacementParticles);
            TestObservingSystems.verifyRegistry(testCase,model,[coefficients tracer replacementParticles],{1:3,4,5:7})
        end

        function outputGroupLifecycleIsIdentityBasedAndAtomic(testCase)
            model = TestObservingSystems.linearModel3D();
            outputGroup = WVModelOutputGroup(model,name="memory-only");
            eulerian = model.eulerianObservingSystem;
            particles = TestObservingSystems.particles(model,"particles",0);
            tracer = WVTracer(model,name="tracer",phi=zeros(model.wvt.spatialMatrixSize));

            outputGroup.addObservingSystem(eulerian);
            outputGroup.addObservingSystem(particles);
            outputGroup.addObservingSystem(tracer);
            testCase.verifyTrue(outputGroup.observingSystemWithName("particles") == particles)
            TestObservingSystems.verifyRegistry(testCase,model,[particles tracer],{1:3,4})

            duplicate = TestObservingSystems.particles(model,"particles",100);
            testCase.verifyError(@()outputGroup.addObservingSystem(duplicate),'')
            testCase.verifyError(@()outputGroup.removeObservingSystem(duplicate),'')
            testCase.verifyEqual(numel(outputGroup.observingSystems),3)
            TestObservingSystems.verifyRegistry(testCase,model,[particles tracer],{1:3,4})

            otherModel = TestObservingSystems.linearModel3D();
            foreign = TestObservingSystems.particles(otherModel,"foreign",0);
            testCase.verifyError(@()outputGroup.addObservingSystem(foreign),'')
            testCase.verifyError(@()outputGroup.removeObservingSystem(foreign),'')
            testCase.verifyError(@()outputGroup.observingSystemWithName("missing"),'')
            testCase.verifyEqual(numel(outputGroup.observingSystems),3)
            TestObservingSystems.verifyRegistry(testCase,model,[particles tracer],{1:3,4})

            outputGroup.removeObservingSystem(particles);
            testCase.verifyEqual(string({outputGroup.observingSystems.name}),["eulerian fields" "tracer"])
            TestObservingSystems.verifyRegistry(testCase,model,tracer,{1})

            outputGroup.removeObservingSystem(eulerian);
            testCase.verifyTrue(outputGroup.observingSystems(1) == tracer)
            TestObservingSystems.verifyRegistry(testCase,model,tracer,{1})
        end

        function outputGroupUsesCanonicalCoefficientObserver(testCase)
            model = TestObservingSystems.linearModel3D();
            canonicalCoefficients = WVCoefficients(model,absTolerance=1e-6);
            model.addFluxedCoefficients(canonicalCoefficients);
            restoredCoefficients = WVCoefficients(model,absTolerance=4e-8);
            outputGroup = WVModelOutputGroup(model,name="memory-only");

            outputGroup.addObservingSystem(restoredCoefficients);

            testCase.verifyTrue(outputGroup.observingSystems(1) == canonicalCoefficients)
            testCase.verifyTrue(model.fluxedObservingSystems(1) == canonicalCoefficients)
            testCase.verifyEqual(canonicalCoefficients.absTolerance,4e-8)
            TestObservingSystems.verifyRegistry(testCase,model,canonicalCoefficients,{1:3})
        end

        function suppliedEulerianAndCoefficientContracts(testCase)
            model3D = TestObservingSystems.linearModel3D();
            eulerian = model3D.eulerianObservingSystem;
            eulerian.setNetCDFOutputVariables('u','A0');
            testCase.verifyEqual(eulerian.fieldNames,["u" "A0"])
            testCase.verifyEqual(eulerian.nOutputVariables,2)
            eulerian.addNetCDFOutputVariables('v','u');
            testCase.verifyEqual(sort(eulerian.fieldNames),sort(["u" "v" "A0"]))
            eulerian.removeNetCDFOutputVariables('A0');
            testCase.verifyFalse(any(eulerian.fieldNames == "A0"))
            fieldsBeforeFailure = eulerian.fieldNames;
            testCase.verifyError(@()eulerian.addNetCDFOutputVariables('not-a-field'),'')
            testCase.verifyEqual(eulerian.fieldNames,fieldsBeforeFailure)

            coefficients3D = WVCoefficients(model3D,absTolerance=2e-7);
            testCase.verifyEqual(double(coefficients3D.nFluxComponents),3)
            testCase.verifyEqual(coefficients3D.lengthOfFluxComponents(),repmat(numel(model3D.wvt.Ap),3,1))
            initial3D = coefficients3D.initialConditions();
            testCase.verifyEqual(initial3D,{model3D.wvt.Ap;model3D.wvt.Am;model3D.wvt.A0})
            tolerances3D = coefficients3D.absErrorTolerance();
            testCase.verifyEqual(cellfun(@numel,tolerances3D),repmat(numel(model3D.wvt.Ap),3,1))
            testCase.verifyTrue(all(cellfun(@(value)all(isfinite(value),"all"),tolerances3D)))

            updated3D = {ones(size(model3D.wvt.Ap));2*ones(size(model3D.wvt.Am));3*ones(size(model3D.wvt.A0))};
            coefficients3D.updateIntegratorValues(17,updated3D);
            testCase.verifyEqual(model3D.wvt.t,17)
            testCase.verifyEqual(model3D.wvt.Ap,updated3D{1})
            testCase.verifyEqual(model3D.wvt.Am,updated3D{2})
            testCase.verifyEqual(model3D.wvt.A0,updated3D{3})

            barotropicModel = TestObservingSystems.linearBarotropicModel();
            testCase.verifyEqual(barotropicModel.eulerianObservingSystem.fieldNames,"A0")
            coefficients2D = WVCoefficients(barotropicModel,absTolerance=2e-7);
            testCase.verifyEqual(double(coefficients2D.nFluxComponents),1)
            testCase.verifyEqual(coefficients2D.lengthOfFluxComponents(),numel(barotropicModel.wvt.A0))
            testCase.verifyEqual(coefficients2D.initialConditions(),{barotropicModel.wvt.A0})
            updatedA0 = {4*ones(size(barotropicModel.wvt.A0))};
            coefficients2D.updateIntegratorValues(19,updatedA0);
            testCase.verifyEqual(barotropicModel.wvt.t,19)
            testCase.verifyEqual(barotropicModel.wvt.A0,updatedA0{1})
        end

        function particleAndTracerFacadesExposeCurrentState(testCase)
            model3D = TestObservingSystems.linearModel3D();
            model3D.wvt.initWithInertialMotions(@(z)0.1*ones(size(z)),@(z)zeros(size(z)));
            x0 = [500 1500];
            y0 = [400 1200];
            z0 = [-250 -750];
            model3D.setFloatPositions(x0,y0,z0,'u',absToleranceXY=3e-4,absToleranceZ=7e-5);
            floats = model3D.fluxedObservingSystemWithName("float");
            testCase.verifyEqual(floats.absToleranceXY,3e-4)
            testCase.verifyEqual(floats.absToleranceZ,7e-5)
            testCase.verifyEqual(floats.lengthOfFluxComponents(),[2;2;2])
            testCase.verifyEqual(floats.initialConditions(),{x0,y0,z0})
            testCase.verifyEqual(floats.absErrorTolerance(),{3e-4,3e-4,7e-5})

            floatFlux = floats.fluxAtTime(0,{x0,y0,z0});
            testCase.verifyEqual(floatFlux{1},0.1*ones(size(x0)),AbsTol=1e-12)
            testCase.verifyEqual(floatFlux{2},zeros(size(y0)),AbsTol=1e-12)
            testCase.verifyEqual(floatFlux{3},zeros(size(z0)),AbsTol=1e-12)

            x1 = x0 + 20;
            y1 = y0 - 10;
            floats.updateIntegratorValues(20,{x1,y1,z0});
            model3D.wvt.t = 20;
            [x,y,z,tracked] = model3D.floatPositions();
            testCase.verifyEqual(x,x1)
            testCase.verifyEqual(y,y1)
            testCase.verifyEqual(z,z0)
            expectedU = model3D.wvt.variableAtPositionWithName(x1,y1,z0,'u');
            testCase.verifyEqual(tracked.u,expectedU)

            barotropicModel = TestObservingSystems.linearBarotropicModel();
            barotropicModel.setDrifterPositions([0 100],[0 200],[],'u',absToleranceXY=9e-4);
            drifters = barotropicModel.fluxedObservingSystemWithName("drifter");
            testCase.verifyTrue(drifters.isXYOnly)
            testCase.verifyEqual(drifters.absToleranceXY,9e-4)
            testCase.verifyEqual(drifters.lengthOfFluxComponents(),[2;2])
            testCase.verifyEqual(drifters.absErrorTolerance(),{9e-4,9e-4})

            phi0 = reshape(1:prod(barotropicModel.wvt.spatialMatrixSize),barotropicModel.wvt.spatialMatrixSize);
            barotropicModel.addTracer(phi0,"dye");
            tracer = barotropicModel.fluxedObservingSystemWithName("dye");
            testCase.verifyTrue(tracer.isXYOnly)
            testCase.verifyEqual(tracer.lengthOfFluxComponents(),numel(phi0))
            testCase.verifyEqual(tracer.initialConditions(),{phi0})
            testCase.verifyEqual(tracer.absErrorTolerance(),{tracer.absTolerance})
            tracer.updateIntegratorValues(0,{2*phi0});
            testCase.verifyEqual(barotropicModel.tracer("dye"),2*phi0)
        end

        function mooringsUsePeriodicOneBasedLowerCellIndices(testCase)
            model = TestObservingSystems.linearModel3D();
            dx = model.wvt.x(2)-model.wvt.x(1);
            dy = model.wvt.y(2)-model.wvt.y(1);
            x = [0 dx/2 dx model.wvt.Lx-dx/2 model.wvt.Lx model.wvt.Lx+dx -dx/2];
            y = [0 dy/2 dy model.wvt.Ly-dy/2 model.wvt.Ly -dy model.wvt.Ly+dy];
            mooring = WVMooring(model,name="test mooring",x=x,y=y,trackedFieldNames={'u'});

            expectedX = mod(x,model.wvt.Lx);
            expectedY = mod(y,model.wvt.Ly);
            testCase.verifyEqual(mooring.x,expectedX)
            testCase.verifyEqual(mooring.y,expectedY)
            testCase.verifyEqual(mooring.x_index,floor(expectedX/dx)+1)
            testCase.verifyEqual(mooring.y_index,floor(expectedY/dy)+1)
            testCase.verifyGreaterThanOrEqual(mooring.x_index,1)
            testCase.verifyLessThanOrEqual(mooring.x_index,model.wvt.Nx)
            testCase.verifyGreaterThanOrEqual(mooring.y_index,1)
            testCase.verifyLessThanOrEqual(mooring.y_index,model.wvt.Ny)
            testCase.verifyEqual(mooring.x_index([1 5]),[1 1])
            testCase.verifyEqual(mooring.x_index([3 6]),[2 2])
            testCase.verifyEqual(mooring.y_index([1 5]),[1 1])
            testCase.verifyEqual(mooring.y_index([3 7]),[2 2])
        end

        function mooringsRejectUnsupportedGeometryAndFields(testCase)
            model3D = TestObservingSystems.linearModel3D();
            testCase.verifyError(@()WVMooring(model3D,x=0,y=0,trackedFieldNames={'not-a-field'}),'')

            barotropicModel = TestObservingSystems.linearBarotropicModel();
            testCase.verifyError(@()WVMooring(barotropicModel,x=0,y=0,trackedFieldNames={'u'}),'')
        end

        function stableIntegratorsAdvanceParticlesAndTracerAnalytically(testCase,integratorType)
            [model,x0,y0,z0,phi0,U] = TestObservingSystems.observerIntegrationModel();
            finalTime = 600;
            if integratorType == "fixed"
                model.setupIntegrator(integratorType="fixed",deltaT=10);
            else
                model.setupIntegrator(integratorType="adaptive",relTolerance=1e-10);
            end

            model.integrateToTime(finalTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            [x,y,z,tracked] = model.floatPositions();
            phi = model.tracer("dye");

            omega = abs(model.wvt.f);
            deltaX = U*sin(omega*finalTime)/omega;
            deltaY = U*(cos(omega*finalTime)-1)/omega;
            expectedPhi = sin(2*pi*(model.wvt.X-deltaX)/model.wvt.Lx);

            testCase.verifyEqual(model.wvt.t,finalTime)
            testCase.verifyEqual(x,x0+deltaX,AbsTol=2e-5)
            testCase.verifyEqual(y,y0+deltaY,AbsTol=2e-5)
            testCase.verifyEqual(z,z0,AbsTol=1e-12)
            testCase.verifyNotEqual(x,x0)
            testCase.verifyNotEqual(y,y0)
            testCase.verifyEqual(phi,expectedPhi,AbsTol=5e-5)
            testCase.verifyNotEqual(phi,phi0)
            testCase.verifyEqual(tracked.u,U*cos(omega*finalTime)*ones(size(x0)),AbsTol=1e-10)
        end
    end

    methods (Static, Access = private)
        function model = linearModel3D()
            wvt = WVTransformConstantStratification([4000 3000 1000],[8 6 5],N0=sqrt(2e-5),latitude=45,isHydrostatic=false,shouldAntialias=false);
            model = WVModel(wvt,shouldUseLinearDynamics=true);
        end

        function model = linearBarotropicModel()
            wvt = WVTransformBarotropicQG([4000 3000],[8 6],latitude=45,shouldAntialias=false);
            model = WVModel(wvt,shouldUseLinearDynamics=true);
        end

        function particles = particles(model,name,xOffset)
            particles = WVLagrangianParticles(model,name=name,x=xOffset,y=0,z=-500,trackedFieldNames={});
        end

        function verifyRegistry(testCase,model,expectedSystems,expectedIndices)
            testCase.verifyEqual(numel(model.fluxedObservingSystems),numel(expectedSystems))
            testCase.verifyEqual(numel(model.indicesForFluxedSystem),numel(expectedIndices))
            expectedNFlux = 0;
            for iSystem = 1:numel(expectedSystems)
                testCase.verifyTrue(model.fluxedObservingSystems(iSystem) == expectedSystems(iSystem))
                testCase.verifyEqual(double(model.indicesForFluxedSystem{iSystem}),expectedIndices{iSystem}(:))
                expectedNFlux = expectedNFlux + double(expectedSystems(iSystem).nFluxComponents);
            end
            testCase.verifyEqual(double(model.nFluxComponents),expectedNFlux)
        end

        function [model,x0,y0,z0,phi0,U] = observerIntegrationModel()
            model = TestObservingSystems.linearModel3D();
            U = 0.1;
            model.wvt.initWithInertialMotions(@(z)U*ones(size(z)),@(z)zeros(size(z)));
            x0 = [500 1500];
            y0 = [400 1200];
            z0 = [-250 -750];
            model.setFloatPositions(x0,y0,z0,'u',absToleranceXY=1e-8,absToleranceZ=1e-8);
            phi0 = sin(2*pi*model.wvt.X/model.wvt.Lx);
            model.addTracer(phi0,"dye");
        end
    end
end
