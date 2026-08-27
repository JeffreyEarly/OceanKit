classdef TestWVModelOutputPersistence < matlab.unittest.TestCase
    properties (TestParameter)
        integratorType = struct(fixed="fixed",adaptive="adaptive")
    end

    properties
        tempFolder string
    end

    methods (TestMethodSetup)
        function createTemporaryFolder(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.tempFolder = string(fixture.Folder);
        end
    end

    methods (Test, TestTags="full")
        function outputRegistrationIsIdentityBasedAndAtomic(testCase)
            model = TestWVModelOutputPersistence.linearModel();
            otherModel = TestWVModelOutputPersistence.linearModel();
            path = fullfile(testCase.tempFolder,"registered.nc");
            outputFile = WVModelOutputFile(model,path);
            model.addOutputFile(outputFile);
            model.addOutputFile(outputFile);
            testCase.verifyTrue(model.outputFileWithName("registered.nc") == outputFile)

            duplicateFile = WVModelOutputFile(model,fullfile(testCase.tempFolder,"other","registered.nc"));
            foreignFile = WVModelOutputFile(otherModel,fullfile(testCase.tempFolder,"foreign.nc"));
            testCase.verifyError(@()model.addOutputFile(duplicateFile),'')
            testCase.verifyError(@()model.addOutputFile(foreignFile),'')
            testCase.verifyEqual(model.outputFileNames,"registered.nc")

            outputGroup = WVModelOutputGroupEvenlySpaced(model,name="group",outputInterval=1);
            outputFile.addOutputGroup(outputGroup);
            outputFile.addOutputGroup(outputGroup);
            testCase.verifyTrue(outputFile.outputGroupWithName("group") == outputGroup)
            duplicateGroup = WVModelOutputGroupEvenlySpaced(model,name="group",outputInterval=2);
            foreignGroup = WVModelOutputGroupEvenlySpaced(otherModel,name="foreign",outputInterval=1);
            testCase.verifyError(@()outputFile.addOutputGroup(duplicateGroup),'')
            testCase.verifyError(@()outputFile.addOutputGroup(foreignGroup),'')
            testCase.verifyEqual(outputFile.outputGroupNames,"group")
            testCase.verifyError(@()model.outputFileWithName("missing"),'')
            testCase.verifyError(@()outputFile.outputGroupWithName("missing"),'')
        end

        function schedulesRemainAnchoredAcrossFilesGroupsAndCalls(testCase)
            model = TestWVModelOutputPersistence.linearModel();
            pathA = fullfile(testCase.tempFolder,"schedule-a.nc");
            pathB = fullfile(testCase.tempFolder,"schedule-b.nc");

            fileA = model.addNewOutputFile(pathA,shouldOverwriteExisting=true);
            groupA = fileA.addNewEvenlySpacedOutputGroup("bounded",initialTime=0.2,finalTime=1.1,outputInterval=0.3);
            groupA.addObservingSystem(WVEulerianFields(model,fieldNames={'u'}));
            groupB = fileA.addNewEvenlySpacedOutputGroup("coincident",initialTime=0,finalTime=0.8,outputInterval=0.2);
            groupB.addObservingSystem(WVEulerianFields(model,fieldNames={'v'}));

            fileB = model.addNewOutputFile(pathB,shouldOverwriteExisting=true);
            groupC = fileB.addNewEvenlySpacedOutputGroup("second-file",initialTime=0.1,finalTime=1.1,outputInterval=0.5);
            groupC.addObservingSystem(WVEulerianFields(model,fieldNames={'u'}));

            model.integrateToTime(0.55,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            model.integrateToTime(1.1,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            model.closeNetCDFFile();

            testCase.verifyEqual(TestWVModelOutputPersistence.readVariable(pathA,'bounded/t'),0.2+(0:3).'*0.3)
            testCase.verifyEqual(TestWVModelOutputPersistence.readVariable(pathA,'coincident/t'),(0:4).'*0.2)
            testCase.verifyEqual(TestWVModelOutputPersistence.readVariable(pathB,'second-file/t'),0.1+(0:2).'*0.5)
            TestWVModelOutputPersistence.verifyTimeVector(testCase,TestWVModelOutputPersistence.readVariable(pathA,'bounded/t'),0.2,1.1)
            TestWVModelOutputPersistence.verifyTimeVector(testCase,TestWVModelOutputPersistence.readVariable(pathA,'coincident/t'),0,0.8)
            TestWVModelOutputPersistence.verifyTimeVector(testCase,TestWVModelOutputPersistence.readVariable(pathB,'second-file/t'),0.1,1.1)
        end

        function stableObserversRoundTripWithSharedIntegratedHandles(testCase)
            path = fullfile(testCase.tempFolder,"all-observers.nc");
            model = TestWVModelOutputPersistence.nonlinearModel();
            outputFile = model.createNetCDFFileForModelOutput(path,outputInterval=1,shouldOverwriteExisting=true);
            model.eulerianObservingSystem.addNetCDFOutputVariables('u');
            model.setFloatPositions([500 1500],[400 1200],[-250 -750],'u',advectionInterpolation="spline",trackedVarInterpolation="linear",absToleranceXY=3e-4,absToleranceZ=7e-5);
            model.addTracer(reshape(1:prod(model.wvt.spatialMatrixSize),model.wvt.spatialMatrixSize),"dye");
            defaultGroup = outputFile.outputGroupWithName(model.defaultOutputGroupName());
            mooring = WVMooring(model,name="mooring",x=[0 1000],y=[0 900],trackedFieldNames={'u'});
            defaultGroup.addObservingSystem(mooring);
            particles = model.fluxedObservingSystemWithName("float");
            tracer = model.fluxedObservingSystemWithName("dye");
            sharedGroup = outputFile.addNewEvenlySpacedOutputGroup("shared",outputInterval=1);
            sharedGroup.addObservingSystem([particles tracer]);

            model.setupIntegrator(integratorType="fixed",deltaT=1);
            model.wvCoefficientFluxedObservingSystem().absTolerance = 3e-8;
            model.integrateToTime(2,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            expectedU = model.wvt.u;
            expectedMooringU = zeros(model.wvt.Nz,length(mooring.x));
            for iMooring = 1:length(mooring.x)
                expectedMooringU(:,iMooring) = expectedU(mooring.x_index(iMooring),mooring.y_index(iMooring),:);
            end
            expectedParticles = particles.initialConditions();
            expectedTracer = tracer.phi;
            model.closeNetCDFFile();

            storedU = TestWVModelOutputPersistence.readVariable(path,'wave-vortex/u');
            storedMooringU = TestWVModelOutputPersistence.readVariable(path,'wave-vortex/mooring_u');
            testCase.verifyEqual(storedU(:,:,:,end),expectedU)
            testCase.verifyEqual(storedMooringU(:,:,end),expectedMooringU)

            TestWVModelOutputPersistence.verifySchema(testCase,path);
            restored = WVModel.modelFromFile(char(path));
            cleanup = onCleanup(@()restored.closeNetCDFFile());

            testCase.verifyFalse(restored.isDynamicsLinear)
            testCase.verifyEqual(restored.wvt.t,2)
            testCase.verifyClass(restored.wvt.forcingWithName("adaptive damping"),"WVAdaptiveDamping")
            testCase.verifyEqual(sort(restored.outputFileNames),"all-observers.nc")
            restoredFile = restored.outputFileWithName("all-observers.nc");
            restoredGroupNames = restoredFile.outputGroupNames();
            testCase.verifyEqual(sort(restoredGroupNames(:)),sort(["wave-vortex"; "shared"]))
            restoredDefault = restoredFile.outputGroupWithName("wave-vortex");
            restoredShared = restoredFile.outputGroupWithName("shared");
            testCase.verifyEqual(restoredDefault.outputInterval,1)
            testCase.verifyEqual(restoredDefault.initialTime,0)
            testCase.verifyEqual(restoredDefault.finalTime,Inf)
            testCase.verifyEqual(restoredDefault.incrementsWrittenToGroup,uint64(3))
            testCase.verifyEqual(restoredDefault.timeOfLastIncrementWrittenToGroup,2)
            testCase.verifyEqual(restoredShared.incrementsWrittenToGroup,uint64(3))

            restoredEulerian = restoredDefault.observingSystemWithName("eulerian fields");
            testCase.verifyEqual(sort(restoredEulerian.fieldNames),sort(["A0" "Am" "Ap" "u"]))
            restoredCoefficients = restoredDefault.observingSystemWithName("wave-vortex coefficient flux");
            testCase.verifyTrue(restoredCoefficients == restored.wvCoefficientFluxedObservingSystem())
            testCase.verifyEqual(restoredCoefficients.absTolerance,3e-8)

            restoredParticles = restoredDefault.observingSystemWithName("float");
            sharedParticles = restoredShared.observingSystemWithName("float");
            testCase.verifyTrue(restoredParticles == sharedParticles)
            testCase.verifyEqual(restoredParticles.initialConditions(),expectedParticles)
            testCase.verifyEqual(restoredParticles.advectionInterpolation,'spline')
            testCase.verifyEqual(restoredParticles.trackedVarInterpolation,'linear')
            testCase.verifyEqual(restoredParticles.absToleranceXY,3e-4)
            testCase.verifyEqual(restoredParticles.absToleranceZ,7e-5)
            testCase.verifyEqual(restoredParticles.trackedFieldNames,"u")

            restoredTracer = restoredDefault.observingSystemWithName("dye");
            testCase.verifyTrue(restoredTracer == restoredShared.observingSystemWithName("dye"))
            testCase.verifyEqual(restoredTracer.phi,expectedTracer)
            restoredMooring = restoredDefault.observingSystemWithName("mooring");
            testCase.verifyEqual(restoredMooring.x,[0 1000])
            testCase.verifyEqual(restoredMooring.y,[0 900])
            testCase.verifyEqual(restoredMooring.x_index,[1 3])
            testCase.verifyEqual(restoredMooring.y_index,[1 2])
            testCase.verifyEqual(restoredMooring.trackedFieldNames,"u")

            restored.closeNetCDFFile();
            clear cleanup
            TestWVModelOutputPersistence.verifyWritableOpen(path);
        end

        function barotropicParticlesAndTracerRoundTrip(testCase)
            path = fullfile(testCase.tempFolder,"barotropic.nc");
            wvt = WVTransformBarotropicQG([4000 3000],[8 6],latitude=45,shouldAntialias=false);
            model = WVModel(wvt,shouldUseLinearDynamics=true);
            model.createNetCDFFileForModelOutput(path,outputInterval=1,shouldOverwriteExisting=true);
            model.setDrifterPositions([0 100],[0 200],[],'u',absToleranceXY=9e-4);
            phi = reshape(1:prod(wvt.spatialMatrixSize),wvt.spatialMatrixSize);
            model.addTracer(phi,"dye");
            model.setupIntegrator(integratorType="fixed",deltaT=1);
            model.integrateToTime(1,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            model.closeNetCDFFile();

            restored = WVModel.modelFromFile(char(path));
            cleanup = onCleanup(@()restored.closeNetCDFFile());
            testCase.verifyTrue(restored.isDynamicsLinear)
            testCase.verifyFalse(any(arrayfun(@(observer)isa(observer,'WVCoefficients'),restored.fluxedObservingSystems)))
            drifter = restored.fluxedObservingSystemWithName("drifter");
            testCase.verifyTrue(drifter.isXYOnly)
            testCase.verifyEmpty(drifter.z)
            testCase.verifyEqual(drifter.x,[0 100])
            testCase.verifyEqual(drifter.y,[0 200])
            testCase.verifyEqual(drifter.absToleranceXY,9e-4)
            testCase.verifyEqual(restored.fluxedObservingSystemWithName("dye").phi,phi)
            restored.closeNetCDFFile();
            clear cleanup
        end

        function eachOutputFileRestoresOnlyItsOwnGraph(testCase)
            model = TestWVModelOutputPersistence.nonlinearModel();
            pathA = fullfile(testCase.tempFolder,"restart-a.nc");
            pathB = fullfile(testCase.tempFolder,"restart-b.nc");
            model.createNetCDFFileForModelOutput(pathA,outputInterval=1,shouldOverwriteExisting=true);
            model.createNetCDFFileForModelOutput(pathB,outputInterval=1,shouldOverwriteExisting=true);
            model.setupIntegrator(integratorType="fixed",deltaT=1);
            model.integrateToTime(1,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            model.closeNetCDFFile();

            restoredA = WVModel.modelFromFile(char(pathA));
            cleanupA = onCleanup(@()restoredA.closeNetCDFFile());
            restoredB = WVModel.modelFromFile(char(pathB));
            cleanupB = onCleanup(@()restoredB.closeNetCDFFile());
            testCase.verifyEqual(restoredA.outputFileNames,"restart-a.nc")
            testCase.verifyEqual(restoredB.outputFileNames,"restart-b.nc")
            testCase.verifyEqual(restoredA.wvt.t,1)
            testCase.verifyEqual(restoredB.wvt.t,1)
            restoredA.closeNetCDFFile();
            restoredB.closeNetCDFFile();
            clear cleanupA cleanupB
        end

        function restartedObserversMatchUninterruptedControl(testCase,integratorType)
            checkpointTime = 300;
            finalTime = 600;
            uninterrupted = TestWVModelOutputPersistence.observerIntegrationModel();
            checkpoint = TestWVModelOutputPersistence.observerIntegrationModel();
            TestWVModelOutputPersistence.configureIntegrator(uninterrupted,integratorType);
            TestWVModelOutputPersistence.configureIntegrator(checkpoint,integratorType);
            path = fullfile(testCase.tempFolder,"observer-restart-"+integratorType+".nc");
            checkpoint.createNetCDFFileForModelOutput(path,outputInterval=checkpointTime,shouldOverwriteExisting=true);

            uninterrupted.integrateToTime(finalTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            checkpoint.integrateToTime(checkpointTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            checkpoint.closeNetCDFFile();

            resumed = WVModel.modelFromFile(char(path));
            cleanup = onCleanup(@()resumed.closeNetCDFFile());
            testCase.verifyTrue(resumed.isDynamicsLinear)
            TestWVModelOutputPersistence.configureIntegrator(resumed,integratorType);
            resumed.integrateToTime(finalTime,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);

            [actualX,actualY,actualZ,actualTracked] = resumed.floatPositions();
            [expectedX,expectedY,expectedZ,expectedTracked] = uninterrupted.floatPositions();
            testCase.verifyEqual(resumed.wvt.t,finalTime)
            testCase.verifyEqual(actualX,expectedX,AbsTol=2e-6)
            testCase.verifyEqual(actualY,expectedY,AbsTol=2e-6)
            testCase.verifyEqual(actualZ,expectedZ,AbsTol=2e-10)
            testCase.verifyEqual(actualTracked.u,expectedTracked.u,AbsTol=2e-10)
            testCase.verifyEqual(resumed.tracer("dye"),uninterrupted.tracer("dye"),AbsTol=5e-6)
            testCase.verifyEqual(resumed.wvt.totalEnergy,uninterrupted.wvt.totalEnergy,RelTol=1e-12)
            resumed.closeNetCDFFile();
            clear cleanup
            testCase.verifyEqual(TestWVModelOutputPersistence.readVariable(path,'wave-vortex/t'),[0; checkpointTime; finalTime])
        end
    end

    methods (Static, Access=private)
        function model = linearModel()
            wvt = WVTransformConstantStratification([4000 3000 1000],[8 6 5],N0=sqrt(2e-5),latitude=45,isHydrostatic=false,shouldAntialias=false);
            wvt.initWithInertialMotions(@(z)0.1*ones(size(z)),@(z)zeros(size(z)));
            model = WVModel(wvt,shouldUseLinearDynamics=true);
        end

        function model = nonlinearModel()
            wvt = WVTransformConstantStratification([4000 3000 1000],[8 6 5],N0=sqrt(2e-5),latitude=45,isHydrostatic=false,shouldAntialias=false);
            wvt.addForcing(WVAdaptiveDamping(wvt));
            model = WVModel(wvt);
        end

        function model = observerIntegrationModel()
            model = TestWVModelOutputPersistence.linearModel();
            model.setFloatPositions([500 1500],[400 1200],[-250 -750],'u',absToleranceXY=1e-8,absToleranceZ=1e-8);
            phi = sin(2*pi*model.wvt.X/model.wvt.Lx);
            model.addTracer(phi,"dye");
        end

        function configureIntegrator(model,integratorType)
            if integratorType == "fixed"
                model.setupIntegrator(integratorType="fixed",deltaT=10);
            else
                model.setupIntegrator(integratorType="adaptive",relTolerance=1e-10);
            end
        end

        function value = readVariable(path,name)
            ncfile = NetCDFFile(path,shouldReadOnly=true);
            cleanup = onCleanup(@()TestWVModelOutputPersistence.closeIfOpen(ncfile));
            value = ncfile.readVariables(name);
            ncfile.close();
            clear cleanup
        end

        function verifyTimeVector(testCase,t,initialTime,finalTime)
            testCase.verifyEqual(t,sort(t))
            testCase.verifyEqual(length(t),length(unique(t)))
            testCase.verifyGreaterThanOrEqual(t,initialTime)
            testCase.verifyLessThanOrEqual(t,finalTime)
        end

        function verifySchema(testCase,path)
            ncfile = NetCDFFile(path,shouldReadOnly=true);
            cleanup = onCleanup(@()TestWVModelOutputPersistence.closeIfOpen(ncfile));
            testCase.verifyEqual(ncfile.attributes('WVModelIsDynamicsLinear'),uint8(0))
            group = ncfile.groupWithName('wave-vortex');
            testCase.verifyEqual(group.attributes('AnnotatedClass'),'WVModelOutputGroupEvenlySpaced')
            timeVariable = group.variableWithName('t');
            testCase.verifyEqual(timeVariable.attributes('axis'),'T')
            testCase.verifyEqual(timeVariable.attributes('calendar'),'standard')
            testCase.verifyEqual(timeVariable.attributes('units'),'seconds since 1970-01-01 00:00:00')
            testCase.verifyTrue(all(group.hasVariableWithName('Ap','Am','A0','u','float_x','float_y','float_z','float_u','dye','mooring_x','mooring_y','mooring_u')))
            testCase.verifyEqual(group.dimensionWithName('t').nPoints,3)
            ncfile.close();
            clear cleanup
        end

        function verifyWritableOpen(path)
            ncfile = NetCDFFile(path,shouldReadOnly=false);
            cleanup = onCleanup(@()TestWVModelOutputPersistence.closeIfOpen(ncfile));
            ncfile.close();
            clear cleanup
        end

        function closeIfOpen(ncfile)
            if ~isempty(ncfile) && isvalid(ncfile) && ~isempty(ncfile.id)
                ncfile.close();
            end
        end
    end
end
