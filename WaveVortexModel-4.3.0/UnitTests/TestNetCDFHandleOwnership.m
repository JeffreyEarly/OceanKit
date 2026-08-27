classdef TestNetCDFHandleOwnership < matlab.unittest.TestCase
    properties
        tempFolder
        transformPath
        modelPath
    end

    methods (TestMethodSetup)
        function createFixtures(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.tempFolder = string(fixture.Folder);
            testCase.transformPath = fullfile(testCase.tempFolder,"transform.nc");
            testCase.modelPath = fullfile(testCase.tempFolder,"model.nc");

            wvt = TestNetCDFHandleOwnership.newTransform();
            ncfile = wvt.writeToFile(testCase.transformPath,shouldOverwriteExisting=true);
            ncfile.close();
        end
    end

    methods (Test, TestTags = "full")
        function oneOutputRestorationClosesFile(testCase)
            wvt = WVTransform.waveVortexTransformFromFile(char(testCase.transformPath));
            testCase.verifyClass(wvt,"WVTransformConstantStratification")

            TestNetCDFHandleOwnership.verifyWritableOpen(testCase,testCase.transformPath);
        end

        function twoOutputRestorationReturnsCallerOwnedReader(testCase)
            [wvt,ncfile] = WVTransform.waveVortexTransformFromFile(char(testCase.transformPath));
            cleanup = onCleanup(@()TestNetCDFHandleOwnership.closeIfOpen(ncfile));

            testCase.verifyClass(wvt,"WVTransformConstantStratification")
            testCase.verifyNotEmpty(ncfile.id)
            testCase.verifyError(@()ncfile.addAttribute("unexpected-write",1),"MATLAB:imagesci:netcdf:libraryFailure")

            ncfile.close();
            TestNetCDFHandleOwnership.verifyWritableOpen(testCase,testCase.transformPath);
            clear cleanup
        end

        function explicitReadWriteRestorationRemainsAvailable(testCase)
            [~,ncfile] = WVTransform.waveVortexTransformFromFile(char(testCase.transformPath),shouldReadOnly=false);
            cleanup = onCleanup(@()TestNetCDFHandleOwnership.closeIfOpen(ncfile));

            testCase.verifyWarningFree(@()ncfile.addAttribute("writable-restoration",1))
            ncfile.close();
            testCase.verifyEqual(ncreadatt(testCase.transformPath,"/","writable-restoration"),1)
            clear cleanup
        end

        function directSubclassOneOutputRestorationClosesFile(testCase)
            [transforms,classNames] = TestNetCDFHandleOwnership.supportedTransforms();
            for iTransform = 1:length(transforms)
                path = fullfile(testCase.tempFolder,classNames(iTransform)+".nc");
                ncfile = transforms{iTransform}.writeToFile(path,shouldOverwriteExisting=true);
                ncfile.close();

                wvt = feval(classNames(iTransform)+".waveVortexTransformFromFile",char(path));
                testCase.verifyClass(wvt,classNames(iTransform))
                TestNetCDFHandleOwnership.verifyWritableOpen(testCase,path);
            end
        end

        function directSubclassRestorationErrorClosesFile(testCase)
            invalidPath = fullfile(testCase.tempFolder,"invalid-direct-transform.nc");
            ncfile = NetCDFFile(invalidPath);
            ncfile.close();

            didThrow = false;
            try
                WVTransformConstantStratification.waveVortexTransformFromFile(char(invalidPath));
            catch
                didThrow = true;
            end
            testCase.verifyTrue(didThrow,"Restoring an invalid transform file must throw an error.")
            TestNetCDFHandleOwnership.verifyWritableOpen(testCase,invalidPath);
        end

        function invalidTransformMetadataClosesProbe(testCase)
            invalidPath = fullfile(testCase.tempFolder,"invalid-transform.nc");
            ncfile = NetCDFFile(invalidPath);
            ncfile.addAttribute("WVTransform","UnavailableWVTransform");
            ncfile.close();

            testCase.verifyError(@()WVTransform.waveVortexTransformFromFile(char(invalidPath)),"WVTransform:InvalidTransformClass")
            TestNetCDFHandleOwnership.verifyWritableOpen(testCase,invalidPath);
        end

        function missingTransformMetadataClosesProbe(testCase)
            missingPath = fullfile(testCase.tempFolder,"missing-transform.nc");
            ncfile = NetCDFFile(missingPath);
            ncfile.close();

            testCase.verifyError(@()WVTransform.waveVortexTransformFromFile(char(missingPath)),"WVTransform:MissingTransformClass")
            TestNetCDFHandleOwnership.verifyWritableOpen(testCase,missingPath);
        end

        function modelRestartOwnsSingleWriter(testCase)
            wvt = TestNetCDFHandleOwnership.newTransform();
            model = WVModel(wvt,shouldUseLinearDynamics=true);
            model.createNetCDFFileForModelOutput(testCase.modelPath,outputInterval=1,shouldOverwriteExisting=true);
            model.integrateToTime(2,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            model.closeNetCDFFile();

            warningState = warning;
            warningCleanup = onCleanup(@()warning(warningState));
            warning("off","all")
            resumedModel = WVModel.modelFromFile(char(testCase.modelPath));
            clear warningCleanup
            cleanup = onCleanup(@()resumedModel.closeNetCDFFile());
            resumedModel.integrateToTime(3,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            resumedModel.closeNetCDFFile();

            ncfile = NetCDFFile(testCase.modelPath,shouldReadOnly=true);
            fileCleanup = onCleanup(@()TestNetCDFHandleOwnership.closeIfOpen(ncfile));
            testCase.verifyEqual(ncfile.readVariables("wave-vortex/t"),(0:3).')
            ncfile.close();
            clear fileCleanup cleanup
        end

        function initializationFailureRollsBackAndAllowsRetry(testCase)
            path = fullfile(testCase.tempFolder,"initialization-failure.nc");
            model = WVModel(TestNetCDFHandleOwnership.newTransform(),shouldUseLinearDynamics=true);
            outputFile = model.createNetCDFFileForModelOutput(path,outputInterval=1,shouldOverwriteExisting=true);
            outputGroup = outputFile.outputGroups(1);
            failingObserver = WVFailingObservingSystem(model,failurePhase="initialize");
            outputGroup.addObservingSystem(failingObserver);

            testCase.verifyError(@()model.integrateToTime(1,shouldShowIntegrationDiagnostics=false,callback=@(~)[]),'')
            testCase.verifyFalse(isfile(path))
            testCase.verifyFalse(outputFile.didInitializeStorage)
            testCase.verifyEmpty(outputFile.ncfile)
            testCase.verifyFalse(outputGroup.didInitializeStorage)
            testCase.verifyEmpty(outputGroup.group)
            testCase.verifyEqual(outputGroup.incrementsWrittenToGroup,uint64(0))
            testCase.verifyEqual(outputGroup.timeOfLastIncrementWrittenToGroup,-Inf)

            outputGroup.removeObservingSystem(failingObserver);
            model.integrateToTime(1,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            model.closeNetCDFFile();
            testCase.verifyTrue(isfile(path))
            TestNetCDFHandleOwnership.verifyWritableOpen(testCase,path);
        end

        function writeFailureClosesAndPreservesFile(testCase)
            path = fullfile(testCase.tempFolder,"write-failure.nc");
            model = WVModel(TestNetCDFHandleOwnership.newTransform(),shouldUseLinearDynamics=true);
            outputFile = model.createNetCDFFileForModelOutput(path,outputInterval=1,shouldOverwriteExisting=true);
            outputFile.outputGroups(1).addObservingSystem(WVFailingObservingSystem(model,failurePhase="write"));

            testCase.verifyError(@()model.integrateToTime(1,shouldShowIntegrationDiagnostics=false,callback=@(~)[]),'')
            testCase.verifyTrue(isfile(path))
            testCase.verifyEmpty(outputFile.ncfile)
            TestNetCDFHandleOwnership.verifyWritableOpen(testCase,path);
        end

        function observerRestorationFailureClosesWriter(testCase)
            path = fullfile(testCase.tempFolder,"restoration-failure.nc");
            model = WVModel(TestNetCDFHandleOwnership.newTransform(),shouldUseLinearDynamics=true);
            model.createNetCDFFileForModelOutput(path,outputInterval=1,shouldOverwriteExisting=true);
            model.integrateToTime(1,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            model.closeNetCDFFile();

            ncfile = NetCDFFile(path,shouldReadOnly=false);
            cleanup = onCleanup(@()TestNetCDFHandleOwnership.closeIfOpen(ncfile));
            observerGroup = TestNetCDFHandleOwnership.groupWithAnnotatedClass(ncfile,'WVEulerianFields');
            testCase.assertNotEmpty(observerGroup)
            observerGroup.addAttribute('AnnotatedClass','UnavailableObservingSystem');
            ncfile.close();
            clear cleanup

            didThrow = false;
            try
                WVModel.modelFromFile(char(path));
            catch
                didThrow = true;
            end
            testCase.verifyTrue(didThrow)
            TestNetCDFHandleOwnership.verifyWritableOpen(testCase,path);
        end
    end

    methods (Static, Access=private)
        function wvt = newTransform()
            N0 = 5.2e-3;
            wvt = WVTransformConstantStratification([1e4,1e4,1000],[8,8,5],N0=N0,latitude=30);
        end

        function [transforms,classNames] = supportedTransforms()
            N0 = 5.2e-3;
            N2 = @(z) N0*N0*ones(size(z));
            transforms = {
                WVTransformConstantStratification([1e4,1e4,1000],[8,8,5],N0=N0,latitude=30)
                WVTransformHydrostatic([1e4,1e4,1000],[8,8,5],N2=N2,latitude=30)
                WVTransformBoussinesq([1e4,1e4,1000],[8,8,5],N2=N2,latitude=30)
                WVTransformStratifiedQG([1e4,1e4,1000],[8,8,5],N2=N2,latitude=30)
                WVTransformBarotropicQG([1e4,1e4],[8,8],latitude=30)
            };
            classNames = [
                "WVTransformConstantStratification"
                "WVTransformHydrostatic"
                "WVTransformBoussinesq"
                "WVTransformStratifiedQG"
                "WVTransformBarotropicQG"
            ];
        end

        function verifyWritableOpen(testCase,path)
            ncfile = NetCDFFile(path,shouldReadOnly=false);
            cleanup = onCleanup(@()TestNetCDFHandleOwnership.closeIfOpen(ncfile));
            testCase.verifyNotEmpty(ncfile.id)
            ncfile.close();
            clear cleanup
        end

        function closeIfOpen(ncfile)
            if ~isempty(ncfile) && isvalid(ncfile) && ~isempty(ncfile.id)
                ncfile.close();
            end
        end

        function matchingGroup = groupWithAnnotatedClass(group,className)
            matchingGroup = NetCDFGroup.empty(0,0);
            for iGroup = 1:length(group.groups)
                candidate = group.groups(iGroup);
                if isKey(candidate.attributes,'AnnotatedClass') && strcmp(candidate.attributes('AnnotatedClass'),className)
                    matchingGroup = candidate;
                    return
                end
                matchingGroup = TestNetCDFHandleOwnership.groupWithAnnotatedClass(candidate,className);
                if ~isempty(matchingGroup)
                    return
                end
            end
        end
    end
end
