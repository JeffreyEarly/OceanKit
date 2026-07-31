classdef TestNetCDFHandleOwnership < matlab.unittest.TestCase
    properties
        tempFolder
        transformPath
        modelPath
    end

    methods (TestMethodSetup)
        function createFixtures(testCase)
            testCase.tempFolder = string(tempname);
            mkdir(testCase.tempFolder);
            testCase.transformPath = fullfile(testCase.tempFolder,"transform.nc");
            testCase.modelPath = fullfile(testCase.tempFolder,"model.nc");

            wvt = TestNetCDFHandleOwnership.newTransform();
            ncfile = wvt.writeToFile(testCase.transformPath,shouldOverwriteExisting=true);
            ncfile.close();
        end
    end

    methods (TestMethodTeardown)
        function removeFixtures(testCase)
            if isfolder(testCase.tempFolder)
                rmdir(testCase.tempFolder,"s");
            end
        end
    end

    methods (Test)
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

            resumedModel = WVModel.modelFromFile(char(testCase.modelPath));
            cleanup = onCleanup(@()resumedModel.closeNetCDFFile());
            resumedModel.integrateToTime(3,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            resumedModel.closeNetCDFFile();

            ncfile = NetCDFFile(testCase.modelPath,shouldReadOnly=true);
            fileCleanup = onCleanup(@()TestNetCDFHandleOwnership.closeIfOpen(ncfile));
            testCase.verifyEqual(ncfile.readVariables("wave-vortex/t"),(0:3).')
            ncfile.close();
            clear fileCleanup cleanup
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
    end
end
