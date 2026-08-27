classdef TestShouldExcludeConjugatesPersistence < matlab.unittest.TestCase
    properties (TestParameter)
        transformType = struct( ...
            constantHydrostatic="constantHydrostatic", ...
            constantNonhydrostatic="constantNonhydrostatic", ...
            hydrostatic="hydrostatic", ...
            boussinesq="boussinesq", ...
            stratifiedQG="stratifiedQG", ...
            barotropicQG="barotropicQG")
    end

    properties
        tempFolder
    end

    methods (TestMethodSetup)
        function createTemporaryFolder(testCase)
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            testCase.tempFolder = string(fixture.Folder);
        end
    end

    methods (Test, TestTags="full")
        function correctedFieldRoundTrips(testCase)
            path = fullfile(testCase.tempFolder,"corrected.nc");
            geometry = testCase.newGeometry(false);
            ncfile = geometry.writeToFile(char(path),shouldOverwriteExisting=true);
            cleanup = onCleanup(@()testCase.closeIfOpen(ncfile));

            testCase.verifyTrue(ncfile.hasVariableWithName('shouldExcludeConjugates'))
            testCase.verifyFalse(ncfile.hasVariableWithName('shouldExludeConjugates'))
            ncfile.close();
            clear cleanup

            restored = testCase.restoreGeometry(path);
            testCase.verifyFalse(restored.shouldExcludeConjugates)
        end

        function legacyShouldExludeConjugatesLoadsThroughout4x(testCase)
            path = fullfile(testCase.tempFolder,"legacy.nc");
            testCase.writeGeometryFileWithoutCanonicalField(path,true);
            ncfile = NetCDFFile(path,shouldReadOnly=false);
            cleanup = onCleanup(@()testCase.closeIfOpen(ncfile));
            ncfile.addVariable('shouldExludeConjugates',{},false);
            ncfile.close();
            clear cleanup

            restored = testCase.verifyWarningFree(@()testCase.restoreGeometry(path));
            testCase.verifyFalse(restored.shouldExcludeConjugates)
        end

        function correctedFieldTakesPrecedence(testCase)
            path = fullfile(testCase.tempFolder,"both.nc");
            geometry = testCase.newGeometry(true);
            ncfile = geometry.writeToFile(char(path),shouldOverwriteExisting=true);
            cleanup = onCleanup(@()testCase.closeIfOpen(ncfile));
            ncfile.addVariable('shouldExludeConjugates',{},false);
            ncfile.close();
            clear cleanup

            restored = testCase.restoreGeometry(path);
            testCase.verifyTrue(restored.shouldExcludeConjugates)
        end

        function missingFieldUsesExistingFailure(testCase)
            path = fullfile(testCase.tempFolder,"missing.nc");
            testCase.writeGeometryFileWithoutCanonicalField(path,true);

            testCase.verifyError(@()testCase.restoreGeometry(path),'')
        end

        function stableTransformRoundTripRemainsGreen(testCase,transformType)
            path = fullfile(testCase.tempFolder,transformType+".nc");
            original = testCase.newTransform(transformType);
            ncfile = original.writeToFile(char(path),shouldOverwriteExisting=true);
            cleanup = onCleanup(@()testCase.closeIfOpen(ncfile));

            testCase.verifyFalse(ncfile.hasVariableWithName('shouldExcludeConjugates'))
            testCase.verifyFalse(ncfile.hasVariableWithName('shouldExludeConjugates'))
            ncfile.close();
            clear cleanup

            restored = WVTransform.waveVortexTransformFromFile(char(path));
            testCase.verifyClass(restored,class(original))
            testCase.verifyTrue(restored.shouldExcludeConjugates)
            testCase.verifyEqual(restored.spatialMatrixSize,original.spatialMatrixSize)
            testCase.verifyEqual(restored.spectralMatrixSize,original.spectralMatrixSize)
        end
    end

    methods (Access=private)
        function writeGeometryFileWithoutCanonicalField(testCase,path,value)
            geometry = testCase.newGeometry(value);
            properties = setdiff(geometry.requiredProperties,{'shouldExcludeConjugates'},'stable');
            ncfile = geometry.writeToFile(char(path),properties{:},shouldOverwriteExisting=true,shouldAddRequiredProperties=false);
            cleanup = onCleanup(@()testCase.closeIfOpen(ncfile));
            ncfile.close();
            clear cleanup
        end
    end

    methods (Static, Access=private)
        function geometry = newGeometry(shouldExcludeConjugates)
            geometry = WVGeometryDoublyPeriodic([4000 3000],[8 6],shouldAntialias=false,shouldExcludeNyquist=false,shouldExcludeConjugates=shouldExcludeConjugates);
        end

        function wvt = newTransform(transformType)
            Lxyz = [4000 3000 1000];
            Nxyz = [8 6 5];
            N2 = @(z) 2e-5*exp(z/4000);
            switch transformType
                case "constantHydrostatic"
                    wvt = WVTransformConstantStratification(Lxyz,Nxyz,N0=sqrt(2e-5),latitude=45,isHydrostatic=true,shouldAntialias=false);
                case "constantNonhydrostatic"
                    wvt = WVTransformConstantStratification(Lxyz,Nxyz,N0=sqrt(2e-5),latitude=45,isHydrostatic=false,shouldAntialias=false);
                case "hydrostatic"
                    wvt = WVTransformHydrostatic(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=false);
                case "boussinesq"
                    wvt = WVTransformBoussinesq(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=false);
                case "stratifiedQG"
                    wvt = WVTransformStratifiedQG(Lxyz,Nxyz,N2=N2,latitude=45,shouldAntialias=false);
                case "barotropicQG"
                    wvt = WVTransformBarotropicQG(Lxyz(1:2),Nxyz(1:2),latitude=45,shouldAntialias=false);
            end
        end

        function closeIfOpen(ncfile)
            if ~isempty(ncfile) && isvalid(ncfile) && ~isempty(ncfile.id)
                ncfile.close();
            end
        end

        function geometry = restoreGeometry(path)
            ncfile = NetCDFFile(path,shouldReadOnly=true);
            cleanup = onCleanup(@()TestShouldExcludeConjugatesPersistence.closeIfOpen(ncfile));
            geometry = WVGeometryDoublyPeriodic.geometryFromGroup(ncfile);
            ncfile.close();
            clear cleanup
        end
    end
end
