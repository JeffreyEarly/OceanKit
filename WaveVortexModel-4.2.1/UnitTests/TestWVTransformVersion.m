classdef TestWVTransformVersion < matlab.unittest.TestCase

    methods (Test, TestTags = "smoke")
        function testVersionMatchesPackageManifest(testCase)
            wvt = TestWVTransformVersion.barotropicTransform();
            expectedVersion = TestWVTransformVersion.packageVersionFromManifest();

            testCase.verifyEqual(wvt.version,expectedVersion);
        end

        function testWriteToFileUsesManifestBackedVersion(testCase)
            expectedVersion = TestWVTransformVersion.packageVersionFromManifest();
            fixture = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            path = fullfile(fixture.Folder,'version.nc');

            wvt = TestWVTransformVersion.barotropicTransform();
            ncfile = wvt.writeToFile(path,shouldOverwriteExisting=true);
            ncfile.close();
            ncfile = NetCDFFile(path,shouldReadOnly=true);
            cleanup = onCleanup(@()TestWVTransformVersion.closeIfOpen(ncfile));

            testCase.verifyEqual(string(ncfile.attributes('model_version')),expectedVersion);
            testCase.verifyTrue(contains(string(ncfile.attributes('source')),expectedVersion));

            ncfile.close();
            clear cleanup
        end
    end

    methods (Static, Access=private)
        function wvt = barotropicTransform()
            wvt = WVTransformBarotropicQG([10e3 10e3],[8 8],latitude=30,shouldAntialias=false);
        end

        function version = packageVersionFromManifest()
            classFilePath = which('WVTransform');
            packageRoot = fileparts(fileparts(classFilePath));
            manifestPath = fullfile(packageRoot,'resources','mpackage.json');
            manifest = jsondecode(fileread(manifestPath));
            version = string(manifest.version);
        end

        function closeIfOpen(ncfile)
            if ~isempty(ncfile) && isvalid(ncfile) && ~isempty(ncfile.id)
                ncfile.close();
            end
        end
    end
end
