classdef TestWVTransformVersion < matlab.unittest.TestCase

    methods (Test)
        function testVersionMatchesPackageManifest(testCase)
            wvt = TestWVTransformVersion.barotropicTransform();
            expectedVersion = TestWVTransformVersion.packageVersionFromManifest();

            testCase.verifyEqual(wvt.version,expectedVersion);
        end

        function testWriteToFileUsesManifestBackedVersion(testCase)
            expectedVersion = TestWVTransformVersion.packageVersionFromManifest();
            path = [tempname,'.nc'];
            cleanup = onCleanup(@() TestWVTransformVersion.deleteIfExists(path));

            wvt = TestWVTransformVersion.barotropicTransform();
            ncfile = wvt.writeToFile(path,shouldOverwriteExisting=true);
            ncfile.close();
            ncfile = NetCDFFile(path,shouldReadOnly=true);

            testCase.verifyEqual(string(ncfile.attributes('model_version')),expectedVersion);
            testCase.verifyTrue(contains(string(ncfile.attributes('source')),expectedVersion));

            ncfile = [];
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

        function deleteIfExists(path)
            if isfile(path)
                try
                    delete(path);
                catch
                end
            end
        end
    end
end
