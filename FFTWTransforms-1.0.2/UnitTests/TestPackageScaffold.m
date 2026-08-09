classdef TestPackageScaffold < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addRepositoryToPath(testCase)
            import matlab.unittest.fixtures.PathFixture
            testCase.applyFixture(PathFixture(TestPackageScaffold.repositoryRoot));
        end
    end

    methods (Test)
        function testPackageManifest(testCase)
            package = matlab.mpm.Package(TestPackageScaffold.repositoryRoot);
            testCase.verifyEqual(package.Name,"FFTWTransforms");
            testCase.verifyEqual(package.DisplayName,"FFTWTransforms");
            versionText = string(package.Version);
            testCase.verifyNotEmpty(regexp(versionText,"^(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)\.(0|[1-9][0-9]*)$",'once'));
            testCase.verifyGreaterThanOrEqual(package.Version.Major,0);
            testCase.verifyGreaterThanOrEqual(package.Version.Minor,0);
            testCase.verifyGreaterThanOrEqual(package.Version.Patch,0);
            testCase.verifyEqual(string(package.ReleaseCompatibility),">=R2024b");
            testCase.verifyEqual(package.Summary,"Reusable FFTW transforms for MATLAB.");
            testCase.verifyEmpty(package.Dependencies);
            testCase.verifyEmpty(package.Folders);
            testCase.verifyNotEmpty(regexp(package.ID,"^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$",'once'));
        end

        function testLicensingAndAuthoringFiles(testCase)
            root = TestPackageScaffold.repositoryRoot;
            header = string(fileread(fullfile(root,"fftw3.h")));
            notices = string(fileread(fullfile(root,"THIRD_PARTY_NOTICES.md")));
            readme = string(fileread(fullfile(root,"README.md")));
            testCase.verifyTrue(contains(header,"license applies *only* to this header file"));
            testCase.verifyTrue(contains(notices,"Copyright (c) 2003, 2007-14 Matteo Frigo"));
            testCase.verifyTrue(contains(notices,"GNU General Public License"));
            testCase.verifyTrue(contains(readme,"does not distribute FFTW libraries"));
            testCase.verifyTrue(isfile(fullfile(root,"LICENSE")));
            testCase.verifyTrue(isfile(fullfile(root,"CHANGELOG.md")));
            testCase.verifyFalse(isfile(fullfile(root,"tools","ci_release.m")));
        end

    end

    methods (Static, Access=private)
        function root = repositoryRoot()
            testDirectory = fileparts(mfilename('fullpath'));
            root = fileparts(testDirectory);
        end
    end
end
