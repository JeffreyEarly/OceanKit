function results = RunAllUnitTests()
testRoot = fileparts(mfilename('fullpath'));
originalFolder = pwd;
cleanup = onCleanup(@() cd(originalFolder));
cd(testRoot);

repoRoot = fileparts(testRoot);
workspaceRoot = fileparts(repoRoot);
dependencyRoot = fullfile(workspaceRoot,'spline-core');
distributionRoot = fullfile(workspaceRoot,'distributions');

addPackageToPath(repoRoot);

if ~isfolder(dependencyRoot)
    error('Expected dependency repo at %s.',dependencyRoot);
end
if ~isfolder(distributionRoot)
    error('Expected dependency repo at %s.',distributionRoot);
end

addPackageToPath(distributionRoot);
addPackageToPath(dependencyRoot);

import matlab.unittest.TestRunner
import matlab.unittest.TestSuite

runner = testrunner("textoutput");
suite = [ ...
    TestSuite.fromClass(?IntegratorUnitTests) ...
    TestSuite.fromClass(?IntegratorWithObstaclesUnitTests) ...
    TestSuite.fromClass(?KernelDensityEstimateUnitTests) ...
    TestSuite.fromClass(?LinearVelocityFieldEstimationUnitTests) ...
    TestSuite.fromClass(?GriddedStreamfunctionUnitTests) ...
    TestSuite.fromClass(?GriddedStreamfunctionBootstrapUnitTests) ...
    TestSuite.fromClass(?OscroftFairComparisonUnitTests) ...
    TestSuite.fromClass(?SecondMomentMethodUnitTests) ...
    ];

results = runner.run(suite);

if nargout == 0
    disp(table(results))
    if any(~[results.Passed])
        error('Some unit tests failed.');
    end
end
end

function addPackageToPath(repoRoot)
metadataPath = fullfile(repoRoot,'resources','mpackage.json');

if isfile(metadataPath)
    metadata = jsondecode(fileread(metadataPath));
    for iFolder = 1:length(metadata.folders)
        addpath(fullfile(repoRoot,metadata.folders(iFolder).path));
    end
end

addpath(repoRoot);
end
