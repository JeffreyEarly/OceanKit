function plan = buildfile
import matlab.buildtool.Task
import matlab.buildtool.TaskGroup

tasks = [
    Task(Actions=@testSmokeTask,Description="Run the fast smoke test category.",DisableIncremental=true)
    Task(Actions=@testFullTask,Description="Run all non-optional, non-exhaustive tests.",DisableIncremental=true)
    Task(Actions=@testExhaustiveTask,Description="Run exhaustive numerical test matrices.",DisableIncremental=true)
    Task(Actions=@testOptionalTask,Description="Run tests that require optional dependencies.",DisableIncremental=true)
    ];

plan = buildplan;
plan("test") = TaskGroup(tasks,TaskNames=["smoke"; "full"; "exhaustive"; "optional"],Description="Run WaveVortexModel test categories.");
documentationTasks = [
    Task(Actions=@documentationBuildTask,Description="Build and transactionally replace committed documentation.",DisableIncremental=true)
    Task(Actions=@documentationCheckTask,Description="Verify committed documentation against a clean generated site.",DisableIncremental=true)
    ];
plan("docs") = TaskGroup(documentationTasks,TaskNames=["build"; "check"],Description="Build or verify website documentation.");
plan("analyze") = Task(Actions=@analyzeTask,Description="Analyze production MATLAB source for correctness findings.",DisableIncremental=true);
plan("kernel:contract") = Task(Actions=@kernelContractTask,Description="Build and run the portable C++ kernel contract tests.",DisableIncremental=true);
plan.DefaultTasks = "test:smoke";
end

function kernelContractTask(~)
repositoryRoot = fileparts(mfilename("fullpath"));
scriptPath = fullfile(repositoryRoot,"tools","compiled-kernel","run_contract_tests.sh");
[status,output] = system(sprintf('"%s"',scriptPath));
fprintf("%s",output);
if status ~= 0
    error("WaveVortexModel:KernelContractTestsFailed","Portable C++ kernel contract tests failed.");
end
end

function documentationBuildTask(~)
withToolsPath(@()build_website_documentation(rootDir=fileparts(mfilename("fullpath"))));
end

function documentationCheckTask(~)
withToolsPath(@()check_website_documentation(rootDir=fileparts(mfilename("fullpath"))));
end

function withToolsPath(action)
repositoryRoot = fileparts(mfilename("fullpath"));
originalPath = path;
pathCleanup = onCleanup(@()path(originalPath));
addpath(fullfile(repositoryRoot,"tools"));
action();
clear pathCleanup
end

function analyzeTask(~)
repositoryRoot = fileparts(mfilename("fullpath"));
toolsFolder = fullfile(repositoryRoot,"tools");
originalPath = path;
pathCleanup = onCleanup(@()path(originalPath));
addpath(toolsFolder);
analyzeProductionCode(repositoryRoot);
clear pathCleanup
end

function testSmokeTask(~)
runTestCategory("smoke","smoke");
end

function testFullTask(~)
runTestCategory("full",["smoke" "full"]);
end

function testExhaustiveTask(~)
runTestCategory("exhaustive","exhaustive");
end

function testOptionalTask(~)
runTestCategory("optional","optional");
end

function runTestCategory(categoryName,selectedTags)
repositoryRoot = fileparts(mfilename("fullpath"));
testFolder = fullfile(repositoryRoot,"UnitTests");
originalPath = path;
pathCleanup = onCleanup(@()path(originalPath));
pathEntries = string(strsplit(path,pathsep));
if ~any(pathEntries == string(testFolder))
    addpath(testFolder);
end

folderSuite = matlab.unittest.TestSuite.fromFolder(testFolder,IncludingSubfolders=true);
rejectAccidentalScriptTests(folderSuite);
suite = formalTestSuite(testFolder);
[expectedPairs,discoveredPairs,selectedSuite] = validateTestDiscovery(testFolder,suite,selectedTags);

missingPairs = setdiff(expectedPairs,discoveredPairs);
unexpectedPairs = setdiff(discoveredPairs,expectedPairs);
if ~isempty(missingPairs) || ~isempty(unexpectedPairs)
    details = discoveryDifferenceDetails(missingPairs,unexpectedPairs);
    error("WaveVortexModel:TestDiscoveryFailed","The %s test category did not match its declared test methods.%s",categoryName,details);
end
if isempty(selectedSuite)
    error("WaveVortexModel:EmptyTestCategory","No tests were discovered for the %s category.",categoryName);
end

runner = testrunner("textoutput");
results = runner.run(selectedSuite);
passedCount = nnz([results.Passed]);
failedCount = nnz([results.Failed]);
incompleteCount = nnz([results.Incomplete]);
duration = sum([results.Duration]);
fprintf("\n%s test summary: total=%d, passed=%d, failed=%d, incomplete=%d, duration=%.3f s\n",categoryName,numel(results),passedCount,failedCount,incompleteCount,duration);

if failedCount > 0 || incompleteCount > 0
    error("WaveVortexModel:UnsuccessfulTestRun","The %s test category was unsuccessful: %d failed and %d incomplete.",categoryName,failedCount,incompleteCount);
end
clear pathCleanup
end

function suite = formalTestSuite(testFolder)
testFiles = dir(fullfile(testFolder,"Test*.m"));
classSuites = cell(numel(testFiles),1);
for iFile = 1:numel(testFiles)
    className = erase(string(testFiles(iFile).name),".m");
    testClass = meta.class.fromName(className);
    if isempty(testClass) || ~isTestCaseClass(testClass)
        error("WaveVortexModel:InvalidFormalTest","%s must define a matlab.unittest.TestCase class named %s.",testFiles(iFile).name,className);
    end
    classSuites{iFile} = matlab.unittest.TestSuite.fromClass(testClass);
end
suite = [classSuites{:}];
end

function rejectAccidentalScriptTests(suite)
discoveredClassNames = string({suite.TestClass});
if any(discoveredClassNames == "")
    accidentalNames = unique(string({suite(discoveredClassNames == "").Name}));
    error("WaveVortexModel:AccidentalScriptTests","UnitTests discovery found script-based tests: %s",strjoin(accidentalNames,", "));
end
end

function [expectedPairs,discoveredPairs,selectedSuite] = validateTestDiscovery(testFolder,suite,selectedTags)
primaryTags = ["smoke" "full" "exhaustive" "optional"];
testFiles = dir(fullfile(testFolder,"Test*.m"));
expectedPairsByFile = cell(numel(testFiles),1);

for iFile = 1:numel(testFiles)
    className = erase(string(testFiles(iFile).name),".m");
    testClass = meta.class.fromName(className);
    if isempty(testClass) || ~isTestCaseClass(testClass)
        error("WaveVortexModel:InvalidFormalTest","%s must define a matlab.unittest.TestCase class named %s.",testFiles(iFile).name,className);
    end

    testMethods = testClass.MethodList.findobj("Test",true);
    if isempty(testMethods)
        error("WaveVortexModel:InvalidFormalTest","%s does not declare any test methods.",className);
    end
    selectedMethodPairs = strings(numel(testMethods),1);
    nSelectedMethods = 0;
    for iMethod = 1:numel(testMethods)
        methodTags = string(testMethods(iMethod).TestTags);
        validatePrimaryTags(className,testMethods(iMethod).Name,methodTags,primaryTags);
        if any(ismember(methodTags,selectedTags))
            nSelectedMethods = nSelectedMethods + 1;
            selectedMethodPairs(nSelectedMethods) = className + "/" + testMethods(iMethod).Name;
        end
    end
    expectedPairsByFile{iFile} = selectedMethodPairs(1:nSelectedMethods);
end
nonemptyExpectedPairs = expectedPairsByFile(~cellfun(@isempty,expectedPairsByFile));
if isempty(nonemptyExpectedPairs)
    expectedPairs = strings(0,1);
else
    expectedPairs = vertcat(nonemptyExpectedPairs{:});
end

discoveredClassNames = string({suite.TestClass});
discoveredMethodNames = string({suite.ProcedureName});

selectedMask = false(size(suite));
for iTest = 1:numel(suite)
    testTags = string(suite(iTest).Tags);
    validatePrimaryTags(discoveredClassNames(iTest),discoveredMethodNames(iTest),testTags,primaryTags);
    selectedMask(iTest) = any(ismember(testTags,selectedTags));
end

selectedSuite = suite(selectedMask);
discoveredPairs = unique(discoveredClassNames(selectedMask) + "/" + discoveredMethodNames(selectedMask));
expectedPairs = unique(expectedPairs);
end

function validatePrimaryTags(className,methodName,testTags,primaryTags)
matchingTags = intersect(testTags,primaryTags);
if numel(matchingTags) ~= 1
    error("WaveVortexModel:InvalidTestClassification","%s/%s must declare exactly one primary test tag from %s; found %s.",className,methodName,strjoin(primaryTags,", "),strjoin(matchingTags,", "));
end
end

function tf = isTestCaseClass(testClass)
tf = testClass.Name == "matlab.unittest.TestCase";
if tf
    return
end
for superclass = testClass.SuperclassList'
    if isTestCaseClass(superclass)
        tf = true;
        return
    end
end
end

function details = discoveryDifferenceDetails(missingPairs,unexpectedPairs)
details = "";
if ~isempty(missingPairs)
    details = details + newline + "Missing: " + strjoin(missingPairs,", ");
end
if ~isempty(unexpectedPairs)
    details = details + newline + "Unexpected: " + strjoin(unexpectedPairs,", ");
end
end
