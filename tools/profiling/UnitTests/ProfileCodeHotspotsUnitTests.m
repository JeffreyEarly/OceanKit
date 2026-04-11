classdef ProfileCodeHotspotsUnitTests < matlab.unittest.TestCase

    properties
        toolRoot
    end

    methods (TestClassSetup)
        function addToolFolderToPath(testCase)
            testCase.toolRoot = fileparts(fileparts(mfilename("fullpath")));
            addpath(testCase.toolRoot)
        end
    end

    methods (Test)
        function analyzeProfileInfoComputesFunctionMetrics(testCase)
            projectRoot = "/tmp/project";
            functionTable = [
                ProfileCodeHotspotsUnitTests.functionEntry( ...
                    completeName="pkg.selfHeavy", functionName="selfHeavy", ...
                    fileName=projectRoot + "/selfHeavy.m", type="M-function", ...
                    totalTime=8, numCalls=2, children=struct("TotalTime", {2, 1}), ...
                    parents=struct("Index", {1}), executedLines=[11 1 2; 12 1 1], isRecursive=false)
                ProfileCodeHotspotsUnitTests.functionEntry( ...
                    completeName="pkg.highCallCount", functionName="highCallCount", ...
                    fileName=projectRoot + "/highCallCount.m", type="M-function", ...
                    totalTime=3, numCalls=100, children=[], parents=struct("Index", {1, 2}), ...
                    executedLines=[20 100 2.5], isRecursive=true)
                ProfileCodeHotspotsUnitTests.functionEntry( ...
                    completeName="toolbox.childHeavy", functionName="childHeavy", ...
                    fileName="/Applications/MATLAB/toolbox/childHeavy.m", type="M-function", ...
                    totalTime=5, numCalls=1, children=struct("Time", {4}), ...
                    parents=[], executedLines=[30 1 1], isRecursive=false)
                ProfileCodeHotspotsUnitTests.functionEntry( ...
                    completeName="pkg.unknownChildren", functionName="unknownChildren", ...
                    fileName=projectRoot + "/unknownChildren.m", type="M-function", ...
                    totalTime=2, numCalls=4, children=struct("Other", 1), parents=[], ...
                    executedLines=[], isRecursive=false)
                ProfileCodeHotspotsUnitTests.functionEntry( ...
                    completeName="builtin.plus", functionName="plus", ...
                    fileName="", type="M-builtin", totalTime=1, numCalls=50, ...
                    children=[], parents=[], executedLines=[], isRecursive=false)
            ];
            profileInfo = struct("FunctionTable", functionTable);

            analysis = analyzeProfileInfo(profileInfo, projectRoots=projectRoot, maxFunctions=10, maxLines=10);

            testCase.verifyEqual(height(analysis.functionMetrics), 5)
            testCase.verifyEqual(analysis.functionMetrics.childTime(1), 3)
            testCase.verifyEqual(analysis.functionMetrics.selfTime(1), 5)
            testCase.verifyEqual(analysis.functionMetrics.selfTimePerCall(1), 2.5)
            testCase.verifyTrue(analysis.functionMetrics.isLeaf(2))
            testCase.verifyTrue(analysis.functionMetrics.isRecursive(2))
            testCase.verifyEqual(analysis.functionMetrics.numParents(2), 2)
            testCase.verifyEqual(analysis.functionMetrics.numChildren(1), 2)
            testCase.verifyEqual(analysis.functionMetrics.numExecutedLines(1), 2)
            testCase.verifyEqual(analysis.functionMetrics.childTime(3), 4)
            testCase.verifyEqual(analysis.functionMetrics.selfTime(3), 1)
            testCase.verifyTrue(isnan(analysis.functionMetrics.childTime(4)))
            testCase.verifyTrue(isnan(analysis.functionMetrics.selfTime(4)))
            testCase.verifyTrue(analysis.functionMetrics.isProjectCode(1))
            testCase.verifyTrue(analysis.functionMetrics.isToolboxCode(3))
            testCase.verifyTrue(analysis.functionMetrics.isBuiltIn(5))
        end

        function analyzeProfileInfoBuildsProjectFilteredRankings(testCase)
            projectRoot = "/tmp/project";
            functionTable = [
                ProfileCodeHotspotsUnitTests.functionEntry( ...
                    completeName="pkg.slow", functionName="slow", ...
                    fileName=projectRoot + "/slow.m", type="M-function", ...
                    totalTime=10, numCalls=2, children=struct("TotalTime", {3}), ...
                    parents=[], executedLines=[10 1 4.5; 11 1 1], isRecursive=false)
                ProfileCodeHotspotsUnitTests.functionEntry( ...
                    completeName="pkg.calls", functionName="calls", ...
                    fileName=projectRoot + "/calls.m", type="M-function", ...
                    totalTime=4, numCalls=40, children=[], parents=[], ...
                    executedLines=[20 40 2.5], isRecursive=false)
                ProfileCodeHotspotsUnitTests.functionEntry( ...
                    completeName="toolbox.helper", functionName="helper", ...
                    fileName="/Applications/MATLAB/toolbox/helper.m", type="M-function", ...
                    totalTime=20, numCalls=1, children=struct("Time", {1}), ...
                    parents=[], executedLines=[30 1 0.5], isRecursive=false)
            ];
            profileInfo = struct("FunctionTable", functionTable);

            analysis = analyzeProfileInfo(profileInfo, projectRoots=projectRoot, maxFunctions=10, maxLines=10);

            testCase.verifyEqual(analysis.topByTotalTime.completeName(1), "toolbox.helper")
            testCase.verifyEqual(analysis.topBySelfTime.completeName(1), "toolbox.helper")
            testCase.verifyEqual(analysis.topByNumCalls.completeName(1), "pkg.calls")
            testCase.verifyEqual(analysis.topByTotalTimePerCall.completeName(1), "toolbox.helper")
            testCase.verifyEqual(analysis.topProjectByTotalTime.completeName(1), "pkg.slow")
            testCase.verifyEqual(analysis.topProjectBySelfTime.completeName(1), "pkg.slow")
            testCase.verifyEqual(analysis.topProjectByNumCalls.completeName(1), "pkg.calls")
            testCase.verifyEqual(analysis.topExternalBySelfTime.completeName(1), "toolbox.helper")
            testCase.verifyEqual(analysis.topLines.lineNumber(1), 10)
            testCase.verifyEqual(analysis.lineMetrics.lineFractionOfFunction(1), 0.45, AbsTol=1e-12)
            testCase.verifyFalse(isempty(analysis.priorityTargets))
            testCase.verifyEqual(analysis.priorityTargets.targetSubtype(1), "selfFunction")
            testCase.verifyNotEmpty(analysis.summary)
        end

        function analyzeProfileInfoClassifiesCallsiteAndComputeLines(testCase)
            projectRoot = string(testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder);
            wrapperPath = fullfile(projectRoot, "wrapper.m");
            computePath = fullfile(projectRoot, "computeKernel.m");
            ProfileCodeHotspotsUnitTests.writeTextFile(wrapperPath, [
                "function y = wrapper(x)"
                "y = slowChild(x);"
                "end"]);
            ProfileCodeHotspotsUnitTests.writeTextFile(computePath, [
                "function y = computeKernel(x, B)"
                "y = B \ x;"
                "end"]);

            functionTable = [
                ProfileCodeHotspotsUnitTests.functionEntry( ...
                    completeName="pkg.wrapper", functionName="wrapper", ...
                    fileName=wrapperPath, type="M-function", ...
                    totalTime=10, numCalls=1, children=struct("TotalTime", {9}), ...
                    parents=[], executedLines=[2 1 8], isRecursive=false)
                ProfileCodeHotspotsUnitTests.functionEntry( ...
                    completeName="pkg.computeKernel", functionName="computeKernel", ...
                    fileName=computePath, type="M-function", ...
                    totalTime=10, numCalls=1, children=[], ...
                    parents=[], executedLines=[2 1 8], isRecursive=false)
            ];
            profileInfo = struct("FunctionTable", functionTable);

            analysis = analyzeProfileInfo(profileInfo, projectRoots=projectRoot, maxFunctions=10, maxLines=10);

            wrapperLine = analysis.lineMetrics(analysis.lineMetrics.fileName == string(wrapperPath), :);
            computeLine = analysis.lineMetrics(analysis.lineMetrics.fileName == string(computePath), :);

            testCase.verifyEqual(wrapperLine.lineKind, "callBoundary")
            testCase.verifyTrue(wrapperLine.isCallBoundary)
            testCase.verifyFalse(wrapperLine.isActionableLine)
            testCase.verifyEqual(wrapperLine.calledTargetName, "slowChild")
            testCase.verifyEqual(computeLine.lineKind, "compute")
            testCase.verifyFalse(computeLine.isCallBoundary)
            testCase.verifyTrue(computeLine.isActionableLine)
            testCase.verifyEqual(analysis.topActionableLines.fileName(1), string(computePath))
            testCase.verifyEqual(analysis.topCallsiteLines.fileName(1), string(wrapperPath))

            computeIndex = find(analysis.priorityTargets.targetSubtype == "computeLine", 1, "first");
            callsiteIndex = find(analysis.priorityTargets.targetSubtype == "callsiteLine", 1, "first");
            testCase.verifyNotEmpty(computeIndex)
            testCase.verifyNotEmpty(callsiteIndex)
            testCase.verifyLessThan(computeIndex, callsiteIndex)
            testCase.verifyTrue(contains(analysis.priorityTargets.reason(callsiteIndex), "call site"))
            testCase.verifyTrue(contains(analysis.priorityTargets.reason(callsiteIndex), "slowChild"))
        end

        function analyzeProfileInfoClassifiesMultilineCallBoundaries(testCase)
            projectRoot = string(testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder);
            wrapperPath = fullfile(projectRoot, "multilineWrapper.m");
            ProfileCodeHotspotsUnitTests.writeTextFile(wrapperPath, [
                "function y = multilineWrapper(x)"
                "y = outer( ..."
                "    x);"
                "z = struct( ..."
                '    "value", inner(x));'
                "y = z;"
                "end"]);

            functionTable = ProfileCodeHotspotsUnitTests.functionEntry( ...
                completeName="pkg.multilineWrapper", functionName="multilineWrapper", ...
                fileName=wrapperPath, type="M-function", ...
                totalTime=10, numCalls=1, children=struct("TotalTime", {9}), ...
                parents=[], executedLines=[3 1 4; 5 1 3], isRecursive=false);
            profileInfo = struct("FunctionTable", functionTable);

            analysis = analyzeProfileInfo(profileInfo, projectRoots=projectRoot, maxFunctions=10, maxLines=10);

            callLines = analysis.lineMetrics(analysis.lineMetrics.fileName == string(wrapperPath), :);
            testCase.verifyEqual(callLines.lineKind, ["callBoundary"; "callBoundary"])
            testCase.verifyEqual(callLines.calledTargetName, ["outer"; "inner"])
            testCase.verifyTrue(all(~callLines.isActionableLine))
            testCase.verifyEqual(callLines.lineText(1), "y = outer( x);")
            testCase.verifyEqual(callLines.lineText(2), "z = struct( ""value"", inner(x));")
        end

        function analyzeProfileInfoClassifiesLowFractionConstructorCallsites(testCase)
            projectRoot = string(testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder);
            wrapperPath = fullfile(projectRoot, "constructorWrapper.m");
            ProfileCodeHotspotsUnitTests.writeTextFile(wrapperPath, [
                "function fit = constructorWrapper(x)"
                "fit = GriddedStreamfunction(x);"
                "end"]);

            functionTable = ProfileCodeHotspotsUnitTests.functionEntry( ...
                completeName="pkg.constructorWrapper", functionName="constructorWrapper", ...
                fileName=wrapperPath, type="M-function", ...
                totalTime=100, numCalls=1, children=struct("TotalTime", {96}), ...
                parents=[], executedLines=[2 1 4], isRecursive=false);
            profileInfo = struct("FunctionTable", functionTable);

            analysis = analyzeProfileInfo(profileInfo, projectRoots=projectRoot, maxFunctions=10, maxLines=10);

            wrapperLine = analysis.lineMetrics(analysis.lineMetrics.fileName == string(wrapperPath), :);
            testCase.verifyEqual(wrapperLine.lineKind, "callBoundary")
            testCase.verifyTrue(wrapperLine.isCallBoundary)
            testCase.verifyEqual(wrapperLine.calledTargetName, "GriddedStreamfunction")
            testCase.verifyTrue(any(analysis.topCallsiteLines.fileName == string(wrapperPath)))
        end

        function compareProfileHotspotsComputesDeltasAndContext(testCase)
            projectRoot = "/tmp/project";
            baseline = struct("FunctionTable", [
                ProfileCodeHotspotsUnitTests.functionEntry( ...
                    completeName="pkg.foo", functionName="foo", ...
                    fileName=projectRoot + "/foo.m", type="M-function", ...
                    totalTime=4, numCalls=2, children=struct("TotalTime", {1}), ...
                    parents=[], executedLines=[10 1 2], isRecursive=false)
                ProfileCodeHotspotsUnitTests.functionEntry( ...
                    completeName="pkg.bar", functionName="bar", ...
                    fileName=projectRoot + "/bar.m", type="M-function", ...
                    totalTime=2, numCalls=10, children=[], parents=[], ...
                    executedLines=[20 10 1], isRecursive=false)
            ]);
            candidate = struct("FunctionTable", [
                ProfileCodeHotspotsUnitTests.functionEntry( ...
                    completeName="pkg.foo", functionName="foo", ...
                    fileName=projectRoot + "/foo.m", type="M-function", ...
                    totalTime=5.5, numCalls=3, children=struct("TotalTime", {1.5}), ...
                    parents=[], executedLines=[10 1 3], isRecursive=false)
                ProfileCodeHotspotsUnitTests.functionEntry( ...
                    completeName="pkg.baz", functionName="baz", ...
                    fileName=projectRoot + "/baz.m", type="M-function", ...
                    totalTime=1, numCalls=1, children=[], parents=[], ...
                    executedLines=[30 1 0.5], isRecursive=false)
            ]);

            comparison = compareProfileHotspots( ...
                analyzeProfileInfo(baseline, projectRoots=projectRoot, label="before"), ...
                analyzeProfileInfo(candidate, projectRoots=projectRoot, label="after"));

            fooRow = comparison.functionDiffs(comparison.functionDiffs.completeName == "pkg.foo", :);
            barRow = comparison.functionDiffs(comparison.functionDiffs.completeName == "pkg.bar", :);
            bazRow = comparison.functionDiffs(comparison.functionDiffs.completeName == "pkg.baz", :);
            newFunctionRows = comparison.regressions( ...
                comparison.regressions.targetType == "function" & comparison.regressions.changeType == "new", :);
            removedFunctionRows = comparison.improvements( ...
                comparison.improvements.targetType == "function" & comparison.improvements.changeType == "removed", :);

            testCase.verifyEqual(comparison.baselineLabel, "before")
            testCase.verifyEqual(comparison.candidateLabel, "after")
            testCase.verifyEqual(fooRow.deltaTotalTime, 1.5, AbsTol=1e-12)
            testCase.verifyEqual(fooRow.deltaSelfTime, 1, AbsTol=1e-12)
            testCase.verifyEqual(fooRow.deltaNumCalls, 1)
            testCase.verifyEqual(barRow.changeType, "removed")
            testCase.verifyEqual(bazRow.changeType, "new")
            testCase.verifyTrue(all(ismember( ...
                ["completeName", "fileName", "functionName", "lineNumber", "changeType"], ...
                string(comparison.regressions.Properties.VariableNames))))
            testCase.verifyFalse(isempty(newFunctionRows))
            testCase.verifyEqual(newFunctionRows.targetName(1), "pkg.baz")
            testCase.verifyTrue(contains(newFunctionRows.reason(1), "newly executed"))
            testCase.verifyFalse(isempty(removedFunctionRows))
            testCase.verifyEqual(removedFunctionRows.targetName(1), "pkg.bar")
            testCase.verifyTrue(contains(removedFunctionRows.reason(1), "no longer executed"))
            testCase.verifyNotEmpty(comparison.summary)
        end

        function profileCodeHotspotsRunsRealProfilerIntegration(testCase)
            analysis = profileCodeHotspots( ...
                @() ProfileCodeHotspotsUnitTests.exerciseToyWorkload(), ...
                projectRoots=string(testCase.toolRoot), ...
                shouldPrintReport=false, ...
                label="toy workload");

            testCase.verifyEqual(analysis.label, "toy workload")
            testCase.verifyGreaterThan(height(analysis.functionMetrics), 0)
            testCase.verifyTrue(any(analysis.functionMetrics.isProjectCode))
            testCase.verifyGreaterThanOrEqual(sum(isfinite(analysis.functionMetrics.selfTime)), 1)
            testCase.verifyGreaterThanOrEqual(height(analysis.lineMetrics), 1)
            testCase.verifyGreaterThan(height(analysis.topProjectBySelfTime), 0)
            testCase.verifyGreaterThan(height(analysis.topActionableLines), 0)
            testCase.verifyNotEqual(analysis.topActionableLines.fileName(1), ...
                fullfile(testCase.toolRoot, "profileCodeHotspots.m"))
            testCase.verifyNotEmpty(analysis.summary)
        end
    end

    methods (Static, Access = private)
        function entry = functionEntry(options)
            arguments
                options.completeName (1,1) string
                options.functionName (1,1) string
                options.fileName (1,1) string
                options.type (1,1) string
                options.totalTime (1,1) double
                options.numCalls (1,1) double
                options.children = []
                options.parents = []
                options.executedLines = []
                options.isRecursive (1,1) logical = false
            end

            entry = struct( ...
                "CompleteName", char(options.completeName), ...
                "FunctionName", char(options.functionName), ...
                "FileName", char(options.fileName), ...
                "Type", char(options.type), ...
                "Children", options.children, ...
                "Parents", options.parents, ...
                "ExecutedLines", options.executedLines, ...
                "IsRecursive", options.isRecursive, ...
                "TotalRecursiveTime", options.totalTime, ...
                "PartialData", [], ...
                "NumCalls", options.numCalls, ...
                "TotalTime", options.totalTime);
        end

        function writeTextFile(filePath, lines)
            fid = fopen(filePath, "w");
            cleanupFile = onCleanup(@() fclose(fid));
            for iLine = 1:numel(lines)
                fprintf(fid, "%s\n", lines(iLine));
            end
        end

        function exerciseToyWorkload()
            for iOuter = 1:8
                ProfileCodeHotspotsUnitTests.parentWorkload(iOuter);
            end
        end

        function parentWorkload(scale)
            ProfileCodeHotspotsUnitTests.selfHeavyLeaf(scale);
            for iInner = 1:4
                ProfileCodeHotspotsUnitTests.childLeaf(scale + iInner);
            end
        end

        function selfHeavyLeaf(scale)
            x = 0;
            for i = 1:(150*scale)
                x = x + sum(sin(1:12));
            end
            if ~isfinite(x)
                error("ProfileCodeHotspotsUnitTests:UnexpectedNegativeValue", ...
                    "Synthetic workload accumulator should stay finite.");
            end
        end

        function childLeaf(scale)
            y = 0;
            for i = 1:(40*scale)
                y = y + sum(cos(1:8));
            end
            if ~isfinite(y)
                error("ProfileCodeHotspotsUnitTests:UnexpectedNegativeValue", ...
                    "Synthetic workload accumulator should stay finite.");
            end
        end
    end
end
