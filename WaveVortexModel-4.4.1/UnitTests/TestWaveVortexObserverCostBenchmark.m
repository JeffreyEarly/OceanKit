classdef TestWaveVortexObserverCostBenchmark < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addBenchmarkPath(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,"Benchmarks")));
        end
    end

    methods (Test,TestTags="smoke")
        function selectionIsValidated(testCase)
            testCase.verifyError(@()runWaveVortexObserverCostBenchmark(caseIds="unknown"),"WaveVortexBenchmark:UnknownObserverCostCase")
            testCase.verifyError(@()runWaveVortexObserverCostBenchmark(caseIds=["coefficient-endpoint" "coefficient-endpoint"]),"WaveVortexBenchmark:DuplicateObserverCostCase")
        end
    end

    methods (Test,TestTags="full")
        function reducedMatrixDeliversMatchedWorkAndGraphs(testCase)
            results = runWaveVortexObserverCostBenchmark(Nxyz=[8 6 5],deltaT=1/1024,integrationStepCount=2,denseOutputPointsPerStep=1,shouldUseFreshProcess=false);
            cases = results.cases;
            work = [cases.work];
            testCase.verifyEqual(string({cases.id}),["coefficient-endpoint" "coefficient-dense-output" "integrated-observer-endpoint" "composite-dense-output"])
            testCase.verifyEqual([cases.integratesObserverState],[false false true true])
            testCase.verifyEqual([cases.usesDenseOutput],[false true false true])
            testCase.verifyEqual([work.acceptedStepCount],2*ones(1,4))
            testCase.verifyEqual([work.rhsEvaluationCount],9*ones(1,4))
            testCase.verifyEqual(vertcat(work.outputRecordCounts),[2 0 0 0; 2 3 0 0; 2 0 0 0; 2 3 3 3])
            testCase.verifyLessThan(cases(1).work.integratedStateElementCount,cases(3).work.integratedStateElementCount)
            testCase.verifyEqual(cases(1).finalState,cases(2).finalState)
            testCase.verifyEqual(cases(3).finalState,cases(4).finalState)
        end
    end
end
