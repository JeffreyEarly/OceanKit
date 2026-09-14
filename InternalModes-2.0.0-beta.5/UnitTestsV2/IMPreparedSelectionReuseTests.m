classdef IMPreparedSelectionReuseTests < matlab.unittest.TestCase
    methods (Test)
        function preparedDiagnosticsPreserveSingularAndRegularPencils(testCase)
            for q = [0 1]
                evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=q,r=1,surfaceBoundary=IMBoundaryCondition.neumann(),bottomBoundary=IMBoundaryCondition.neumann());
                solver = IMSolverSpectral(nEVP=32).configuredForEVP(evp);
                [A,~,samples] = evp.assembleConfigured(solver);
                expected = evp.modeSelectionDiagnostics(solver,A);
                profile clear
                profile on
                cleanup = onCleanup(@() profile('off'));
                actual = evp.preparedModeSelectionDiagnostics(samples,A);
                info = profile('info');
                profile off
                clear cleanup
                testCase.verifyEqual(actual,expected);
                functions = info.FunctionTable;
                index = contains(string({functions.FunctionName}),'zeroModeAssessment');
                testCase.verifyEqual(sum([functions(index).NumCalls]),1);
            end
        end
    end
end
