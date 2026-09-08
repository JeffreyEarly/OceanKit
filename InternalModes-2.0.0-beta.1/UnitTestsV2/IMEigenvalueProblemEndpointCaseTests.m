classdef IMEigenvalueProblemEndpointCaseTests < matlab.unittest.TestCase

    properties (Access = private)
        originalPath
    end

    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            repoRoot = fileparts(fileparts(mfilename('fullpath')));
            testCase.originalPath = path;
            addpath(repoRoot);
        end
    end

    methods (TestClassTeardown)
        function restoreRepositoryPath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function endpointAlgebraUsesYassinIndexing(testCase)
            surface = IMBoundaryCondition(a=0, b=1, c=1, d=0);
            bottom = IMBoundaryCondition(a=0, b=1, c=1, d=0);

            testCase.verifyTrue(surface.isEigenvalueDependent())
            testCase.verifyEqual(surface.determinant("surface"), 1, AbsTol=0)
            testCase.verifyEqual(surface.endpointWeightCoefficient("surface"), 1, AbsTol=0)
            testCase.verifyEqual(bottom.determinant("bottom"), -1, AbsTol=0)
            testCase.verifyEqual(bottom.endpointWeightCoefficient("bottom"), -1, AbsTol=0)
        end

        function boundaryConditionValidationUsesBuiltInArgumentValidation(testCase)
            boundary = IMBoundaryCondition(a=1, b=1, c=1, d=0);
            badCoefficientThrows = false;
            badDeterminantLocationThrows = false;
            badWeightLocationThrows = false;
            badRobinLocationThrows = false;
            badToleranceThrows = false;
            badRobinCoefficientThrows = false;

            try
                IMBoundaryCondition(a=Inf);
            catch
                badCoefficientThrows = true;
            end
            try
                boundary.determinant("middle");
            catch
                badDeterminantLocationThrows = true;
            end
            try
                boundary.endpointWeightCoefficient("middle");
            catch
                badWeightLocationThrows = true;
            end
            try
                boundary.robinEnergyCoefficient("middle");
            catch
                badRobinLocationThrows = true;
            end
            try
                boundary.isEigenvalueDependent(Inf);
            catch
                badToleranceThrows = true;
            end
            try
                IMBoundaryCondition.robin(1, Inf);
            catch
                badRobinCoefficientThrows = true;
            end

            testCase.verifyTrue(badCoefficientThrows)
            testCase.verifyTrue(badDeterminantLocationThrows)
            testCase.verifyTrue(badWeightLocationThrows)
            testCase.verifyTrue(badRobinLocationThrows)
            testCase.verifyTrue(badToleranceThrows)
            testCase.verifyTrue(badRobinCoefficientThrows)
            testCase.verifyError(@() IMBoundaryCondition(a=0, b=0, c=0, d=0), ...
                "IMBoundaryCondition:DegenerateCondition")
        end

        function endpointDiagnosticsUseBuiltInArgumentValidation(testCase)
            evp = IMEigenvalueProblem();
            badLocationThrows = false;
            badToleranceThrows = false;

            try
                evp.endpointWeights("middle");
            catch
                badLocationThrows = true;
            end
            try
                evp.negativeEndpointWeightCount(tolerance=Inf);
            catch
                badToleranceThrows = true;
            end

            testCase.verifyTrue(badLocationThrows)
            testCase.verifyTrue(badToleranceThrows)
            testCase.verifyEmpty(evp.endpointWeights("surface"))
            testCase.verifyEmpty(evp.endpointWeights("bottom"))
        end

        function tinyEigenvalueSideCoefficientsAreInactive(testCase)
            boundary = IMBoundaryCondition(a=1, b=0, c=eps, d=0);
            evp = IMEigenvalueProblem(surfaceBoundary=boundary);

            testCase.verifyFalse(boundary.isEigenvalueDependent())
            testCase.verifyEmpty(evp.endpointWeights("surface"))
        end

        function inactiveRobinEnergyUsesYassinEndpointSigns(testCase)
            boundary = IMBoundaryCondition(a=2, b=3);

            testCase.verifyEqual(boundary.robinEnergyCoefficient("surface"), -2/3, AbsTol=0)
            testCase.verifyEqual(boundary.robinEnergyCoefficient("bottom"), 2/3, AbsTol=0)
        end

        function negativeEndpointWeightCountIsEndpointBounded(testCase)
            zeroCountEVP = IMEigenvalueProblem(p=1, q=0, r=1);
            oneCountEVP = IMEigenvalueProblem(p=1, q=0, r=1, ...
                bottomBoundary=IMBoundaryCondition(a=0, b=1, c=1, d=0));
            twoCountEVP = IMEigenvalueProblem(p=1, q=0, r=1, ...
                surfaceBoundary=IMBoundaryCondition(a=0, b=1, c=-1, d=0), ...
                bottomBoundary=IMBoundaryCondition(a=0, b=1, c=1, d=0));

            testCase.verifyEqual(zeroCountEVP.negativeEndpointWeightCount(), 0)
            testCase.verifyEqual(oneCountEVP.negativeEndpointWeightCount(), 1)
            testCase.verifyEqual(twoCountEVP.negativeEndpointWeightCount(), 2)
        end

        function case1PositiveMetricAndNonnegativeQSuppressesNegativeCandidates(testCase)
            solver = testCase.canonicalSolver();
            evp = IMEigenvalueProblem(p=1, q=0, r=1);
            [A, ~] = evp.assemble(solver);

            diagnostics = evp.definitenessDiagnostics(solver);
            bounds = evp.negativeEigenvalueBounds(solver, A);
            selection = evp.selectModes([-100; 1; 2; 3], 2, solver, A);

            testCase.verifyTrue(diagnostics.metricPositive)
            testCase.verifyTrue(diagnostics.quadraticFormNonnegative)
            testCase.verifyEqual(bounds.assessmentLevel, diagnostics.assessmentLevel)
            testCase.verifyEqual(bounds.minNegativeEigenvalueCount, 0)
            testCase.verifyEqual(bounds.maxNegativeEigenvalueCount, 0)
            testCase.verifyEqual(selection.modeSelectionDiagnostics.zeroModeStatus, "absent")
            testCase.verifyEqual(selection.sortIndex(:), [2; 3])
            testCase.verifyEqual(selection.modeNumber, [1 2])
        end

        function case2NegativeEndpointWeightCountWithZeroAbsentHasExactNegativeCount(testCase)
            solver = testCase.canonicalSolver();
            bottom = IMBoundaryCondition(a=0, b=1, c=1, d=0);
            evp = IMEigenvalueProblem(p=1, q=1, r=1, bottomBoundary=bottom);
            [A, ~] = evp.assemble(solver);

            diagnostics = evp.definitenessDiagnostics(solver);
            bounds = evp.negativeEigenvalueBounds(solver, A);
            removedZeroStatusField = "zero" + "EigenvalueStatus";

            testCase.verifyEqual(diagnostics.negativeEndpointWeightCount, 1)
            testCase.verifyEqual(bounds.negativeEndpointWeightCount, 1)
            testCase.verifyTrue(diagnostics.quadraticFormNonnegative)
            testCase.verifyEqual(bounds.zeroModeStatus, "absent")
            testCase.verifyFalse(isfield(bounds, removedZeroStatusField))
            testCase.verifyEqual(bounds.minNegativeEigenvalueCount, 1)
            testCase.verifyEqual(bounds.maxNegativeEigenvalueCount, 1)
        end

        function case2WithoutLeftMatrixKeepsNegativeEndpointWeightCountAsUpperBound(testCase)
            solver = testCase.canonicalSolver();
            bottom = IMBoundaryCondition(a=0, b=1, c=1, d=0);
            evp = IMEigenvalueProblem(p=1, q=1, r=1, bottomBoundary=bottom);

            bounds = evp.negativeEigenvalueBounds(solver);

            testCase.verifyEqual(bounds.negativeEndpointWeightCount, 1)
            testCase.verifyEqual(bounds.zeroModeStatus, "unchecked")
            testCase.verifyEqual(bounds.minNegativeEigenvalueCount, 0)
            testCase.verifyEqual(bounds.maxNegativeEigenvalueCount, 1)
        end

        function modeSelectionDiagnosticsInferZeroModeFromLeftMatrix(testCase)
            solver = testCase.canonicalSolver();
            zeroEVP = IMEigenvalueProblem(p=1, q=0, r=1, ...
                surfaceBoundary=IMBoundaryCondition.neumann(), ...
                bottomBoundary=IMBoundaryCondition.neumann());
            constrainedEVP = IMEigenvalueProblem(p=1, q=0, r=1);
            [zeroA, ~] = zeroEVP.assemble(solver);
            [constrainedA, ~] = constrainedEVP.assemble(solver);

            present = zeroEVP.modeSelectionDiagnostics(solver, zeroA);
            absent = constrainedEVP.modeSelectionDiagnostics(solver, constrainedA);

            testCase.verifyEqual(present.zeroModeStatus, "present")
            testCase.verifyEqual(present.zeroModeCount, 1)
            testCase.verifyLessThanOrEqual(present.zeroModeSingularValue, present.zeroModeTolerance)
            testCase.verifyEqual(absent.zeroModeStatus, "absent")
            testCase.verifyEqual(absent.zeroModeCount, 0)
            testCase.verifyGreaterThan(absent.zeroModeSingularValue, absent.zeroModeTolerance)
        end

        function selectModesIncludesZeroOnlyWhenDiagnosticsReportPresent(testCase)
            solver = testCase.canonicalSolver();
            zeroEVP = IMEigenvalueProblem(p=1, q=0, r=1, ...
                surfaceBoundary=IMBoundaryCondition.neumann(), ...
                bottomBoundary=IMBoundaryCondition.neumann());
            constrainedEVP = IMEigenvalueProblem(p=1, q=0, r=1);
            [zeroA, ~] = zeroEVP.assemble(solver);
            [constrainedA, ~] = constrainedEVP.assemble(solver);

            zeroSelection = zeroEVP.selectModes([-100; 0; 1; 2], 3, solver, zeroA);
            constrainedSelection = constrainedEVP.selectModes([0; 1; 2], 2, solver, constrainedA);

            testCase.verifyEqual(zeroSelection.sortIndex(:), [2; 3; 4])
            testCase.verifyEqual(zeroSelection.modeNumber, [0 1 2])
            testCase.verifyEqual(zeroSelection.modeSelectionDiagnostics.zeroModeStatus, "present")
            testCase.verifyEqual(constrainedSelection.sortIndex(:), [2; 3])
            testCase.verifyEqual(constrainedSelection.modeNumber, [1 2])
            testCase.verifyEqual(constrainedSelection.modeSelectionDiagnostics.zeroModeStatus, "absent")
        end

        function case3InactiveRobinDirectionsSetFiniteSearchCap(testCase)
            solver = testCase.canonicalSolver();
            surface = IMBoundaryCondition(a=1, b=1);
            bottom = IMBoundaryCondition(a=-1, b=1);
            evp = IMEigenvalueProblem(p=1, q=0, r=1, surfaceBoundary=surface, bottomBoundary=bottom);

            diagnostics = evp.definitenessDiagnostics(solver);
            bounds = evp.negativeEigenvalueBounds(solver);

            testCase.verifyTrue(diagnostics.metricPositive)
            testCase.verifyEqual(diagnostics.endpointNumeratorNegativeDirections, 2)
            testCase.verifyEqual(bounds.minNegativeEigenvalueCount, 0)
            testCase.verifyEqual(bounds.maxNegativeEigenvalueCount, 2)
        end

        function case3ActiveEndpointInertiaSetsFiniteSearchCap(testCase)
            solver = testCase.canonicalSolver();
            surface = IMBoundaryCondition(a=0, b=1, c=1, d=1);
            evp = IMEigenvalueProblem(p=1, q=0, r=1, surfaceBoundary=surface);

            diagnostics = evp.definitenessDiagnostics(solver);
            bounds = evp.negativeEigenvalueBounds(solver);

            testCase.verifyTrue(diagnostics.metricPositive)
            testCase.verifyEqual(diagnostics.endpointNumeratorNegativeDirections, 1)
            testCase.verifyEqual(bounds.minNegativeEigenvalueCount, 0)
            testCase.verifyEqual(bounds.maxNegativeEigenvalueCount, 1)
        end

        function case3GeostrophicWKBRetainsAssessedNegativeMode(testCase)
            [solver, evp] = testCase.geostrophicWKBProblem();
            [A, B] = evp.assemble(solver);
            eigenvalues = eig(A, B);
            valid = isfinite(real(eigenvalues)) & isfinite(imag(eigenvalues)) ...
                & abs(imag(eigenvalues)) < 1e-8*max(1,abs(real(eigenvalues)));
            finiteRealEigenvalues = real(eigenvalues(valid));

            diagnostics = evp.definitenessDiagnostics(solver);
            bounds = evp.negativeEigenvalueBounds(solver, A);
            basisSet = solver.solveEVP(evp, nModes=4);

            testCase.verifyGreaterThanOrEqual(nnz(finiteRealEigenvalues < 0), 1)
            testCase.verifyEqual(diagnostics.endpointNumeratorNegativeDirections, 2)
            testCase.verifyEqual(bounds.maxNegativeEigenvalueCount, 2)
            testCase.verifyEqual(bounds.zeroModeStatus, "absent")
            testCase.verifyLessThan(basisSet.eigenvalues(1), 0)
            testCase.verifyEqual(basisSet.modeNumber(1), -1)
        end

        function hydrostaticFModesInferBarotropicZeroMode(testCase)
            zDomain = [-1000 0];
            N0 = 5.2e-3;
            N2 = @(z) N0*N0*(1 + 0*z);
            solver = IMSolverSpectral(nEVP=32);
            evp = IMInternalModes.hydrostaticFModes(N2=N2, zDomain=zDomain);
            removedZeroProperty = "has" + "ZeroMode";

            basisSet = solver.solveEVP(evp, nModes=3);

            testCase.verifyFalse(isprop(evp, removedZeroProperty))
            testCase.verifyEqual(basisSet.modeNumber(1), 0)
            testCase.verifyEqual(basisSet.modeSelectionDiagnostics.zeroModeStatus, "present")
            testCase.verifyEqual(basisSet.modeSelectionDiagnostics.zeroModeCount, 1)
        end

        function robinNeumannGeostrophicWKBDoesNotInventZeroMode(testCase)
            D = 4000;
            N0 = 5.2e-3;
            b = 1300;
            g = 9.81;
            zDomain = [-D 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            p = @(z,ctx) 1 ./ ctx.N2(z);
            q = @(z,~) zeros(size(z));
            r = @(z,ctx) ones(size(z))/ctx.g;
            surfaceBoundary = IMBoundaryCondition(a=-1/g, b=1);
            bottomBoundary = IMBoundaryCondition.neumann();
            evp = IMInternalModes(name="unforced-APV-modes", formulation="F", N2=N2, zDomain=zDomain, ...
                p=p, q=q, r=r, g=g, surfaceBoundary=surfaceBoundary, bottomBoundary=bottomBoundary);
            solver = IMSolverSpectral(nEVP=128, coordinateKind="wkb");

            basisSet = solver.solveEVP(evp, nModes=4);

            testCase.verifyEqual(basisSet.modeSelectionDiagnostics.zeroModeStatus, "absent")
            testCase.verifyFalse(any(basisSet.modeNumber == 0))
            testCase.verifyEqual(basisSet.modeNumber, 1:4)
            testCase.verifyGreaterThan(min(diff(basisSet.eigenvalues)), 0)
        end

        function case4FailedCoefficientSignsReturnUnknown(testCase)
            solver = testCase.canonicalSolver();
            evp = IMEigenvalueProblem(p=-1, q=0, r=1);

            bounds = evp.negativeEigenvalueBounds(solver);

            testCase.verifyEqual(bounds.minNegativeEigenvalueCount, 0)
            testCase.verifyEqual(bounds.maxNegativeEigenvalueCount, "unknown")
            testCase.verifyTrue(contains(bounds.reason, "fail"))
        end

        function case4NonfiniteSamplesReturnUnknown(testCase)
            solver = testCase.canonicalSolver();
            evp = IMEigenvalueProblem(p=@(z,~) NaN(size(z)), q=0, r=1);

            bounds = evp.negativeEigenvalueBounds(solver);

            testCase.verifyEqual(bounds.minNegativeEigenvalueCount, 0)
            testCase.verifyEqual(bounds.maxNegativeEigenvalueCount, "unknown")
            testCase.verifyTrue(contains(bounds.reason, "nonfinite"))
        end

        function case4DegenerateActiveEndpointReturnsUnknown(testCase)
            solver = testCase.canonicalSolver();
            surface = IMBoundaryCondition(a=0, b=1, c=0, d=1);
            evp = IMEigenvalueProblem(p=1, q=0, r=1, surfaceBoundary=surface);

            diagnostics = evp.definitenessDiagnostics(solver);
            bounds = evp.negativeEigenvalueBounds(solver);

            testCase.verifyTrue(diagnostics.hasDegenerateEndpointMetric)
            testCase.verifyEmpty(evp.endpointWeights("surface"))
            testCase.verifyEqual(bounds.minNegativeEigenvalueCount, 0)
            testCase.verifyEqual(bounds.maxNegativeEigenvalueCount, "unknown")
            testCase.verifyTrue(contains(bounds.reason, "degenerate"))
        end
    end

    methods (Static, Access = private)
        function solver = canonicalSolver()
            solver = IMSolverSpectral(nEVP=24);
        end

        function [solver, evp] = geostrophicWKBProblem()
            D = 4000;
            N0 = 5.2e-3;
            b = 1300;
            g = 9.81;
            g0 = -N0*N0*b;
            gd = -N0*N0*b/10;
            zDomain = [-D 0];
            N2 = @(z) N0*N0*exp(2*z/b);

            p = @(z,ctx) 1 ./ ctx.N2(z);
            q = @(z,~) zeros(size(z));
            r = @(z,ctx) ones(size(z))/ctx.g;
            surfaceBoundary = IMBoundaryCondition(a=-(1/g + 1/g0), b=1);
            bottomBoundary = IMBoundaryCondition(a=1/gd, b=1);

            evp = IMInternalModes(name="unforced-APV-modes", formulation="F", N2=N2, zDomain=zDomain, ...
                p=p, q=q, r=r, g=g, surfaceBoundary=surfaceBoundary, bottomBoundary=bottomBoundary);
            solver = IMSolverSpectral(nEVP=128, coordinateKind="wkb");
        end
    end
end
