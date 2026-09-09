classdef IMHydrostaticBoundaryConditionTests < matlab.unittest.TestCase

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
        function constructionStoresPhysicalCoefficients(testCase)
            law = IMHydrostaticBoundaryCondition(a=1, b=2, c=3, d=4, e=5);

            testCase.verifyEqual(law.a, 1)
            testCase.verifyEqual(law.b, 2)
            testCase.verifyEqual(law.c, 3)
            testCase.verifyEqual(law.d, 4)
            testCase.verifyEqual(law.e, 5)
        end

        function validationRejectsInvalidCoefficientSets(testCase)
            badCoefficientThrows = false;
            try
                IMHydrostaticBoundaryCondition(a=Inf);
            catch
                badCoefficientThrows = true;
            end

            testCase.verifyTrue(badCoefficientThrows)
            testCase.verifyError(@() IMHydrostaticBoundaryCondition(), "IMHydrostaticBoundaryCondition:DegenerateLaw")
        end

        function fFormulationConversionScalesCanonicalCoefficients(testCase)
            g = 10;
            law = IMHydrostaticBoundaryCondition(a=2, b=3, c=4, d=5);
            boundary = law.canonicalBoundary(formulation="F", g=g);

            testCase.verifyBoundaryCoefficients(boundary, -2, 3, -4/g, 5/g)
        end

        function gFormulationConversionScalesCanonicalCoefficients(testCase)
            g = 10;
            law = IMHydrostaticBoundaryCondition(a=2, b=3, c=4, e=5);
            boundary = law.canonicalBoundary(formulation="G", g=g);

            testCase.verifyBoundaryCoefficients(boundary, 5, 2, 3/g, 4/g)
        end

        function nonlinearSlicesThrowClearError(testCase)
            fLaw = IMHydrostaticBoundaryCondition(a=1, e=2);
            gLaw = IMHydrostaticBoundaryCondition(a=1, d=2);

            testCase.verifyError(@() fLaw.canonicalBoundary(formulation="F", g=9.81), "IMHydrostaticBoundaryCondition:NonlinearBoundaryLaw")
            testCase.verifyError(@() gLaw.canonicalBoundary(formulation="G", g=9.81), "IMHydrostaticBoundaryCondition:NonlinearBoundaryLaw")
        end

        function convertedFBoundariesMatchManualCanonicalBoundaries(testCase)
            g = 9.81;
            testCase.verifyConvertedBoundaryMatchesManual("F", IMHydrostaticBoundaryCondition(a=2, b=3), IMBoundaryCondition(a=-2, b=3), g)
            testCase.verifyConvertedBoundaryMatchesManual("F", IMHydrostaticBoundaryCondition(b=2, c=3), IMBoundaryCondition(a=0, b=2, c=-3/g, d=0), g)
            testCase.verifyConvertedBoundaryMatchesManual("F", IMHydrostaticBoundaryCondition(b=2, c=3, d=4), IMBoundaryCondition(a=0, b=2, c=-3/g, d=4/g), g)
        end

        function convertedGBoundariesMatchManualCanonicalBoundaries(testCase)
            g = 9.81;
            testCase.verifyConvertedBoundaryMatchesManual("G", IMHydrostaticBoundaryCondition(a=2, b=3), IMBoundaryCondition(a=0, b=2, c=3/g, d=0), g)
            testCase.verifyConvertedBoundaryMatchesManual("G", IMHydrostaticBoundaryCondition(a=2, b=3, c=4), IMBoundaryCondition(a=0, b=2, c=3/g, d=4/g), g)
            testCase.verifyConvertedBoundaryMatchesManual("G", IMHydrostaticBoundaryCondition(a=2, e=5), IMBoundaryCondition(a=5, b=2, c=0, d=0), g)
        end
    end

    methods
        function verifyBoundaryCoefficients(testCase, boundary, a, b, c, d)
            testCase.verifyEqual(boundary.a, a, AbsTol=0)
            testCase.verifyEqual(boundary.b, b, AbsTol=0)
            testCase.verifyEqual(boundary.c, c, RelTol=1e-14, AbsTol=1e-14)
            testCase.verifyEqual(boundary.d, d, RelTol=1e-14, AbsTol=1e-14)
        end

        function verifyConvertedBoundaryMatchesManual(testCase, formulation, law, manualBoundary, g)
            [N2, zDomain] = testCase.profile();
            convertedBoundary = law.canonicalBoundary(formulation=formulation, g=g);
            convertedEVP = testCase.evpFor(formulation, convertedBoundary, N2, zDomain, g);
            manualEVP = testCase.evpFor(formulation, manualBoundary, N2, zDomain, g);

            testCase.verifyBoundaryCoefficients(convertedBoundary, manualBoundary.a, manualBoundary.b, manualBoundary.c, manualBoundary.d)
            testCase.verifyInnerProductSpecsEqual(convertedEVP.innerProduct("F"), manualEVP.innerProduct("F"))
            testCase.verifyInnerProductSpecsEqual(convertedEVP.innerProduct("G"), manualEVP.innerProduct("G"))
        end

        function verifyInnerProductSpecsEqual(testCase, actual, expected)
            testCase.verifyEqual(actual.variable, expected.variable)
            testCase.verifyEqual(actual.hasInnerProduct, expected.hasInnerProduct)
            testCase.verifyEndpointWeightsEqual(actual.surfaceWeights, expected.surfaceWeights)
            testCase.verifyEndpointWeightsEqual(actual.bottomWeights, expected.bottomWeights)
            testCase.verifyEndpointTermsEqual(actual.endpointInnerProductTerms, expected.endpointInnerProductTerms)
        end

        function verifyEndpointWeightsEqual(testCase, actual, expected)
            testCase.verifyEqual(numel(actual), numel(expected))
            for iWeight = 1:numel(actual)
                testCase.verifyEqual(actual(iWeight).location, expected(iWeight).location)
                testCase.verifyEqual(actual(iWeight).coefficient, expected(iWeight).coefficient, RelTol=1e-14, AbsTol=1e-14)
                testCase.verifyEqual(actual(iWeight).c, expected(iWeight).c, RelTol=1e-14, AbsTol=1e-14)
                testCase.verifyEqual(actual(iWeight).d, expected(iWeight).d, RelTol=1e-14, AbsTol=1e-14)
            end
        end

        function verifyEndpointTermsEqual(testCase, actual, expected)
            testCase.verifyEqual(numel(actual), numel(expected))
            for iTerm = 1:numel(actual)
                testCase.verifyEqual(actual(iTerm).location, expected(iTerm).location)
                testCase.verifyEqual(actual(iTerm).coefficient, expected(iTerm).coefficient, RelTol=1e-14, AbsTol=1e-14)
                testCase.verifyEqual(actual(iTerm).variable, expected(iTerm).variable)
                testCase.verifyEqual(actual(iTerm).catalogCase, expected(iTerm).catalogCase)
            end
        end
    end

    methods (Static, Access = private)
        function evp = evpFor(formulation, boundary, N2, zDomain, g)
            if string(formulation) == "F"
                evp = IMInternalModes.hydrostaticFModes(N2=N2, zDomain=zDomain, g=g, surfaceBoundary=boundary);
            else
                evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, g=g, surfaceBoundary=boundary);
            end
        end

        function [N2, zDomain] = profile()
            zDomain = [-1000 0];
            N0 = 5.2e-3;
            N2 = @(z) N0*N0*(1 + 0.1*z/abs(zDomain(1)));
        end
    end
end
