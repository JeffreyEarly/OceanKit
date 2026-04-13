classdef DistributionsUnitTests < matlab.unittest.TestCase

    methods (Test)
        function normalDistributionDefaultsToUnitSigma(testCase)
            distribution = NormalDistribution();

            testCase.verifyEqual(distribution.sigma, 1);
            testCase.verifyEqual(distribution.variance, 1);
        end

        function rayleighDistributionDefaultsToUnitSigma(testCase)
            distribution = RayleighDistribution();

            testCase.verifyEqual(distribution.sigma, 1);
            testCase.verifyEqual(distribution.variance, (4-pi)/2, AbsTol=1e-12);
        end

        function uniformDistributionDefaultsToCenteredUnitInterval(testCase)
            distribution = UniformDistribution();

            testCase.verifyEqual(distribution.a, -0.5);
            testCase.verifyEqual(distribution.b, 0.5);
        end

        function uniformDistributionSupportsPartialEndpointDefaults(testCase)
            leftDefault = UniformDistribution(a=-2);
            rightDefault = UniformDistribution(b=3);

            testCase.verifyEqual(leftDefault.a, -2);
            testCase.verifyEqual(leftDefault.b, 0.5);
            testCase.verifyEqual(rightDefault.a, -0.5);
            testCase.verifyEqual(rightDefault.b, 3);
        end

        function normalDistributionRoundTripsThroughAnnotatedPersistence(testCase)
            distribution = NormalDistribution(sigma=1.75);

            testCase.roundTripDistribution(distribution, @(path) NormalDistribution.annotatedClassFromFile(path), [-3; 0; 2]);
        end

        function rayleighDistributionRoundTripsThroughAnnotatedPersistence(testCase)
            distribution = RayleighDistribution(sigma=2.25);

            testCase.roundTripDistribution(distribution, @(path) RayleighDistribution.annotatedClassFromFile(path), [0; 1; 4]);
        end

        function chi2DistributionRoundTripsThroughAnnotatedPersistence(testCase)
            distribution = Chi2Distribution(k=5);

            testCase.roundTripDistribution(distribution, @(path) Chi2Distribution.annotatedClassFromFile(path), [0.5; 3; 8]);
        end

        function uniformDistributionRoundTripsThroughAnnotatedPersistence(testCase)
            distribution = UniformDistribution(a=-2.5, b=0.75);

            testCase.roundTripDistribution(distribution, @(path) UniformDistribution.annotatedClassFromFile(path), [-2; -0.5; 0.5]);
        end

        function studentTDistributionRoundTripsWithSigmaConstructor(testCase)
            distribution = StudentTDistribution(nu=6, sigma=1.25);

            testCase.roundTripDistribution(distribution, @(path) StudentTDistribution.annotatedClassFromFile(path), [-2; 0; 2]);
        end

        function studentTDistributionPersistenceCanonicalizesSigma(testCase)
            distribution = StudentTDistribution(nu=8, variance=3.5);
            distributionRead = testCase.roundTripDistribution(distribution, @(path) StudentTDistribution.annotatedClassFromFile(path), [-3; 0; 3]);

            testCase.verifyEqual(distributionRead.sigma, distribution.sigma, AbsTol=1e-12);
        end

        function twoDimDistanceDistributionRoundTripPreservesNestedComponent(testCase)
            distribution = TwoDimDistanceDistribution(distribution1d=NormalDistribution(sigma=0.85));
            distributionRead = testCase.roundTripDistribution(distribution, @(path) TwoDimDistanceDistribution.annotatedClassFromFile(path), [0; 0.5; 1.5; 3]);

            testCase.verifyClass(distributionRead.distribution1d, class(distribution.distribution1d));
        end

        function addedDistributionRoundTripPreservesWeightsAndOrder(testCase)
            distribution = AddedDistribution(scalings=[0.2 0.5 0.3], distributions=[NormalDistribution(sigma=0.75); StudentTDistribution(nu=4.5, sigma=1.1); UniformDistribution(a=-1.0, b=1.5)]);
            distributionRead = testCase.roundTripDistribution(distribution, @(path) AddedDistribution.annotatedClassFromFile(path), [-2; 0; 2]);

            testCase.verifySize(distributionRead.distributions, [3 1]);
            testCase.verifyEqual(distributionRead.scalings, distribution.scalings);
            testCase.verifyClass(distributionRead.distributions(1), class(distribution.distributions(1)));
            testCase.verifyClass(distributionRead.distributions(2), class(distribution.distributions(2)));
            testCase.verifyClass(distributionRead.distributions(3), class(distribution.distributions(3)));
        end

        function rhoRoundTripsWhenPresent(testCase)
            distribution = NormalDistribution(sigma=0.9);
            distribution.rho = @(t) exp(-(t/3).^2);
            distributionRead = testCase.roundTripDistribution(distribution, @(path) NormalDistribution.annotatedClassFromFile(path), [-1.5; 0; 1.5]);

            testCase.verifyEqual(string(func2str(distributionRead.rho)), string(func2str(distribution.rho)));
        end

        function warningOnlyRhoHandlesStillRoundTrip(testCase)
            distribution = NormalDistribution(sigma=1.0);
            cleanupToken = onCleanup(@() []);
            distribution.rho = @(t) cleanupToken;
            path = strcat(tempname, '.nc');
            cleanup = onCleanup(@() DistributionsUnitTests.deleteIfPresent(path)); %#ok<NASGU>
            lastwarn('', '');
            distribution.writeToFile(path, shouldOverwriteExisting=true);
            [~, warningIdentifier] = lastwarn;
            distributionRead = NormalDistribution.annotatedClassFromFile(path);
            distributionReadAbstract = Distribution.distributionFromFile(path);

            testCase.verifyTrue(isempty(warningIdentifier) || strcmp(warningIdentifier, 'NetCDFGroup:PotentiallyUnserializableFunctionHandle'));
            testCase.verifyClass(distributionRead, 'NormalDistribution');
            testCase.verifyClass(distributionReadAbstract, 'NormalDistribution');
            testCase.assertDistributionEquivalent(distribution, distributionRead, [-1; 0; 1]);
            testCase.assertDistributionEquivalent(distribution, distributionReadAbstract, [-1; 0; 1]);
            testCase.verifyEqual(string(func2str(distributionRead.rho)), string(func2str(distribution.rho)));
        end

        function normalTwoDimDistanceMatchesRayleighWithinNumericalTolerance(testCase)
            sigma = 1;
            normal = NormalDistribution(sigma=sigma);
            radialDistribution = TwoDimDistanceDistribution(distribution1d=normal);
            rayleigh = RayleighDistribution(sigma=sigma);
            z = linspace(0, 5, 1000)';

            pdfError = max(abs(radialDistribution.pdf(z) - rayleigh.pdf(z)));
            cdfError = max(abs(radialDistribution.cdf(z) - rayleigh.cdf(z)));

            testCase.verifyLessThanOrEqual(pdfError, 0.012);
            testCase.verifyLessThanOrEqual(cdfError, 4e-4);
        end
    end

    methods (Access = private)
        function distributionRead = roundTripDistribution(testCase, distribution, concreteReader, queryPoints)
            arguments
                testCase (1,1) matlab.unittest.TestCase
                distribution (1,1) Distribution
                concreteReader (1,1) function_handle
                queryPoints {mustBeNumeric, mustBeReal}
            end
            path = strcat(tempname, '.nc');
            cleanup = onCleanup(@() DistributionsUnitTests.deleteIfPresent(path)); %#ok<NASGU>
            distribution.writeToFile(path, shouldOverwriteExisting=true);

            distributionRead = concreteReader(path);
            distributionReadAbstract = Distribution.distributionFromFile(path);

            testCase.verifyClass(distributionRead, class(distribution));
            testCase.verifyClass(distributionReadAbstract, class(distribution));
            testCase.assertDistributionEquivalent(distribution, distributionRead, queryPoints);
            testCase.assertDistributionEquivalent(distribution, distributionReadAbstract, queryPoints);
        end

        function assertDistributionEquivalent(testCase, original, restored, queryPoints)
            pdfOriginal = arrayfun(@(z) original.pdf(z), queryPoints);
            pdfRestored = arrayfun(@(z) restored.pdf(z), queryPoints);
            cdfOriginal = arrayfun(@(z) original.cdf(z), queryPoints);
            cdfRestored = arrayfun(@(z) restored.cdf(z), queryPoints);

            pdfError = max(abs(pdfOriginal - pdfRestored));
            cdfError = max(abs(cdfOriginal - cdfRestored));
            varianceScale = max([1; abs(original.variance); abs(restored.variance)]);

            testCase.verifyTrue(original.isequal(restored));
            testCase.verifyLessThan(pdfError, 1e-10);
            testCase.verifyLessThan(cdfError, 1e-10);
            testCase.verifyLessThan(abs(original.variance - restored.variance), 1e-10*varianceScale);

            if isempty(original.rho)
                testCase.verifyEmpty(restored.rho);
            else
                testCase.verifyEqual(string(func2str(original.rho)), string(func2str(restored.rho)));
            end
        end
    end

    methods (Static, Access = private)
        function deleteIfPresent(path)
            if isfile(path)
                delete(path);
            end
        end
    end

end
