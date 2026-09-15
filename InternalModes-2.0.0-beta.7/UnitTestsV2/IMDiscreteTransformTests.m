classdef IMDiscreteTransformTests < matlab.unittest.TestCase

    properties (Access = private)
        originalPath
    end

    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
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
        function scalarBuilderStoresExpectedGalerkinMatrices(testCase)
            [basisSet, z, dz] = testCase.regularBasis(4);
            transform = basisSet.discreteTransform(z=z, weights=dz, nModes=3);

            expectedBasis = basisSet.u(z);
            expectedBasis = expectedBasis(:,1:3);
            expectedMetric = diag(dz);
            expectedGram = expectedBasis.'*expectedMetric*expectedBasis;
            expectedGram = 0.5*(expectedGram + expectedGram.');
            continuousGram = basisSet.gramMatrix();
            expectedTarget = diag(diag(continuousGram(1:3,1:3)));
            scale = 1./sqrt(abs(diag(expectedTarget)));
            expectedRelativeError = norm(scale.*(expectedGram - expectedTarget).*scale.',2);

            testCase.verifyClass(transform, "IMDiscreteTransform")
            testCase.verifyEqual(transform.z, z, AbsTol=0)
            testCase.verifyEqual(transform.weights, dz, AbsTol=0)
            testCase.verifyEqual(transform.modeNumber, basisSet.modeNumber(1:3), AbsTol=0)
            testCase.verifyEqual(transform.normalization, "unity")
            testCase.verifyEqual(transform.inverseMatrix, expectedBasis, RelTol=1e-13, AbsTol=1e-13)
            testCase.verifyEqual(transform.metricMatrix, expectedMetric, AbsTol=0)
            testCase.verifyEqual(transform.gramMatrix, expectedGram, RelTol=1e-13, AbsTol=1e-13)
            testCase.verifyEqual(transform.targetGramMatrix, expectedTarget, RelTol=1e-12, AbsTol=1e-12)
            testCase.verifyEqual(transform.relativeGramOperatorError, expectedRelativeError, RelTol=1e-12, AbsTol=1e-14)
            testCase.verifyEqual(transform.inverseMatrixConditionNumber, cond(expectedBasis), RelTol=1e-12)
            testCase.verifyEqual(transform.gramConditionNumber, cond(expectedGram), RelTol=1e-12)
            testCase.verifyTrue(transform.targetGramIsPositiveDefinite)
            testCase.verifyFalse(transform.hasNegativeWeights)
        end

        function forwardAndBackTransformsRoundTripRetainedCoefficients(testCase)
            [basisSet, z, dz] = testCase.regularBasis(4);
            transform = basisSet.discreteTransform(z=z, weights=dz, nModes=3);
            coefficients = [1 -2; 0.5 3; -4 0.25];
            values = transform.inverseMatrix*coefficients;

            testCase.verifyEqual(transform.forwardMatrix*transform.inverseMatrix, eye(3), RelTol=1e-11, AbsTol=1e-12)
            testCase.verifyEqual(transform.transformForward(values), coefficients, RelTol=1e-11, AbsTol=1e-12)
            testCase.verifyEqual(transform.transformForward(values), transform.forwardMatrix*values, RelTol=1e-13, AbsTol=1e-13)
            testCase.verifyEqual(transform.transformBack(coefficients), values, RelTol=1e-13, AbsTol=1e-13)
            testCase.verifyEqual(transform.transformBack(coefficients), transform.inverseMatrix*coefficients, RelTol=1e-13, AbsTol=1e-13)
            sampledProjector = transform.inverseMatrix*transform.forwardMatrix;
            testCase.verifyEqual(sampledProjector*sampledProjector, sampledProjector, RelTol=1e-11, AbsTol=1e-12)
            testCase.verifyLessThan(transform.roundTripError, 1e-11)
            testCase.verifyFalse(isprop(transform,"basisMatrix"))
            testCase.verifyFalse(isprop(transform,"basisConditionNumber"))
            testCase.verifyFalse(isprop(transform,"increments"))
            testCase.verifyFalse(isprop(transform,"hasNegativeIncrements"))
            testCase.verifyFalse(isprop(transform,"relativeGramError"))
            testCase.verifyFalse(ismethod(transform,"project"))
            testCase.verifyFalse(ismethod(transform,"reconstruct"))
        end

        function galerkinResidualIsOrthogonalInSampleMetric(testCase)
            [basisSet, z, dz] = testCase.regularBasis(4);
            transform = basisSet.discreteTransform(z=z, weights=dz, nModes=3);
            values = [z.^3 + 0.2*cos(7*z), exp(z) - z];
            coefficients = transform.transformForward(values);
            residual = values - transform.transformBack(coefficients);
            normalResidual = transform.inverseMatrix.'*transform.metricMatrix*residual;

            testCase.verifyLessThan(norm(normalResidual,2), 1e-11*max(1,norm(values,2)))
        end

        function scalarBuilderSelectsLeadingModePrefix(testCase)
            [basisSet, z, dz] = testCase.regularBasis(4);
            transform = basisSet.discreteTransform(z=z, weights=dz, nModes=2);

            testCase.verifySize(transform.inverseMatrix, [length(z) 2])
            testCase.verifySize(transform.forwardMatrix, [2 length(z)])
            testCase.verifyEqual(transform.modeNumber, basisSet.modeNumber(1:2), AbsTol=0)
            testCase.verifyError(@() basisSet.discreteTransform(z=z, weights=dz, nModes=5), "IMBasisSet:InvalidDiscreteModeCount")
            testCase.verifyError(@() basisSet.discreteTransform(z=z(1:2), weights=dz(1:2), nModes=3), "IMBasisSet:InsufficientDiscreteSamples")
        end

        function valueEndpointTermEntersSampleMetric(testCase)
            zDomain = [-1 0];
            solver = IMSolverSpectral(nEVP=96);
            surfaceBoundary = IMBoundaryCondition(a=0, b=1, c=1, d=0);
            evp = IMEigenvalueProblem(zDomain=zDomain, p=1, q=0, r=1, surfaceBoundary=surfaceBoundary, bottomBoundary=IMBoundaryCondition.dirichlet());
            basisSet = solver.solveEVP(evp, nModes=3);
            z = linspace(zDomain(1),zDomain(2),65).';
            dz = testCase.trapezoidalIncrements(z);
            transform = basisSet.discreteTransform(z=z, weights=dz, nModes=2, gramTolerance=1e6);
            endpointWeight = evp.innerProduct().surfaceWeights(1);
            expectedMetric = diag(dz);
            expectedMetric(end,end) = expectedMetric(end,end) + endpointWeight.coefficient*endpointWeight.c^2;

            testCase.verifyEqual(transform.metricMatrix, expectedMetric, AbsTol=1e-14)
            testCase.verifyEqual(transform.gramMatrix, transform.inverseMatrix.'*expectedMetric*transform.inverseMatrix, RelTol=1e-12, AbsTol=1e-12)
            testCase.verifyError(@() basisSet.discreteTransform(z=z(1:end-1), weights=dz(1:end-1), nModes=2), "IMBasisSet:MissingDiscreteEndpointSample")
        end

        function derivativeEndpointTermIsRejected(testCase)
            zDomain = [-1 0];
            solver = IMSolverSpectral(nEVP=96);
            surfaceBoundary = IMBoundaryCondition(a=1, b=0, c=0, d=1);
            evp = IMEigenvalueProblem(zDomain=zDomain, p=1, q=0, r=1, surfaceBoundary=surfaceBoundary, bottomBoundary=IMBoundaryCondition.dirichlet());
            basisSet = solver.solveEVP(evp, nModes=2);
            z = linspace(zDomain(1),zDomain(2),65).';
            dz = testCase.trapezoidalIncrements(z);

            testCase.verifyError(@() basisSet.discreteTransform(z=z, weights=dz), "IMBasisSet:UnsupportedDiscreteEndpointMetric")
        end

        function signedWeightsAreRetainedAndFlagged(testCase)
            [basisSet, z, dz] = testCase.regularBasis(3);
            dz(2) = -dz(2);
            transform = basisSet.discreteTransform(z=z, weights=dz, nModes=2, gramTolerance=1e6);

            testCase.verifyEqual(transform.weights, dz, AbsTol=0)
            testCase.verifyTrue(transform.hasNegativeWeights)
        end

        function indefiniteTargetsUseAbsoluteNormScaling(testCase)
            targetGram = diag([-1 2]);
            transform = IMDiscreteTransform(z=[-1; 0], weights=[-1; 1], modeNumber=[-1 1], normalization="unity", ...
                inverseMatrix=eye(2), metricMatrix=diag([-2 3]), targetGramMatrix=targetGram);
            scale = 1./sqrt(abs(diag(targetGram)));
            expected = norm(scale.*(transform.gramMatrix - targetGram).*scale.',2);

            testCase.verifyFalse(transform.targetGramIsPositiveDefinite)
            testCase.verifyTrue(transform.hasNegativeWeights)
            testCase.verifyEqual(transform.relativeGramOperatorError, expected, AbsTol=1e-14)
        end

        function invalidDiscreteInputsThrowStructuredErrors(testCase)
            [basisSet, z, dz] = testCase.regularBasis(3);
            outside = z;
            outside(1) = basisSet.zDomain(1) - 0.1;
            duplicate = z;
            duplicate(2) = duplicate(1);

            testCase.verifyError(@() basisSet.discreteTransform(z=flipud(z), weights=dz), "IMBasisSet:InvalidDiscreteGrid")
            testCase.verifyError(@() basisSet.discreteTransform(z=duplicate, weights=dz), "IMBasisSet:InvalidDiscreteGrid")
            testCase.verifyError(@() basisSet.discreteTransform(z=outside, weights=dz), "IMBasisSet:InvalidDiscreteGrid")
            testCase.verifyError(@() basisSet.discreteTransform(z=z, weights=dz(1:end-1)), "IMBasisSet:InvalidDiscreteWeights")
            testCase.verifyError(@() basisSet.discreteTransform(z=z, weights=zeros(size(dz))), "IMBasisSet:InvalidDiscreteWeights")
        end

        function matrixObjectDiagnosesSingularAndRejectsMalformedOperations(testCase)
            singular = IMDiscreteTransform(z=[-1; -0.5; 0],weights=ones(3,1),modeNumber=[1 2],normalization="unity", ...
                inverseMatrix=ones(3,2),metricMatrix=eye(3),targetGramMatrix=eye(2));
            testCase.verifyEqual(singular.sampledGramRank,1)
            testCase.verifyEqual(singular.relativeGramOperatorError,Inf)
            testCase.verifyGreaterThan(singular.roundTripError,0.9)

            transform = IMDiscreteTransform(z=[-1; 0], weights=ones(2,1), modeNumber=[1 2], normalization="unity", ...
                inverseMatrix=eye(2), metricMatrix=eye(2), targetGramMatrix=eye(2));
            testCase.verifyError(@() transform.transformForward(ones(3,1)), "IMDiscreteTransform:InvalidSampleCount")
            testCase.verifyError(@() transform.transformBack(ones(3,1)), "IMDiscreteTransform:InvalidCoefficientCount")
        end
    end

    methods (Access = private)
        function [basisSet, z, dz] = regularBasis(testCase, nModes)
            zDomain = [-1 0];
            solver = IMSolverSpectral(nEVP=96);
            evp = IMEigenvalueProblem(zDomain=zDomain, p=1, q=0, r=1, surfaceBoundary=IMBoundaryCondition.dirichlet(), bottomBoundary=IMBoundaryCondition.dirichlet());
            basisSet = solver.solveEVP(evp, nModes=nModes);
            z = linspace(zDomain(1),zDomain(2),65).';
            dz = testCase.trapezoidalIncrements(z);
        end

        function dz = trapezoidalIncrements(~, z)
            dz = [(z(2)-z(1))/2; (z(3:end)-z(1:end-2))/2; (z(end)-z(end-1))/2];
        end
    end
end
