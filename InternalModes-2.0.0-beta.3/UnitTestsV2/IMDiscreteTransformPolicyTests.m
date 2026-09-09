classdef IMDiscreteTransformPolicyTests < matlab.unittest.TestCase

    properties (Access = private)
        originalPath
        cosineBasis
        z
        weights
    end

    methods (TestClassSetup)
        function configureCosineReference(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.originalPath = path;
            addpath(repoRoot);

            zDomain = [-4000 0];
            solver = IMSolverSpectral(nEVP=160);
            evp = IMEigenvalueProblem(zDomain=zDomain,p=1,q=0,r=1, ...
                surfaceBoundary=IMBoundaryCondition.neumann(),bottomBoundary=IMBoundaryCondition.neumann());
            testCase.cosineBasis = solver.solveEVP(evp,nModes=30);
            testCase.z = linspace(zDomain(1),zDomain(2),13).';
            testCase.weights = [0.5; ones(11,1); 0.5]*(diff(zDomain)/12);
        end
    end

    methods (TestClassTeardown)
        function restoreRepositoryPath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function gramPolicyReturnsLargestConsecutivePrefix(testCase)
            [transform,assessment] = testCase.cosineBasis.discreteTransform(z=testCase.z,weights=testCase.weights);

            testCase.verifyEqual(assessment.candidateModeCount,13)
            testCase.verifyEqual(assessment.retainedModeCount,12)
            testCase.verifyEqual(transform.modeNumber,0:11,AbsTol=0)
            testCase.verifyTrue(all(assessment.prefixDiagnostics.gramAccepted(1:12)))
            testCase.verifyFalse(assessment.prefixDiagnostics.gramAccepted(13))
            testCase.verifyEqual(assessment.gramPolicy.maximumAcceptedModeCount,12)
            testCase.verifyEqual(assessment.gramPolicy.limitingValue,assessment.prefixDiagnostics.gramError(13),AbsTol=0)
            testCase.verifyEqual(assessment.limitingPolicy,"gram")
            testCase.verifyEqual(transform.weights,assessment.candidateTransform.weights,AbsTol=0)
        end

        function explicitBandFailureReportsPolicyValueToleranceAndCount(testCase)
            try
                testCase.cosineBasis.discreteTransform(z=testCase.z,weights=testCase.weights,nModes=13);
                testCase.assertFail("Expected the strict Gram policy to reject the requested band.")
            catch exception
                testCase.verifyEqual(string(exception.identifier),"IMBasisSet:DiscreteTransformPolicyFailed")
                testCase.verifySubstring(exception.message,"gram")
                testCase.verifySubstring(exception.message,"accepted 12")
                testCase.verifySubstring(exception.message,"tolerance 0.01")
                testCase.verifySubstring(exception.message,"first failing value")
            end
        end

        function rankDeficientCandidateReducesToFullRankPrefix(testCase)
            solver = IMSolverSpectral(nEVP=64);
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1,surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
            basis = solver.solveEVP(evp,nModes=6);
            sampleZ = linspace(-1,0,5).';
            sampleWeights = [0.5;ones(3,1);0.5]/4;

            [transform,assessment] = basis.discreteTransform(z=sampleZ,weights=sampleWeights);

            testCase.verifyEqual(assessment.candidateModeCount,5)
            testCase.verifyEqual(assessment.candidateTransform.sampledGramRank,3)
            testCase.verifyEqual(assessment.prefixDiagnostics.sampledGramRank,[1;2;3;3;3])
            testCase.verifyEqual(assessment.gramPolicy.maximumAcceptedModeCount,3)
            testCase.verifyEqual(assessment.retainedModeCount,3)
            testCase.verifyEqual(transform.sampledGramRank,3)
            testCase.verifyError(@() basis.discreteTransform(z=sampleZ,weights=sampleWeights,nModes=5),"IMBasisSet:DiscreteTransformPolicyFailed")
        end

        function leakageMatchesIndependentRejectedModeProjection(testCase)
            [~,assessment] = testCase.cosineBasis.discreteTransform(z=testCase.z,weights=testCase.weights, ...
                leakageTolerance=1e-8,nCheckModes=20);
            prefixCount = 6;
            prefixTransform = testCase.cosineBasis.discreteTransform(z=testCase.z,weights=testCase.weights,nModes=prefixCount);
            sourceValues = testCase.cosineBasis.u(testCase.z);
            sourceGram = testCase.cosineBasis.gramMatrix();
            rejectedIndices = (prefixCount+1):20;
            coefficients = prefixTransform.forwardMatrix*sourceValues(:,rejectedIndices);
            retainedNorms = abs(diag(prefixTransform.targetGramMatrix));
            rejectedNorms = abs(diag(sourceGram));
            independentLeakage = sqrt(sum(retainedNorms.*abs(coefficients).^2,1)./rejectedNorms(rejectedIndices).');
            [expectedMaximum,iLimiting] = max(independentLeakage);

            testCase.verifyEqual(assessment.prefixDiagnostics.leakageError(prefixCount),expectedMaximum,RelTol=1e-12,AbsTol=1e-13)
            testCase.verifyEqual(assessment.prefixDiagnostics.leakageLimitingModeNumber(prefixCount), ...
                testCase.cosineBasis.modeNumber(rejectedIndices(iLimiting)),AbsTol=0)
            testCase.verifyEqual(assessment.leakagePolicy.maximumAcceptedModeCount,5)
            testCase.verifyEqual(assessment.retainedModeCount,5)
        end

        function leakageDefaultChecksTwiceTheCandidateBand(testCase)
            [~,assessment] = testCase.cosineBasis.discreteTransform(z=testCase.z,weights=testCase.weights,leakageTolerance=2);

            testCase.verifyEqual(assessment.candidateModeCount,13)
            testCase.verifyEqual(assessment.leakagePolicy.nCheckModes,26)
            testCase.verifyEqual(length(assessment.leakagePolicy.error),13)
        end

        function enabledPoliciesImposeTheirMinimumLimit(testCase)
            [transform,assessment] = testCase.cosineBasis.discreteTransform(z=testCase.z,weights=testCase.weights, ...
                leakageTolerance=1e-8,nCheckModes=20,quadraticAliasingTolerance=1e-8);

            testCase.verifyEqual(assessment.gramPolicy.maximumAcceptedModeCount,12)
            testCase.verifyEqual(assessment.leakagePolicy.maximumAcceptedModeCount,5)
            testCase.verifyEqual(assessment.quadraticAliasingPolicy.maximumAcceptedModeCount,8)
            testCase.verifyEqual(assessment.retainedModeCount,5)
            testCase.verifyEqual(length(transform.modeNumber),5)
            testCase.verifyEqual(assessment.limitingPolicy,"leakage")
            testCase.verifyEqual(transform.weights,assessment.candidateTransform.weights,AbsTol=0)
        end

        function cosineQuadraticPolicyHasExactDCTICutoff(testCase)
            Nz = length(testCase.z);
            M = Nz-1;
            Kmax = floor((2*M-1)/3);
            [transform,assessment] = testCase.cosineBasis.discreteTransform(z=testCase.z,weights=testCase.weights, ...
                quadraticAliasingTolerance=1e-8);
            firstRejectedRow = Kmax+2;

            testCase.verifyEqual(Kmax,7)
            testCase.verifyEqual(assessment.retainedModeCount,Kmax+1)
            testCase.verifyEqual(transform.modeNumber,0:Kmax,AbsTol=0)
            testCase.verifyLessThan(assessment.prefixDiagnostics.quadraticAliasingError(Kmax+1),1e-8)
            testCase.verifyEqual(assessment.prefixDiagnostics.quadraticAliasingError(firstRejectedRow),1/sqrt(3),RelTol=1e-8)
            testCase.verifyEqual(assessment.prefixDiagnostics.quadraticLimitingModeNumberI(firstRejectedRow),Kmax+1,AbsTol=0)
            testCase.verifyEqual(assessment.prefixDiagnostics.quadraticLimitingModeNumberJ(firstRejectedRow),Kmax+1,AbsTol=0)
            testCase.verifyFalse(assessment.prefixDiagnostics.quadraticAccepted(firstRejectedRow))
        end

        function strictQuadraticBandDoesNotSilentlyReduce(testCase)
            try
                testCase.cosineBasis.discreteTransform(z=testCase.z,weights=testCase.weights,nModes=9,quadraticAliasingTolerance=1e-8);
                testCase.assertFail("Expected the strict quadratic policy to reject the requested band.")
            catch exception
                testCase.verifyEqual(string(exception.identifier),"IMBasisSet:DiscreteTransformPolicyFailed")
                testCase.verifySubstring(exception.message,"quadraticAliasing")
                testCase.verifySubstring(exception.message,"accepted 8")
                testCase.verifySubstring(exception.message,"tolerance 1e-08")
            end
        end

        function leakageCheckCountAndUnavailableModesAreActionable(testCase)
            testCase.verifyError(@() testCase.cosineBasis.discreteTransform(z=testCase.z,weights=testCase.weights, ...
                leakageTolerance=1e-8,nCheckModes=13),"IMBasisSet:InvalidLeakageCheckModeCount")

            solver = IMSolverSpectral(nEVP=8);
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1, ...
                surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
            basisSet = solver.solveEVP(evp,nModes=3);
            samplePoints = linspace(-1,0,5).';
            sampleWeights = [0.5; ones(3,1); 0.5]/4;
            testCase.verifyError(@() basisSet.discreteTransform(z=samplePoints,weights=sampleWeights,nModes=2, ...
                gramTolerance=1,leakageTolerance=1,nCheckModes=20),"IMBasisSet:AuxiliaryModeUnavailable")
        end

        function normPoliciesRejectSignedTargetsButGramRemainsAvailable(testCase)
            solver = IMSolverSpectral(nEVP=96);
            signedSurface = IMBoundaryCondition(a=0,b=1,c=-1,d=0);
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1,surfaceBoundary=signedSurface, ...
                bottomBoundary=IMBoundaryCondition.dirichlet());
            basisSet = solver.solveEVP(evp,nModes=5);
            samplePoints = linspace(-1,0,33).';
            sampleWeights = [0.5; ones(31,1); 0.5]/32;
            [~,assessment] = basisSet.discreteTransform(z=samplePoints,weights=sampleWeights,nModes=3,gramTolerance=10);

            testCase.verifyFalse(assessment.candidateTransform.targetGramIsPositiveDefinite)
            testCase.verifyFalse(assessment.leakagePolicy.enabled)
            testCase.verifyFalse(assessment.quadraticAliasingPolicy.enabled)
            testCase.verifyError(@() basisSet.discreteTransform(z=samplePoints,weights=sampleWeights,nModes=3,gramTolerance=10, ...
                leakageTolerance=1),"IMBasisSet:UnavailableDiscreteTransformPolicy")
            testCase.verifyError(@() basisSet.discreteTransform(z=samplePoints,weights=sampleWeights,nModes=3,gramTolerance=10, ...
                quadraticAliasingTolerance=1),"IMBasisSet:UnavailableDiscreteTransformPolicy")
        end
    end
end
