classdef IMQuadratureWeightFitTests < matlab.unittest.TestCase

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
        function omittedWeightsUseDefaultPhysicalFit(testCase)
            basisSet = testCase.regularBasis(4);
            z = -1 + linspace(0,1,17).'.^1.4;
            [weights, fit] = basisSet.quadratureWeightsForPoints(z=z,nModes=3);
            [transform, assessment] = basisSet.discreteTransform(z=z,nModes=3);
            builderFit = assessment.weightFit;

            testCase.verifyClass(fit,"IMQuadratureWeightFit")
            testCase.verifyClass(builderFit,"IMQuadratureWeightFit")
            testCase.verifyEqual(weights,fit.weights,AbsTol=0)
            testCase.verifyEqual(transform.weights,builderFit.weights,AbsTol=0)
            testCase.verifyEqual(transform.forwardMatrix,builderFit.transform.forwardMatrix,AbsTol=0)
            testCase.verifyEqual(transform.inverseMatrix,builderFit.transform.inverseMatrix,AbsTol=0)
            testCase.verifyEqual(transform.weights,fit.weights,AbsTol=1e-12)
            testCase.verifyEqual(transform.forwardMatrix,fit.transform.forwardMatrix,RelTol=1e-12,AbsTol=1e-12)
            testCase.verifyEqual(fit.objectiveName,"normalizedGramFrobenius")
            testCase.verifyTrue(fit.nonnegativeConstraint)
            testCase.verifyTrue(fit.depthConstraint)
            testCase.verifyGreaterThan(fit.exitFlag,0)
            testCase.verifyGreaterThanOrEqual(min(fit.weights),-1e-12)
            testCase.verifyEqual(sum(fit.weights),1,AbsTol=1e-10)
            testCase.verifyLessThanOrEqual(fit.residualNorm,fit.geometricResidualNorm + 1e-10)
            targetNorms = diag(fit.transform.targetGramMatrix);
            scale = 1./sqrt(abs(targetNorms));
            fittedGramMismatch = scale.*(fit.transform.gramMatrix-fit.transform.targetGramMatrix).*scale.';
            geometricGramMismatch = scale.*(fit.geometricTransform.gramMatrix-fit.geometricTransform.targetGramMatrix).*scale.';
            testCase.verifyEqual(fit.residualNorm,norm(fittedGramMismatch,"fro"),RelTol=1e-12,AbsTol=1e-13)
            testCase.verifyEqual(fit.geometricResidualNorm,norm(geometricGramMismatch,"fro"),RelTol=1e-12,AbsTol=1e-13)
            testCase.verifyEqual(fit.transform.relativeGramOperatorError,norm(fittedGramMismatch,2),RelTol=1e-12,AbsTol=1e-13)
            testCase.verifyFalse(ismethod(basisSet,"fitQuadrature"))
            testCase.verifyFalse(isprop(fit,"fittedTransform"))
            testCase.verifyFalse(isprop(fit,"fittedIncrements"))
            testCase.verifyEmpty(meta.class.fromName("IMQuadratureFit"))
        end

        function suppliedWeightsBypassFitting(testCase)
            basisSet = testCase.regularBasis(3);
            z = linspace(-1,0,21).';
            dz = testCase.geometricWeights(z,[-1 0]);
            [transform, assessment] = basisSet.discreteTransform(z=z,weights=dz,nModes=2);

            testCase.verifyEqual(transform.weights,dz,AbsTol=0)
            testCase.verifyEqual(transform.metricMatrix,diag(dz),AbsTol=0)
            testCase.verifyClass(assessment,"IMDiscreteTransformAssessment")
            testCase.verifyEmpty(assessment.weightFit)
        end

        function geometricWeightsUseFullDomainControlVolumes(testCase)
            basisSet = testCase.regularBasis(3);
            zInterior = [-0.9; -0.63; -0.31; -0.08];
            [~, fitInterior] = basisSet.quadratureWeightsForPoints(z=zInterior,nModes=2);
            expectedInterior = testCase.geometricWeights(zInterior,[-1 0]);

            zEndpoints = [-1; -0.72; -0.4; -0.15; 0];
            [~, fitEndpoints] = basisSet.quadratureWeightsForPoints(z=zEndpoints,nModes=2);
            expectedEndpoints = testCase.geometricWeights(zEndpoints,[-1 0]);

            testCase.verifyEqual(fitInterior.geometricWeights,expectedInterior,AbsTol=1e-14)
            testCase.verifyEqual(fitEndpoints.geometricWeights,expectedEndpoints,AbsTol=1e-14)
            testCase.verifyEqual(sum(fitInterior.geometricWeights),1,AbsTol=1e-14)
            testCase.verifyEqual(sum(fitEndpoints.geometricWeights),1,AbsTol=1e-14)
        end

        function normalizedGramSystemIncludesFixedEndpointContribution(testCase)
            zDomain = [-1 0];
            solver = IMSolverSpectral(nEVP=96);
            surfaceBoundary = IMBoundaryCondition(a=0,b=1,c=1,d=0);
            evp = IMEigenvalueProblem(zDomain=zDomain,p=1,q=0,r=1,surfaceBoundary=surfaceBoundary,bottomBoundary=IMBoundaryCondition.dirichlet());
            basisSet = solver.solveEVP(evp,nModes=3);
            z = linspace(-1,0,17).';
            [~, fit] = basisSet.quadratureWeightsForPoints(z=z,nModes=2);

            Phi = fit.geometricTransform.inverseMatrix;
            target = fit.geometricTransform.targetGramMatrix;
            targetNorms = diag(target);
            endpointWeight = evp.innerProduct().surfaceWeights(1);
            endpointGram = endpointWeight.coefficient*endpointWeight.c^2*(Phi(end,:).'*Phi(end,:));
            expectedPairs = [1 1; 1 2; 2 2];
            expectedA = zeros(3,length(z));
            expectedB = zeros(3,1);
            iRow = 0;
            for iMode = 1:2
                for jMode = iMode:2
                    iRow = iRow + 1;
                    rowFactor = 1;
                    if iMode ~= jMode
                        rowFactor = sqrt(2);
                    end
                    scale = sqrt(abs(targetNorms(iMode)*targetNorms(jMode)));
                    expectedA(iRow,:) = rowFactor*(Phi(:,iMode).*Phi(:,jMode)).'/scale;
                    expectedB(iRow) = rowFactor*(target(iMode,jMode) - endpointGram(iMode,jMode))/scale;
                end
            end

            testCase.verifySize(fit.objectiveMatrix,[size(expectedPairs,1) length(z)])
            testCase.verifyEqual(fit.objectiveMatrix,expectedA,RelTol=1e-12,AbsTol=1e-12)
            testCase.verifyEqual(fit.objectiveTarget,expectedB,RelTol=1e-12,AbsTol=1e-12)
            testCase.verifyEqual(fit.residual,expectedA*fit.weights - expectedB,AbsTol=1e-13)
            testCase.verifyEqual(fit.geometricResidual,expectedA*fit.geometricWeights - expectedB,AbsTol=1e-13)
            testCase.verifyError(@() basisSet.quadratureWeightsForPoints(z=z(1:end-1),nModes=2),"IMBasisSet:MissingDiscreteEndpointSample")
        end

        function derivativeEndpointMetricRemainsUnsupported(testCase)
            solver = IMSolverSpectral(nEVP=96);
            surfaceBoundary = IMBoundaryCondition(a=1,b=0,c=0,d=1);
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1,surfaceBoundary=surfaceBoundary,bottomBoundary=IMBoundaryCondition.dirichlet());
            basisSet = solver.solveEVP(evp,nModes=2);
            z = linspace(-1,0,17).';

            testCase.verifyError(@() basisSet.quadratureWeightsForPoints(z=z),"IMBasisSet:UnsupportedDiscreteEndpointMetric")
        end

        function fittedRuleImprovesVariableWeightGramDiagnostic(testCase)
            zDomain = [-1 0];
            solver = IMSolverSpectral(nEVP=128);
            evp = IMEigenvalueProblem(zDomain=zDomain,p=1,q=0,r=@(z) exp(2*z),surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
            basisSet = solver.solveEVP(evp,nModes=6);
            z = -1 + linspace(0,1,15).'.^1.7;
            [~, fit] = basisSet.quadratureWeightsForPoints(z=z,nModes=5);

            testCase.verifyLessThan(fit.transform.relativeGramOperatorError,fit.geometricTransform.relativeGramOperatorError)
            testCase.verifyLessThanOrEqual(fit.residualNorm,fit.geometricResidualNorm + 1e-10)
        end

        function physicalDepthScalingConvergesForHydrostaticModes(testCase)
            D = 4000;
            N0 = 5.2e-3;
            b = 1300;
            g = 9.81;
            zDomain = [-D 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=zDomain,g=g, ...
                surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
            solver = IMSolverSpectral(nEVP=160,coordinateKind="wkb");
            basisSet = solver.solveEVP(evp,nModes=8);
            sigma = linspace(0,1,24).';
            z = zDomain(1) + D*(1 - (1 - sigma).^2);
            [~, fit] = basisSet.quadratureWeightsForPoints(z=z,nModes=4,variables="G");

            testCase.verifyGreaterThan(fit.exitFlag,0)
            testCase.verifyGreaterThanOrEqual(min(fit.weights),-1e-10)
            testCase.verifyEqual(sum(fit.weights),D,AbsTol=1e-8)
            testCase.verifyLessThan(fit.residualNorm,fit.geometricResidualNorm)
            fittedDiagnostics = fit.transform.channelDiagnostics(variable="G");
            geometricDiagnostics = fit.geometricTransform.channelDiagnostics(variable="G");
            testCase.verifyLessThan(fittedDiagnostics.relativeGramOperatorError,geometricDiagnostics.relativeGramOperatorError)
        end

        function dimensionNormalizedRegularizedObjectiveMatchesDefinition(testCase)
            basisSet = testCase.regularBasis(4);
            z = -1 + linspace(0,1,13).'.^1.3;
            nModes = 3;
            lambda = 1e-2;
            [~, baseline] = basisSet.quadratureWeightsForPoints(z=z,nModes=nModes);
            objective = @(context) testCase.dimensionNormalizedRegularizedObjective(context,lambda,nModes);
            [~, fit] = basisSet.quadratureWeightsForPoints(z=z,nModes=nModes,objective=objective);
            expectedA = [baseline.objectiveMatrix/nModes; sqrt(lambda/length(z))*diag(1./baseline.geometricWeights)];
            expectedB = [baseline.objectiveTarget/nModes; sqrt(lambda/length(z))*ones(length(z),1)];

            testCase.verifyEqual(fit.objectiveName,"dimensionNormalizedRegularizedGram")
            testCase.verifyEqual(fit.objectiveMatrix,expectedA,RelTol=1e-13,AbsTol=1e-13)
            testCase.verifyEqual(fit.objectiveTarget,expectedB,RelTol=1e-13,AbsTol=1e-13)
        end

        function regularizationPreservesAccuracyOnClusteredExponentialGrid(testCase)
            [basisSet,z] = testCase.exponentialBasisAndClusteredGrid();
            nModes = 8;
            lambda = 1e-6;
            [~, pureFit] = basisSet.quadratureWeightsForPoints(z=z,nModes=nModes,variables="G");
            objective = @(context) testCase.dimensionNormalizedRegularizedObjective(context,lambda,nModes);
            [~, regularizedFit] = basisSet.quadratureWeightsForPoints(z=z,nModes=nModes,variables="G",objective=objective);
            pureDisplacement = testCase.relativeWeightDisplacement(pureFit.weights,pureFit.geometricWeights);
            regularizedDisplacement = testCase.relativeWeightDisplacement(regularizedFit.weights,pureFit.geometricWeights);

            regularizedDiagnostics = regularizedFit.transform.channelDiagnostics(variable="G");
            pureDiagnostics = pureFit.transform.channelDiagnostics(variable="G");
            testCase.verifyLessThanOrEqual(regularizedDiagnostics.relativeGramOperatorError,max(10*pureDiagnostics.relativeGramOperatorError+1e-12,1e-8))
            testCase.verifyLessThan(regularizedDisplacement,0.01*pureDisplacement)
            testCase.verifyGreaterThan(min(regularizedFit.weights),0)
        end

        function regularizationPreservesAccuracyOnModeRootGrid(testCase)
            [basisSet,~] = testCase.exponentialBasisAndClusteredGrid();
            nModes = 8;
            lambda = 1e-6;
            z = basisSet.pointsFromModeRoots(nModes=nModes);
            [~, pureFit] = basisSet.quadratureWeightsForPoints(z=z,nModes=nModes,variables="G");
            objective = @(context) testCase.dimensionNormalizedRegularizedObjective(context,lambda,nModes);
            [~, regularizedFit] = basisSet.quadratureWeightsForPoints(z=z,nModes=nModes,variables="G",objective=objective);
            pureDisplacement = testCase.relativeWeightDisplacement(pureFit.weights,pureFit.geometricWeights);
            regularizedDisplacement = testCase.relativeWeightDisplacement(regularizedFit.weights,pureFit.geometricWeights);

            regularizedDiagnostics = regularizedFit.transform.channelDiagnostics(variable="G");
            pureDiagnostics = pureFit.transform.channelDiagnostics(variable="G");
            testCase.verifyLessThanOrEqual(regularizedDiagnostics.relativeGramOperatorError,1.10*pureDiagnostics.relativeGramOperatorError + 1e-12)
            testCase.verifyLessThan(regularizedDisplacement,pureDisplacement)
            testCase.verifyGreaterThan(min(regularizedFit.weights),0)
        end

        function regularizationPreservesEndpointAndZeroModeConstraints(testCase)
            lambda = 1e-6;
            endpointSolver = IMSolverSpectral(nEVP=96);
            endpointEVP = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1,surfaceBoundary=IMBoundaryCondition(a=0,b=1,c=1,d=0),bottomBoundary=IMBoundaryCondition.dirichlet());
            endpointBasis = endpointSolver.solveEVP(endpointEVP,nModes=4);
            endpointZ = linspace(-1,0,17).';
            endpointObjective = @(context) testCase.dimensionNormalizedRegularizedObjective(context,lambda,3);
            [endpointWeights,endpointFit] = endpointBasis.quadratureWeightsForPoints(z=endpointZ,nModes=3,objective=endpointObjective);

            D = 4000;
            N0 = 5.2e-3;
            b = 1300;
            g = 9.81;
            zDomain = [-D 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            fSolver = IMSolverSpectral(nEVP=160,coordinateKind="wkb");
            fBasis = fSolver.solveEVP(IMInternalModes.hydrostaticFModes(N2=N2,zDomain=zDomain,g=g),nModes=6);
            sigma = linspace(0,1,18).';
            fZ = zDomain(1) + D*(1-(1-sigma).^2);
            fObjective = @(context) testCase.dimensionNormalizedRegularizedObjective(context,lambda,4);
            [fWeights,fFit] = fBasis.quadratureWeightsForPoints(z=fZ,nModes=4,variables="F",objective=fObjective);

            testCase.verifyGreaterThanOrEqual(min(endpointWeights),-1e-12)
            testCase.verifyEqual(sum(endpointWeights),1,AbsTol=1e-10)
            testCase.verifyTrue(endpointFit.transform.targetGramIsPositiveDefinite)
            testCase.verifyEqual(fBasis.modeNumber(1),0)
            testCase.verifyGreaterThanOrEqual(min(fWeights),-1e-10)
            testCase.verifyEqual(sum(fWeights),D,AbsTol=1e-8)
            fDiagnostics = fFit.transform.channelDiagnostics(variable="F");
            testCase.verifyTrue(fDiagnostics.targetGramIsPositiveDefinite)
        end

        function customObjectiveCanReplaceNormalizedGramSystem(testCase)
            basisSet = testCase.regularBasis(3);
            z = [-1; -0.73; -0.44; -0.18; 0];
            desiredWeights = [0.08; 0.17; 0.24; 0.29; 0.22];
            objective = @(context) struct("A",eye(size(context.inverseMatrix,1)),"b",desiredWeights,"name","prescribedWeights");
            [~, fit] = basisSet.quadratureWeightsForPoints(z=z,nModes=2,objective=objective);

            testCase.verifyEqual(fit.objectiveName,"prescribedWeights")
            testCase.verifyEqual(fit.objectiveMatrix,eye(length(z)),AbsTol=0)
            testCase.verifyEqual(fit.objectiveTarget,desiredWeights,AbsTol=0)
            testCase.verifyEqual(fit.weights,desiredWeights,AbsTol=1e-10)
            testCase.verifyEqual(fit.residualNorm,norm(fit.objectiveMatrix*fit.weights-fit.objectiveTarget,2),AbsTol=1e-14)
        end

        function customObjectiveCanReweightAndAppendRows(testCase)
            basisSet = testCase.regularBasis(3);
            z = -1 + linspace(0,1,9).'.^1.3;
            [~, baseline] = basisSet.quadratureWeightsForPoints(z=z,nModes=2);
            expectedA = [2*baseline.objectiveMatrix; eye(length(z))];
            expectedB = [2*baseline.objectiveTarget; baseline.geometricWeights];
            objective = @(context) struct("A",[2*context.normalizedGramA; eye(length(context.z))], ...
                "b",[2*context.normalizedGramB; context.geometricWeights],"name","weightedRegularizedGram");
            [~, fit] = basisSet.quadratureWeightsForPoints(z=z,nModes=2,objective=objective);

            testCase.verifyEqual(fit.objectiveName,"weightedRegularizedGram")
            testCase.verifyEqual(fit.objectiveMatrix,expectedA,RelTol=1e-13,AbsTol=1e-13)
            testCase.verifyEqual(fit.objectiveTarget,expectedB,RelTol=1e-13,AbsTol=1e-13)
        end

        function customObjectiveReceivesNormalizedGramModePairs(testCase)
            basisSet = testCase.regularBasis(4);
            z = -1 + linspace(0,1,11).'.^1.2;
            expectedPairs = [1 1; 1 2; 1 3; 2 2; 2 3; 3 3];
            objective = @(context) testCase.modePairAwareObjective(context,expectedPairs);
            [~, fit] = basisSet.quadratureWeightsForPoints(z=z,nModes=3,objective=objective);

            testCase.verifyEqual(fit.objectiveName,"modePairAwareNormalizedGram")
            testCase.verifySize(fit.objectiveMatrix,[size(expectedPairs,1) length(z)])
            testCase.verifyEqual(fit.residual,fit.objectiveMatrix*fit.weights-fit.objectiveTarget,AbsTol=1e-13)
        end

        function constraintsCanBeDisabledIndependently(testCase)
            basisSet = testCase.regularBasis(3);
            z = [-1; -0.75; -0.5; -0.25; 0];
            signedWeights = [-0.02; 0.2; 0.3; 0.3; 0.22];
            signedObjective = @(context) struct("A",eye(length(context.z)),"b",signedWeights,"name","signed");
            [~, signedFit] = basisSet.quadratureWeightsForPoints(z=z,nModes=2,objective=signedObjective,nonnegative=false);

            shallowWeights = 0.75*testCase.geometricWeights(z,[-1 0]);
            shallowObjective = @(context) struct("A",eye(length(context.z)),"b",shallowWeights,"name","shallow");
            [~, shallowFit] = basisSet.quadratureWeightsForPoints(z=z,nModes=2,objective=shallowObjective,constrainDepth=false);

            testCase.verifyFalse(signedFit.nonnegativeConstraint)
            testCase.verifyTrue(signedFit.depthConstraint)
            testCase.verifyEqual(signedFit.weights,signedWeights,AbsTol=1e-10)
            testCase.verifyTrue(signedFit.transform.hasNegativeWeights)
            testCase.verifyEqual(sum(signedFit.weights),1,AbsTol=1e-10)
            testCase.verifyTrue(shallowFit.nonnegativeConstraint)
            testCase.verifyFalse(shallowFit.depthConstraint)
            testCase.verifyEqual(shallowFit.weights,shallowWeights,AbsTol=1e-10)
            testCase.verifyEqual(sum(shallowFit.weights),0.75,AbsTol=1e-10)
        end

        function malformedObjectivesThrowStructuredErrors(testCase)
            basisSet = testCase.regularBasis(3);
            z = linspace(-1,0,7).';
            missingTarget = @(context) struct("A",eye(length(context.z)));
            wrongColumns = @(context) struct("A",ones(2,length(context.z) + 1),"b",ones(2,1));
            nonfinite = @(context) struct("A",NaN(1,length(context.z)),"b",1);
            failing = @(context) error("Test:ObjectiveFailure","failed for %d points",length(context.z));

            testCase.verifyError(@() basisSet.quadratureWeightsForPoints(z=z,nModes=2,objective="notAnObjective"),"IMBasisSet:UnknownQuadratureObjective")
            testCase.verifyError(@() basisSet.quadratureWeightsForPoints(z=z,nModes=2,objective="normalizedGram"),"IMBasisSet:UnknownQuadratureObjective")
            testCase.verifyError(@() basisSet.quadratureWeightsForPoints(z=z,nModes=2,objective=missingTarget),"IMBasisSet:InvalidQuadratureObjective")
            testCase.verifyError(@() basisSet.quadratureWeightsForPoints(z=z,nModes=2,objective=wrongColumns),"IMBasisSet:InvalidQuadratureObjective")
            testCase.verifyError(@() basisSet.quadratureWeightsForPoints(z=z,nModes=2,objective=nonfinite),"IMBasisSet:InvalidQuadratureObjective")
            testCase.verifyError(@() basisSet.quadratureWeightsForPoints(z=z,nModes=2,objective=failing),"IMBasisSet:QuadratureObjectiveFailed")
        end
    end

    methods (Access = private)
        function basisSet = regularBasis(~,nModes)
            solver = IMSolverSpectral(nEVP=96);
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1,surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
            basisSet = solver.solveEVP(evp,nModes=nModes);
        end

        function dz = geometricWeights(~,z,zDomain)
            edges = [zDomain(1); 0.5*(z(1:end-1) + z(2:end)); zDomain(2)];
            dz = diff(edges);
        end

        function specification = dimensionNormalizedRegularizedObjective(~,context,lambda,nModes)
            nSamples = length(context.z);
            regularizationMatrix = sqrt(lambda/nSamples)*diag(1./context.geometricWeights);
            regularizationTarget = sqrt(lambda/nSamples)*ones(nSamples,1);
            specification = struct("A",[context.normalizedGramA/nModes; regularizationMatrix],"b",[context.normalizedGramB/nModes; regularizationTarget],"name","dimensionNormalizedRegularizedGram");
        end

        function specification = modePairAwareObjective(~,context,expectedPairs)
            if ~isequal(context.normalizedGramModePairs,expectedPairs)
                error("Test:UnexpectedNormalizedGramModePairs","The normalized Gram rows do not use the expected upper-triangle mode-pair ordering.");
            end
            specification = struct("A",context.normalizedGramA,"b",context.normalizedGramB,"name","modePairAwareNormalizedGram");
        end

        function [basisSet,z] = exponentialBasisAndClusteredGrid(~)
            D = 4000;
            N0 = 5.2e-3;
            b = 1300;
            g = 9.81;
            zDomain = [-D 0];
            N2 = @(z) N0*N0*exp(2*z/b);
            evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=zDomain,g=g,surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
            solver = IMSolverSpectral(nEVP=160,coordinateKind="wkb");
            basisSet = solver.solveEVP(evp,nModes=12);
            basisSet.normalization = "geostrophic";
            sigma = linspace(0,1,24).';
            z = zDomain(1) + D*(1-(1-sigma).^2);
        end

        function displacement = relativeWeightDisplacement(~,weights,geometricWeights)
            displacement = norm((weights-geometricWeights)./geometricWeights)/sqrt(length(weights));
        end
    end
end
