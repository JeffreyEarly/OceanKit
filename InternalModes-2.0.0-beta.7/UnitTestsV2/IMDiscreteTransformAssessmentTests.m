classdef IMDiscreteTransformAssessmentTests < matlab.unittest.TestCase

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
        function exactPointConstructionCoversStratificationsAndCoordinates(testCase)
            N0 = 5.2e-3;
            zDomain = [-4000 0];
            coordinateKinds = ["z" "wkb" "density"];
            constantN2 = @(z) N0*N0*ones(size(z));
            exponentialN2 = @(z) N0*N0*exp(2*z/1300);
            evps = [IMInternalModes.hydrostaticFModes(N2=constantN2,zDomain=zDomain), ...
                IMInternalModes.hydrostaticFModes(N2=exponentialN2,zDomain=zDomain), ...
                IMInternalModes.hydrostaticGModes(N2=exponentialN2,zDomain=zDomain)];

            for iEVP = 1:length(evps)
                for iCoordinate = 1:length(coordinateKinds)
                    solver = IMSolverSpectral(nEVP=128,coordinateKind=coordinateKinds(iCoordinate));
                    basisSet = solver.solveEVP(evps(iEVP),nModes=12);
                    [transform,assessment] = basisSet.discreteTransform(nPoints=8);

                    testCase.verifyClass(assessment,"IMInternalModesDiscreteTransformAssessment")
                    testCase.verifyEqual(assessment.requestedPointCount,8)
                    testCase.verifyEqual(assessment.actualPointCount,8)
                    testCase.verifyEqual(length(transform.z),8)
                    testCase.verifyEqual(transform.z(1),zDomain(1),AbsTol=0)
                    testCase.verifyEqual(transform.z(end),zDomain(2),AbsTol=0)
                    testCase.verifyEqual(max(assessment.certificationSearch.modeCount),6)
                    testCase.verifyEqual(assessment.candidateModeCount,assessment.retainedModeCount)
                    testCase.verifyEqual(assessment.weightFitModeCount,assessment.retainedModeCount)
                    testCase.verifyGreaterThanOrEqual(assessment.retainedModeCount,1)
                    testCase.verifyLessThanOrEqual(assessment.retainedModeCount,6)
                    testCase.verifyEqual(assessment.gramPolicy.tolerance,1e-2,AbsTol=0)
                    testCase.verifyLessThanOrEqual(assessment.prefixDiagnostics.gramError(assessment.retainedModeCount),1e-2)
                    testCase.verifyEqual(assessment.gridDesign.kind,"modeRoot")
                    testCase.verifyEqual(assessment.gridDesign.pointCount,8)
                end
            end
        end

        function assessmentReproducesFixedRuleAndNestedFit(testCase)
            basisSet = testCase.constantFBasis(10);
            z = -4000 + 4000*linspace(0,1,17).'.^1.25;
            [transform,assessment] = basisSet.discreteTransform(z=z,nModes=6);
            candidate = assessment.candidateTransform;
            fit = assessment.weightFit;
            finalRow = assessment.prefixDiagnostics(end,:);

            testCase.verifyClass(fit,"IMInternalModesQuadratureWeightFit")
            testCase.verifyEqual(fit.transform.weights,candidate.weights,AbsTol=0)
            testCase.verifyEqual(transform.weights,candidate.weights,AbsTol=0)
            testCase.verifyEqual(assessment.candidateModeNumber,basisSet.modeNumber(1:6),AbsTol=0)
            testCase.verifyEqual(assessment.retainedModeNumber,transform.modeNumber,AbsTol=0)
            channelErrors = zeros(size(candidate.availableVariables));
            for iVariable = 1:length(candidate.availableVariables)
                variable = candidate.availableVariables(iVariable);
                testCase.verifyEqual(fit.transform.forwardMatrix(variable=variable),candidate.forwardMatrix(variable=variable),AbsTol=0)
                testCase.verifyEqual(fit.transform.inverseMatrix(variable=variable),candidate.inverseMatrix(variable=variable),AbsTol=0)
                testCase.verifyEqual(assessment.transform.forwardMatrix(variable=variable),transform.forwardMatrix(variable=variable),AbsTol=0)
                channelErrors(iVariable) = candidate.relativeGramOperatorError(variable=variable);
                variableDiagnostics = assessment.variablePrefixDiagnostics(variable=variable);
                testCase.verifyEqual(variableDiagnostics.roundTripError(end),candidate.roundTripError(variable=variable),AbsTol=0)
                testCase.verifyEqual(variableDiagnostics.inverseMatrixConditionNumber(end),candidate.inverseMatrixConditionNumber(variable=variable),AbsTol=0)
                testCase.verifyEqual(variableDiagnostics.gramConditionNumber(end),candidate.gramConditionNumber(variable=variable),AbsTol=0)
            end
            testCase.verifyEqual(finalRow.gramError,max(channelErrors),AbsTol=0)
            testCase.verifyEqual(height(assessment.prefixDiagnostics),assessment.candidateModeCount)
            testCase.verifyEqual(assessment.limitingPolicy,"none")
        end

        function namedGridAndExactFitMakeProvenanceAndBandExplicit(testCase)
            basisSet = testCase.constantFBasis(10);
            [z,gridDesign] = basisSet.modeRootGrid(modeCount=6);
            [transform,assessment] = basisSet.fitDiscreteTransform(z=z,modeCount=6,variables=["F","G"],gridDesign=gridDesign,gramTolerance=1);

            testCase.verifyEqual(length(transform.modeNumber),6)
            testCase.verifyEqual(assessment.retainedModeCount,6)
            testCase.verifyEqual(assessment.weightFitModeCount,6)
            testCase.verifyEmpty(assessment.certificationSearch)
            testCase.verifyEqual(assessment.gridDesign.kind,"modeRoot")
            testCase.verifyEqual(assessment.gridDesign.sourceFamily,"hydrostatic")
            testCase.verifyEqual(assessment.gridDesign.generatingVariable,"F")
            testCase.verifyEqual(assessment.gridDesign.generatingModeIndex,7)
            testCase.verifySubstring(assessment.gridDesign.interpretationForG,"G-extrema-like")
        end

        function certifiedSelectionIsUnaffectedByAdditionalSolvedModes(testCase)
            shortBasis = testCase.constantFBasis(10);
            longBasis = testCase.constantFBasis(16);
            [z,gridDesign] = shortBasis.modeRootGrid(modeCount=6);
            [shortTransform,shortAssessment] = shortBasis.certifiedDiscreteTransform( ...
                z=z,variables=["F","G"],gridDesign=gridDesign,gramTolerance=1);
            [longTransform,longAssessment] = longBasis.certifiedDiscreteTransform( ...
                z=z,variables=["F","G"],gridDesign=gridDesign,gramTolerance=1);

            testCase.verifyEqual(longAssessment.retainedModeCount,shortAssessment.retainedModeCount)
            testCase.verifyEqual(longTransform.modeNumber,shortTransform.modeNumber,AbsTol=0)
            testCase.verifyEqual(longTransform.weights,shortTransform.weights,RelTol=1e-12,AbsTol=1e-12)
        end

        function certifiedSearchProbeMatchesCompletePrefixAssessment(testCase)
            basisSet = testCase.constantFBasis(10);
            z = linspace(-4000,0,9).';
            variables = ["F","G"];
            gramTolerance = 1e-2;
            [transform,assessment] = basisSet.certifiedDiscreteTransform(z=z,variables=variables,gramTolerance=gramTolerance);
            search = assessment.certificationSearch;

            testCase.verifyGreaterThan(height(search),1)
            for iSearch = 1:height(search)
                modeCount = search.modeCount(iSearch);
                [referenceTransform,referenceAssessment] = basisSet.discreteTransform( ...
                    z=z,nModes=modeCount,variables=variables,gramTolerance=gramTolerance,allowRetainedPrefix=true);
                referenceAccepted = referenceAssessment.gramPolicy.accepted(end);
                referenceError = referenceAssessment.gramPolicy.error(end);
                fullBandAccepted = referenceError <= gramTolerance;

                testCase.verifyTrue(search.fitSucceeded(iSearch))
                testCase.verifyEqual(referenceAccepted,fullBandAccepted)
                testCase.verifyEqual(search.accepted(iSearch),referenceAccepted)
                if isfinite(referenceError)
                    testCase.verifyEqual(search.gramError(iSearch),referenceError,RelTol=2e-12,AbsTol=2e-14)
                else
                    testCase.verifyEqual(search.gramError(iSearch),referenceError)
                end
                if referenceAccepted
                    testCase.verifyEqual(search.retainedModeCount(iSearch),modeCount)
                    testCase.verifyEqual(transform.weights,referenceTransform.weights,RelTol=2e-12,AbsTol=2e-12)
                else
                    testCase.verifyTrue(isnan(search.retainedModeCount(iSearch)))
                end
            end
        end

        function signedAPVFullBandAcceptanceImpliesPrefixAcceptance(testCase)
            D = 4000;
            N0 = 5.2e-3;
            N2 = @(z) N0*N0*exp(2*z/1300);
            zDomain = [-D 0];
            g0 = -integral(N2,zDomain(1),zDomain(2));
            solver = IMSolverSpectral(nEVP=128,coordinateKind="wkb");
            evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=zDomain,g0=g0,gd=Inf,surfaceBoundary="freeSurface");
            basisSet = solver.solveEVP(evp,nModes=16);
            [z,gridDesign] = basisSet.modeRootGrid(nPoints=16);
            [~,assessment] = basisSet.certifiedDiscreteTransform(z=z,gridDesign=gridDesign,variables=["F","G"],gramTolerance=1e-2);

            search = assessment.certificationSearch;
            acceptedRow = find(search.accepted,1);
            modeCount = search.modeCount(acceptedRow);
            [referenceTransform,referenceAssessment] = basisSet.discreteTransform(z=z,nModes=modeCount,variables=["F","G"],gramTolerance=1e-2,allowRetainedPrefix=true);
            signedG = referenceAssessment.candidateTransform.targetGramMatrix(variable="G");

            testCase.verifyTrue(any(diag(signedG) < 0))
            testCase.verifyFalse(referenceTransform.hasNegativeWeights)
            testCase.verifyTrue(all(referenceAssessment.gramPolicy.accepted))
            testCase.verifyEqual(referenceAssessment.retainedModeCount,modeCount)
        end

        function fixedWavenumberWavePagesCertifyIndependently(testCase)
            N0 = 5.2e-3;
            zDomain = [-4000 0];
            N2 = @(z) N0*N0*exp(2*z/1300);
            z = linspace(zDomain(1),zDomain(2),17).';
            solver = IMSolverSpectral(nEVP=128,coordinateKind="wkb");
            surfaceBoundary = IMBoundaryCondition(a=0,b=1,c=1,d=0);
            k = 2*pi./[200e3 20e3];

            for iK = 1:length(k)
                evp = IMInternalModes.waveModesAtWavenumber(N2=N2,zDomain=zDomain,k=k(iK),f0=1e-4,surfaceBoundary=surfaceBoundary);
                basisSet = solver.solveEVP(evp,nModes=16);
                [transform,assessment] = basisSet.certifiedDiscreteTransform(z=z,variables="G",gramTolerance=1e-2);

                testCase.verifyEqual(length(transform.modeNumber),assessment.retainedModeCount)
                testCase.verifyEqual(assessment.weightFitModeCount,assessment.retainedModeCount)
                testCase.verifyLessThanOrEqual(assessment.gramPolicy.error(end),1e-2)
                testCase.verifyTrue(assessment.certificationSearch.accepted(end))
            end
        end

        function suppliedWeightsBypassOnlyTheFit(testCase)
            basisSet = testCase.constantFBasis(8);
            z = linspace(-4000,0,17).';
            weights = testCase.trapezoidalWeights(z);
            [transform,assessment] = basisSet.discreteTransform(z=z,weights=weights,nModes=6);

            testCase.verifyEmpty(assessment.weightFit)
            testCase.verifyEqual(transform.weights,weights,AbsTol=0)
            testCase.verifyEqual(assessment.candidateTransform.weights,weights,AbsTol=0)
            testCase.verifyEqual(assessment.requestedPointCount,length(z))
            testCase.verifyEqual(assessment.actualPointCount,length(z))
        end

        function unattainablePointCountsReportNearbyCounts(testCase)
            basisSet = testCase.constantFBasis(12);
            try
                basisSet.discreteTransform(nPoints=15);
                testCase.assertFail("Expected an unattainable exact point-count error.")
            catch exception
                testCase.verifyEqual(string(exception.identifier),"IMBasisSet:UnattainableDiscretePointCount")
                testCase.verifySubstring(exception.message,"exactly 15 points")
                testCase.verifySubstring(exception.message,"14")
            end
        end

        function pointAndBandSpecificationsAreUnambiguous(testCase)
            basisSet = testCase.constantFBasis(8);
            z = linspace(-4000,0,9).';
            weights = testCase.trapezoidalWeights(z);

            testCase.verifyError(@() basisSet.discreteTransform(),"IMBasisSet:InvalidDiscretePointSpecification")
            testCase.verifyError(@() basisSet.discreteTransform(nPoints=8,z=z),"IMBasisSet:InvalidDiscretePointSpecification")
            testCase.verifyError(@() basisSet.discreteTransform(weights=weights),"IMBasisSet:InvalidDiscretePointSpecification")
            testCase.verifyError(@() basisSet.discreteTransform(nPoints=8,nModes=6),"IMBasisSet:InvalidDiscretePointSpecification")
            testCase.verifyError(@() basisSet.discreteTransform(nPoints=[7 8]),"IMBasisSet:InvalidDiscreteTransformOption")
        end
    end

    methods (Access = private)
        function basisSet = constantFBasis(~,nModes)
            N0 = 5.2e-3;
            zDomain = [-4000 0];
            N2 = @(z) N0*N0*ones(size(z));
            solver = IMSolverSpectral(nEVP=128);
            evp = IMInternalModes.hydrostaticFModes(N2=N2,zDomain=zDomain);
            basisSet = solver.solveEVP(evp,nModes=nModes);
        end

        function weights = trapezoidalWeights(~,z)
            weights = [(z(2)-z(1))/2; (z(3:end)-z(1:end-2))/2; (z(end)-z(end-1))/2];
        end
    end
end
