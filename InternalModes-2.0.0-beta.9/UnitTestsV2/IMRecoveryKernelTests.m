classdef IMRecoveryKernelTests < matlab.unittest.TestCase
    % Mathematical controls for shared recovery, independent of mode solving.
    methods (Test)
        function directWeightedRankSurvivesRoundedNormalMatrix(testCase)
            A = [1 1;0 1e-8];
            prepared = prepareIMRecovery(A,observationFactor=eye(2));
            testCase.verifyEqual(rank(A'*A),1);
            testCase.verifyEqual(prepared.diagnostics.observation.rank,2);
            testCase.verifyEqual(prepared.diagnostics.fit.rank,2);
            testCase.verifyEqual(prepared.forwardMatrix*A,eye(2),AbsTol=1e-14);
        end

        function fittingAndObservationNormsRemainIndependent(testCase)
            A = [1;1];
            fitted = prepareIMRecovery(A,observationFactor=eye(2),fitFactor=[1 0]);
            optimal = prepareIMRecovery(A,observationFactor=eye(2));
            testCase.verifyEqual(fitted.forwardMatrix,[1 0],AbsTol=1e-14);
            testCase.verifyEqual(optimal.forwardMatrix,[0.5 0.5],AbsTol=1e-14);
            testCase.verifyFalse(fitted.diagnostics.noiseOptimalForObservationNorm);
            testCase.verifyTrue(optimal.diagnostics.noiseOptimalForObservationNorm);
            testCase.verifyEqual(fitted.diagnostics.observationNoiseGain,1,AbsTol=1e-14);
            testCase.verifyEqual(optimal.diagnostics.observationNoiseGain,1/sqrt(2),AbsTol=1e-14);
        end

        function zeroFitWeightsDoNotDeletePhysicalObservations(testCase)
            B = [2 3;4 5];
            prepared = prepareIMRecovery(eye(2),observationFactor=eye(2),fitFactor=[1 0],physicalDual=B,failureMode="diagnostic");
            testCase.verifyEqual(prepared.diagnostics.observation.rank,2);
            testCase.verifyEqual(prepared.diagnostics.fit.rank,1);
            testCase.verifyEqual(prepared.status,"fitRankDeficient");
            testCase.verifyEmpty(prepared.forwardMatrix);
            testCase.verifyEqual(prepared.physicalDual,B);
            testCase.verifyError(@() prepareIMRecovery(eye(2),observationFactor=eye(2),fitFactor=[1 0]),"IMRecovery:FitRankDeficient");
        end

        function algebraicAdaptersDoNotInventNoiseEvidence(testCase)
            A = [1 0;0 2;1 1];
            ls = prepareIMRecovery(A,fitFactor=eye(3));
            weak = prepareIMRecovery(A,policy="correctedPhysicalDual",physicalDual=[1 0 1;0 1 2]);
            testCase.verifyLessThan(norm(ls.forwardMatrix*A-eye(2)),1e-14);
            testCase.verifyLessThan(norm(weak.forwardMatrix*A-eye(2)),1e-14);
            testCase.verifyFalse(ls.diagnostics.observation.available);
            testCase.verifyFalse(weak.diagnostics.observation.available);
            testCase.verifyFalse(weak.diagnostics.supportsPrincipalGeometry);
            testCase.verifyTrue(isnan(ls.diagnostics.observationNoiseGain));
            testCase.verifyTrue(isnan(weak.diagnostics.observationNoiseGain));
        end

        function anExplicitEmptyFittingSpaceDoesNotDefaultToObservations(testCase)
            prepared = prepareIMRecovery(eye(2),observationFactor=eye(2),fitFactor=zeros(0,2),failureMode="diagnostic");
            testCase.verifyEqual(prepared.status,"fitRankDeficient");
            testCase.verifyEqual(prepared.diagnostics.fit.rank,0);
            testCase.verifyEqual(prepared.diagnostics.observation.rank,2);
            testCase.verifyEmpty(prepared.forwardMatrix);
        end

        function rowPreparationPreservesWeakMapAndRepairsRankDiagnosis(testCase)
            A = [eye(2);zeros(1,2)];
            B = [1e12 0 2e12;0 1 3];
            prepared = prepareIMRecovery(A,policy="correctedPhysicalDual",observationFactor=eye(3),physicalDual=B);
            testCase.verifyEqual(prepared.diagnostics.overlap.raw.rank,1);
            testCase.verifyEqual(prepared.diagnostics.test.rank,2);
            testCase.verifyEqual(prepared.diagnostics.overlap.rank,2);
            testCase.verifyEqual(prepared.forwardMatrix,[1 0 2;0 1 3],AbsTol=2e-14);
            testCase.verifyEqual(prepared.physicalDual,B);
            testCase.verifyLessThan(prepared.diagnostics.scaledDefiningEquationResidual,1e-14);
            transformedB = diag([1e-6*exp(0.8i),1e6*exp(-0.2i)])*B;
            transformedB = transformedB([2 1],:);
            transformed = prepareIMRecovery(A,policy="correctedPhysicalDual",observationFactor=eye(3),physicalDual=transformedB);
            testCase.verifyEqual(transformed.forwardMatrix,prepared.forwardMatrix,AbsTol=3e-14);
            testCase.verifyEqual(transformed.diagnostics.principalCosines,prepared.diagnostics.principalCosines,AbsTol=1e-14);
        end

        function sampleTestAndOverlapDeficienciesAreDistinct(testCase)
            A = [eye(2);zeros(1,2)];
            sample = prepareIMRecovery([ones(3,1),ones(3,1)],observationFactor=eye(3),failureMode="diagnostic");
            test = prepareIMRecovery(A,policy="correctedPhysicalDual",observationFactor=eye(3),physicalDual=[1 0 0;2 0 0],failureMode="diagnostic");
            overlap = prepareIMRecovery(A,policy="correctedPhysicalDual",observationFactor=eye(3),physicalDual=[1 0 0;0 0 1],failureMode="diagnostic");
            testCase.verifyEqual(sample.status,"sampleRankDeficient");
            testCase.verifyEqual(test.status,"testRankDeficient");
            testCase.verifyEqual(overlap.status,"overlapRankDeficient");
            testCase.verifyEqual(overlap.diagnostics.test.rank,2);
            testCase.verifyEqual(overlap.diagnostics.overlap.rank,1);
            testCase.verifyError(@() prepareIMRecovery(A,policy="correctedPhysicalDual",observationFactor=eye(3),physicalDual=[1 0 0;0 0 1]),"IMRecovery:OverlapRankDeficient");
        end

        function complexDenseCoordinatesAndUnitsAreCovariant(testCase)
            A = [1 0;0 2;1i 1;2 -1i];
            B = [1 0 1i 1;0 1 2 1i];
            C = [2 0.2i;0.1 0.5];
            D = diag([1 2 0.5 1.5]);
            T = [1+0.3i 0.2;0.1i 2-0.2i];
            U = diag([1e-4 20 0.2 3]);
            for policy = ["sampledLeastSquares","correctedPhysicalDual"]
                original = prepareIMRecovery(A,policy=policy,coefficientScale=C,observationFactor=D,physicalDual=B);
                coordinates = prepareIMRecovery(A*T,policy=policy,coefficientScale=T\C,observationFactor=D,physicalDual=T\B);
                units = prepareIMRecovery(U*A,policy=policy,coefficientScale=C,observationFactor=D/U,physicalDual=B/U);
                testCase.verifyEqual(coordinates.forwardMatrix,T\original.forwardMatrix,AbsTol=2e-13);
                testCase.verifyEqual(units.forwardMatrix*U,original.forwardMatrix,AbsTol=2e-13);
                testCase.verifyEqual(coordinates.diagnostics.observation.singularValues,original.diagnostics.observation.singularValues,AbsTol=1e-13);
                testCase.verifyEqual(units.diagnostics.observationNoiseGain,original.diagnostics.observationNoiseGain,AbsTol=1e-13);
            end
        end

        function observationAndColumnOrderingPreserveTheMap(testCase)
            A = [1 0;0 2;1i 1;2 -1i];
            B = [1 0 1i 1;0 1 2 1i];
            labels = ["second" "first"];
            p = [3 1 4 2]; c = [2 1];
            original = prepareIMRecovery(A,policy="correctedPhysicalDual",physicalDual=B,observationFactor=eye(4),columnLabels=labels);
            permuted = prepareIMRecovery(A(p,c),policy="correctedPhysicalDual",physicalDual=B(c,p),observationFactor=eye(4),columnLabels=labels(c));
            testCase.verifyEqual(permuted.forwardMatrix,original.forwardMatrix(c,p),AbsTol=1e-13);
        end

        function inactiveAlignmentAndEmptySpansAreExact(testCase)
            A = [zeros(3,1),[1;2;3]];
            prepared = prepareIMRecovery(A,activeColumnMask=[false true],coefficientScale=2,observationFactor=eye(3));
            testCase.verifyEqual(prepared.forwardMatrix(1,:),zeros(1,3));
            testCase.verifyEqual(prepared.forwardMatrix*A,diag([0 1]),AbsTol=1e-14);
            empty = prepareIMRecovery(zeros(3,0),observationFactor=eye(3));
            inactive = prepareIMRecovery(zeros(3,2),activeColumnMask=[false false],observationFactor=eye(3));
            testCase.verifySize(empty.forwardMatrix,[0 3]);
            testCase.verifyEqual(inactive.forwardMatrix,zeros(2,3));
            testCase.verifyEqual(inactive.diagnostics.fit.rank,0);
            testCase.verifyEqual(inactive.diagnostics.scaledRoundTripError,0);
            activeZero = prepareIMRecovery(zeros(3,1),observationFactor=eye(3),failureMode="diagnostic");
            testCase.verifyEqual(activeZero.status,"sampleRankDeficient");
        end

        function extendedSpanRejectsExactlyItsRepresentedTail(testCase)
            extended = [zeros(4,1),eye(4)];
            labels = ["G0" "one" "two" "three" "four"];
            retained = [1 3 2];
            prepared = prepareIMRecovery(extended(:,retained),policy="extendedSpan",columnLabels=labels(retained),activeColumnMask=[false true true],observationFactor=eye(4),extendedSynthesis=extended,extendedColumnLabels=labels,extendedActiveColumnMask=[false true true true true],retainedColumnIndices=retained);
            testCase.verifyEqual(prepared.forwardMatrix*extended(:,retained),diag([0 1 1]),AbsTol=1e-14);
            testCase.verifyEqual(prepared.forwardMatrix*extended(:,[4 5]),zeros(3,2),AbsTol=1e-14);
            testCase.verifyEqual(prepared.diagnostics.extendedRejectionError,0,AbsTol=1e-14);
            testCase.verifyEqual(prepared.diagnostics.extended.fit.rank,4);
            testCase.verifyFalse(prepared.diagnostics.noiseOptimalForObservationNorm);
        end

        function extendedFailureCannotBeApplied(testCase)
            extended = [1 1;0 0];
            prepared = prepareIMRecovery(extended(:,1),policy="extendedSpan",columnLabels="one",observationFactor=eye(2),extendedSynthesis=extended,extendedColumnLabels=["one" "two"],extendedActiveColumnMask=[true true],retainedColumnIndices=1,failureMode="diagnostic");
            testCase.verifyEqual(prepared.status,"extendedRecoveryFailed");
            testCase.verifyEqual(prepared.diagnostics.extended.observation.rank,1);
            testCase.verifyEmpty(prepared.forwardMatrix);
        end

        function leastSquaresNoiseOrderingDoesNotOrderOmittedResponses(testCase)
            A = [eye(2);zeros(1,2)];
            B = [1 0 2;0 1 3];
            ls = prepareIMRecovery(A,observationFactor=eye(3));
            weak = prepareIMRecovery(A,policy="correctedPhysicalDual",observationFactor=eye(3),physicalDual=B);
            outputFactors = {[2 1i;0 0],[0 3]};
            for j = 1:numel(outputFactors)
                H = outputFactors{j};
                covarianceDifference = H*(weak.forwardMatrix*weak.forwardMatrix'-ls.forwardMatrix*ls.forwardMatrix')*H';
                testCase.verifyGreaterThanOrEqual(min(eig((covarianceDifference+covarianceDifference')/2)),-1e-13);
                testCase.verifyGreaterThanOrEqual(norm(H*weak.forwardMatrix,2),norm(H*ls.forwardMatrix,2)-1e-13);
            end
            omitted = [-2;-3;1];
            testCase.verifyLessThan(norm(weak.forwardMatrix*omitted),1e-14);
            testCase.verifyGreaterThan(norm(ls.forwardMatrix*omitted),1);
        end

        function nonmodalPhysicalDualAndRecoveryRemainDifferent(testCase)
            A = [eye(2);zeros(1,2)];
            B = [2 0 4;0 3 9];
            prepared = prepareIMRecovery(A,policy="correctedPhysicalDual",observationFactor=eye(3),physicalDual=B);
            source = [1;2;3i];
            testCase.verifyEqual(prepared.physicalDual,B);
            testCase.verifyGreaterThan(norm(B*source-prepared.forwardMatrix*source),1);
            samples = A*[1 2i;3-1i 4];
            testCase.verifyEqual(prepared.forwardMatrix*samples,[1 2i;3-1i 4],AbsTol=2e-14);
        end

        function invalidContractsFailBeforeFactorization(testCase)
            testCase.verifyError(@() prepareIMRecovery(eye(2),columnLabels=["same" "same"],observationFactor=eye(2)),"IMRecovery:InvalidColumnLabels");
            testCase.verifyError(@() prepareIMRecovery(eye(2),activeColumnMask=[true false],observationFactor=eye(2)),"IMRecovery:InvalidInactiveColumn");
            testCase.verifyError(@() prepareIMRecovery(eye(2),observationFactor=[1 0]),"IMRecovery:SingularCoordinateFactor");
            testCase.verifyError(@() prepareIMRecovery(eye(2),coefficientScale=eye(3),observationFactor=eye(2)),"IMRecovery:InvalidCoordinateFactor");
            testCase.verifyError(@() prepareIMRecovery(eye(2),fitFactor=ones(2,3)),"IMRecovery:InvalidFitFactor");
            testCase.verifyError(@() prepareIMRecovery(eye(2),observationFactor=eye(2),rankTolerance=1),"IMRecovery:InvalidRankTolerance");
        end

        function finiteExtremeRowsDoNotOverflowTestNormalization(testCase)
            A = [1;0;0;0];
            B = 1e308*ones(1,4);
            prepared = prepareIMRecovery(A,policy="correctedPhysicalDual",physicalDual=B,observationFactor=eye(4));
            testCase.verifyTrue(prepared.success);
            testCase.verifyEqual(prepared.forwardMatrix,ones(1,4),AbsTol=2e-14);
            testCase.verifyEqual(prepared.physicalDual,B);
            testCase.verifyEqual(prepared.diagnostics.test.rank,1);
            testCase.verifyTrue(isinf(prepared.diagnostics.test.rowNorms));
            testCase.verifyTrue(all(isfinite(prepared.diagnostics.test.rowNormalizedLengths)));
        end

        function finiteExtremeSynthesisIsScaledBeforeSvdAndSolve(testCase)
            A = 1e308*ones(4,1);
            for policy = ["sampledLeastSquares","correctedPhysicalDual"]
                prepared = prepareIMRecovery(A,policy=policy,observationFactor=eye(4),physicalDual=ones(1,4));
                testCase.verifyTrue(prepared.success);
                testCase.verifyEqual(1e308*prepared.forwardMatrix,0.25*ones(1,4),AbsTol=2e-15);
                testCase.verifyLessThan(prepared.diagnostics.roundTripError,3e-15);
                testCase.verifyEqual(prepared.diagnostics.observation.rank,1);
                testCase.verifyFalse(prepared.diagnostics.observation.singularValuesRepresentable);
                testCase.verifyTrue(all(isfinite(prepared.diagnostics.observation.scaledSingularValues)));
                testCase.verifyEqual(prepared.diagnostics.observation.status,"scaledSpectrum");
            end
        end

        function extremeWeakProductsRetainScaledEvidenceAndResiduals(testCase)
            for alpha = [1e-200 1e200]
                A = alpha*eye(2);
                prepared = prepareIMRecovery(A,policy="correctedPhysicalDual",physicalDual=A,observationFactor=eye(2));
                testCase.verifyEqual(alpha*prepared.forwardMatrix,eye(2),AbsTol=2e-15);
                testCase.verifyEqual(prepared.diagnostics.overlap.raw.rank,2);
                testCase.verifyFalse(prepared.diagnostics.overlap.raw.singularValuesRepresentable);
                testCase.verifyLessThan(prepared.diagnostics.scaledDefiningEquationResidual,2e-15);
                testCase.verifyEqual(prepared.diagnostics.definingEquationResidual,0,AbsTol=0);
                testCase.verifyEqual(prepared.diagnostics.roundTripError,0,AbsTol=0);
            end
        end

        function diagnosticRangeLimitsAndResolvedSpecificationAreExplicit(testCase)
            for alpha = [1e-200 1e200]
                prepared = prepareIMRecovery(eye(2),coefficientScale=alpha*eye(2),observationFactor=alpha*eye(2));
                testCase.verifyEqual(prepared.forwardMatrix,eye(2),AbsTol=2e-15);
                testCase.verifyEqual(prepared.diagnostics.observation.rank,2);
                testCase.verifyEqual(prepared.diagnostics.observationNoiseGainStatus,"outsideFloatingPointRange");
                testCase.verifyEqual(prepared.specification.coefficientScale,alpha*eye(2));
                testCase.verifyEqual(prepared.specification.observationFactor,alpha*eye(2));
                testCase.verifyEqual(prepared.specification.fitFactor,alpha*eye(2));
            end
            algebraic = prepareIMRecovery(eye(2),fitFactor=eye(2));
            testCase.verifyEmpty(algebraic.specification.observationFactor);
            failed = prepareIMRecovery(realmin*eps,observationFactor=1,failureMode="diagnostic");
            testCase.verifyEqual(failed.status,"unrepresentableRecovery");
            testCase.verifyEmpty(failed.forwardMatrix);
            testCase.verifyEqual(failed.specification.observationFactor,1);
        end

        function weakDiagnosticsRemainSmallForManyObservations(testCase)
            count = 20000;
            A = [ones(count,1),linspace(-1,1,count).'];
            B = A'/count;
            prepared = prepareIMRecovery(A,policy="correctedPhysicalDual",physicalDual=B);
            testCase.verifySize(prepared.forwardMatrix,[2 count]);
            testCase.verifyEqual(prepared.forwardMatrix*A,eye(2),AbsTol=3e-13);
            testCase.verifyLessThan(prepared.diagnostics.scaledDefiningEquationResidual,3e-13);
            testCase.verifyLessThan(prepared.diagnostics.definingEquationResidual,3e-13);
            testCase.verifyEmpty(prepared.specification.observationFactor);
        end

        function unrepresentableWhiteningReturnsAnExplicitFailedDiagnostic(testCase)
            A = diag([1e-308 1e10]);
            D = diag([1e308 1e-10]);
            B = eye(2);
            prepared = prepareIMRecovery(A,policy="correctedPhysicalDual",physicalDual=B,observationFactor=D,failureMode="diagnostic");
            testCase.verifyFalse(prepared.success);
            testCase.verifyEqual(prepared.status,"unrepresentableIntermediate");
            testCase.verifyEmpty(prepared.forwardMatrix);
            testCase.verifyEqual(prepared.physicalDual,B);
            testCase.verifyEqual(prepared.specification.observationFactor,D);
            testCase.verifyError(@() prepareIMRecovery(A,policy="correctedPhysicalDual",physicalDual=B,observationFactor=D),"IMRecovery:UnrepresentableIntermediate");
        end
    end
end
