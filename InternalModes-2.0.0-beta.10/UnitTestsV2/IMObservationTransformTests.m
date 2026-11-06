classdef IMObservationTransformTests < matlab.unittest.TestCase
    % Public generic recovery API, independent of scalar point/Gram metadata.
    methods (Test)
        function complexObservationsPreserveBatchDimensions(testCase)
            A = [1 1i; 2-1i 1; -1 2+3i; 0.5i -2];
            B = [1 0 2i 0; 0 1 0 -1i];
            metadata = struct(observations=["u" "v" "surface" "derivative"]);
            transform = IMDiscreteTransform.fromObservationMatrices(A,observationFactor=diag([1 2 3 4]),columnLabels=["wave:+" "wave:-"],physicalDual=B,metadata=metadata);
            coefficients = reshape(complex(1:12,12:-1:1),2,3,2);
            values = transform.transformBack(coefficients);
            testCase.verifySize(values,[4 3 2])
            testCase.verifyEqual(transform.transformForward(values),coefficients,AbsTol=1e-12)
            testCase.verifyEqual(transform.physicalDualMatrix,B)
            testCase.verifyEqual(transform.metadata,metadata)
            testCase.verifyEqual(transform.sampleCount,4)
            testCase.verifyEqual(transform.columnCount,2)
            testCase.verifyTrue(transform.success)
            testCase.verifyFalse(isprop(transform,"z"))
            testCase.verifyFalse(isprop(transform,"gramMatrix"))
            testCase.verifyFalse(isprop(transform,"inverseMatrix"))
            testCase.verifyFalse(isprop(transform,"sampledGramRank"))
        end

        function diagnosticFailurePreservesSynthesisAndPhysicalDual(testCase)
            A = ones(3,2); B = [1 0 0;0 1 0];
            transform = IMDiscreteTransform.fromObservationMatrices(A,observationFactor=eye(3),physicalDual=B,failureMode="diagnostic");
            testCase.verifyFalse(transform.success)
            testCase.verifyEmpty(transform.forwardMatrix)
            testCase.verifyEqual(transform.synthesisMatrix,A)
            testCase.verifyEqual(transform.physicalDualMatrix,B)
            testCase.verifyEqual(transform.transformBack([1;2]),3*ones(3,1))
            testCase.verifyError(@() transform.transformForward(ones(3,1)),"IMDiscreteTransform:UnavailableRecovery")
        end

        function declaredInactiveAndEmptyCoordinatesRemainAligned(testCase)
            A = [ones(3,1) zeros(3,1) (1:3).'];
            transform = IMDiscreteTransform.fromObservationMatrices(A,observationFactor=eye(3),activeColumnMask=[true false true],columnLabels=["F:0" "G:0" "internal:1"]);
            testCase.verifyEqual(transform.forwardMatrix(2,:),zeros(1,3))
            testCase.verifyEqual(transform.forwardMatrix*A,diag([1 0 1]),AbsTol=1e-12)
            empty = IMDiscreteTransform.fromObservationMatrices(zeros(3,0),observationFactor=eye(3));
            testCase.verifySize(empty.transformForward(zeros(3,2,4)),[0 2 4])
            testCase.verifyEqual(empty.transformBack(zeros(0,2,4)),zeros(3,2,4))
        end

        function outputSelectionKeepsThePreparedExtension(testCase)
            A = [1 1;0 1;1 0]; labels=["first" "second"];
            full = IMDiscreteTransform.fromObservationMatrices(A,observationFactor=eye(3),columnLabels=labels);
            selected = full.selectOutputs(1);
            refitted = IMDiscreteTransform.fromObservationMatrices(A(:,1),observationFactor=eye(3),columnLabels="first");
            testCase.verifyEqual(selected.synthesisMatrix,A(:,1))
            testCase.verifyEqual(selected.forwardMatrix,full.forwardMatrix(1,:))
            testCase.verifyLessThan(abs(selected.forwardMatrix*A(:,2)),1e-12)
            testCase.verifyGreaterThan(abs(refitted.forwardMatrix*A(:,2)),0.4)
            testCase.verifyEqual(selected.recoveryPolicy,"extendedSpan")
            testCase.verifyEqual(selected.recoverySpecification.extendedColumnLabels,labels)
            testCase.verifyLessThan(selected.diagnostics.roundTripError,1e-12)
            testCase.verifyLessThan(selected.diagnostics.extendedRejectionError,1e-12)
            testCase.verifyFalse(selected.diagnostics.fit.available)
            testCase.verifyTrue(selected.diagnostics.extended.fit.available)
            testCase.verifyTrue(isnan(selected.diagnostics.observationNoiseGain))
            testCase.verifyError(@() full.selectOutputs([1 1]),"IMDiscreteTransform:InvalidColumns")
        end

        function leadingDimensionsAreValidatedByTheGenericContract(testCase)
            transform = IMDiscreteTransform.fromObservationMatrices(eye(2),observationFactor=eye(2));
            testCase.verifyError(@() transform.transformForward(ones(3,1)),"IMDiscreteTransform:InvalidSampleCount")
            testCase.verifyError(@() transform.transformBack(ones(3,1)),"IMDiscreteTransform:InvalidCoefficientCount")
        end

        function numericalSpecificationPreservesFactorsAndComposesViews(testCase)
            A = [1 0 1;0 1 1;1 1 0;2 -1 1];
            D = diag([1 2 3 4]); F = [1 0 0 0;0 1 0 0;0 0 1 0];
            metadata = struct(outputSelection="caller-owned",experiment="factor fixture");
            fitted = IMDiscreteTransform.fromObservationMatrices(A(:,1:2),coefficientScale=[2;3],observationFactor=D,fitFactor=F,rankTolerance=1e-10,metadata=metadata);
            testCase.verifyEqual(fitted.recoverySpecification.coefficientScale,diag([2 3]))
            testCase.verifyEqual(fitted.recoverySpecification.observationFactor,D)
            testCase.verifyEqual(fitted.recoverySpecification.fitFactor,F)
            testCase.verifyEqual(fitted.recoverySpecification.rankTolerance,1e-10)
            testCase.verifyEqual(fitted.metadata,metadata)
            labels = ["first" "second" "third"]; retained = [3 1];
            extended = IMDiscreteTransform.fromObservationMatrices(A(:,retained),policy="extendedSpan",columnLabels=labels(retained),coefficientScale=[5;7],observationFactor=D,fitFactor=F,extendedSynthesis=A,extendedColumnLabels=labels,extendedActiveColumnMask=true(1,3),extendedCoefficientScale=[2;3;4],retainedColumnIndices=retained,metadata=metadata);
            testCase.verifyEqual(extended.recoverySpecification.extendedSynthesis,A)
            testCase.verifyEqual(extended.recoverySpecification.extendedCoefficientScale,diag([2 3 4]))
            testCase.verifyEqual(extended.recoverySpecification.extendedRecoverySpecification.fitFactor,F)
            firstView = extended.selectOutputs([2 1]);
            secondView = firstView.selectOutputs(2);
            testCase.verifyEqual(firstView.recoverySpecification.retainedColumnIndices,[1 3])
            testCase.verifyEqual(secondView.recoverySpecification.retainedColumnIndices,3)
            testCase.verifyEqual(secondView.recoverySpecification.extendedSynthesis,A)
            testCase.verifyEqual(secondView.recoverySpecification.extendedColumnLabels,labels)
            testCase.verifyEqual(secondView.metadata,metadata)
            testCase.verifyEmpty(secondView.recoverySpecification.coefficientScale)
            testCase.verifyFalse(secondView.recoverySpecification.retainedCoefficientScaleAvailable)
            testCase.verifyTrue(isnan(secondView.diagnostics.observationNoiseGain))
            testCase.verifyLessThan(secondView.diagnostics.extendedRejectionError,1e-12)
        end
    end
end
