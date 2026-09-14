classdef IMProjectionTests < matlab.unittest.TestCase
    properties (Access = private)
        originalPath
    end

    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            testCase.originalPath = path;
            addpath(fileparts(fileparts(mfilename("fullpath"))));
        end
    end

    methods (TestClassTeardown)
        function restoreRepositoryPath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function signedPairingAndPositiveErrorNormStayDistinct(testCase)
            projection = IMProjection(ones(2,1),diag([2 -1]),1,majorantGramMatrix=3,columnLabels="external");
            testCase.verifyEqual(projection.project([1;0]),2,AbsTol=1e-14)
            testCase.verifyEqual(projection.productError([1;0],0,1),2*sqrt(3),AbsTol=1e-14)
            testCase.verifyEqual(projection.leakage([1;0],1),2*sqrt(3),AbsTol=1e-14)
            testCase.verifyNotEqual(projection.project([1;0]),2/3)
        end

        function complexPhaseLeavesPositiveNormsInvariant(testCase)
            projection = IMProjection(ones(2,1),diag([2 -1]),1,majorantGramMatrix=3);
            phase = exp(1i*0.7);
            testCase.verifyEqual(projection.project(phase*[1;0]),2*phase,AbsTol=1e-14)
            testCase.verifyEqual(projection.leakage(phase*[1;0],1),2*sqrt(3),AbsTol=1e-14)
            testCase.verifyEqual(projection.productError(phase*[1;0],phase,1),sqrt(3),AbsTol=1e-14)
        end

        function rankDefectAndInactiveColumnsRemainExplicit(testCase)
            projection = IMProjection([1 1;1 1],eye(2),eye(2));
            testCase.verifyEqual(projection.sampledGramRank,1)
            testCase.verifyEqual(projection.gramError,Inf)
            testCase.verifyEqual(projection.project([2;2]),[1;1],AbsTol=1e-14)
            testCase.verifyEqual(projection.roundTripError,1,AbsTol=1e-14)
            inactive = IMProjection([0 1;0 2],eye(2),diag([0 5]),activeColumnMask=[false true],columnLabels=["barotropic" "1"]);
            testCase.verifyEqual(inactive.forwardMatrix(1,:),[0 0],AbsTol=0)
            testCase.verifyLessThan(inactive.roundTripError,1e-14)
            testCase.verifyEqual(inactive.prefix(1).activeColumnMask,false)
            testCase.verifyEqual(inactive.prefix(1).gramError,0)
            testCase.verifyEqual(inactive.prefix(0).sampledBasis,zeros(2,0))
            testCase.verifyEqual(inactive.prefix(2).columnLabels,["barotropic" "1"])
        end

        function denseSignedTargetsRequireAnExplicitPositiveMajorant(testCase)
            target = [1 2;2 -1];
            testCase.verifyError(@() IMProjection(eye(2),target,target),"IMProjection:MajorantRequired")
            projection = IMProjection(eye(2),target,target,majorantGramMatrix=3*eye(2));
            testCase.verifyEqual(projection.forwardMatrix,eye(2),AbsTol=1e-14)
            testCase.verifyEqual(projection.gramError,0,AbsTol=0)
            testCase.verifyFalse(projection.targetGramIsPositiveDefinite)
            testCase.verifyEqual(projection.productError([1;2],zeros(2,1),5),sqrt(3),AbsTol=1e-14)
            testCase.verifyError(@() IMProjection(eye(2),target,target,majorantGramMatrix=target),"IMProjection:InvalidMajorant")
        end

        function zeroProductsDistinguishZerosFromInconsistentReferenceData(testCase)
            projection = IMProjection([1;0],eye(2),1);
            sampled = [0 0 0 1e-100;0 1 0 0];
            reference = [0 0 1 0];
            testCase.verifyEqual(projection.productError(sampled,reference,[0;0;0;1e-200]),[0 Inf Inf 1],AbsTol=1e-14)
            testCase.verifyEqual(projection.leakage(sampled(:,1:2),[0;0]),[0 Inf])
        end

        function cosineAndSineProductsMatchExactIdentitiesAndDetectAliasing(testCase)
            z = linspace(0,pi,13).';
            weights = [0.5;ones(11,1);0.5]*(pi/12);
            labels = 0:8;
            cosines = cos(z*labels);
            sines = sin(z*labels);
            f = IMProjection(cosines,diag(weights),diag([pi repmat(pi/2,1,8)]));
            g = IMProjection(sines,diag(weights),diag([0 repmat(pi/2,1,8)]),activeColumnMask=[false true(1,8)]);
            ff = zeros(9,1);
            ff([2 6]) = 0.5;
            gg = ff;
            gg(6) = -0.5;
            fg = ff;
            testCase.verifyEqual(f.project(cosines(:,3).*cosines(:,4)),ff,AbsTol=2e-14)
            testCase.verifyEqual(f.project(sines(:,3).*sines(:,4)),gg,AbsTol=2e-14)
            testCase.verifyEqual(g.project(cosines(:,3).*sines(:,4)),fg,AbsTol=2e-14)
            testCase.verifyLessThan(f.productError(cosines(:,3).*cosines(:,4),ff,pi/4),2e-14)
            exact = zeros(9,1);
            exact(1) = 0.5;
            testCase.verifyGreaterThan(f.productError(cosines(:,9).^2,exact,3*pi/8),0.5)
        end

        function numericalScalarAndAlignedTransformsAgree(testCase)
            z = linspace(-1000,0,17).';
            weights = [0.5;ones(15,1);0.5]*(1000/16);
            solver = IMSolverSpectral(nEVP=64);
            evp = IMEigenvalueProblem(zDomain=[-1000 0],p=1,q=0,r=1,surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
            scalar = solver.solveEVP(evp,nModes=4);
            oldScalar = scalar.discreteTransform(z=z,weights=weights,nModes=4);
            projection = IMProjection(oldScalar.inverseMatrix,oldScalar.metricMatrix,oldScalar.targetGramMatrix);
            testCase.verifyEqual(projection.forwardMatrix,oldScalar.forwardMatrix,AbsTol=2e-13)
            testCase.verifyEqual(projection.gramError,oldScalar.relativeGramOperatorError,AbsTol=2e-13)
            testCase.verifyEqual(projection.roundTripError,oldScalar.roundTripError,AbsTol=2e-13)
            family = solver.solveEVP(IMInternalModes.hydrostaticFModes(N2=@(z) 1e-4*ones(size(z)),zDomain=[-1000 0]),nModes=4);
            oldAligned = family.discreteTransform(z=z,weights=weights,nModes=4,variables=["F" "G"]);
            for variable = ["F" "G"]
                data = oldAligned.channelDiagnostics(variable=variable);
                projection = IMProjection(oldAligned.inverseMatrix(variable=variable),oldAligned.metricMatrix(variable=variable),oldAligned.targetGramMatrix(variable=variable),majorantGramMatrix=oldAligned.targetMajorantGramMatrix(variable=variable),activeColumnMask=oldAligned.activeModeMask(variable=variable));
                testCase.verifyEqual(projection.forwardMatrix,oldAligned.forwardMatrix(variable=variable),AbsTol=2e-13)
                testCase.verifyEqual(projection.gramError,data.relativeGramOperatorError,AbsTol=2e-13)
                testCase.verifyEqual(projection.roundTripError,data.roundTripError,AbsTol=2e-13)
            end
        end

        function prescribedComplexDualUsesItsSuppliedCoefficientSystem(testCase)
            pairing = [1 2i 3;2 -1i 1];
            sampleGram = [2 0.5;0.5 -3];
            targetGram = [3 0.25;0.25 -2];
            majorant = [4 0.2;0.2 3];
            provenance = struct(source="physical-energy-dual",component="u");
            projection = IMProjection.fromPrescribedDual(pairing,sampleGram,targetGram,majorantGramMatrix=majorant,columnLabels=["mode1+" "mode1-"],provenance=provenance);
            values = [1+2i;3-1i;-2];
            expected = sampleGram \ (pairing*values);
            testCase.verifyEqual(projection.project(values),expected,AbsTol=1e-14)
            testCase.verifyEqual(projection.sampleCount,3)
            testCase.verifyEqual(projection.columnCount,2)
            testCase.verifyEmpty(projection.sampledBasis)
            testCase.verifyEmpty(projection.metricMatrix)
            testCase.verifyEqual(projection.projectionKind,"prescribedDual")
            testCase.verifyFalse(projection.supportsGramAssessment)
            testCase.verifyFalse(projection.supportsRoundTripAssessment)
            testCase.verifyTrue(isnan(projection.gramError))
            testCase.verifyTrue(isnan(projection.roundTripError))
            testCase.verifyTrue(isnan(projection.inverseMatrixConditionNumber))
            testCase.verifyEqual(projection.gramMatrix,sampleGram,AbsTol=0)
            reference = targetGram \ [1+1i;2-1i];
            expectedError = sqrt(real((expected-reference)'*majorant*(expected-reference))/7);
            testCase.verifyEqual(projection.productError(values,reference,7),expectedError,AbsTol=1e-14)
            prefix = projection.prefix(1);
            testCase.verifyEqual(prefix.project(values),(pairing(1,:)*values)/sampleGram(1,1),AbsTol=1e-14)
            testCase.verifyEqual(prefix.provenance,provenance)
            testCase.verifyEqual(prefix.columnLabels,"mode1+")
            testCase.verifyTrue(isnan(prefix.gramError))
            testCase.verifyEqual(projection.prefix(0).forwardMatrix,zeros(0,3))
        end

        function prescribedDualZeroInactiveAndRankDefectRemainExplicit(testCase)
            provenance = struct(source="test-dual");
            projection = IMProjection.fromPrescribedDual([0 0;1 1i],diag([0 -2]),diag([0 -2]),activeColumnMask=[false true],provenance=provenance);
            testCase.verifyEqual(projection.project([1;1]),[0;-(1+1i)/2],AbsTol=1e-14)
            testCase.verifyEqual(projection.productError(zeros(2,1),zeros(2,1),0),0)
            testCase.verifyEqual(projection.productError([0;1],zeros(2,1),0),Inf)
            defect = IMProjection.fromPrescribedDual(eye(2),ones(2),eye(2),provenance=provenance);
            testCase.verifyEqual(defect.sampledGramRank,1)
            testCase.verifyTrue(isnan(defect.gramError))
            testCase.verifyEqual(defect.forwardMatrix,pinv(ones(2)),AbsTol=1e-14)
            offDiagonal = IMProjection.fromPrescribedDual(eye(2),[0 1;1 0],[0 1;1 0],majorantGramMatrix=eye(2),provenance=provenance);
            testCase.verifyEqual(offDiagonal.project([1;2]),[2;1],AbsTol=0)
        end

        function prescribedDualRejectsInvalidContracts(testCase)
            provenance = struct(source="test-dual");
            testCase.verifyError(@() IMProjection.fromPrescribedDual(eye(2),1,eye(2),provenance=provenance),"IMProjection:InvalidShape")
            testCase.verifyError(@() IMProjection.fromPrescribedDual(eye(2),[1 1;0 1],eye(2),provenance=provenance),"IMProjection:NonSymmetricMatrix")
            testCase.verifyError(@() IMProjection.fromPrescribedDual(eye(2),eye(2),eye(2),majorantGramMatrix=diag([1 -1]),provenance=provenance),"IMProjection:InvalidMajorant")
            testCase.verifyError(@() IMProjection.fromPrescribedDual(eye(2),eye(2),eye(2),activeColumnMask=[false true],provenance=provenance),"IMProjection:InvalidInactiveColumn")
            testCase.verifyError(@() IMProjection.fromPrescribedDual(eye(2),eye(2),eye(2),provenance=struct()),"IMProjection:MissingProvenance")
        end

        function publicInputsRejectInvalidMathematicalData(testCase)
            testCase.verifyError(@() IMProjection(eye(2),1,eye(2)),"IMProjection:InvalidShape")
            testCase.verifyError(@() IMProjection(eye(2),[1 1;0 1],eye(2)),"IMProjection:NonSymmetricMatrix")
            testCase.verifyError(@() IMProjection(eye(2),eye(2),diag([0 1])),"IMProjection:InvalidTargetNorm")
            testCase.verifyError(@() IMProjection(eye(2),eye(2),eye(2),activeColumnMask=[false true]),"IMProjection:InvalidInactiveColumn")
            projection = IMProjection(eye(2),eye(2),eye(2));
            testCase.verifyError(@() projection.project(ones(3,1)),"IMProjection:InvalidShape")
            testCase.verifyError(@() projection.productError(ones(2,1),ones(3,1),1),"IMProjection:InvalidShape")
            testCase.verifyError(@() projection.leakage(ones(2,1),[1;2]),"IMProjection:InvalidShape")
            testCase.verifyError(@() projection.prefix(3),"IMProjection:InvalidColumnCount")
        end
    end
end
