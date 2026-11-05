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
        function unrepresentablePhysicalDualDoesNotVetoPositiveRecovery(testCase)
            projection = IMProjection(1e-100,1e100,1e-310,majorantGramMatrix=1e-310);
            testCase.verifyEqual(projection.physicalDualStatus,"unrepresentableCoefficientSystem")
            testCase.verifyEmpty(projection.physicalDualMatrix)
            testCase.verifyTrue(projection.recovery.success)
            testCase.verifyEqual(projection.recovery.transformForward(1e-100),1,AbsTol=1e-14)
            testCase.verifyError(@() projection.projectPhysical(1),"IMProjection:UnavailablePhysicalDual")

            source = IMProjection.fromPrescribedDual(1,1e-310,1,majorantGramMatrix=1,provenance=struct(kind="rangeControl"));
            testCase.verifyEqual(source.physicalDualStatus,"unrepresentableCoefficientSystem")
            testCase.verifyEmpty(source.recovery)
            testCase.verifyEmpty(source.physicalDualMatrix)
            testCase.verifyError(@() source.projectPhysical(1),"IMProjection:UnavailablePhysicalDual")
            testCase.verifyError(@() source.projectPhysical(1),"IMProjection:UnavailablePhysicalDual")
            testCase.verifyEqual(source.samplePairingMatrix,1)

            tinySource = IMProjection.fromPrescribedDual(1e-310,1e100,1,majorantGramMatrix=1,provenance=struct(kind="underflowControl"));
            testCase.verifyEqual(tinySource.physicalDualStatus,"unrepresentableCoefficientSystem")
            testCase.verifyEmpty(tinySource.physicalDualMatrix)
            testCase.verifyError(@() tinySource.projectPhysical(1e100),"IMProjection:UnavailablePhysicalDual")
            tinyPhysical = IMProjection(1,1e-310,1e100,majorantGramMatrix=1e100);
            testCase.verifyEqual(tinyPhysical.physicalDualStatus,"unrepresentableCoefficientSystem")
            testCase.verifyTrue(tinyPhysical.recovery.success)
            testCase.verifyEqual(tinyPhysical.recovery.transformForward(1),1,AbsTol=1e-14)

            physicalOnly = IMProjection(1e-100,-1e100,-1e-310,majorantGramMatrix=1e-310,prepareRecovery=false);
            testCase.verifyEqual(physicalOnly.physicalDualStatus,"unrepresentableCoefficientSystem")
            testCase.verifyEmpty(physicalOnly.recovery)
            testCase.verifyError(@() IMProjection(1e-100,-1e100,-1e-310,majorantGramMatrix=1e-310),"IMProjection:UnavailablePhysicalDual")
        end

        function signedPairingAndPositiveErrorNormStayDistinct(testCase)
            projection = IMProjection(ones(2,1),diag([2 -1]),1,majorantGramMatrix=3,columnLabels="external");
            testCase.verifyEqual(projection.recovery.transformForward([1;0]),2,AbsTol=1e-14)
            testCase.verifyEqual(projection.productError([1;0],0,1),2*sqrt(3),AbsTol=1e-14)
            testCase.verifyEqual(projection.leakage([1;0],1),2*sqrt(3),AbsTol=1e-14)
            testCase.verifyNotEqual(projection.recovery.transformForward([1;0]),2/3)
        end

        function complexPhaseLeavesPositiveNormsInvariant(testCase)
            projection = IMProjection(ones(2,1),diag([2 -1]),1,majorantGramMatrix=3);
            phase = exp(1i*0.7);
            testCase.verifyEqual(projection.recovery.transformForward(phase*[1;0]),2*phase,AbsTol=1e-14)
            testCase.verifyEqual(projection.leakage(phase*[1;0],1),2*sqrt(3),AbsTol=1e-14)
            testCase.verifyEqual(projection.productError(phase*[1;0],phase,1),sqrt(3),AbsTol=1e-14)
        end

        function rankDefectAndInactiveColumnsRemainExplicit(testCase)
            projection = IMProjection([1 1;1 1],eye(2),eye(2));
            testCase.verifyEqual(projection.sampledGramRank,1)
            testCase.verifyEqual(projection.gramError,3,AbsTol=1e-14)
            testCase.verifyFalse(projection.recovery.success)
            testCase.verifyEmpty(projection.recovery.forwardMatrix)
            testCase.verifyError(@() projection.recovery.transformForward([2;2]),"IMDiscreteTransform:UnavailableRecovery")
            testCase.verifyEqual(projection.projectPhysical([2;2]),[4;4],AbsTol=1e-14)
            testCase.verifyEqual(projection.roundTripError,Inf)
            inactive = IMProjection([0 1;0 2],eye(2),diag([0 5]),activeColumnMask=[false true],columnLabels=["barotropic" "1"]);
            testCase.verifyEqual(inactive.recovery.forwardMatrix(1,:),[0 0],AbsTol=0)
            testCase.verifyLessThan(inactive.roundTripError,1e-14)
            testCase.verifyEqual(inactive.prefix(1).activeColumnMask,false)
            testCase.verifyEqual(inactive.prefix(1).gramError,0)
            testCase.verifyEqual(inactive.prefix(0).sampledBasis,zeros(2,0))
            testCase.verifyEqual(inactive.prefix(2).columnLabels,["barotropic" "1"])
        end

        function weakStratificationKeepsTheAnalyticalDCTRank(testCase)
            % Physical F norms span many orders of magnitude even though
            % the energy-normalized columns are orthogonal on this grid.
            D = 1000; g = 9.81; nPoints = 65;
            z = linspace(0,D,nPoints).';
            w = D/(nPoints-1)*[0.5;ones(nPoints-2,1);0.5];
            j = 1:nPoints-2;
            for N0 = [1e-6 1e-7 1e-10]
                h = N0^2*D^2./(g*(j*pi).^2);
                H = [D h];
                Phi = [ones(nPoints,1) cos(z*(j*pi/D)).*sqrt(2*h/D)];
                projection = IMProjection(Phi,diag(w),diag(H));
                testCase.verifyEqual(projection.sampledGramRank,nPoints-1)
                testCase.verifyTrue(projection.targetGramIsPositiveDefinite)
                testCase.verifyLessThan(projection.gramError,1e-12)
                coefficients = (1+1i*sin(0:nPoints-2)).'/sqrt(H).';
                recovered = projection.recovery.transformForward(Phi*coefficients);
                testCase.verifyLessThan(norm(sqrt(H).'.*(recovered-coefficients))/norm(sqrt(H).'.*coefficients),1e-12)
            end
        end

        function signedDensePairingsAreInvariantToColumnUnits(testCase)
            metric = [2 0.3 0.2;0.3 -1 0.4;0.2 0.4 3];
            majorant = [4 0.2 0.1;0.2 2 0.3;0.1 0.3 5];
            columnScales = [1e-12 1 1e8];
            scaling = diag(columnScales);
            projection = IMProjection(scaling,metric,scaling*metric*scaling,majorantGramMatrix=scaling*majorant*scaling);
            testCase.verifyEqual(projection.sampledGramRank,3)
            testCase.verifyFalse(projection.targetGramIsPositiveDefinite)
            values = [1;2i;3-1i];
            testCase.verifyEqual(columnScales.'.*projection.recovery.transformForward(values),values,AbsTol=1e-12)
            pairing = [1 2i;3 -1;2i 4];
            prescribed = IMProjection.fromPrescribedDual(scaling*pairing,scaling*metric*scaling,scaling*metric*scaling,majorantGramMatrix=scaling*majorant*scaling,provenance=struct(source="scaled signed dual"));
            testCase.verifyEqual(prescribed.sampledGramRank,3)
            testCase.verifyFalse(prescribed.targetGramIsPositiveDefinite)
            testCase.verifyEqual(columnScales.'.*prescribed.projectPhysical([1;2i]),metric\(pairing*[1;2i]),AbsTol=1e-12)
        end

        function zeroWeightsIgnoreSamplesWithDenseCoefficientMetric(testCase)
            z = linspace(-1,1,7).';
            scales = [1e-10 1 1e8];
            Phi = [ones(size(z)) z z.^2].*scales;
            w = [0;ones(5,1);0];
            H = Phi.'*(w.*Phi);
            projection = IMProjection(Phi,diag(w),H,majorantGramMatrix=H);
            coefficients = [1+1i;2;3]./scales.';
            values = Phi*coefficients;
            values([1 end]) = 1e12;
            testCase.verifyEqual(projection.sampledGramRank,3)
            testCase.verifyTrue(projection.targetGramIsPositiveDefinite)
            testCase.verifyEqual(projection.recovery.forwardMatrix(:,[1 end]),zeros(3,2),AbsTol=0)
            testCase.verifyEqual(scales.'.*projection.recovery.transformForward(values),scales.'.*coefficients,AbsTol=1e-12)
        end

        function deficientPairingPreservesUnitsAndInactiveZeros(testCase)
            Phi = [0 1 1;0 2 2;0 -1 -1];
            H = diag([0 1 1]);
            scales = [1 1e-10 1e5];
            base = IMProjection(Phi,eye(3),H,activeColumnMask=[false true true]);
            scaled = IMProjection(Phi.*scales,eye(3),diag(scales)*H*diag(scales),activeColumnMask=[false true true]);
            values = [1+1i;2;3-1i];
            testCase.verifyEqual(scaled.sampledGramRank,1)
            testCase.verifyEqual(scaled.gramError,base.gramError,AbsTol=1e-12)
            testCase.verifyFalse(scaled.recovery.success)
            testCase.verifyEmpty(scaled.recovery.forwardMatrix)
            testCase.verifyError(@() scaled.recovery.transformForward(values),"IMDiscreteTransform:UnavailableRecovery")
            testCase.verifyEqual(scales.'.*scaled.projectPhysical(values),base.projectPhysical(values),AbsTol=1e-12)
            testCase.verifyEqual(scaled.prefix(1).recovery.forwardMatrix,zeros(1,3),AbsTol=0)
        end

        function denseSignedTargetsRequireAnExplicitPositiveMajorant(testCase)
            target = [1 2;2 -1];
            testCase.verifyError(@() IMProjection(eye(2),target,target),"IMProjection:MajorantRequired")
            projection = IMProjection(eye(2),target,target,majorantGramMatrix=3*eye(2));
            testCase.verifyEqual(projection.recovery.forwardMatrix,eye(2),AbsTol=1e-14)
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
            testCase.verifyEqual(f.recovery.transformForward(cosines(:,3).*cosines(:,4)),ff,AbsTol=2e-14)
            testCase.verifyEqual(f.recovery.transformForward(sines(:,3).*sines(:,4)),gg,AbsTol=2e-14)
            testCase.verifyEqual(g.recovery.transformForward(cosines(:,3).*sines(:,4)),fg,AbsTol=2e-14)
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
            projection = IMProjection(oldScalar.inverseMatrix,oldScalar.metricMatrix,oldScalar.targetGramMatrix,majorantGramMatrix=oldScalar.projection.majorantGramMatrix);
            testCase.verifyEqual(projection.recovery.forwardMatrix,oldScalar.forwardMatrix,AbsTol=2e-13)
            testCase.verifyEqual(projection.gramError,oldScalar.relativeGramOperatorError,AbsTol=2e-13)
            testCase.verifyEqual(projection.roundTripError,oldScalar.roundTripError,AbsTol=2e-13)
            family = solver.solveEVP(IMInternalModes.hydrostaticFModes(rotationRate=0,N2=@(z) 1e-4*ones(size(z)),zDomain=[-1000 0],surfaceBoundary=IMBoundaryCondition.neumann()),nModes=4);
            oldAligned = family.discreteTransform(z=z,weights=weights,nModes=4,variables=["F" "G"]);
            for variable = ["F" "G"]
                data = oldAligned.channelDiagnostics(variable=variable);
                projection = IMProjection(oldAligned.inverseMatrix(variable=variable),oldAligned.metricMatrix(variable=variable),oldAligned.targetGramMatrix(variable=variable),majorantGramMatrix=oldAligned.targetMajorantGramMatrix(variable=variable),activeColumnMask=oldAligned.activeModeMask(variable=variable));
                testCase.verifyEqual(projection.recovery.forwardMatrix,oldAligned.forwardMatrix(variable=variable),AbsTol=2e-13)
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
            testCase.verifyEqual(projection.projectPhysical(values),expected,AbsTol=1e-14)
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
            testCase.verifyEqual(prefix.projectPhysical(values),(pairing(1,:)*values)/sampleGram(1,1),AbsTol=1e-14)
            testCase.verifyEqual(prefix.provenance,provenance)
            testCase.verifyEqual(prefix.columnLabels,"mode1+")
            testCase.verifyTrue(isnan(prefix.gramError))
            testCase.verifyEqual(projection.prefix(0).physicalDualMatrix,zeros(0,3))
        end

        function prescribedDualZeroInactiveAndRankDefectRemainExplicit(testCase)
            provenance = struct(source="test-dual");
            projection = IMProjection.fromPrescribedDual([0 0;1 1i],diag([0 -2]),diag([0 -2]),activeColumnMask=[false true],provenance=provenance);
            testCase.verifyEqual(projection.projectPhysical([1;1]),[0;-(1+1i)/2],AbsTol=1e-14)
            testCase.verifyEqual(projection.productError(zeros(2,1),zeros(2,1),0),0)
            testCase.verifyEqual(projection.productError([0;1],zeros(2,1),0),Inf)
            defect = IMProjection.fromPrescribedDual(eye(2),ones(2),eye(2),provenance=provenance);
            testCase.verifyEqual(defect.sampledGramRank,1)
            testCase.verifyTrue(isnan(defect.gramError))
            testCase.verifyEmpty(defect.recovery)
            testCase.verifyEmpty(defect.physicalDualMatrix)
            testCase.verifyEqual(defect.samplePairingMatrix,eye(2))
            testCase.verifyError(@() defect.projectPhysical([1;2]),"IMProjection:UnavailablePhysicalDual")
            offDiagonal = IMProjection.fromPrescribedDual(eye(2),[0 1;1 0],[0 1;1 0],majorantGramMatrix=eye(2),provenance=provenance);
            testCase.verifyEqual(offDiagonal.projectPhysical([1;2]),[2;1],AbsTol=0)
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
            testCase.verifyError(@() IMProjection([1e-15 1;0 2],eye(2),diag([0 5]),activeColumnMask=[false true]),"IMProjection:InvalidInactiveColumn")
            projection = IMProjection(eye(2),eye(2),eye(2));
            testCase.verifyError(@() projection.recovery.transformForward(ones(3,1)),"IMDiscreteTransform:InvalidSampleCount")
            testCase.verifyError(@() projection.productError(ones(2,1),ones(3,1),1),"IMProjection:InvalidShape")
            testCase.verifyError(@() projection.leakage(ones(2,1),[1;2]),"IMProjection:InvalidShape")
            testCase.verifyError(@() projection.prefix(3),"IMProjection:InvalidColumnCount")
        end

        function recoveryDoesNotDependOnRoundedGramRank(testCase)
            Phi = [1 1;0 1e-8];
            projection = IMProjection(Phi,eye(2),eye(2));
            testCase.verifyEqual(projection.sampledGramRank,1)
            testCase.verifyEqual(projection.gramError,norm(Phi.'*Phi-eye(2),2),AbsTol=1e-14)
            testCase.verifyTrue(projection.recovery.success)
            testCase.verifyEqual(projection.recovery.diagnostics.fit.rank,2)
            testCase.verifyEqual(projection.recovery.transformForward(Phi),eye(2),AbsTol=2e-8)
            testCase.verifyEqual(projection.physicalDualMatrix,Phi.',AbsTol=0)
            testCase.verifyGreaterThan(norm(projection.projectPhysical(Phi)-eye(2),2),0.9)
        end

        function singularSignedRecoveryRetainsPhysicalProjection(testCase)
            projection = IMProjection(ones(2,1),diag([1 -1]),1);
            testCase.verifyFalse(projection.recovery.success)
            testCase.verifyEqual(projection.recovery.status,"overlapRankDeficient")
            testCase.verifyEqual(projection.sampledGramRank,0)
            testCase.verifyEqual(projection.gramError,1,AbsTol=0)
            testCase.verifyEqual(projection.projectPhysical([3;1]),2,AbsTol=0)
            testCase.verifyError(@() projection.recovery.transformForward([3;1]),"IMDiscreteTransform:UnavailableRecovery")
            testCase.verifyError(@() IMProjection(ones(2,1),diag([1 -1]),1,failureMode="error"),"IMRecovery:OverlapRankDeficient")
        end

        function physicalRecordDoesNotPrepareOrApplyAnInverse(testCase)
            projection = IMProjection([1;2],eye(2),3,prepareRecovery=false);
            testCase.verifyEmpty(projection.recovery)
            testCase.verifyFalse(isprop(projection,"forwardMatrix"))
            testCase.verifyFalse(projection.supportsRoundTripAssessment)
            testCase.verifyTrue(isnan(projection.roundTripError))
            testCase.verifyEqual(projection.projectPhysical([3;0]),1,AbsTol=1e-14)
            testCase.verifyFalse(ismethod(projection,"project"))
            testCase.verifyEmpty(projection.prefix(1).recovery)
        end
    end
end
