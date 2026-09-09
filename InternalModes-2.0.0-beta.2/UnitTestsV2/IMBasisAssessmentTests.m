classdef IMBasisAssessmentTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addPackage(testCase)
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fileparts(fileparts(mfilename("fullpath")))));
        end
    end
    methods (Test)
        function scalarControlsUseTheCommonContract(testCase)
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1,surfaceBoundary=IMBoundaryCondition.neumann(),bottomBoundary=IMBoundaryCondition.neumann());
            basis = IMSolverSpectral(nEVP=64).solveEVP(evp,nModes=18);
            collection = IMBasisCollection({basis});
            z = linspace(-1,0,13).';
            w = [0.5;ones(11,1);0.5]/12;
            [~,legacy] = basis.discreteTransform(z=z,weights=w,leakageTolerance=1e-8,nCheckModes=18,quadraticAliasingTolerance=1e-8);
            count = 13;
            p = collection.projectionOnGrid(z,w,columns=1:count);
            checkGram = basis.gramMatrix();
            leakage = struct(sampledValues=basis.u(z),normSquared=diag(checkGram),columnLabels=string(basis.modeNumber));
            zr = basis.solver.innerProductGrid(basis.zDomain);
            cardinal = eye(numel(zr));
            wr = zeros(size(zr));
            for iPoint = 1:numel(zr)
                wr(iPoint) = basis.solver.integrateInnerProduct(zr,cardinal(:,iPoint),basis.zDomain);
            end
            [zr,order] = sort(zr);
            wr = wr(order);
            referenceProjection = collection.projectionOnGrid(zr,wr,columns=1:count);
            U = referenceProjection.sampledBasis;
            Mr = referenceProjection.metricMatrix;
            pairs = nchoosek(1:count,2).';
            pairs = [pairs [1:count;1:count]];
            samples = p.sampledBasis(:,pairs(1,:)).*p.sampledBasis(:,pairs(2,:));
            referenceProducts = U(:,pairs(1,:)).*U(:,pairs(2,:));
            products = struct(sampledValues=samples,referencePairings=U.'*Mr*referenceProducts,normSquared=sum(referenceProducts.*(Mr*referenceProducts),1),inputColumns=pairs,inputLabels=string(basis.modeNumber(pairs)),referenceStatus="qualified",referenceProvenance="same explicit reference rule as the established trigonometric regression");
            result = collection.assess(z,w,columns=1:count,prefixColumnCounts=1:count,leakage=leakage,products=products);
            gram = result.measurements(result.measurements.quantity == "gram",:);
            leak = result.measurements(result.measurements.quantity == "leakage",:);
            quadratic = result.measurements(result.measurements.quantity == "quadraticAliasing",:);
            testCase.verifyEqual(gram.value,legacy.prefixDiagnostics.gramError,AbsTol=1e-11);
            testCase.verifyEqual(leak.value,legacy.prefixDiagnostics.leakageError,AbsTol=1e-11);
            testCase.verifyEqual(quadratic.value,legacy.prefixDiagnostics.quadraticAliasingError,AbsTol=1e-10);
            decision = checkBasisAssessment(result,gramTolerance=1e-2,quadraticAliasingTolerance=1e-8);
            testCase.verifyEqual(decision.status,"rejected");
            testCase.verifyEqual(decision.requestedColumnCount,count);
            testCase.verifyEqual(decision.largestExaminedAcceptedPrefix,legacy.quadraticAliasingPolicy.maximumAcceptedModeCount);
            testCase.verifySize(result.projection.sampledBasis,[13 13]);
            testCase.verifyFalse(result.coverage.superpositionGuarantee);
        end

        function inconclusiveReferenceCannotAcceptOrChangeCounts(testCase)
            p = IMProjection(eye(2),eye(2),eye(2),columnLabels=["-1" "0"]);
            products = struct(sampledValues=eye(2),referencePairings=eye(2),normSquared=[1 1],inputColumns=[1 2;1 2],inputLabels=["-1" "0";"-1" "0"],referenceStatus="inconclusive");
            result = p.assess(products=products);
            decision = checkBasisAssessment(result,gramTolerance=0,quadraticAliasingTolerance=0);
            testCase.verifyEqual(decision.status,"inconclusive");
            testCase.verifyEqual(decision.requestedColumnCount,2);
            testCase.verifyEqual(result.measurements.value(3),0);
            testCase.verifyEqual(result.projection.columnLabels,["-1" "0"]);
            testCase.verifyError(@() checkBasisAssessment(result),"checkBasisAssessment:InvalidPolicy");
            gramOnly = p.assess();
            testCase.verifyEqual(checkBasisAssessment(gramOnly,quadraticAliasingTolerance=1).status,"inconclusive");
        end

        function denseSignedTargetsResolveEachPrefixIndependently(testCase)
            target = [1 0.4;0.4 -1];
            p = IMProjection(eye(2),target,target,majorantGramMatrix=eye(2),columnLabels=["surfaceMode" "interiorMode"]);
            values = [1;2];
            products = struct(sampledValues=values,referencePairings=target*values,normSquared=5,inputColumns=[0;0],inputLabels=["mixedA";"mixedB"],referenceStatus="qualified");
            result = p.assess(prefixColumnCounts=[1 2],products=products);
            rows = result.measurements(result.measurements.quantity == "quadraticAliasing",:);
            testCase.verifyEqual(rows.value,[0;0],AbsTol=1e-14);
            testCase.verifyEqual(rows.limitingInputI,["mixedA";"mixedA"]);
        end

        function prescribedPhysicalDualDoesNotClaimScalarGramQuality(testCase)
            pairing = [1 1i 0;0 1 -1i];
            gram = diag([2 5]);
            p = IMProjection.fromPrescribedDual(pairing,gram,gram,majorantGramMatrix=gram,columnLabels=["wavePlus" "waveMinus"],provenance=struct(kind="supplied full-state energy dual"));
            samples = [1;2;3];
            products = struct(sampledValues=samples,referencePairings=pairing*samples,normSquared=14,inputColumns=[0;0],inputLabels=["sourceA";"sourceB"],referenceStatus="qualified");
            result = p.assess(products=products,identity=struct(modeNumber=[1 1],frequencySigns=[1 -1]));
            testCase.verifyEqual(result.measurements.status(1),"unsupported");
            testCase.verifyTrue(isnan(result.measurements.value(1)));
            testCase.verifyEqual(checkBasisAssessment(result,quadraticAliasingTolerance=1e-12).status,"accepted");
            testCase.verifyEqual(checkBasisAssessment(result,gramTolerance=1e-2).status,"inconclusive");
            testCase.verifyEmpty(p.sampledBasis);
            testCase.verifyEqual(result.projection.columnCount,2);
        end

        function selectedColumnOrdinalsAreExplicit(testCase)
            evp = IMEigenvalueProblem(zDomain=[-1 0],p=1,q=0,r=1,surfaceBoundary=IMBoundaryCondition.dirichlet(),bottomBoundary=IMBoundaryCondition.dirichlet());
            basis = IMSolverSpectral(nEVP=32).solveEVP(evp,nModes=3);
            collection = IMBasisCollection({basis});
            z = linspace(-1,0,17).';
            w = [0.5;ones(15,1);0.5]/16;
            products = struct(sampledValues=zeros(17,1),referencePairings=zeros(2,1),normSquared=0,inputColumns=[2;2],inputLabels=["1";"1"],referenceStatus="qualified");
            result = collection.assess(z,w,columns=[3 1],prefixColumnCounts=[1 2],products=products);
            rows = result.measurements(result.measurements.quantity == "quadraticAliasing",:);
            testCase.verifyEqual(result.projection.columnLabels,["3" "1"]);
            testCase.verifyEqual(rows.examinedCount,[0;1]);
            products.inputColumns = [3;3];
            testCase.verifyError(@() collection.assess(z,w,columns=[3 1],products=products),"IMProjection:InvalidProducts");
        end

        function inactiveReferenceCannotMasqueradeAsZeroProduct(testCase)
            p = IMProjection(zeros(2,1),eye(2),0,activeColumnMask=false);
            products = struct(sampledValues=zeros(2,1),referencePairings=1,normSquared=0,inputColumns=[1;1],inputLabels=["zero";"zero"],referenceStatus="qualified");
            testCase.verifyError(@() p.assess(products=products),"IMProjection:InactiveReferencePairing");
        end

        function endpointCoordinatesAreNeverAutomaticallyPrefixed(testCase)
            p = IMProjection(eye(2),eye(2),eye(2),columnLabels=["surface" "bottom"]);
            result = p.assess(columnKind="endpoint");
            testCase.verifyEqual(result.prefixColumnCounts,2);
            testCase.verifyError(@() p.assess(columnKind="endpoint",prefixColumnCounts=[1 2]),"IMProjection:EndpointPrefix");
        end

        function structuredResultsValidateAtTheDecisionBoundary(testCase)
            testCase.verifyError(@() checkBasisAssessment(struct(),gramTolerance=1),"checkBasisAssessment:InvalidAssessment");
            result = IMProjection(eye(2),eye(2),eye(2)).assess();
            result.measurements(1,:) = [];
            testCase.verifyError(@() checkBasisAssessment(result,gramTolerance=1),"checkBasisAssessment:InvalidAssessment");
            result = IMProjection(eye(2),eye(2),eye(2)).assess();
            result.prefixColumnCounts = [2 1];
            testCase.verifyError(@() checkBasisAssessment(result,gramTolerance=1),"checkBasisAssessment:InvalidAssessment");
        end

        function cumulativePolicyPreservesEarlierFailure(testCase)
            p = IMProjection([1 0;0 1],eye(2),eye(2));
            leakage = struct(sampledValues=[2;0],normSquared=1,columnLabels="2");
            result = p.assess(prefixColumnCounts=[1 2],leakage=leakage);
            decision = checkBasisAssessment(result,leakageTolerance=0.5);
            testCase.verifyEqual(decision.status,"rejected");
            testCase.verifyEqual(decision.largestExaminedAcceptedPrefix,0);
        end
    end
end
