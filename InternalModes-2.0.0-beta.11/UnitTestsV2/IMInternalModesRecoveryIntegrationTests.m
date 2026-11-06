classdef IMInternalModesRecoveryIntegrationTests < matlab.unittest.TestCase
    % Verify physical F/G records and shared numerical recovery remain distinct.

    methods (Test)
        function positiveFitDoesNotUseRoundedGramRank(testCase)
            A = [1 1;0 1e-8];
            transform = preparedF(A);
            recovery = transform.recovery(variable="F");
            physical = transform.physicalProjection(variable="F");
            testCase.verifyTrue(recovery.success)
            testCase.verifyEqual(recovery.recoveryPolicy,"sampledLeastSquares")
            testCase.verifyEqual(recovery.diagnostics.fit.rank,2)
            testCase.verifyEqual(physical.sampledGramRank,1)
            testCase.verifyEqual(transform.forwardMatrix(variable="F")*A,eye(2),AbsTol=1e-12)
            testCase.verifyEqual(physical.samplePairingMatrix,A',AbsTol=0)
            testCase.verifyEqual(physical.physicalDualMatrix,A',AbsTol=0)
            testCase.verifyEmpty(physical.recovery)
            testCase.verifyGreaterThan(norm(recovery.forwardMatrix*[0;1]-physical.physicalDualMatrix*[0;1]),1e8)
            testCase.verifyFalse(recovery.diagnostics.observation.available)
        end

        function failedRecoveryLeavesPhysicalDualAvailable(testCase)
            transform = preparedF([1 0;1 0]);
            recovery = transform.recovery(variable="F");
            physical = transform.physicalProjection(variable="F");
            testCase.verifyFalse(recovery.success)
            testCase.verifyEqual(recovery.diagnostics.fit.rank,1)
            testCase.verifyFalse(transform.hasForwardTransform(variable="F"))
            testCase.verifyEmpty(recovery.forwardMatrix)
            testCase.verifyEqual(physical.physicalDualMatrix,[1 1;0 0],AbsTol=0)
            testCase.verifyEqual(transform.modeProjectionFunctional([2;3],variable="F"),[5;0],AbsTol=0)
            testCase.verifyError(@() transform.forwardMatrix(variable="F"),"IMInternalModesDiscreteTransform:UnavailableRecovery")
        end

        function densePositiveAndSignedMetricsKeepTheirPhysicalPairings(testCase)
            A = [1 0;1 1;0 1];
            metrics = {[1 1 0;1 2 0;0 0 1],diag([1 -1 1])};
            policies = ["sampledLeastSquares" "correctedPhysicalDual"];
            for j = 1:numel(metrics)
                W = metrics{j};
                transform = preparedF(A,metric=W);
                recovery = transform.recovery(variable="F");
                physical = transform.physicalProjection(variable="F");
                testCase.verifyTrue(recovery.success)
                testCase.verifyEqual(recovery.recoveryPolicy,policies(j))
                testCase.verifyEqual(physical.samplePairingMatrix,A'*W,AbsTol=0)
                testCase.verifyEqual(physical.physicalDualMatrix,A'*W,AbsTol=0)
                testCase.verifyEqual(recovery.forwardMatrix,(A'*W*A)\(A'*W),AbsTol=2e-14)
                testCase.verifyFalse(recovery.diagnostics.observation.available)
            end
        end

        function ordinaryPrefixRefitsItsDeclaredSmallerSpan(testCase)
            A = [1 0;1 1;0 1];
            transform = preparedF(A);
            prefix = transform.prefixTransform(1);
            testCase.verifyEqual(prefix.forwardMatrix(variable="F"),[.5 .5 0],AbsTol=1e-14)
            fullForward = transform.forwardMatrix(variable="F");
            physical = prefix.physicalProjection(variable="F");
            testCase.verifyGreaterThan(norm(prefix.forwardMatrix(variable="F")-fullForward(1,:)),.1)
            testCase.verifyEqual(physical.physicalDualMatrix,[1 1 0],AbsTol=0)
        end

        function alignedApplicationPreservesComplexTrailingBatches(testCase)
            A = [1 0;1 1;0 1];
            transform = preparedF(A);
            coefficients = reshape(complex(1:12,12:-1:1),[2 2 3]);
            samples = transform.transformBack(coefficients,variable="F");
            testCase.verifySize(samples,[3 2 3])
            testCase.verifyEqual(transform.transformForward(samples,variable="F"),coefficients,AbsTol=2e-14)
            pairings = transform.modeProjectionFunctional(samples,variable="F");
            expected = reshape(A'*reshape(samples,3,[]),[2 2 3]);
            testCase.verifyEqual(pairings,expected,AbsTol=0)
        end

        function emptyAlignedFamilyPreservesBatchShapes(testCase)
            transform = preparedF(zeros(2,0));
            samples = zeros(2,2,3);
            testCase.verifySize(transform.transformForward(samples,variable="F"),[0 2 3])
            testCase.verifySize(transform.modeProjectionFunctional(samples,variable="F"),[0 2 3])
            testCase.verifySize(transform.transformBack(zeros(0,2,3),variable="F"),[2 2 3])
        end

        function preparedRecoveryMustKeepTheDeclaredScientificCoordinates(testCase)
            A = [1 0;1 1;0 1];
            valid = IMDiscreteTransform.fromObservationMatrices(A,columnLabels=["1" "2"],fitFactor=eye(3));
            extended = preparedF(A,recovery=valid.selectOutputs(1:2));
            testCase.verifyEqual(extended.forwardOperator(variable="F"),"extendedSpan")
            wrongLabels = IMDiscreteTransform.fromObservationMatrices(A,columnLabels=["other:1" "other:2"],fitFactor=eye(3));
            testCase.verifyError(@() preparedF(A,recovery=wrongLabels),"IMInternalModesDiscreteTransform:InvalidRecovery")
            wrongMetadata = IMDiscreteTransform.fromObservationMatrices(A,columnLabels=["1" "2"],fitFactor=eye(3),metadata=struct(variable="G"));
            testCase.verifyError(@() preparedF(A,recovery=wrongMetadata),"IMInternalModesDiscreteTransform:InvalidRecovery")
        end

        function nativeAndAnalyticalFixedGridBridgesReusePreparedRecovery(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            grid = IMVerticalGrid.physicalDepth(background.zDomain,17);
            for solver = {IMSolverSpectral(nEVP=48),IMSolverAnalytical()}
                basis = solver{1}.solveEVP(IMInternalModes.hydrostaticFModes(background),nModes=4);
                F = basis.F(grid.z); G = basis.G(grid.z);
                metadata = struct(normalization=string(basis.normalization),z=grid.z);
                recoveryF = IMDiscreteTransform.fromObservationMatrices(F,columnLabels=basis.columnLabels,fitFactor=diag(sqrt(grid.weights)),metadata=metadata);
                recoveryG = IMDiscreteTransform.fromObservationMatrices(G,columnLabels=basis.columnLabels,activeColumnMask=[false true true true],fitFactor=diag(sqrt(grid.weights)),metadata=metadata);
                transform = basis.discreteTransform(grid,recoveryF=recoveryF,recoveryG=recoveryG,forwardOperator="galerkin");
                testCase.verifyEqual(transform.forwardMatrix(variable="F"),recoveryF.forwardMatrix,AbsTol=0);
                testCase.verifyEqual(transform.forwardMatrix(variable="G"),recoveryG.forwardMatrix,AbsTol=0);
                testCase.verifyEqual(transform.inverseMatrix(variable="F"),F,AbsTol=0);
                testCase.verifyEqual(transform.inverseMatrix(variable="G"),G,AbsTol=0);
                testCase.verifyEqual(transform.transformForward(F,variable="F"),eye(4),AbsTol=3e-13);
                testCase.verifyEqual(transform.transformForward(G,variable="G"),diag([0 1 1 1]),AbsTol=3e-13);
                wrong = IMDiscreteTransform.fromObservationMatrices(1.1*F,columnLabels=basis.columnLabels,fitFactor=diag(sqrt(grid.weights)),metadata=metadata);
                testCase.verifyError(@() basis.discreteTransform(grid,recoveryF=wrong),"IMInternalModesDiscreteTransform:InvalidRecovery");
            end
        end

        function rootExtendedSpanSurvivesPrefixAndReweighting(testCase)
            background = IMBackground(N2=@(z) 1e-4*ones(size(z)),zDomain=[-1000 0],surfaceCondition="rigidLid");
            transform = IMInternalModesDiscreteTransform.hydrostaticRoot(background,9,nInternalModes=3);
            prefix = transform.prefixTransform(3);
            weights = transform.weights.*linspace(.9,1.1,9).';
            weights = 1000*weights/sum(weights);
            reweighted = transform.transformWithWeights(weights);
            for variable = ["F","G"]
                recovery = transform.recovery(variable=variable);
                testCase.verifyEqual(recovery.recoveryPolicy,"extendedSpan")
                fullForward = transform.forwardMatrix(variable=variable);
                originalPhysical = transform.physicalProjection(variable=variable);
                newPhysical = reweighted.physicalProjection(variable=variable);
                testCase.verifyEqual(prefix.forwardMatrix(variable=variable),fullForward(1:3,:),AbsTol=0)
                testCase.verifyEqual(reweighted.forwardMatrix(variable=variable),fullForward,AbsTol=0)
                testCase.verifyNotEqual(newPhysical.physicalDualMatrix,originalPhysical.physicalDualMatrix)
                testCase.verifyEmpty(newPhysical.recovery)
            end
            recoveryG = transform.recovery(variable="G");
            forwardG = transform.forwardMatrix(variable="G");
            testCase.verifyEqual(recoveryG.activeColumnMask,[false true true true])
            testCase.verifyEqual(forwardG(:,[1 end]),zeros(4,2),AbsTol=0)
            testCase.verifyEqual(forwardG(1,:),zeros(1,9),AbsTol=0)
        end
    end
end

function transform = preparedF(A,options)
arguments (Input)
    A (:,:) double
    options.recovery = []
    options.metric = []
end
nSamples = size(A,1);
nColumns = size(A,2);
z = linspace(0,1,nSamples).';
if isempty(options.metric), metric = eye(nSamples); else, metric = options.metric; end
channelF = struct(available=true,reason="",activeModeMask=true(1,nColumns),metricMatrix=metric,targetGramMatrix=eye(nColumns),targetMajorantGramMatrix=eye(nColumns));
if ~isempty(options.recovery), channelF.recovery = options.recovery; end
channelG = struct(available=false,reason="G is not part of this numerical F-channel control.",activeModeMask=false(1,nColumns),metricMatrix=[],targetGramMatrix=[],targetMajorantGramMatrix=[]);
transform = IMInternalModesDiscreteTransform(z=z,weights=ones(nSamples,1)/nSamples,modeNumber=1:nColumns,h=ones(1,nColumns),normalization="numerical control",inverseF=A,inverseG=zeros(size(A)),endpointF=A([end 1],:),endpointG=zeros(2,nColumns),channelData=struct(F=channelF,G=channelG),primaryVariable="F",zDomain=[0 1],g=9.81,modeFamily="numerical F-channel control",N2Values=ones(nSamples,1));
end
