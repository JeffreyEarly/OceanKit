classdef IMNormalizationCovarianceTests < matlab.unittest.TestCase
    % Check physical invariance under immutable coefficient-coordinate changes.

    methods (TestClassSetup)
        function addProvider(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,"UnitTestsV2","TestSupport")));
        end
    end

    methods (Test)
        function complexObservationsCovaryForEveryRecoveryPolicy(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            basis = createIMTestBasis("wave",background,k=1e-3,solver=IMSolverAnalytical(),nModes=5);
            normalized = basis.withNormalization("maxAbsG");
            T = normalizationMap(basis,normalized);
            z = linspace(-970,-30,19).';
            A = [basis.F(z);(1+2i)*basis.G(z)];
            normalizedA = [normalized.F(z);(1+2i)*normalized.G(z)];
            testCase.verifyEqual(normalizedA,A*T,AbsTol=3e-14);
            nRows = size(A,1); nColumns = size(A,2);
            D = diag(linspace(.7,1.8,nRows));
            F = diag(linspace(1.6,.6,nRows));
            C = diag(linspace(.8,1.3,nColumns));
            complement = null(A');
            B = A\eye(nRows)+.1*(1:nColumns).'*complement(:,1)';
            retained = [3 1 2];
            Tr = T(retained,retained);
            Cr = C(retained,retained);
            retainedLabels = basis.columnLabels(retained);
            coefficients = [1+2i;-.3+.5i;2-1i];
            source = A(:,retained)*coefficients+.4*A(:,4)-.2i*A(:,5)+(1-.7i)*complement(:,1);
            testCase.verifyGreaterThan(norm(B*complement(:,1)),.1);
            maps = cell(1,3);
            policies = ["sampledLeastSquares","correctedPhysicalDual","extendedSpan"];
            for j = 1:numel(policies)
                options = struct(policy=policies(j),columnLabels=retainedLabels,coefficientScale=Cr,observationFactor=D,fitFactor=F,physicalDual=B(retained,:));
                newOptions = options;
                newOptions.coefficientScale = Tr\Cr;
                newOptions.physicalDual = Tr\B(retained,:);
                if policies(j) == "extendedSpan"
                    options.extendedSynthesis = A;
                    options.extendedColumnLabels = basis.columnLabels;
                    options.extendedActiveColumnMask = true(1,nColumns);
                    options.extendedCoefficientScale = C;
                    options.retainedColumnIndices = retained;
                    newOptions.extendedSynthesis = normalizedA;
                    newOptions.extendedColumnLabels = normalized.columnLabels;
                    newOptions.extendedActiveColumnMask = true(1,nColumns);
                    newOptions.extendedCoefficientScale = T\C;
                    newOptions.retainedColumnIndices = retained;
                end
                pairs = namedargs2cell(options); newPairs = namedargs2cell(newOptions);
                original = IMDiscreteTransform.fromObservationMatrices(A(:,retained),pairs{:});
                changed = IMDiscreteTransform.fromObservationMatrices(normalizedA(:,retained),newPairs{:});
                testCase.verifyTrue(original.success && changed.success);
                testCase.verifyEqual(changed.forwardMatrix,Tr\original.forwardMatrix,RelTol=3e-12,AbsTol=2e-13);
                testCase.verifyEqual(changed.physicalDualMatrix,Tr\original.physicalDualMatrix,AbsTol=0);
                testCase.verifyEqual(changed.transformForward(source),Tr\original.transformForward(source),RelTol=3e-12,AbsTol=3e-12);
                testCase.verifyEqual(changed.transformBack(changed.transformForward(source)),original.transformBack(original.transformForward(source)),RelTol=3e-12,AbsTol=3e-12);
                testCase.verifyEqual(changed.physicalDualMatrix*source,Tr\(original.physicalDualMatrix*source),RelTol=5e-14,AbsTol=3e-14);
                testCase.verifyEqual(changed.recoverySpecification.observationFactor,original.recoverySpecification.observationFactor,AbsTol=0);
                testCase.verifyEqual(changed.recoverySpecification.fitFactor,original.recoverySpecification.fitFactor,AbsTol=0);
                testCase.verifyEqual(changed.diagnostics.observationNoiseGain,original.diagnostics.observationNoiseGain,RelTol=3e-12,AbsTol=3e-12);
                maps{j} = original.forwardMatrix;
                if policies(j) == "extendedSpan"
                    testCase.verifyEqual(changed.forwardMatrix*normalizedA(:,[4 5]),zeros(3,2),AbsTol=3e-13);
                    testCase.verifyEqual(changed.recoverySpecification.extendedCoefficientScale,T\C,AbsTol=0);
                else
                    testCase.verifyEqual(changed.recoverySpecification.coefficientScale,Tr\Cr,AbsTol=0);
                end
            end
            % This control includes observations outside the retained and
            % extended spans; covariance must preserve policy differences.
            testCase.verifyGreaterThan(norm((maps{1}-maps{2})*source),.01);
            testCase.verifyGreaterThan(norm((maps{1}-maps{3})*source),.01);
        end

        function continuousWaveEnergyChangesByCongruence(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            basis = createIMTestBasis("wave",background,k=1e-3,solver=IMSolverAnalytical(),nModes=3);
            normalized = basis.withNormalization("maxAbsG");
            T = normalizationMap(basis,normalized);
            [z,w] = legpts(80,background.zDomain); w = w(:);
            kD = basis.k*background.D;
            F = basis.F(z); G = basis.G(z);
            newF = normalized.F(z); newG = normalized.G(z);
            H = F'*(w.*F)+kD^2*(G'*(w.*G));
            newH = newF'*(w.*newF)+kD^2*(newG'*(w.*newG));
            testCase.verifyEqual(newH,T'*H*T,RelTol=5e-13,AbsTol=2e-10);
            coefficients = [1+.2i;-.3+.4i;2-.1i];
            changedCoefficients = T\coefficients;
            testCase.verifyEqual(changedCoefficients'*newH*changedCoefficients,coefficients'*H*coefficients,RelTol=5e-13,AbsTol=2e-10);
        end

        function signedPhysicalPairingAndMajorantCovaryTogether(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            energy = IMGeneralizedEnergyBoundary(g0=-.04,gd=.03);
            basis = createIMTestBasis("MDA",background,energyBoundary=energy,nModes=4,solver=IMSolverSpectral(nEVP=64,coordinateKind="z"));
            normalized = basis.withNormalization("rmsG");
            T = normalizationMap(basis,normalized);
            [z,w] = legpts(100,background.zDomain); w = w(:);
            observationZ = [z;background.zDomain([2 1]).'];
            A = basis.G(observationZ); changedA = normalized.G(observationZ);
            % The signed generalized-energy G pairing is defined by the
            % interior N2/g integral and the independent g0/g, gd/g traces.
            metric = [w.*background.N2(z)/background.g;energy.g0/background.g;energy.gd/background.g];
            W = diag(metric);
            H = A'*W*A;
            changedH = changedA'*W*changedA;
            Hplus = A'*(abs(metric).*A);
            changedHplus = changedA'*(abs(metric).*changedA);
            testCase.verifyEqual(H,diag(basis.normSignatures),AbsTol=3e-11);
            testCase.verifyTrue(any(eig(H) < 0));
            testCase.verifyEqual(changedH,T'*H*T,RelTol=8e-13,AbsTol=3e-12);
            testCase.verifyEqual(changedHplus,T'*Hplus*T,RelTol=8e-13,AbsTol=3e-12);
            testCase.verifyEqual(nnz(eig(changedH) < 0),nnz(eig(H) < 0));
            nColumns = basis.nColumns;
            coefficients = complex((1:nColumns).'/nColumns,linspace(.1,.7,nColumns).');
            changedCoefficients = T\coefficients;
            testCase.verifyEqual(changedCoefficients'*changedH*changedCoefficients,coefficients'*H*coefficients,RelTol=8e-13,AbsTol=3e-12);
            B = H\(A'*W);
            changedB = changedH\(changedA'*W);
            C = chol(Hplus,"lower")'\eye(nColumns);
            D = diag(sqrt(abs(metric)));
            original = IMDiscreteTransform.fromObservationMatrices(A,policy="correctedPhysicalDual",columnLabels=basis.columnLabels,physicalDual=B,coefficientScale=C,observationFactor=D);
            changed = IMDiscreteTransform.fromObservationMatrices(changedA,policy="correctedPhysicalDual",columnLabels=normalized.columnLabels,physicalDual=changedB,coefficientScale=T\C,observationFactor=D);
            source = complex(sin(observationZ/37),cos(observationZ/53));
            testCase.verifyEqual(changedB,T\B,RelTol=2e-12,AbsTol=3e-12);
            testCase.verifyEqual(changed.forwardMatrix,T\original.forwardMatrix,RelTol=2e-12,AbsTol=3e-12);
            testCase.verifyEqual(changed.physicalDualMatrix*source,T\(original.physicalDualMatrix*source),RelTol=2e-12,AbsTol=3e-12);
        end

        function physicalDepthFGAdaptersPreserveActiveAndInactiveCoordinates(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            full = createIMTestBasis("hydrostaticF",background,solver=IMSolverAnalytical(),nModes=4);
            grid = IMVerticalGrid.physicalDepth(background,33);
            for selected = {[1 2 3 4],[4 2 1],[]}
                basis = full.selectColumns(selected{1});
                normalized = basis.withNormalization("maxAbsF");
                original = basis.discreteTransform(grid);
                changed = normalized.discreteTransform(grid);
                verifyAlignedCovariance(testCase,original,changed,normalizationMap(basis,normalized));
                testCase.verifyEqual(original.z,grid.z,AbsTol=0);
                testCase.verifyEqual(changed.weights,grid.weights,AbsTol=0);
                testCase.verifyEqual(changed.basis.columnLabels,basis.columnLabels);
            end
        end

        function anExistingRootRuleSupportsTheSameNormalizationCovariance(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            [original,grid,~,basis] = IMInternalModesDiscreteTransform.hydrostaticRoot(background,9,nInternalModes=3,solver=IMSolverSpectral(nEVP=64,coordinateKind="z"),forwardOperator="galerkin");
            normalized = basis.withNormalization("maxAbsF");
            changed = normalized.discreteTransform(grid);
            testCase.verifyEqual(grid.kind,"modeRoots");
            testCase.verifyEqual(changed.z,original.z,AbsTol=0);
            testCase.verifyEqual(changed.weights,original.weights,AbsTol=0);
            verifyAlignedCovariance(testCase,original,changed,normalizationMap(basis,normalized));
        end

        function defaultRootRecoveryPreservesItsExtendedSpanUnderNormalization(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            [original,grid,~,basis] = IMInternalModesDiscreteTransform.hydrostaticRoot(background,9,nInternalModes=3,solver=IMSolverSpectral(nEVP=64,coordinateKind="z"));
            normalized = basis.withNormalization("maxAbsF");
            T = normalizationMap(basis,normalized);
            recoveries = struct();
            for variable = ["F","G"]
                oldRecovery = original.recovery(variable=variable);
                specification = oldRecovery.recoverySpecification;
                testCase.verifyEqual(oldRecovery.recoveryPolicy,"extendedSpan");
                testCase.verifyFalse(specification.retainedCoefficientScaleAvailable);
                retained = specification.retainedColumnIndices;
                active = specification.extendedActiveColumnMask;
                extendedT = eye(size(specification.extendedSynthesis,2));
                extendedT(retained,retained) = T;
                extendedA = specification.extendedSynthesis*extendedT;
                normalizedA = normalized.evaluate(grid.z,variable=variable);
                testCase.verifyEqual(extendedA(:,retained),normalizedA,RelTol=5e-13,AbsTol=3e-13);
                % Bind the retained entries to their actual evaluated
                % values, preserving exact public synthesis alignment.
                extendedA(:,retained) = normalizedA;
                extendedB = specification.extendedPhysicalDual;
                if ~isempty(extendedB), extendedB = extendedT\extendedB; end
                activeT = extendedT(active,active);
                metadata = oldRecovery.metadata;
                metadata.sourceNormalization = metadata.normalization;
                metadata = rmfield(metadata,"normalization");
                metadata.retainedNormalization = normalized.normalization;
                metadata.coefficientCoordinateChange = diag(extendedT).';
                metadata.retainedColumnIndices = retained;
                fullRecovery = IMDiscreteTransform.fromObservationMatrices(extendedA, ...
                    policy=specification.extendedPolicy,columnLabels=specification.extendedColumnLabels, ...
                    activeColumnMask=active,coefficientScale=activeT\specification.extendedCoefficientScale, ...
                    observationFactor=specification.observationFactor,fitFactor=specification.fitFactor, ...
                    physicalDual=extendedB,rankTolerance=specification.rankTolerance,failureMode=specification.failureMode,metadata=metadata);
                testCase.verifyTrue(fullRecovery.success);
                recoveries.(variable) = fullRecovery.selectOutputs(retained);
            end
            changed = normalized.discreteTransform(grid,recoveryF=recoveries.F,recoveryG=recoveries.G,forwardOperator="truncatedCollocation");
            testCase.verifyEqual(changed.z,original.z,AbsTol=0);
            testCase.verifyEqual(changed.weights,original.weights,AbsTol=0);
            testCase.verifyEqual(normalized.constructionRecord,basis.constructionRecord);
            physicalSource = complex(sin((1:grid.nPoints).'/3),cos((1:grid.nPoints).'/5));
            for variable = ["F","G"]
                oldRecovery = original.recovery(variable=variable);
                newRecovery = changed.recovery(variable=variable);
                oldSpecification = oldRecovery.recoverySpecification;
                newSpecification = newRecovery.recoverySpecification;
                oldPhysical = original.physicalProjection(variable=variable);
                newPhysical = changed.physicalProjection(variable=variable);
                testCase.verifyEqual(changed.forwardOperator(variable=variable),"truncatedCollocation");
                testCase.verifyEqual(newRecovery.recoveryPolicy,"extendedSpan");
                testCase.verifyFalse(isfield(newRecovery.metadata,"normalization"));
                testCase.verifyEqual(newRecovery.metadata.retainedNormalization,normalized.normalization);
                testCase.verifyEqual(newSpecification.extendedPolicy,oldSpecification.extendedPolicy);
                testCase.verifyEqual(newSpecification.extendedColumnLabels,oldSpecification.extendedColumnLabels);
                testCase.verifyEqual(newSpecification.extendedActiveColumnMask,oldSpecification.extendedActiveColumnMask);
                testCase.verifyEqual(newSpecification.retainedColumnIndices,oldSpecification.retainedColumnIndices);
                extendedT = eye(numel(oldSpecification.extendedColumnLabels));
                extendedT(oldSpecification.retainedColumnIndices,oldSpecification.retainedColumnIndices) = T;
                activeT = extendedT(oldSpecification.extendedActiveColumnMask,oldSpecification.extendedActiveColumnMask);
                testCase.verifyEqual(newSpecification.extendedCoefficientScale,activeT\oldSpecification.extendedCoefficientScale,AbsTol=0);
                testCase.verifyEqual(newSpecification.extendedSynthesis,oldSpecification.extendedSynthesis*extendedT,RelTol=5e-13,AbsTol=3e-13);
                testCase.verifyEqual(newSpecification.fitFactor,oldSpecification.fitFactor,AbsTol=0);
                testCase.verifyEqual(newSpecification.observationFactor,oldSpecification.observationFactor,AbsTol=0);
                testCase.verifyEqual(newRecovery.activeColumnMask,oldRecovery.activeColumnMask);
                testCase.verifyFalse(newSpecification.retainedCoefficientScaleAvailable);
                testCase.verifyEmpty(newSpecification.coefficientScale);
                testCase.verifyEqual(newRecovery.forwardMatrix,T\oldRecovery.forwardMatrix,RelTol=3e-12,AbsTol=3e-12);
                testCase.verifyEqual(newPhysical.physicalDualMatrix,T\oldPhysical.physicalDualMatrix,RelTol=3e-12,AbsTol=3e-12);
                testCase.verifyEqual(newPhysical.targetGramMatrix,T'*oldPhysical.targetGramMatrix*T,RelTol=3e-12,AbsTol=3e-9);
                testCase.verifyEqual(newPhysical.metricMatrix,oldPhysical.metricMatrix,AbsTol=0);
                testCase.verifyEqual(newRecovery.transformForward(physicalSource),T\oldRecovery.transformForward(physicalSource),RelTol=3e-12,AbsTol=3e-12);
                testCase.verifyEqual(newRecovery.transformBack(newRecovery.transformForward(physicalSource)),oldRecovery.transformBack(oldRecovery.transformForward(physicalSource)),RelTol=3e-12,AbsTol=3e-12);
                omitted = newSpecification.extendedActiveColumnMask & ~ismember(1:numel(newSpecification.extendedColumnLabels),newSpecification.retainedColumnIndices);
                testCase.verifyEqual(newRecovery.forwardMatrix*newSpecification.extendedSynthesis(:,omitted),zeros(basis.nColumns,nnz(omitted)),AbsTol=3e-12);
                testCase.verifyEqual(newRecovery.forwardMatrix(~newRecovery.activeColumnMask,:),zeros(nnz(~newRecovery.activeColumnMask),grid.nPoints),AbsTol=0);
                testCase.verifyEqual(newSpecification.extendedPhysicalDual,oldSpecification.extendedPhysicalDual,AbsTol=0);
            end
        end

        function neutralAPVPhysicalDualAndMajorantCovaryTogether(testCase)
            intervals = [-400 -200;-100 0];
            N2 = @(z) 1e-4*double(~((z >= -400 & z <= -200) | (z >= -100 & z <= 0)));
            background = IMBackground(N2=N2,zDomain=[-1000 0],surfaceCondition="freeSurface");
            full = IMNeutralLabelModesBasis(background,neutralIntervals=intervals,nModesPerInterval=[2 2]);
            [first,wFirst] = legpts(80,intervals(1,:));
            [second,wSecond] = legpts(80,intervals(2,:));
            z = [first;second]; weights = [wFirst(:);wSecond(:)];
            physicalSource = 1e-6*complex(sin(z/23),cos(z/17));
            for selected = {[1 2 3 4],[4 1 3],[]}
                basis = full.selectColumns(selected{1});
                normalized = basis.withNormalization("maxAbsG");
                T = normalizationMap(basis,normalized);
                original = basis.projectionOnGrid(z,weights,variable="APV");
                changed = normalized.projectionOnGrid(z,weights,variable="APV");
                testCase.verifyEqual(normalized.APV(z),basis.APV(z)*T,RelTol=5e-14,AbsTol=1e-19);
                testCase.verifyEqual(changed.physicalDualMatrix,T\original.physicalDualMatrix,RelTol=5e-13,AbsTol=3e-9);
                testCase.verifyEqual(changed.projectPhysical(physicalSource),T\original.projectPhysical(physicalSource),RelTol=5e-13,AbsTol=3e-13);
                testCase.verifyEqual(changed.majorantGramMatrix,T'*original.majorantGramMatrix*T,RelTol=5e-14,AbsTol=1e-23);
                testCase.verifyEqual(changed.provenance.continuousAPVGramMatrix,changed.majorantGramMatrix,AbsTol=0);
                testCase.verifyEqual(changed.targetGramMatrix,eye(basis.nColumns),AbsTol=0);
                coefficients = complex((1:basis.nColumns).',.2*ones(basis.nColumns,1));
                source = basis.APV(z)*coefficients;
                testCase.verifyEqual(changed.projectPhysical(source),T\original.projectPhysical(source),RelTol=5e-13,AbsTol=3e-13);
                testCase.verifyEqual(normalized.APV(z)*(T\coefficients),source,RelTol=5e-14,AbsTol=1e-19);
                testCase.verifyEmpty(changed.recovery);
            end
        end
    end
end

function T = normalizationMap(original,changed)
oldScales = original.normalizationFactors();
newScales = changed.normalizationFactors();
T = diag(oldScales./newScales);
end

function verifyAlignedCovariance(testCase,original,changed,T)
source = complex(sin((1:numel(original.z)).'/3),cos((1:numel(original.z)).'/5));
for variable = ["F","G"]
    active = original.activeModeMask(variable=variable);
    oldA = original.inverseMatrix(variable=variable);
    newA = changed.inverseMatrix(variable=variable);
    oldR = original.forwardMatrix(variable=variable);
    newR = changed.forwardMatrix(variable=variable);
    oldPhysical = original.physicalProjection(variable=variable);
    newPhysical = changed.physicalProjection(variable=variable);
    oldRecovery = original.recovery(variable=variable);
    newRecovery = changed.recovery(variable=variable);
    testCase.verifyEqual(changed.activeModeMask(variable=variable),active);
    testCase.verifyEqual(newA,oldA*T,RelTol=5e-13,AbsTol=3e-13);
    testCase.verifyEqual(newR,T\oldR,RelTol=3e-12,AbsTol=3e-12);
    testCase.verifyEqual(newPhysical.physicalDualMatrix,T\oldPhysical.physicalDualMatrix,RelTol=3e-12,AbsTol=3e-12);
    testCase.verifyEqual(newPhysical.targetGramMatrix,T'*oldPhysical.targetGramMatrix*T,RelTol=3e-12,AbsTol=3e-9);
    testCase.verifyEqual(newPhysical.majorantGramMatrix,T'*oldPhysical.majorantGramMatrix*T,RelTol=3e-12,AbsTol=3e-9);
    testCase.verifyEqual(newPhysical.metricMatrix,oldPhysical.metricMatrix,AbsTol=0);
    testCase.verifyEqual(newRecovery.recoverySpecification.fitFactor,oldRecovery.recoverySpecification.fitFactor,AbsTol=0);
    testCase.verifyEqual(newRecovery.recoverySpecification.observationFactor,oldRecovery.recoverySpecification.observationFactor,AbsTol=0);
    activeT = T(active,active);
    testCase.verifyEqual(newRecovery.recoverySpecification.coefficientScale,activeT\oldRecovery.recoverySpecification.coefficientScale,RelTol=3e-12,AbsTol=3e-12);
    testCase.verifyEqual(changed.transformForward(source,variable=variable),T\original.transformForward(source,variable=variable),RelTol=3e-12,AbsTol=3e-12);
    testCase.verifyEqual(newPhysical.projectPhysical(source),T\oldPhysical.projectPhysical(source),RelTol=3e-12,AbsTol=3e-12);
    testCase.verifyEqual(newA*(newR*source),oldA*(oldR*source),RelTol=3e-12,AbsTol=3e-12);
    testCase.verifyEqual(newR(~active,:),zeros(nnz(~active),size(newR,2)),AbsTol=0);
    testCase.verifyEqual(newA(:,~active),zeros(size(newA,1),nnz(~active)),AbsTol=0);
end
end
