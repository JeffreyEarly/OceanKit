classdef IMInternalModesDiscreteTransformPolicyTests < matlab.unittest.TestCase

    properties (Access = private)
        originalPath
        cosineBasis
        z
        weights
    end

    methods (TestClassSetup)
        function configureCosineSineFamily(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.originalPath = path;
            addpath(repoRoot);
            addpath(fullfile(repoRoot,"UnitTestsV2","TestSupport"));
            D = 4000;
            N2 = @(z) (5.2e-3)^2*ones(size(z));
            solver = IMSolverSpectral(nEVP=160);
            testCase.cosineBasis = solver.solveEVP(IMInternalModes.hydrostaticFModes(N2=N2,zDomain=[-D 0]),nModes=30);
            testCase.z = linspace(-D,0,13).';
            testCase.weights = [0.5;ones(11,1);0.5]*(D/12);
        end
    end

    methods (TestClassTeardown)
        function restorePath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function coupledDCTIDSTIHasExactTwoThirdsCutoff(testCase)
            M = length(testCase.z)-1;
            Kmax = floor((2*M-1)/3);
            [transform,assessment] = testCase.cosineBasis.discreteTransform(z=testCase.z,weights=testCase.weights,variables=["F","G"],quadraticAliasingTolerance=1e-8);
            firstRejected = Kmax+2;

            testCase.verifyEqual(Kmax,7)
            testCase.verifyEqual(assessment.retainedModeCount,Kmax+1)
            testCase.verifyEqual(transform.modeNumber,0:Kmax,AbsTol=0)
            testCase.verifyLessThan(assessment.prefixDiagnostics.quadraticAliasingError(Kmax+1),1e-8)
            testCase.verifyGreaterThan(assessment.prefixDiagnostics.quadraticAliasingError(firstRejected),0.9)
            testCase.verifyEqual(assessment.prefixDiagnostics.quadraticLimitingChannel(firstRejected),"FG->G")
            testCase.verifyEqual(assessment.prefixDiagnostics.quadraticLimitingModeNumberI(firstRejected),Kmax+1,AbsTol=0)
            testCase.verifyEqual(assessment.prefixDiagnostics.quadraticLimitingModeNumberJ(firstRejected),Kmax+1,AbsTol=0)
        end

        function eachCoupledProductChannelMatchesIndependentTrigExpansion(testCase)
            nModes = 6;
            transform = testCase.cosineBasis.discreteTransform(z=testCase.z,weights=testCase.weights,nModes=nModes,variables=["F","G"],gramTolerance=1);
            F = transform.inverseMatrix(variable="F");
            G = transform.inverseMatrix(variable="G");
            for pair = [2 3;2 4;3 4].'
                iMode = pair(1);
                jMode = pair(2);
                ff = transform.transformForward(F(:,iMode).*F(:,jMode),variable="F");
                fg = transform.transformForward(F(:,iMode).*G(:,jMode),variable="G");
                gg = transform.transformForward(G(:,iMode).*G(:,jMode),variable="F");
                testCase.verifyGreaterThan(norm(ff),0)
                testCase.verifyGreaterThan(norm(fg),0)
                testCase.verifyGreaterThan(norm(gg),0)
                testCase.verifyEqual(ff,transform.forwardMatrix(variable="F")*(F(:,iMode).*F(:,jMode)),AbsTol=1e-13)
                testCase.verifyEqual(fg,transform.forwardMatrix(variable="G")*(F(:,iMode).*G(:,jMode)),AbsTol=1e-13)
                testCase.verifyEqual(gg,transform.forwardMatrix(variable="F")*(G(:,iMode).*G(:,jMode)),AbsTol=1e-13)
            end
        end

        function gramPolicyUsesWorstRequestedChannel(testCase)
            [~,assessment] = testCase.cosineBasis.discreteTransform(z=testCase.z,weights=testCase.weights,variables=["F","G"]);
            for nPrefix = 1:assessment.candidateModeCount
                f = assessment.variablePrefixDiagnostics(variable="F");
                g = assessment.variablePrefixDiagnostics(variable="G");
                expected = max(f.gramError(nPrefix),g.gramError(nPrefix));
                testCase.verifyEqual(assessment.prefixDiagnostics.gramError(nPrefix),expected,AbsTol=0)
            end
        end

        function zeroNormClassificationUsesPerModeMajorantScale(testCase)
            basis = IMInternalModesBasisTestAccess(testCase.cosineBasis);
            sampled = [zeros(4,1) ones(4,1) 1e30*ones(4,1)];
            targetGram = diag([0 1 1e60]);
            targetMajorantGram = diag([0 1 1e60]);

            testCase.verifyEqual(basis.classifyTransformColumns(sampled,targetGram,targetMajorantGram), ...
                [false true true])
        end

        function signedAPVUsesMajorantForQuadraticErrorsButStillRejectsLeakage(testCase)
            N2 = @(z) 1e-4*ones(size(z));
            evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=[-4000 0],g0=-1,gd=-1,surfaceBoundary="rigidLid");
            solver = IMSolverSpectral(nEVP=160);
            basis = solver.solveEVP(evp,nModes=6);
            sampleZ = linspace(-4000,0,25).';
            sampleWeights = [0.5;ones(23,1);0.5]*(4000/24);
            transform = basis.discreteTransform(z=sampleZ,weights=sampleWeights,nModes=4,variables="G",gramTolerance=10);
            [quadraticTransform,assessment] = basis.discreteTransform(z=sampleZ,weights=sampleWeights,nModes=4, ...
                variables=["F","G"],gramTolerance=10,quadraticAliasingTolerance=1);
            testCase.verifyFalse(transform.targetGramIsPositiveDefinite(variable="G"))
            testCase.verifyEqual(quadraticTransform.modeNumber(1),-1,AbsTol=0)
            testCase.verifyGreaterThan(min(eig(quadraticTransform.targetMajorantGramMatrix(variable="G"))),0)
            testCase.verifyTrue(all(isfinite(assessment.quadraticAliasingPolicy.error)))
            testCase.verifyEqual(assessment.prefixDiagnostics.quadraticLimitingModeNumberI(1),-1,AbsTol=0)
            testCase.verifyEqual(assessment.prefixDiagnostics.quadraticLimitingModeNumberJ(1),-1,AbsTol=0)
            testCase.verifyEqual(assessment.quadraticAliasingPolicy.errorNorm,"inducedHilbertMajorant")
            testCase.verifyEqual(assessment.quadraticAliasingPolicy.projectionPairing,"signedPontryagin")
            testCase.verifyError(@() basis.discreteTransform(z=sampleZ,weights=sampleWeights,nModes=4,variables="G",gramTolerance=10,leakageTolerance=1), ...
                "IMBasisSet:UnavailableDiscreteTransformPolicy")
        end

        function exponentialCoupledPolicyUsesRefinedContinuousReference(testCase)
            D = 4000;
            N0 = 5.2e-3;
            b = 1300;
            N2 = @(z) N0*N0*exp(2*z/b);
            solver = IMSolverSpectral(nEVP=192,coordinateKind="wkb");
            basis = solver.solveEVP(IMInternalModes.hydrostaticFModes(N2=N2,zDomain=[-D 0]),nModes=12);
            sampleZ = linspace(-D,0,25).';
            sampleWeights = [0.5;ones(23,1);0.5]*(D/24);
            [~,assessment] = basis.discreteTransform(z=sampleZ,weights=sampleWeights,nModes=6,variables=["F","G"],gramTolerance=1,quadraticAliasingTolerance=1);
            testCase.verifyTrue(all(isfinite(assessment.quadraticAliasingPolicy.error)))
            testCase.verifyTrue(all(ismember(assessment.quadraticAliasingPolicy.limitingChannel,["FF->F","FG->G","GG->F"])))
        end

        function batchedQuadraticAssessmentMatchesScalarIntegrals(testCase)
            D = 4000;
            N0 = 5.2e-3;
            b = 1300;
            N2 = @(z) N0*N0*exp(2*z/b);
            solver = IMSolverSpectral(nEVP=128,coordinateKind="wkb");
            basis = solver.solveEVP(IMInternalModes.hydrostaticFModes(N2=N2,zDomain=[-D 0]),nModes=10);
            sampleZ = linspace(-D,0,19).';
            sampleWeights = [0.5;ones(17,1);0.5]*(D/18);
            nModes = 5;
            [~,assessment] = basis.discreteTransform(z=sampleZ,weights=sampleWeights,nModes=nModes,variables=["F","G"], ...
                gramTolerance=10,quadraticAliasingTolerance=10);
            [expectedError,expectedChannel,expectedLabelI,expectedLabelJ] = scalarQuadraticReference(basis,sampleZ,sampleWeights,nModes);

            testCase.verifyEqual(assessment.quadraticAliasingPolicy.error,expectedError,AbsTol=2e-11)
            testCase.verifyEqual(assessment.quadraticAliasingPolicy.limitingChannel,expectedChannel)
            testCase.verifyEqual(assessment.quadraticAliasingPolicy.limitingModeNumberI,expectedLabelI,AbsTol=0)
            testCase.verifyEqual(assessment.quadraticAliasingPolicy.limitingModeNumberJ,expectedLabelJ,AbsTol=0)
        end

        function exactPointCountAndStrictFamilyBandKeepIssueSevenContracts(testCase)
            [transform,assessment] = testCase.cosineBasis.discreteTransform(nPoints=6,variables=["F","G"],gramTolerance=1);
            testCase.verifyEqual(length(transform.z),6)
            testCase.verifyEqual(assessment.requestedPointCount,6)
            testCase.verifyEqual(assessment.actualPointCount,6)
            testCase.verifyEqual(assessment.candidateModeCount,4)
            testCase.verifyEqual(assessment.weightFit.transform.modeNumber,0:3,AbsTol=0)
            testCase.verifyError(@() testCase.cosineBasis.discreteTransform(nPoints=2,variables=["F","G"]), ...
                "IMBasisSet:UnattainableDiscretePointCount")
            testCase.verifyError(@() testCase.cosineBasis.discreteTransform(z=testCase.z,weights=testCase.weights,nModes=13,variables=["F","G"],gramTolerance=10), ...
                "IMBasisSet:StrictDiscreteModeCountRejected")
        end
    end
end

function [errors,limitingChannel,limitingLabelI,limitingLabelJ] = scalarQuadraticReference(basis,z,weights,nModes)
variables = ["F","G"];
transform = basis.discreteTransform(z=z,weights=weights,nModes=nModes,variables=variables,gramTolerance=10);
integrationGrid = basis.solver.innerProductGrid(basis.zDomain);
continuousValues = struct(F=basis.F(integrationGrid),G=basis.G(integrationGrid));
sampleValues = struct(F=transform.inverseMatrix(variable="F"),G=transform.inverseMatrix(variable="G"));
channels = ["F" "F" "F";"G" "G" "F";"F" "G" "G"];
errors = zeros(nModes,1);
limitingChannel = strings(nModes,1);
limitingLabelI = nan(nModes,1);
limitingLabelJ = nan(nModes,1);
for iChannel = 1:size(channels,1)
    sourceA = channels(iChannel,1);
    sourceB = channels(iChannel,2);
    target = channels(iChannel,3);
    channelName = sourceA+sourceB+"->"+target;
    spec = basis.evp.innerProduct(target);
    coefficientContext = basis.evp.contextForSolver(basis.solver);
    interiorWeight = IMEigenvalueProblem.evaluateCoefficient(spec.interiorWeight,integrationGrid,coefficientContext);
    if isscalar(interiorWeight)
        interiorWeight = interiorWeight*ones(size(integrationGrid));
    end
    pairIndices = referenceProductPairs(nModes,sourceA == sourceB);
    nPairs = size(pairIndices,1);
    sampleProducts = zeros(length(z),nPairs);
    productNorms = zeros(nPairs,1);
    pairings = zeros(nModes,nPairs);
    usable = false(nPairs,1);
    for iPair = 1:nPairs
        iMode = pairIndices(iPair,1);
        jMode = pairIndices(iPair,2);
        sampleProduct = sampleValues.(char(sourceA))(:,iMode).*sampleValues.(char(sourceB))(:,jMode);
        product = continuousValues.(char(sourceA))(:,iMode).*continuousValues.(char(sourceB))(:,jMode);
        sampleProducts(:,iPair) = sampleProduct;
        productNorm = basis.solver.integrateInnerProduct(integrationGrid,interiorWeight.*product.*product,basis.zDomain);
        productScale = max(abs(product))^2*diff(basis.zDomain);
        zeroTolerance = 1e3*eps(max(1,productScale));
        usable(iPair) = ~(abs(productNorm) <= zeroTolerance && norm(sampleProduct) <= sqrt(zeroTolerance));
        productNorms(iPair) = productNorm;
        for kMode = 1:nModes
            pairings(kMode,iPair) = basis.solver.integrateInnerProduct(integrationGrid, ...
                interiorWeight.*continuousValues.(char(target))(:,kMode).*product,basis.zDomain);
        end
    end
    metric = transform.metricMatrix(variable=target);
    inverseFull = transform.inverseMatrix(variable=target);
    targetGramFull = transform.targetGramMatrix(variable=target);
    targetMajorantGramFull = transform.targetMajorantGramMatrix(variable=target);
    activeFull = transform.activeModeMask(variable=target);
    for nPrefix = 1:nModes
        selected = usable & max(pairIndices,[],2) <= nPrefix;
        if ~any(selected)
            continue
        end
        active = activeFull(1:nPrefix);
        inverse = inverseFull(:,1:nPrefix);
        sampledGram = inverse.'*metric*inverse;
        targetGram = targetGramFull(1:nPrefix,1:nPrefix);
        targetMajorantGram = targetMajorantGramFull(1:nPrefix,1:nPrefix);
        sampledCoefficients = zeros(nPrefix,nnz(selected));
        continuousCoefficients = zeros(nPrefix,nnz(selected));
        sampledCoefficients(active,:) = sampledGram(active,active)\(inverse(:,active).'*metric*sampleProducts(:,selected));
        continuousCoefficients(active,:) = targetGram(active,active)\pairings(active,selected);
        difference = sampledCoefficients-continuousCoefficients;
        numerator = real(sum(conj(difference(active,:)).*(targetMajorantGram(active,active)*difference(active,:)),1));
        values = sqrt(max(0,numerator)./productNorms(selected).');
        [value,iLimiting] = max(values);
        if value >= errors(nPrefix)
            selectedPairs = find(selected);
            pair = pairIndices(selectedPairs(iLimiting),:);
            errors(nPrefix) = value;
            limitingChannel(nPrefix) = channelName;
            limitingLabelI(nPrefix) = transform.modeNumber(pair(1));
            limitingLabelJ(nPrefix) = transform.modeNumber(pair(2));
        end
    end
end
end

function pairIndices = referenceProductPairs(nModes,isSymmetric)
pairIndices = zeros(0,2);
for iMode = 1:nModes
    if isSymmetric
        jModes = iMode:nModes;
    else
        jModes = 1:nModes;
    end
    pairIndices = [pairIndices;repmat(iMode,length(jModes),1) jModes(:)]; %#ok<AGROW>
end
end
