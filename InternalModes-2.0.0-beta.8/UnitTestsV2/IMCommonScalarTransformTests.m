classdef IMCommonScalarTransformTests < matlab.unittest.TestCase
    % Check physical source pairings independently of sample recovery.

    methods (TestClassSetup)
        function configurePath(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,"UnitTestsV2","TestSupport")));
        end
    end

    methods (Test)
        function scalarAdaptersRequireAnExplicitObservationVariable(testCase)
            background = signedBackground("rigidLid");
            basis = IMSolverAnalytical().solveEVP(IMInternalModes.hydrostaticFModes(background),nModes=3);
            [z,w] = endpointRule(background.zDomain,16);
            testCase.verifyError(@() basis.projectionOnGrid(z,w),"IMScalarProjection:MissingVariable");
            testCase.verifyError(@() basis.scalarTransform(z,w),"IMScalarProjection:MissingVariable");
        end

        function nativeHydrostaticActivityAndExplicitSelectionsRemainScientific(testCase)
            background = signedBackground("rigidLid");
            [z,w] = endpointRule(background.zDomain,64);
            for solver = {IMSolverSpectral(nEVP=48),IMSolverFiniteDifference(z=linspace(-1,0,65).'),IMSolverAnalytical()}
                basis = solver{1}.solveEVP(IMInternalModes.hydrostaticFModes(background),nModes=4);
                transform = basis.scalarTransform(z,w,variable="G");
                testCase.verifyEqual(transform.activeColumnMask,[false true true true]);
                testCase.verifyEqual(transform.transformForward(basis.G(z)),diag([0 1 1 1]),AbsTol=3e-11);
                collection = IMBasisCollection({basis},basisIndex=[1 1]);
                selected = basis.selectColumns(2:4).withNormalization("rmsG");
                values = collection.evaluate(z,variable="G",columns=2:4,normalization="rmsG");
                testCase.verifyEqual(values,repmat(selected.G(z),1,1,2),AbsTol=3e-12);
                testCase.verifySize(collection.evaluate(z,variable="G",columns=[]),[numel(z) 0 2]);
                empty = basis.projectionOnGrid(z,w,variable="G",columns=[],prepareRecovery=false);
                testCase.verifySize(empty.physicalDualMatrix,[0 numel(z)]);
            end
        end

        function apvPairingsMatchIndependentFormsAndNonmodalSources(testCase)
            background = signedBackground();
            energy = IMGeneralizedEnergyBoundary(g0=-.4,gd=.3);
            solvers = {IMSolverAnalytical(),IMSolverSpectral(nEVP=64,coordinateKind="z")};
            [z,w] = endpointRule(background.zDomain,96);
            for j = 1:numel(solvers)
                basis = solvers{j}.solveEVP(IMInternalModes.geostrophicAPVModes(background,g0=energy.g0,gd=energy.gd),nModes=4);
                for variable = ["F","G"]
                    physical = basis.projectionOnGrid(z,w,variable=variable,prepareRecovery=false);
                    [H,Hplus] = independentForms(basis,variable,energy,192);
                    [Href,~] = independentForms(basis,variable,energy,256);
                    testCase.verifyEqual(H,Href,RelTol=3e-12,AbsTol=3e-12);
                    tolerance = 3e-9;
                    if j == 1, tolerance = 3e-5; end
                    C = chol(Hplus)\eye(basis.nColumns);
                    targetError = norm(C'*(physical.targetGramMatrix-H)*C,2);
                    majorantError = norm(C'*(physical.majorantGramMatrix-Hplus)*C,2);
                    fprintf("Common APV %s/%s physical operator errors: H %.3g, Hplus %.3g.\n",string(class(solvers{j})),variable,targetError,majorantError);
                    testCase.verifyLessThanOrEqual(targetError,tolerance);
                    testCase.verifyLessThanOrEqual(majorantError,tolerance);
                    A = basis.evaluate(z,variable=variable);
                    metric = independentMetric(background,energy,"geostrophicAPV",variable,z,w);
                    P = A'*diag(metric);
                    testCase.verifyEqual(physical.samplePairingMatrix,P,AbsTol=2e-13);
                    testCase.verifyEqual(physical.columnLabels,basis.columnLabels);
                    testCase.verifyTrue(all(physical.activeColumnMask));
                    input = independentSource(z,variable);
                    expected = physical.targetGramMatrix\(P*input);
                    testCase.verifyEqual(physical.projectPhysical(input),expected,RelTol=2e-12,AbsTol=2e-12);
                    [zr,wr] = endpointRule(background.zDomain,256);
                    refinedMetric = independentMetric(background,energy,"geostrophicAPV",variable,zr,wr);
                    continuous = Href\(basis.evaluate(zr,variable=variable)'*(refinedMetric.*independentSource(zr,variable)));
                    testCase.verifyEqual(physical.projectPhysical(input),continuous,RelTol=2*tolerance,AbsTol=2*tolerance);
                    testCase.verifyEmpty(physical.recovery);
                    if variable == "G"
                        testCase.verifyTrue(any(eig(H) < -.01));
                        testCase.verifyGreaterThan(min(eig(Hplus)),0);
                    end
                end
            end
        end

        function literalEndpointCoefficientsSurvivePhysicalUnitScaling(testCase)
            [z,w] = endpointRule([-1 0],16);
            for g = [2 1e-16]
                background = IMBackground.constantStratification(N0=sqrt(g),zDomain=[-1 0],g=g,rotationRate=.5,latitude=90,surfaceCondition="freeSurface");
                for g0 = [.3 1e-15]
                    boundary = IMGeneralizedEnergyBoundary(g0=g0,gd=Inf);
                    basis = IMSolverAnalytical().solveEVP(IMInternalModes.geostrophicAPVModes(background,g0=boundary.g0,gd=boundary.gd),nModes=0);
                    physical = basis.projectionOnGrid(z,w,variable="G",prepareRecovery=false);
                    expected = g0/(g+g0);
                    testCase.verifyEqual(physical.metricMatrix(end,end),expected,RelTol=1e-14);
                    testCase.verifyError(@() basis.projectionOnGrid(z(1:end-1),w(1:end-1),variable="G",prepareRecovery=false),"IMProjection:MissingEndpointSample");
                    mda = IMSolverSpectral(nEVP=32,coordinateKind="z").solveEVP(IMInternalModes.meanDensityAnomalyModes(background,g0=boundary.g0,gd=boundary.gd),nModes=1);
                    physical = mda.projectionOnGrid(z,w,variable="G",prepareRecovery=false);
                    testCase.verifyEqual(physical.metricMatrix(end,end),g0/g,RelTol=1e-14);
                end
            end
        end

        function mdaUsesItsOwnSignedGFormAndExplicitRecoveryPolicy(testCase)
            background = signedBackground();
            energy = IMGeneralizedEnergyBoundary(g0=-.4,gd=.3);
            basis = IMSolverSpectral(nEVP=64,coordinateKind="z").solveEVP(IMInternalModes.meanDensityAnomalyModes(background,g0=energy.g0,gd=energy.gd),nModes=5);
            [z,w] = endpointRule(background.zDomain,96);
            [H,Hplus] = independentForms(basis,"G",energy,256);
            transform = basis.scalarTransform(z,w,variable="G");
            testCase.verifyClass(transform,"IMScalarDiscreteTransform");
            testCase.verifyTrue(transform.success);
            testCase.verifyEqual(transform.targetGramMatrix,H,RelTol=3e-9,AbsTol=3e-9);
            testCase.verifyEqual(transform.projection.majorantGramMatrix,Hplus,RelTol=3e-9,AbsTol=3e-9);
            testCase.verifyTrue(any(eig(H) < -.01));
            testCase.verifyEqual(transform.columnLabels,basis.columnLabels);
            testCase.verifyEqual(transform.modeNumber,basis.modeNumber);
            testCase.verifyEqual(transform.z,z,AbsTol=0);
            testCase.verifyEqual(transform.weights,w,AbsTol=0);
            A = basis.G(z);
            metric = independentMetric(background,energy,"meanDensityAnomaly","G",z,w);
            B = transform.targetGramMatrix\(A'*diag(metric));
            source = independentSource(z,"G");
            testCase.verifyEqual(transform.physicalDualMatrix,B,RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(transform.projection.projectPhysical(source),B*source,RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(transform.transformForward(A),eye(basis.nColumns),AbsTol=3e-11);
            testCase.verifyError(@() basis.scalarTransform(z,w,variable="F"),"IMScalarProjection:UnavailableProjection");
            testCase.verifyError(@() basis.scalarTransform(z,w,variable="G",policy="sampledLeastSquares"),"IMScalarProjection:MissingFitFactor");
            fit = diag(linspace(.6,1.4,numel(z)));
            leastSquares = basis.scalarTransform(z,w,variable="G",policy="sampledLeastSquares",fitFactor=fit);
            testCase.verifyTrue(leastSquares.success);
            testCase.verifyEqual(leastSquares.forwardMatrix,(fit*A)\fit,RelTol=2e-12,AbsTol=2e-12);
            testCase.verifyEqual(leastSquares.physicalDualMatrix,B,RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyGreaterThan(norm((leastSquares.forwardMatrix-transform.forwardMatrix)*source),1e-5);
        end

        function requiredEndpointTracesAndInactiveLimitsAreLiteral(testCase)
            background = signedBackground();
            energy = IMGeneralizedEnergyBoundary(g0=-.4,gd=.3);
            basis = IMSolverSpectral(nEVP=48,coordinateKind="z").solveEVP(IMInternalModes.geostrophicAPVModes(background,g0=energy.g0,gd=energy.gd),nModes=3);
            [z,w] = endpointRule(background.zDomain,48);
            testCase.verifyError(@() basis.scalarTransform(z(1:end-1),w(1:end-1),variable="G"),"IMProjection:MissingEndpointSample");
            testCase.verifyError(@() basis.scalarTransform(z(2:end),w(2:end),variable="G"),"IMProjection:MissingEndpointSample");
            interior = basis.scalarTransform(z(2:end-1),w(2:end-1),variable="F");
            testCase.verifyTrue(interior.success);
            rigid = signedBackground("rigidLid");
            noTraces = IMSolverSpectral(nEVP=48,coordinateKind="z").solveEVP(IMInternalModes.geostrophicAPVModes(rigid,g0=Inf,gd=Inf),nModes=3);
            unaugmented = noTraces.scalarTransform(z(2:end-1),w(2:end-1),variable="G");
            testCase.verifyTrue(unaugmented.success);
            testCase.verifyEqual(unaugmented.metricMatrix,diag(w(2:end-1).*rigid.N2(z(2:end-1))/rigid.g),AbsTol=0);
            freeInfinite = IMSolverSpectral(nEVP=48,coordinateKind="z").solveEVP(IMInternalModes.geostrophicAPVModes(background,g0=Inf,gd=Inf),nModes=2);
            testCase.verifyError(@() freeInfinite.scalarTransform(z(2:end-1),w(2:end-1),variable="G"),"IMProjection:MissingEndpointSample");
            freePhysical = freeInfinite.projectionOnGrid(z,w,variable="G",prepareRecovery=false);
            testCase.verifyEqual(freePhysical.metricMatrix(end,end),1,AbsTol=0);
            singular = IMSolverSpectral(nEVP=48,coordinateKind="z").solveEVP(IMInternalModes.geostrophicAPVModes(background,g0=-background.g,gd=Inf),nModes=3);
            singularPhysical = singular.projectionOnGrid(z(2:end-1),w(2:end-1),variable="G",prepareRecovery=false);
            testCase.verifyEqual(singularPhysical.metricMatrix,diag(w(2:end-1).*background.N2(z(2:end-1))/background.g),AbsTol=0);
        end

        function unsupportedDerivativeAndCompanionObservationsStayUnsupported(testCase)
            background = signedBackground();
            solver = IMSolverSpectral(nEVP=48,coordinateKind="z");
            evp = IMInternalModes(zDomain=background.zDomain,N2=@(z) ones(size(z)),formulation="G",modeFamily="hydrostatic",r=.5,rotationRate=background.rotationRate,latitude=background.latitude,g=background.g, ...
                surfaceBoundary=IMBoundaryCondition(a=1,b=0,c=0,d=1),bottomBoundary=IMBoundaryCondition.dirichlet());
            source = solver.solveEVP(evp,nModes=3);
            wrapper = source.withNormalization("rmsG");
            [z,w] = endpointRule(background.zDomain,48);
            testCase.verifySize(wrapper.G(z),[numel(z) 3]);
            testCase.verifyError(@() wrapper.scalarTransform(z,w,variable="G"),"IMProjection:UnsupportedDerivativeTrace");
            testCase.verifyError(@() wrapper.selectColumns([]).scalarTransform(z,w,variable="G"),"IMProjection:UnsupportedDerivativeTrace");
            coupled = solver.solveEVP(IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes(background,k=.8,alpha0=.7,alphaD=.4),nModes=2);
            testCase.verifySize(coupled.G(z),[numel(z) coupled.nColumns]);
            testCase.verifyError(@() coupled.scalarTransform(z,w,variable="G"),"IMScalarProjection:UnavailableProjection");
            waveBackground = IMBackground.constantStratification(N0=2,zDomain=[-1 0],g=2,rotationRate=background.rotationRate,latitude=33,surfaceCondition="freeSurface");
            wave = solver.solveEVP(IMInternalModes.waveModesAtWavenumber(waveBackground,k=.3),nModes=3);
            testCase.verifyError(@() wave.scalarTransform(z,w,variable="F"),"IMScalarProjection:UnavailableProjection");
            zero = solver.solveGeostrophicZeroAPVModes(IMGeostrophicZeroAPVModes.atWavenumber(background,.3,endpoints=["surface","bottom"]));
            testCase.verifyError(@() zero.scalarTransform(z,w,variable="G"),"IMScalarProjection:UnavailableProjection");
        end

        function selectedNormalizationCovariesWithFixedRecoveryFactors(testCase)
            background = signedBackground();
            energy = IMGeneralizedEnergyBoundary(g0=-.4,gd=.3);
            full = IMSolverSpectral(nEVP=64,coordinateKind="z").solveEVP(IMInternalModes.geostrophicAPVModes(background,g0=energy.g0,gd=energy.gd),nModes=4);
            columns = [full.nColumns 1 2];
            basis = full.selectColumns(columns);
            changed = basis.withNormalization("maxAbsG");
            T = normalizationMap(basis,changed);
            [z,w] = endpointRule(background.zDomain,64);
            physical = basis.projectionOnGrid(z,w,variable="G",prepareRecovery=false);
            changedPhysical = changed.projectionOnGrid(z,w,variable="G",prepareRecovery=false);
            fullPhysical = full.projectionOnGrid(z,w,variable="G",prepareRecovery=false);
            % Retained native normalizations may accumulate roundoff in a
            % different matrix shape. Compare the actual physical forms in
            % dimensionless coefficient units, rather than requiring bits.
            selectedScale = chol(physical.majorantGramMatrix)\eye(basis.nColumns);
            targetDifference = selectedScale'*(physical.targetGramMatrix-fullPhysical.targetGramMatrix(columns,columns))*selectedScale;
            pairingScale = max(abs(fullPhysical.samplePairingMatrix(columns,:)),[],2);
            pairingDifference = max(abs(physical.samplePairingMatrix-fullPhysical.samplePairingMatrix(columns,:)),[],2)./max(pairingScale,realmin);
            testCase.verifyLessThan(norm(targetDifference,2),3e-13);
            testCase.verifyLessThan(max(pairingDifference),3e-13);
            testCase.verifyEqual(physical.physicalDualMatrix,physical.targetGramMatrix\physical.samplePairingMatrix,RelTol=3e-12,AbsTol=3e-12);
            C = chol(physical.majorantGramMatrix)\eye(basis.nColumns);
            D = diag(linspace(.8,1.3,numel(z)));
            F = diag(linspace(1.4,.6,numel(z)));
            original = basis.scalarTransform(z,w,variable="G",coefficientScale=C,observationFactor=D,fitFactor=F);
            normalized = changed.scalarTransform(z,w,variable="G",coefficientScale=T\C,observationFactor=D,fitFactor=F);
            testCase.verifyTrue(original.success && normalized.success);
            testCase.verifyEqual(changed.G(z),basis.G(z)*T,RelTol=3e-14,AbsTol=3e-14);
            testCase.verifyEqual(changedPhysical.targetGramMatrix,T'*physical.targetGramMatrix*T,RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(changedPhysical.majorantGramMatrix,T'*physical.majorantGramMatrix*T,RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(changedPhysical.samplePairingMatrix,T'*physical.samplePairingMatrix,RelTol=3e-14,AbsTol=3e-14);
            testCase.verifyEqual(normalized.physicalDualMatrix,T\original.physicalDualMatrix,RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(normalized.forwardMatrix,T\original.forwardMatrix,RelTol=3e-12,AbsTol=3e-12);
            input = independentSource(z,"G");
            testCase.verifyEqual(normalized.transformForward(input),T\original.transformForward(input),RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(normalized.transformBack(normalized.transformForward(input)),original.transformBack(original.transformForward(input)),RelTol=3e-12,AbsTol=3e-12);
            testCase.verifyEqual(normalized.recoverySpecification.observationFactor,D,AbsTol=0);
            testCase.verifyEqual(normalized.recoverySpecification.fitFactor,F,AbsTol=0);
            testCase.verifyEqual(normalized.recoverySpecification.coefficientScale,T\C,AbsTol=0);
            testCase.verifyEqual(original.columnLabels,full.columnLabels(columns));
            supplied = IMDiscreteTransform.fromObservationMatrices(basis.G(z),policy="sampledLeastSquares",columnLabels=basis.columnLabels,fitFactor=F,physicalDual=physical.physicalDualMatrix);
            rebound = basis.scalarTransform(z,w,variable="G",recovery=supplied);
            testCase.verifyEqual(rebound.forwardMatrix,supplied.forwardMatrix,AbsTol=0);
            testCase.verifyEqual(rebound.recoverySpecification,supplied.recoverySpecification);
            collection = IMBasisCollection({basis,changed},basisIndex=[2 1 2]);
            for page = [3 2 1]
                expected = collection.bases{collection.basisIndex(page)}.projectionOnGrid(z,w,variable="G",prepareRecovery=false);
                actual = collection.projectionOnGrid(z,w,page=page,variable="G",prepareRecovery=false);
                testCase.verifyEqual(actual.targetGramMatrix,expected.targetGramMatrix,AbsTol=0);
                testCase.verifyEqual(actual.physicalDualMatrix,expected.physicalDualMatrix,AbsTol=0);
                testCase.verifyEqual(actual.columnLabels,expected.columnLabels);
            end
        end

        function hydrostaticActivityEmptyAndFailedRecoveryRemainScientific(testCase)
            background = signedBackground("rigidLid");
            full = IMSolverAnalytical().solveEVP(IMInternalModes.hydrostaticFModes(background),nModes=4);
            [z,w] = endpointRule(background.zDomain,64);
            for variable = ["F","G"]
                transform = full.scalarTransform(z,w,variable=variable);
                testCase.verifyTrue(transform.success);
                active = true(1,full.nColumns);
                if variable == "G", active(full.modeInfo.branch == "barotropic") = false; end
                testCase.verifyEqual(transform.activeColumnMask,active);
                testCase.verifyEqual(transform.forwardMatrix*transform.inverseMatrix,diag(double(active)),AbsTol=3e-12);
                empty = full.selectColumns([]).scalarTransform(z,w,variable=variable);
                explicitEmpty = full.projectionOnGrid(z,w,variable=variable,columns=[],prepareRecovery=false);
                testCase.verifySize(empty.forwardMatrix,[0 numel(z)]);
                testCase.verifySize(empty.inverseMatrix,[numel(z) 0]);
                testCase.verifySize(empty.targetGramMatrix,[0 0]);
                testCase.verifyEmpty(empty.columnLabels);
                testCase.verifySize(explicitEmpty.physicalDualMatrix,[0 numel(z)]);
            end
            endpoints = background.zDomain.';
            failed = full.scalarTransform(endpoints,[.5;.5],variable="G");
            testCase.verifyFalse(failed.success);
            testCase.verifyEqual(failed.activeColumnMask,[false true true true]);
            testCase.verifyEqual(failed.projection.physicalDualStatus,"available");
            testCase.verifySize(failed.physicalDualMatrix,[full.nColumns 2]);
            input = [1+2i;-.4+.7i];
            testCase.verifyEqual(failed.projection.projectPhysical(input),failed.physicalDualMatrix*input,AbsTol=0);
            testCase.verifyEmpty(failed.forwardMatrix);
            mda = IMSolverSpectral(nEVP=48,coordinateKind="z").solveEVP(IMInternalModes.meanDensityAnomalyModes(background,g0=0,gd=0),nModes=3);
            constant = find(isinf(mda.h));
            testCase.verifyEqual(numel(constant),1);
            mdaTransform = mda.scalarTransform(z,w,variable="G");
            testCase.verifyTrue(mdaTransform.success);
            testCase.verifyTrue(mdaTransform.activeColumnMask(constant));
            testCase.verifyGreaterThan(norm(mdaTransform.inverseMatrix(:,constant)),1);
            free = signedBackground();
            for factory = {@IMInternalModes.hydrostaticFModes,@IMInternalModes.hydrostaticGModes}
                basis = IMSolverAnalytical().solveEVP(factory{1}(free),nModes=3);
                for variable = ["F","G"]
                    transform = basis.scalarTransform(z,w,variable=variable);
                    testCase.verifyTrue(transform.success);
                    testCase.verifyEqual(transform.columnLabels,basis.columnLabels);
                    testCase.verifyEqual(transform.transformForward(basis.evaluate(z,variable=variable)),eye(basis.nColumns),AbsTol=3e-11);
                end
            end
        end

        function neutralDualAndFixedRuleReuseDoNotResolveTheEVP(testCase)
            background = signedBackground();
            counter = IMScalarTransformTestCounter();
            solver = IMScalarTransformCountingSolver(counter);
            energy = IMGeneralizedEnergyBoundary(g0=-.4,gd=.3);
            full = solver.solveEVP(IMInternalModes.geostrophicAPVModes(background,g0=energy.g0,gd=energy.gd),nModes=4);
            initialCount = counter.solveCount;
            testCase.verifyGreaterThan(initialCount,0);
            [z,w] = endpointRule(background.zDomain,64);
            selected = full.selectColumns([3 1]).withNormalization("maxAbsG");
            projection = selected.projectionOnGrid(z,w,variable="G",prepareRecovery=false);
            transform = selected.scalarTransform(z,w,variable="G");
            testCase.verifyTrue(transform.success);
            testCase.verifyEqual(counter.solveCount,initialCount);
            testCase.verifyEqual(transform.z,z,AbsTol=0);
            testCase.verifyEqual(transform.weights,w,AbsTol=0);
            testCase.verifyEqual(transform.physicalDualMatrix,projection.physicalDualMatrix,AbsTol=0);
            neutralBackground = IMBackground(N2=@(z) zeros(size(z)),zDomain=[-1 0],rotationRate=background.rotationRate,latitude=33,g=2,surfaceCondition="freeSurface");
            neutral = IMNeutralLabelModesBasis(neutralBackground,neutralIntervals=[-1 0],g0=0,gd=0,nModesPerInterval=3);
            physical = neutral.projectionOnGrid(z,w,variable="APV",prepareRecovery=false);
            adapter = neutral.scalarTransform(z,w,variable="APV");
            testCase.verifyTrue(adapter.success);
            testCase.verifyEqual(adapter.physicalDualMatrix,physical.physicalDualMatrix,AbsTol=0);
            testCase.verifyEqual(adapter.projection.projectionKind,"prescribedDual");
            testCase.verifyFalse(adapter.projection.supportsGramAssessment);
            testCase.verifyEqual(adapter.inverseMatrix,neutral.APV(z),AbsTol=0);
            testCase.verifyEqual(adapter.transformForward(neutral.APV(z)),eye(neutral.nColumns),AbsTol=3e-12);
            input = independentSource(z,"F");
            testCase.verifyEqual(adapter.projection.projectPhysical(input),physical.physicalDualMatrix*input,AbsTol=0);
            empty = neutral.selectColumns([]).scalarTransform(z,w,variable="APV");
            testCase.verifySize(empty.forwardMatrix,[0 numel(z)]);
            testCase.verifySize(empty.physicalDualMatrix,[0 numel(z)]);
        end
    end
end

function background = signedBackground(surfaceCondition)
arguments
    surfaceCondition (1,1) string = "freeSurface"
end
background = IMBackground.constantStratification(N0=1,zDomain=[-1 0],g=2,rotationRate=1/(2*sind(33)),latitude=33,surfaceCondition=surfaceCondition);
end

function [z,w] = endpointRule(domain,n)
[interior,weights] = legpts(n,domain);
z = [domain(1);interior;domain(2)];
w = [0;weights(:);0];
end

function values = independentSource(z,variable)
if variable == "F"
    values = [exp(z)+.2*z.^3+1i*cos(2.3*z),sin(1.7*z)+1i*(1+z.^2)];
else
    values = [cos(1.9*z)+.4*z+1i*exp(.7*z),1+z.^3+1i*sin(2.7*z)];
end
end

function [H,Hplus] = independentForms(basis,variable,energy,n)
[z,w] = endpointRule(basis.zDomain,n);
A = basis.evaluate(z,variable=variable);
metric = independentMetric(basis.evp.background,energy,string(basis.evp.name),variable,z,w);
H = A'*(metric.*A);
Hplus = A'*(abs(metric).*A);
end

function metric = independentMetric(background,energy,family,variable,z,w)
if variable == "F"
    metric = w;
    return
end
metric = w.*background.N2(z)/background.g;
if ismember(family,["meanDensityAnomaly","meanDensityAnomalyModes"])
    surface = energy.g0/background.g;
elseif background.surfaceCondition == "freeSurface"
    if isinf(energy.g0)
        surface = 1;
    elseif energy.g0 == -background.g
        surface = 0;
    else
        surface = energy.g0/(background.g+energy.g0);
    end
else
    surface = energy.g0/background.g;
end
bottom = energy.gd/background.g;
if isfinite(surface), metric(end) = metric(end)+surface; end
if isfinite(bottom), metric(1) = metric(1)+bottom; end
end

function T = normalizationMap(basis,changed)
old = basis.normalizationFactors();
new = changed.normalizationFactors();
T = diag(old./new);
end
