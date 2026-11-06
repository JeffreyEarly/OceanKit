classdef IMBasisCapabilityTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function configurePath(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,"UnitTestsV2","TestSupport")));
        end
    end
    methods (Test)
        function waveDerivativesAgreeAcrossBackendsAndRepeatedPages(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            for solver = {IMSolverSpectral(nEVP=48),IMSolverAnalytical()}
                basis = createIMTestBasis("wave",bg,solver=solver{1});
                collection = IMBasisCollection({basis},basisIndex=[1 1 1]);
                testCase.verifyEqual(basis.capabilities.variables,["F","G"]);
                testCase.verifyEqual(basis.capabilities.derivativeOrders,{0:1,0:2});
                testCase.verifyEqual(collection.metadata.capabilities,basis.capabilities);
                z = linspace(-1000,0,21).';
                Fz = collection.evaluate(z,variable="F",derivativeOrder=1,pages=[3 1 3],columns=[3 1 3]);
                Gzz = collection.evaluate(z,variable="G",derivativeOrder=2,pages=[3 1 3],columns=[3 1 3]);
                testCase.verifyEqual(Fz,bg.D*Gzz,AbsTol=3e-14);
                scalar = basis.Fz(z);
                testCase.verifyEqual(Fz,repmat(scalar(:,[3 1 3]),1,1,3),AbsTol=3e-17);
            end
        end
        function projectionAndRecoveryCapabilitiesRemainSeparate(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            hydro = createIMTestBasis("hydrostaticF",bg,solver=IMSolverAnalytical());
            wave = createIMTestBasis("wave",bg,solver=IMSolverAnalytical());
            testCase.verifyEqual([hydro.capabilities.projectionOperations.variable],["F","G"]);
            testCase.verifyEmpty(wave.capabilities.projectionOperations);
            for item = {hydro,wave}
                recovery = item{1}.capabilities.sampledRecovery;
                testCase.verifyEqual(recovery.availability,"requiresPreparedSystem");
                testCase.verifyFalse(isfield(recovery,"rank"));
                testCase.verifyEqual(recovery.noiseEvidence,"requiresDeclaredObservationFactor");
            end
            transform = IMDiscreteTransform.fromObservationMatrices(wave.F(-500),columnLabels=wave.columnLabels,observationFactor=1,failureMode="diagnostic");
            testCase.verifyFalse(transform.success);
        end
        function neutralTracesAndAPVDualRequirementsAreExplicit(testCase)
            bg = IMBackground(N2=@(z) 1e-4*double(z < -200),zDomain=[-1000 0],surfaceCondition="freeSurface");
            basis = IMNeutralLabelModesBasis(bg,neutralIntervals=[-200 0],nModesPerInterval=2);
            testCase.verifyEqual(basis.capabilities.derivativeConvention,"classicalIntervalInterior");
            testCase.verifyEqual(basis.capabilities.traceConvention,"insideIntervalEndpoints");
            testCase.verifyEqual(basis.capabilities.projectionOperations.kind,"prescribedContinuousDual");
            Gz = basis.Gz([-201;-200;-100;0]);
            testCase.verifyEqual(Gz(1,:),[0 0]);
            testCase.verifyGreaterThan(abs(Gz(2,:)),[0 0]);
            testCase.verifyError(@() basis.evaluate(-100,variable="APV",derivativeOrder=1),"IMNeutralLabelModesBasis:UnsupportedDerivative");
        end
        function emptySelectionsKeepTheirDeclaredOperations(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            full = createIMTestBasis("wave",bg,solver=IMSolverAnalytical()); empty = full.selectColumns([]);
            testCase.verifyEqual(empty.capabilities,full.capabilities);
            testCase.verifySize(empty.evaluate([-1000;0],variable="G",derivativeOrder=2),[2 0]);
            testCase.verifyError(@() empty.evaluate(-500,variable="F",derivativeOrder=2),"IMAnalyticalInternalModesBasis:UnsupportedDerivative");
        end
        function unsupportedRequestsDoNotEvaluateIntrinsicProfiles(testCase)
            calls = 0;
            problem = IMInternalModes(N2=@(z) ones(size(z)),zDomain=[-1 0],formulation="G",p=1,q=0,r=1,FfromGz=@diagnostic);
            basis = IMSolverSpectral(nEVP=32).solveEVP(problem,nModes=2);
            calls = 0;
            description = basis.capabilities; testCase.verifyEqual(description.derivativeOrders{1},0);
            testCase.verifyError(@() basis.evaluate(-.5,variable="F",derivativeOrder=1),"IMInternalModesBasis:UnsupportedDerivative");
            testCase.verifyEqual(calls,0);
            function values = diagnostic(~,dG,~,~)
                calls = calls+1; values = dG.^2+1;
            end
        end
        function productFactorsPreserveTheSameDerivativeContract(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            basis = createIMTestBasis("wave",bg,solver=IMSolverAnalytical()); c = IMBasisCollection({basis});
            f = c.productFactor(id="Fz",variable="F",derivativeOrder=1);
            g = c.productFactor(id="Gzz",variable="G",derivativeOrder=2);
            z = [-900;-600;-100]; testCase.verifyEqual(f.evaluate(z,[3 1],"ref"),bg.D*g.evaluate(z,[3 1],"ref"),AbsTol=2e-14);
            testCase.verifyError(@() c.productFactor(id="bad",variable="G",derivativeOrder=3),"IMBasisCollection:UnsupportedFactor");
        end
        function fixedFrequencyDoesNotClaimAFixedWavenumberEnergyNorm(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            basis = createIMTestBasis("frequency",bg,solver=IMSolverAnalytical());
            testCase.verifyFalse(ismember("waveEnergy",basis.capabilities.normalizationNames));
            testCase.verifyTrue(ismember("kSquared",basis.capabilities.spectrumNames));
            testCase.verifyEqual(basis.Fz([-900;-100]),bg.D*basis.Gzz([-900;-100]),AbsTol=3e-14);
        end
        function endpointResponsesExposeTheirRegularPairDerivatives(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            problem = IMGeostrophicZeroAPVModes.atWavenumber(bg,1e-3,endpoints=["surface","bottom"]);
            basis = IMSolverAnalytical().solveGeostrophicZeroAPVModes(problem); z = [-1000;-600;0];
            testCase.verifyEqual(basis.Fz(z),-bg.N2(z)/bg.g.*basis.G(z),AbsTol=2e-15);
            testCase.verifyEqual(basis.Gz(z),-bg.g*1e-6/bg.f0^2.*basis.F(z),AbsTol=2e-13);
            testCase.verifyEmpty(basis.capabilities.projectionOperations);
        end
        function requestedEmptyAndSelectedEmptyPreserveEVPNormalizations(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            full = createIMTestBasis("enstrophy",bg,nModes=2); empty = createIMTestBasis("enstrophy",bg,nModes=0);
            testCase.verifyEqual(empty.capabilities,full.selectColumns([]).capabilities);
            testCase.verifyEqual(empty.normalization,"generalizedPotentialEnstrophy");
            testCase.verifySize(empty.Fz([-1000;-500;0]),[3 0]);
        end
        function knownIntrinsicSectorsDeclareTheirFirstFDerivative(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            basis = createIMTestBasis("hydrostaticF",bg,nModes=1);
            testCase.verifyEqual(basis.Fz([-1000;0]),zeros(2,1),AbsTol=0);
            testCase.verifyEqual(basis.capabilities.derivativeOrders,{0:1,0:2});
            testCase.verifyError(@() basis.evaluate(-500,variable="Fz"),"IMInternalModesBasis:UnsupportedVariable");
        end
    end
end
