classdef IMModeBasisProvenanceTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function configurePath(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,"UnitTestsV2","TestSupport")));
        end
    end
    methods (Test)
        function generalizedEnergyBoundaryRemainsIndependentOfPhysicalSurface(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            energy = IMGeneralizedEnergyBoundary(g0=-.04,gd=.03);
            b = createIMTestBasis("APV",bg,energyBoundary=energy);
            r = b.scientificSpecification;
            testCase.verifyEqual(r.background.surfaceCondition,"rigidLid");
            testCase.verifyEqual(r.resolved.energyBoundary,struct(g0=-.04,gd=.03));
            testCase.verifyEqual(r.requested.energyBoundary,struct(g0=-.04,gd=.03));
            testCase.verifyEqual(r.resolved.horizontal.kind,"independent");
        end
        function resolvedDefaultAlphaIsDistinctFromRequestedDefault(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            b = createIMTestBasis("enstrophy",bg); r = b.scientificSpecification;
            testCase.verifyEmpty(r.requested.alpha0);
            testCase.verifyEqual(r.resolved.evp.parameters.alpha0Source,"stratification");
            testCase.verifyEqual(r.resolved.evp.parameters.bEffective,250,AbsTol=5e-12);
            testCase.verifyEqual(r.resolved.evp.parameters.alpha0,bg.f0^2/250,RelTol=2e-14);
        end
        function constructionUsesConfiguredCoordinateWithoutClaimingConvergence(testCase)
            bg = IMBackground.exponentialStratification(N0=.01,b=500,zDomain=[-1000 0],surfaceCondition="freeSurface");
            for coordinate = ["z","wkb","bandwidth"]
                b = createIMTestBasis("wave",bg,solver=IMSolverSpectral(nEVP=48,coordinateKind=coordinate)); r = b.constructionRecord;
                testCase.verifyEqual(r.requested.solver.coordinateKind,coordinate);
                testCase.verifyEqual(r.resolved.solver.coordinateKind,coordinate);
                testCase.verifySize(r.resolved.solver.zNative,[48 1]);
                testCase.verifyEqual(r.selection.request.nModes,3);
                testCase.verifyNumElements(r.selection.sourceEigenpairIndices,3);
                testCase.verifyEqual(r.evidence.solveAccuracy.status,"unverified");
                testCase.verifyEqual(r.evidence.referenceConvergence.status,"unverified");
            end
        end
        function finiteDifferenceRequestedGridIsReproducible(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface"); z = linspace(-1000,0,81).';
            b = createIMTestBasis("wave",bg,solver=IMSolverFiniteDifference(z=z)); r = b.constructionRecord;
            testCase.verifyEqual(sort(r.requested.solver.zNative),z);
            testCase.verifyEqual(r.resolved.solver.zNative,r.requested.solver.zNative);
        end
        function emptyEVPKeepsItsScientificAndConstructionRequests(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            b = createIMTestBasis("wave",bg,nModes=0,solver=IMSolverSpectral(nEVP=48,coordinateKind="bandwidth"));
            testCase.verifyEqual(b.scientificSpecification.requested.k,1e-3);
            testCase.verifyEqual(b.scientificSpecification.resolved.evp.name,"waveModesAtWavenumber");
            testCase.verifyEqual(b.constructionRecord.selection.request.nModes,0);
            testCase.verifyEqual(b.constructionRecord.requested.solver.nEVP,48);
            testCase.verifyEmpty(b.constructionRecord.selection.sourceEigenpairIndices);
        end
        function analyticalAndBarotropicBypassesDescribeActualConstruction(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            b = createIMTestBasis("wave",bg,solver=IMSolverAnalytical());
            testCase.verifyEqual(b.constructionRecord.resolved.backendKind,"analytical");
            testCase.verifyEqual(b.constructionRecord.resolved.sourceClass,"IMAnalyticalInternalModesBasis");
            barotropic = createIMTestBasis("hydrostaticF",bg,solver=IMSolverAnalytical(),nModes=1);
            testCase.verifyTrue(barotropic.scientificSpecification.resolved.includesAnalyticalBarotropic);
            testCase.verifyEqual(barotropic.h,Inf);
            testCase.verifyEqual(barotropic.constructionRecord.selection.request.nModes,1);
        end
        function DDIntervalsPreserveFixedEndpointAndCoefficientRoles(testCase)
            bg = IMBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0],surfaceCondition="freeSurface"); intervals = [-1000 -800;-200 0];
            b = IMNeutralLabelModesBasis(bg,neutralIntervals=intervals,nModesPerInterval=[1 2],g0=.03,gd=Inf); r = b.scientificSpecification;
            testCase.verifyEqual(r.requested.nModesPerInterval,[1 2]); testCase.verifyEqual(r.resolved.intervals,intervals);
            testCase.verifyEqual(r.resolved.intervalBoundaryTypes,repmat("Dirichlet",2,2));
            testCase.verifyEqual(r.resolved.energyBoundaryRole,"provenanceOnly");
            testCase.verifyEqual(r.resolved.energyBoundary,struct(g0=.03,gd=Inf));
            testCase.verifyEqual(b.G(intervals(:)),zeros(4,3));
        end
        function criticalIntervalsAndEndpointSelectionsRemainExplicit(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            critical = IMBackground(N2=@(z) bg.f0^2*ones(size(z)),zDomain=bg.zDomain,surfaceCondition="freeSurface");
            b = IMZeroDepthWaveModesBasis(critical,1e-3,criticalIntervals=[-900 -700],nModesPerInterval=2);
            testCase.verifyEqual(b.scientificSpecification.requested.criticalIntervals,[-900 -700]);
            testCase.verifyEqual(b.scientificSpecification.resolved.intervalKind,"critical");
            testCase.verifyEqual(b.constructionRecord.resolved.backendKind,"analyticalDDKernel");
            endpoint = IMSolverAnalytical().solveGeostrophicZeroAPVModes(IMGeostrophicZeroAPVModes.atWavenumber(bg,1e-3,endpoints="surface"));
            testCase.verifyEqual(endpoint.scientificSpecification.resolved.endpoints,"surface");
            testCase.verifyEqual(endpoint.scientificSpecification.background.profileKind,bg.profileKind);
            testCase.verifyEqual(endpoint.scientificSpecification.background.profileParameters,bg.profileParameters);
            testCase.verifyEqual(endpoint.scientificSpecification.background.rotationRate,bg.rotationRate);
            testCase.verifyEqual(endpoint.scientificSpecification.background.latitude,bg.latitude);
        end
        function selectedRecordsComposeWithoutLosingOriginalRequest(testCase)
            bg = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            original = createIMTestBasis("frequency",bg,nModes=4,solver=IMSolverAnalytical());
            selected = original.selectColumns([4 2 1]).selectColumns([2 1]);
            testCase.verifyEqual(selected.scientificSpecification.requested,original.scientificSpecification.requested);
            testCase.verifyEqual(selected.constructionRecord.requested,original.constructionRecord.requested);
            testCase.verifyEqual(selected.constructionRecord.selection.initialColumnIndices,original.basisContractState.initialColumnIndices([2 4]));
            testCase.verifyEqual(selected.columnLabels,original.columnLabels([2 4]));
            testCase.verifyEqual(selected.kSquared,original.kSquared([2 4]));
            record = selected.scientificSpecification; record.resolved.columnLabels(:) = "changed";
            testCase.verifyEqual(selected.columnLabels,original.columnLabels([2 4]));
        end
    end
end
