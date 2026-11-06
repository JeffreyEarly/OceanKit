classdef IMModeCreationTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addProvider(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
        end
    end

    methods (Test)
        function bulkRequestsPreservePagesAndIndependentCounts(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            collection = IMSolverAnalytical().solveWaveModesAtWavenumbers(background,[2e-4 1e-4 2e-4],nModes=[3 4 3]);
            testCase.verifyEqual(collection.kappa,[2e-4 1e-4 2e-4])
            testCase.verifyEqual(collection.basisIndex,[1 2 1])
            testCase.verifyEqual(numel(collection.bases),2)
            testCase.verifyEqual(collection.bases{1}.nColumns,3)
            testCase.verifyEqual(collection.bases{2}.nColumns,4)
            z = linspace(-1000,0,21).';
            F = collection.evaluate(z,variable="F",pages=[3 1]);
            testCase.verifyEqual(F(:,:,1),F(:,:,2))
            distinct = IMSolverAnalytical().solveWaveModesAtWavenumbers(background,[1e-4 1e-4],nModes=[2 3]);
            testCase.verifyEqual(cellfun(@(basis) basis.nColumns,distinct.bases),[2 3])
        end

        function intervalDeclarationsAndNormalizationContractsAreValidated(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            energy = IMGeneralizedEnergyBoundary(g0=0,gd=0);
            testCase.verifyError(@() IMNeutralLabelModesBasis(background,neutralIntervals=[-200 0],g0=energy.g0,gd=energy.gd,nModesPerInterval=0),"IMNeutralLabelModesBasis:InconsistentInterval")
            weak = IMBackground(N2=@(z) 1e-30*ones(size(z)),zDomain=[-1000 0],surfaceCondition="freeSurface");
            testCase.verifyError(@() IMNeutralLabelModesBasis(weak,neutralIntervals=[-1000 0],nModesPerInterval=0),"IMNeutralLabelModesBasis:InconsistentInterval")
            testCase.verifyError(@() IMZeroDepthWaveModesBasis(background,1e-4,criticalIntervals=[-200 0],nModesPerInterval=0),"IMZeroDepthWaveModesBasis:InconsistentInterval")
            frequency = IMSolverAnalytical().solveEVP(IMInternalModes.waveModesAtFrequency(background,omega=.002),nModes=2);
            testCase.verifyFalse(ismember("waveEnergy",frequency.normalizationNames()))
            testCase.verifyError(@() IMInternalModes.waveModesAtFrequency(background,omega=abs(background.f0)),"IMInternalModes:InertialFrequency")
            testCase.verifyError(@() IMSolverAnalytical().solveEVP(IMInternalModes.meanDensityAnomalyModes(background,g0=energy.g0,gd=energy.gd),nModes=2),"IMSolverAnalytical:UnsupportedFamily")
        end

        function criticalRigidColumnHasAnEmptyRegularWaveInventory(testCase)
            f0 = 2*7.2921e-5*sind(33);
            background = IMBackground(N2=@(z) f0^2*ones(size(z)),zDomain=[-1000 0],surfaceCondition="rigidLid");
            basis = IMSolverSpectral().solveEVP(IMInternalModes.waveModesAtWavenumber(background,k=1e-4),nModes=0);
            testCase.verifyEqual(basis.nColumns,0)
            testCase.verifyEqual(basis.h,zeros(1,0))
            testCase.verifyError(@() IMSolverSpectral().solveEVP(IMInternalModes.waveModesAtWavenumber(background,k=1e-4),nModes=1),"IMSolver:InsufficientModes")
        end

        function wavesKeepExplicitCountedLadders(testCase)
            N2 = @(z) 1e-4*max(min((-z-100)/100,1),0).^2;
            background = IMBackground(N2=N2,zDomain=[-1000 0],surfaceCondition="freeSurface");
            k = 2e-4;
            basis = IMSolverSpectral(nEVP=160,coordinateKind="z").solveEVP(IMInternalModes.waveModesAtWavenumber(background,k=k),branch=["surface","internal","subinertial"],nModes=[1 3 2]);
            testCase.verifyClass(basis,"IMInternalModesBasis")
            testCase.verifyEqual(basis.nColumns,6)
            testCase.verifyEqual(basis.modeInfo.branch,["surface";"internal";"internal";"internal";"subinertial";"subinertial"])
            testCase.verifyEqual(basis.modeInfo.branchIndex,[1;1;2;3;1;2])
            testCase.verifyEqual(string(basis.modeInfo.columnLabel),basis.columnLabels(:))
            testCase.verifyEqual(basis.normalization,"waveEnergy")
            testCase.verifyGreaterThan(basis.h(1:4),zeros(1,4))
            testCase.verifyLessThan(basis.h(5:6),zeros(1,2))
            testCase.verifyLessThan(diff(basis.frequency(2:4)),zeros(1,2))
            testCase.verifyGreaterThan(diff(basis.frequency(5:6)),0)
            testCase.verifyGreaterThan(basis.frequency(2:4),abs(background.f0)*ones(1,3))
            testCase.verifyLessThan(basis.frequency(5:6),abs(background.f0)*ones(1,2))
            z = linspace(-1000,0,137).';
            testCase.verifyEqual(basis.F(z),background.D*basis.Gz(z),RelTol=2e-11,AbsTol=2e-11)
            for j = 1:basis.nColumns
                norm = integral(@(z) waveNormIntegrand(basis,z,j,k,background.D),-1000,0,AbsTol=1e-8,RelTol=1e-10)/background.D;
                testCase.verifyEqual(norm,1,AbsTol=2e-8)
            end
        end

        function zeroWaveCountsReturnNoImplicitSpecialModes(testCase)
            free = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            rigid = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            surface = IMSolverAnalytical().solveEVP(IMInternalModes.waveModesAtWavenumber(free,k=2e-4),nModes=0);
            empty = IMSolverAnalytical().solveEVP(IMInternalModes.waveModesAtWavenumber(rigid,k=2e-4),nModes=0);
            testCase.verifyEqual(surface.nColumns,0)
            testCase.verifyEmpty(surface.modeInfo.branch)
            explicit = IMSolverAnalytical().solveEVP(IMInternalModes.waveModesAtWavenumber(free,k=2e-4),nModes=1,branch="surface");
            testCase.verifyEqual(explicit.nColumns,1)
            testCase.verifyEqual(explicit.modeInfo.branch,"surface")
            testCase.verifyEqual(empty.nColumns,0)
            testCase.verifySize(empty.G([-1000;0]),[2 0])
        end

        function analyticalAndNumericalWavesShareIdentityAndConvention(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],latitude=-30,surfaceCondition="freeSurface");
            k = 2e-4;
            exact = IMSolverAnalytical().solveEVP(IMInternalModes.waveModesAtWavenumber(background,k=k),nModes=5);
            numerical = IMSolverSpectral(nEVP=64,coordinateKind="z").solveEVP(IMInternalModes.waveModesAtWavenumber(background,k=k),nModes=5);
            testCase.verifyEqual(numerical.columnLabels,exact.columnLabels)
            testCase.verifyEqual(numerical.modeInfo.branch,exact.modeInfo.branch)
            testCase.verifyEqual(numerical.modeInfo.branchIndex,exact.modeInfo.branchIndex)
            testCase.verifyEqual(numerical.h,exact.h,RelTol=2e-7)
            z = linspace(-1000,0,101).';
            testCase.verifyEqual(exact.F(z),1000*exact.Gz(z),AbsTol=2e-12)
            testCase.verifyEqual(numerical.F(z),exact.F(z),RelTol=2e-6,AbsTol=2e-7)
            testCase.verifyEqual(numerical.G(z),exact.G(z),RelTol=2e-6,AbsTol=2e-7)
        end

        function waveNormalizationsShareNamesAndScalePairedProfiles(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            k = 2e-4;
            z = linspace(-1000,0,20001).';
            names = ["rmsF","rmsG","maxAbsF","maxAbsG"];
            for name = names
                basis = IMSolverAnalytical().solveEVP(IMInternalModes.waveModesAtWavenumber(background,k=k),nModes=3,normalization=name);
                testCase.verifyEqual(basis.normalization,name)
                F = basis.F(z);
                G = basis.G(z);
                testCase.verifyEqual(F,1000*basis.Gz(z),RelTol=2e-12,AbsTol=2e-12)
                switch name
                    case "rmsF"
                        norms = trapz(z,abs(F).^2)/1000;
                    case "rmsG"
                        norms = trapz(z,abs(G).^2)/1000;
                    case "maxAbsF"
                        norms = max(abs(F),[],1);
                    case "maxAbsG"
                        norms = max(abs(G),[],1);
                end
                testCase.verifyEqual(norms,ones(1,basis.nColumns),AbsTol=2e-7)
            end
        end

        function apvCountsFiniteNegativeAndZeroSectorsWithinTotal(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            negativeBoundary = IMGeneralizedEnergyBoundary(g0=-.04,gd=0);
            zeroBoundary = IMGeneralizedEnergyBoundary(g0=Inf,gd=Inf);
            for solver = {IMSolverAnalytical(),IMSolverSpectral(nEVP=64,coordinateKind="z")}
                negative = solver{1}.solveEVP(IMInternalModes.geostrophicAPVModes(background,g0=negativeBoundary.g0,gd=negativeBoundary.gd),nModes=1,branch="negative");
                zero = solver{1}.solveEVP(IMInternalModes.geostrophicAPVModes(background,g0=zeroBoundary.g0,gd=zeroBoundary.gd),nModes=1,branch="zeroEigenvalue");
                testCase.verifyEqual(negative.nColumns,1)
                testCase.verifyLessThan(negative.eigenvalues,0)
                testCase.verifyEqual(negative.modeInfo.branch,"negative")
                testCase.verifyEqual(zero.nColumns,1)
                testCase.verifyEqual(zero.eigenvalues,0)
                testCase.verifyEqual(zero.h,Inf)
                testCase.verifyEqual(zero.modeInfo.branch,"zeroEigenvalue")
                positive = solver{1}.solveEVP(IMInternalModes.geostrophicAPVModes(background,g0=negativeBoundary.g0,gd=negativeBoundary.gd),nModes=4);
                testCase.verifyEqual(positive.nColumns,4)
                testCase.verifyEqual(nnz(positive.eigenvalues > 0),3)
                testCase.verifyEqual(positive.normalization,"rmsF")
                z = linspace(-1000,0,10001).';
                testCase.verifyEqual(trapz(z,abs(positive.F(z)).^2)/1000,ones(1,4),AbsTol=1e-7)
            end
        end

        function meanDensityConstantConsumesRequestedInventory(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            energyBoundary = IMGeneralizedEnergyBoundary(g0=.02,gd=-.01);
            solver = IMSolverSpectral(nEVP=64,coordinateKind="z");
            basis = solver.solveEVP(IMInternalModes.meanDensityAnomalyModes(background,g0=energyBoundary.g0,gd=energyBoundary.gd),branch=["negative","constant"],nModes=[1 1],normalization="meanDensityAnomalyEnergy");
            testCase.verifyGreaterThanOrEqual(basis.nColumns,1)
            zero = basis.eigenvalues == 0;
            testCase.verifyEqual(nnz(zero),1)
            testCase.verifyEqual(nnz(basis.eigenvalues > 0),0)
            testCase.verifyEqual(basis.normalization,"meanDensityAnomalyEnergy")
            z = linspace(-1000,0,101).';
            G = basis.G(z);
            testCase.verifyEqual(G(:,zero),repmat(G(1,zero),numel(z),1),AbsTol=2e-10)
            regular = solver.solveEVP(IMInternalModes.meanDensityAnomalyModes(background,g0=energyBoundary.g0,gd=energyBoundary.gd),nModes=5,normalization="meanDensityAnomalyEnergy");
            testCase.verifyEqual(nnz(regular.eigenvalues > 0),3)
            testCase.verifyEqual(nnz(regular.eigenvalues <= 0),basis.nColumns)
        end

        function zeroAPVReturnsExplicitlyRequestedEndpoints(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            solver = IMSolverSpectral(nEVP=64,coordinateKind="z");
            both = solver.solveGeostrophicZeroAPVModes(IMGeostrophicZeroAPVModes.atWavenumber(background,2e-4,endpoints=["surface","bottom"]));
            surface = solver.solveGeostrophicZeroAPVModes(IMGeostrophicZeroAPVModes.atWavenumber(background,2e-4,endpoints="surface"));
            bottom = solver.solveGeostrophicZeroAPVModes(IMGeostrophicZeroAPVModes.atWavenumber(background,2e-4,endpoints="bottom"));
            testCase.verifyEqual(both.nColumns,2)
            testCase.verifyEqual(surface.nColumns,1)
            testCase.verifyEqual(bottom.nColumns,1)
            testCase.verifyEqual(both.normalization,"endpointResponse")
            testCase.verifyEqual(both.columnLabels,[surface.columnLabels bottom.columnLabels])
            z = linspace(-1000,0,37).';
            testCase.verifyEqual(both.F(z),[surface.F(z) bottom.F(z)],RelTol=2e-10,AbsTol=2e-10)
        end

        function neutralLabelsKeepFixedEndpointsAndExactCounts(testCase)
            background = IMBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0],surfaceCondition="freeSurface");
            z = linspace(-1000,0,101).';
            empty = IMNeutralLabelModesBasis(background,neutralIntervals=[-1000 0],nModesPerInterval=0);
            testCase.verifyEqual(empty.nColumns,0)
            testCase.verifyEqual(empty.normalization,"rmsG")
            testCase.verifySize(empty.G(z),[numel(z) 0])
            testCase.verifyEmpty(empty.columnLabels)
            for acceleration = [0 .1 -.1 Inf]
                energyBoundary = IMGeneralizedEnergyBoundary(g0=acceleration,gd=acceleration);
                basis = IMNeutralLabelModesBasis(background,neutralIntervals=[-1000 0],g0=energyBoundary.g0,gd=energyBoundary.gd,nModesPerInterval=3);
                testCase.verifyEqual(basis.nColumns,3)
                testCase.verifyEqual(basis.G([-1000;0]),zeros(2,3))
                testCase.verifyEqual(basis.modeInfo.branch,repmat("neutralLabel",3,1))
                testCase.verifyEqual(basis.modeInfo.localModeNumber,(1:3).')
                testCase.verifyEqual(basis.G(z),sqrt(2)*sin((z+1000)*(pi*(1:3)/1000)),AbsTol=3e-15)
                weights = 10*ones(size(z)); weights([1 end]) = 5;
                projection = basis.projectionOnGrid(z,weights,variable="APV");
                testCase.verifyEqual(projection.projectPhysical(basis.APV(z)),eye(3),AbsTol=2e-14)
                testCase.verifyEqual(projection.provenance.g0,acceleration)
                testCase.verifyEqual(projection.provenance.gd,acceleration)
            end
        end

        function hydrostaticFactoriesCountExternalModeWithinTotal(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            for solver = {IMSolverAnalytical(),IMSolverSpectral(nEVP=64,coordinateKind="z")}
                F = solver{1}.solveEVP(IMInternalModes.hydrostaticFModes(background),nModes=3);
                G = solver{1}.solveEVP(IMInternalModes.hydrostaticGModes(background),nModes=3);
                testCase.verifyEqual(F.nColumns,3)
                testCase.verifyEqual(G.nColumns,3)
                testCase.verifyEqual(F.h,G.h,RelTol=2e-7)
                testCase.verifyEqual(solver{1}.solveEVP(IMInternalModes.hydrostaticGModes(background),nModes=0).nColumns,0)
            end
        end

        function inertialModesExposeAuxiliaryLadders(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            basis = IMSolverAnalytical().solveEVP(IMInternalModes.waveModesAtWavenumber(background,k=0),nModes=4,normalization="rmsF");
            testCase.verifyEqual(basis.nColumns,4)
            testCase.verifyEqual(basis.normalization,"rmsF")
            testCase.verifyEqual(basis.frequency,abs(background.f0)*ones(1,4),AbsTol=2e-15)
            z = linspace(-1000,0,101).';
            testCase.verifyEqual(basis.F(z),1000*basis.Gz(z),AbsTol=2e-12)
            testCase.verifyEqual(basis.evaluate(z,variable="F",derivativeOrder=1),1000*basis.Gzz(z),AbsTol=2e-12)
        end

        function generalizedEnstrophyFactoryReturnsRequestedLadder(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            basis = IMSolverSpectral(nEVP=64,coordinateKind="z").solveEVP(IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes(background,k=2e-4),nModes=3);
            testCase.verifyEqual(basis.nColumns,3)
            testCase.verifyEqual(basis.normalization,"generalizedPotentialEnstrophy")
            testCase.verifyGreaterThan(basis.eigenvalues,zeros(1,3))
        end

        function fixedFrequencyWavesUseTheSameProfileConvention(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            omega = .002;
            basis = IMSolverAnalytical().solveEVP(IMInternalModes.waveModesAtFrequency(background,omega=omega),nModes=4);
            testCase.verifyEqual(basis.nColumns,4)
            testCase.verifyEqual(basis.normalization,"rmsF")
            testCase.verifyEqual(basis.frequency,omega*ones(1,4),AbsTol=2e-15)
            z = linspace(-1000,0,10001).';
            F = basis.F(z);
            testCase.verifyEqual(F,1000*basis.Gz(z),AbsTol=2e-12)
            testCase.verifyEqual(trapz(z,abs(F).^2)/1000,ones(1,4),AbsTol=1e-7)
        end

        function exactZeroDepthWavesHaveFiniteRequestedAuxiliaryBasis(testCase)
            f0 = 2*7.2921e-5*sind(33);
            background = IMBackground(N2=@(z) f0^2*ones(size(z)),zDomain=[-1000 0],surfaceCondition="freeSurface");
            basis = IMZeroDepthWaveModesBasis(background,2e-4,criticalIntervals=[-1000 0],nModesPerInterval=2);
            testCase.verifyEqual(basis.nColumns,2)
            testCase.verifyEqual(basis.h,zeros(1,2))
            testCase.verifyEqual(basis.frequency,abs(f0)*ones(1,2),AbsTol=2e-15)
            z = linspace(-1000,0,101).';
            testCase.verifyEqual(basis.F(z),1000*basis.Gz(z),AbsTol=2e-12)
            testCase.verifyEqual(basis.evaluate(z,variable="F",derivativeOrder=1),1000*basis.Gzz(z),AbsTol=2e-12)
            empty = IMZeroDepthWaveModesBasis(background,2e-4,criticalIntervals=[-1000 0],nModesPerInterval=0);
            testCase.verifyEqual(empty.nColumns,0)
            testCase.verifyEqual(empty.h,zeros(1,0))
        end

        function numericalRequestsAreExactAndRejectUnavailableSector(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            small = IMSolverSpectral(nEVP=16,coordinateKind="z");
            testCase.verifyError(@() small.solveEVP(IMInternalModes.waveModesAtWavenumber(background,k=2e-4),nModes=30),"IMSolver:InsufficientModes")
            testCase.verifyError(@() IMSolverAnalytical().solveEVP(IMInternalModes.waveModesAtWavenumber(background,k=2e-4),nModes=1,branch="subinertial"),"IMSolver:InsufficientModes")
            testCase.verifyError(@() small.solveEVP(IMInternalModes.waveModesAtWavenumber(background,k=2e-4),nModes=[1 2]),"IMSolver:InvalidModeRequest")
            testCase.verifyError(@() small.solveEVP(IMInternalModes.waveModesAtWavenumber(background,k=2e-4),nModes=2,branch=["surface","internal"]),"IMSolver:InvalidModeRequest")
        end
    end
end

function values = waveNormIntegrand(basis,z,column,k,D)
F = basis.F(z(:));
G = basis.G(z(:));
values = reshape(abs(F(:,column)).^2+(k*D)^2*abs(G(:,column)).^2,size(z));
end
