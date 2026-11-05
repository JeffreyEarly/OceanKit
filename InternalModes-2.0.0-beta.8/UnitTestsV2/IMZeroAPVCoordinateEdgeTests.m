classdef IMZeroAPVCoordinateEdgeTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function configurePath(testCase)
            root = fileparts(fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root,"UnitTestsV2","TestSupport")));
        end
    end

    methods (Test)
        function originallyEmptyFamiliesRemainEmptyWithoutRestoringDirections(testCase)
            background=IMBackground.constantStratification(N0=1,zDomain=[-1 0],rotationRate=1/(2*sind(33)),latitude=33,g=2,surfaceCondition="freeSurface");
            problem=IMGeostrophicZeroAPVModes.atWavenumber(background,[1.2 .7],endpoints=strings(1,0));
            for solver={IMSolverAnalytical(),IMSolverSpectral(nEVP=32),IMSolverFiniteDifference(z=linspace(-1,0,33).')}
                b=solver{1}.solveGeostrophicZeroAPVModes(problem);
                o=b.orthogonalizeGeneralizedEnergy(g0=Inf,gd=Inf);
                testCase.verifySize(o.F([-1;-.3;0]),[3 0 2]);
                testCase.verifySize(o.Gzz([-1;-.3;0]),[3 0 2]);
                testCase.verifySize(o.boundaryResponseMatrix,[0 0 2]);
                testCase.verifySize(o.energyMatrix,[0 0 2]);
                testCase.verifySize(o.boundaryDepth,[0 2]);
                testCase.verifyEqual(o.nColumns,0);
            end
        end

        function infinitiesCannotDeactivateRetainedEndpoints(testCase)
            background=IMBackground.constantStratification(N0=1,zDomain=[-1 0],rotationRate=1/(2*sind(33)),latitude=33,g=2,surfaceCondition="freeSurface");
            for solver={IMSolverAnalytical(),IMSolverSpectral(nEVP=48)}
                full=solver{1}.solveGeostrophicZeroAPVModes(IMGeostrophicZeroAPVModes.atWavenumber(background,1.2));
                prefix=string(class(full));
                testCase.verifyError(@() full.orthogonalizeGeneralizedEnergy(g0=Inf,gd=0),prefix+":InactiveEndpointRequired");
                testCase.verifyError(@() full.generalizedEnergyMatrix(g0=NaN,gd=0),prefix+":InvalidEnergyCoefficient");
                testCase.verifyError(@() full.generalizedEnergyMatrix(g0=0,gd=-Inf),prefix+":InvalidEnergyCoefficient");
                bottom=solver{1}.solveGeostrophicZeroAPVModes(IMGeostrophicZeroAPVModes.atWavenumber(background,1.2,endpoints="bottom"));
                o=bottom.orthogonalizeGeneralizedEnergy(g0=Inf,gd=.4);
                testCase.verifyTrue(all(isfinite(o.boundaryDepth)));
                testCase.verifyEqual(bottom.generalizedEnergyMatrix(g0=Inf,gd=.4),bottom.generalizedEnergyMatrix(g0=0,gd=.4),AbsTol=0);
            end
        end

        function schurTunedNullEnergyRemainsAStateWithoutAUnitEnergyNorm(testCase)
            background=IMBackground.constantStratification(N0=1,zDomain=[-1 0],rotationRate=1/(2*sind(33)),latitude=33,g=2,surfaceCondition="freeSurface");
            b=IMSolverAnalytical().solveGeostrophicZeroAPVModes(IMGeostrophicZeroAPVModes.atWavenumber(background,1.2));
            H=b.energyMatrix; B0=b.surfaceBuoyancyMatrix; Bd=b.bottomBuoyancyMatrix;
            K=H+.4*Bd; g0=-(K(1,1)-K(1,2)^2/K(2,2))/B0(1,1);
            o=b.orthogonalizeGeneralizedEnergy(g0=g0,gd=.4);
            [~,index]=min(abs(o.boundaryDepth)); nullState=o.selectColumns(index);
            testCase.verifyLessThan(abs(nullState.boundaryDepth),1e-12);
            testCase.verifyGreaterThan(norm(nullState.F([-1;-.4;0])),.01);
            testCase.verifyEqual(nullState.normalization,"endpointResponse");
            testCase.verifyEqual(nullState.normSignatures,1);
            testCase.verifyFalse(ismember("generalizedEnergy",nullState.normalizationNames()));
            testCase.verifyError(@() nullState.withNormalization("generalizedEnergy"),"IMAnalyticalGeostrophicZeroAPVModesBasis:UnsupportedNormalization");
            testCase.verifyTrue(all(isfinite(nullState.F([-1;-.4;0],normalization="maxAbsF"))));
        end

        function rigidNeutralNullStateCannotReceiveUnitEnergyFromRoundoff(testCase)
            background=IMBackground(N2=@(z)zeros(size(z)),zDomain=[-1 0],rotationRate=1/(2*sind(33)),latitude=33,g=2,surfaceCondition="rigidLid");
            b=IMSolverSpectral(nEVP=48).solveGeostrophicZeroAPVModes(IMGeostrophicZeroAPVModes.atWavenumber(background,.3));
            o=b.orthogonalizeGeneralizedEnergy(g0=0,gd=0);
            [~,index]=min(abs(o.boundaryDepth)); n=o.selectColumns(index);
            z=linspace(-1,0,17).';
            testCase.verifyLessThan(norm(n.F(z)),1e-11);
            testCase.verifyGreaterThan(norm(n.G(z)),1);
            testCase.verifyEqual(n.normalization,"endpointResponse");
            testCase.verifyEqual(n.normSignatures,1);
            testCase.verifyTrue(all(isfinite(n.G(z,normalization="rmsG"))));
            testCase.verifyError(@() n.withNormalization("generalizedEnergy"),"IMGeostrophicZeroAPVModesBasis:UnsupportedNormalization");
        end

        function nullEnergyPagesRetainFiniteEndpointCoordinates(testCase)
            background=IMBackground.constantStratification(N0=1,zDomain=[-1 0],rotationRate=1/(2*sind(33)),latitude=33,g=2,surfaceCondition="freeSurface");
            b=IMSolverAnalytical().solveGeostrophicZeroAPVModes(IMGeostrophicZeroAPVModes.atWavenumber(background,[1.2 .7]));
            H=b.energyMatrix; B0=b.surfaceBuoyancyMatrix; Bd=b.bottomBuoyancyMatrix;
            K=H(:,:,1)+.4*Bd(:,:,1); g0=-(K(1,1)-K(1,2)^2/K(2,2))/B0(1,1,1);
            o=b.orthogonalizeGeneralizedEnergy(g0=g0,gd=.4);
            values=o.evaluate([-1;-.3;0],variable="G",normalization="endpointResponse",pages=2);
            testCase.verifySize(values,[3 2]);
            testCase.verifyTrue(all(isfinite(o.evaluate([-1;-.3;0],variable="G",normalization="endpointResponse",pages=1)),"all"));
            testCase.verifyEqual(o.normSignatures,ones(2,2));
        end
    end
end
