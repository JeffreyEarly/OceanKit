classdef IMModeBasisTests < matlab.unittest.TestCase
    methods (TestClassSetup)
        function configurePath(testCase)
            originalPath = path;
            addpath(fileparts(fileparts(mfilename("fullpath"))));
            testCase.addTeardown(@() path(originalPath));
        end
    end

    methods (Test)
        function numericalAndAnalyticalWavesShareTheDepthConvention(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            k = 1e-3;
            evp = IMInternalModes.waveModesAtWavenumber(N2=background.N2,zDomain=background.zDomain,k=k,rotationRate=background.rotationRate,latitude=background.latitude,surfaceBoundary=IMBoundaryCondition.dirichlet());
            numerical = IMSolverSpectral(nEVP=64,coordinateKind="z").solveEVP(evp,nModes=4);
            solution = background.analyticalSolution();
            analytical = solution.internalModes(evp,nModes=4);
            numericBasis = testCase.wrap(background,numerical,"waveEnergy",true,k);
            exactBasis = testCase.wrap(background,analytical,"waveEnergy",true,k);
            z = linspace(-1000,0,131).';
            for current = {numericBasis,exactBasis}
                basis = current{1};
                testCase.verifyEqual(basis.F(z),1000*basis.Gz(z),AbsTol=1e-13)
                testCase.verifyEqual(basis.evaluate(z,variable="F",derivativeOrder=1),1000*basis.Gzz(z),AbsTol=1e-12)
                for j = 1:4
                    value = integral(@(zz) waveEnergyDensity(basis,zz,j,k),-1000,0,AbsTol=1e-10)/1000;
                    testCase.verifyEqual(value,1,AbsTol=2e-10)
                end
            end
            testCase.verifyEqual(numericBasis.F(z),exactBasis.F(z),AbsTol=2e-9)
            testCase.verifyEqual(numericBasis.G(z),exactBasis.G(z),AbsTol=2e-9)
            testCase.verifyEqual(numericBasis.h,exactBasis.h,RelTol=2e-9)
        end

        function normalizationsScaleBothWaveProfilesTogether(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            evp = IMInternalModes.waveModesAtWavenumber(N2=background.N2,zDomain=background.zDomain,k=1e-3,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solution = background.analyticalSolution();
            source = solution.internalModes(evp,nModes=3);
            rmsF = testCase.wrap(background,source,"rmsF",true,1e-3);
            rmsG = testCase.wrap(background,source,"rmsG",true,1e-3);
            z = linspace(-1000,0,201).';
            testCase.verifyEqual(trapz(z,rmsF.F(z).^2)/1000,ones(1,3),AbsTol=2e-13)
            testCase.verifyEqual(trapz(z,rmsG.G(z).^2)/1000,ones(1,3),AbsTol=2e-13)
            testCase.verifyEqual(rmsF.F(z).*rmsG.G(z),rmsG.F(z).*rmsF.G(z),AbsTol=1e-14)
            % Changing the active normalization of a value copy preserves the source.
            expected = rmsF.F(z);
            source.normalization = "maxAbsF";
            source.F(z);
            testCase.verifyEqual(rmsF.F(z),expected)
            testCase.verifyEqual(rmsF.normalization,"rmsF")
        end

        function analyticalSecondDerivativesRespectRotatingWaveEquation(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            k = 1e-3;
            evp = IMInternalModes.waveModesAtWavenumber(N2=background.N2,zDomain=background.zDomain,k=k,rotationRate=background.rotationRate,latitude=background.latitude,surfaceBoundary=IMBoundaryCondition.dirichlet());
            solution = background.analyticalSolution();
            source = solution.internalModes(evp,nModes=4);
            basis = testCase.wrap(background,source,"waveEnergy",true,k);
            z = linspace(-1000,0,131).';
            expected = (k^2-(background.N2(z)-background.f0^2)/background.g.*basis.eigenvalues).*basis.G(z);
            testCase.verifyEqual(basis.Gzz(z),expected,AbsTol=1e-16)
            testCase.verifyEqual(basis.evaluate(z,variable="G",derivativeOrder=1),basis.Gz(z))
            testCase.verifyEqual(basis.evaluate(z,variable="G",derivativeOrder=2),basis.Gzz(z))
        end

        function signedMDAEnergyRetainsItsMetricSignature(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            g0 = -.04;
            gd = .03;
            evp = IMInternalModes.meanDensityAnomalyModes(N2=background.N2,zDomain=background.zDomain,g0=g0,gd=gd);
            source = IMSolverSpectral(nEVP=96,coordinateKind="z").solveEVP(evp,nModes=5);
            basis = testCase.wrap(background,source,"meanDensityAnomalyEnergy",false,0);
            [z,weights] = source.solver.nativeQuadratureRule(background.zDomain);
            G = basis.G(z);
            endpoints = basis.G(background.zDomain.');
            gram = G'*(weights.*(background.N2(z)/background.g).*G)+(g0/background.g)*(endpoints(2,:)'*endpoints(2,:))+(gd/background.g)*(endpoints(1,:)'*endpoints(1,:));
            testCase.verifyEqual(gram,diag(basis.normSignatures),AbsTol=2e-9)
            testCase.verifyTrue(any(basis.normSignatures < 0))
            testCase.verifyEqual(basis.F(0),zeros(1,5),AbsTol=1e-13)
            testCase.verifyEqual(basis.Gz(z),basis.evaluate(z,variable="G",derivativeOrder=1))
        end

        function generalizedEnstrophyUsesTheSameNamedPairing(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            evp = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes(N2=background.N2,zDomain=background.zDomain,k=1e-3,rotationRate=background.rotationRate,latitude=background.latitude);
            source = IMSolverSpectral(nEVP=64,coordinateKind="z").solveEVP(evp,nModes=4);
            basis = testCase.wrap(background,source,"generalizedPotentialEnstrophy",false,1e-3);
            z = linspace(-1000,0,31).';
            testCase.verifyEqual(basis.F(z),source.F(z,normalization="generalizedPotentialEnstrophy"),AbsTol=2e-13,RelTol=5e-14)
            testCase.verifyEqual(basis.G(z),source.G(z,normalization="generalizedPotentialEnstrophy"),AbsTol=2e-13,RelTol=5e-14)
            % Native companion derivatives are differentiated independently.
            % Use the per-column residual criterion of the geostrophic
            % equation tests rather than requiring algebraic equality.
            Gzz = basis.Gzz(z);
            expectedGzz = -background.N2(z)/background.g.*basis.G(z)./basis.h;
            scale = max(max(abs(Gzz),[],1),max(abs(expectedGzz),[],1));
            testCase.verifyLessThanOrEqual(max(abs(Gzz-expectedGzz),[],1)./max(scale,realmin),1e-7)
        end

        function neutralLabelsExposeOnlyTheirScientificChannels(testCase)
            background = IMBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0],surfaceCondition="freeSurface");
            source = IMNeutralLabelModesBasis([-1000 0],zDomain=background.zDomain,nModesPerInterval=2);
            basis = testCase.wrap(background,source,"rmsG",false,0);
            z = [-1000;-700;-100;0];
            testCase.verifyEqual(basis.nColumns,2)
            testCase.verifyEqual(basis.G(z),source.G(z))
            testCase.verifyEqual(basis.Gz(z),source.Gz(z))
            testCase.verifyEqual(basis.Gzz(z),source.Gzz(z))
            testCase.verifyEqual(basis.APV(z),source.APV(z))
            testCase.verifyFalse(ismember("F",basis.capabilities.variables))
            testCase.verifyFalse(ismethod(basis,"F"))
            testCase.verifyFalse(isprop(basis,"h"))
            testCase.verifyFalse(isprop(basis,"eigenvalues"))
            testCase.verifyError(@() basis.G(z,normalization="rmsF"),"IMNeutralLabelModesBasis:UnsupportedNormalization")
        end

        function scaledSelectedNeutralProjectionRecoversItsOwnCoefficients(testCase)
            background = IMBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0],surfaceCondition="freeSurface");
            source = IMNeutralLabelModesBasis([-800 -600],zDomain=background.zDomain,nModesPerInterval=3);
            basis = source.selectColumns([3 1]);
            basis.normalization = "maxAbsG";
            z = linspace(-800,-600,4001).';
            weights = .05*ones(size(z));
            weights([1 end]) = .025;
            projection = basis.projectionOnGrid(z,weights,variable="APV");
            coefficients = [2.5;-.75];
            testCase.verifyEqual(projection.projectPhysical(basis.APV(z)*coefficients),coefficients,AbsTol=2e-12)
            testCase.verifyEqual(projection.columnLabels,basis.columnLabels)
            originalG = source.G(z);
            testCase.verifyEqual(basis.G(z),originalG(:,[3 1])./sqrt(10),AbsTol=1e-14)
        end

        function fullyNeutralLabelsAllHaveObservableAPV(testCase)
            background = IMBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0],surfaceCondition="freeSurface");
            source = IMNeutralLabelModesBasis([-1000 0],zDomain=background.zDomain,nModesPerInterval=2);
            basis = testCase.wrap(background,source,"rmsG",false,0);
            z = linspace(-1000,0,1001).';
            weights = ones(size(z));
            weights([1 end]) = .5;
            projection = basis.projectionOnGrid(z,weights,variable="APV");
            testCase.verifyEqual(projection.physicalDualMatrix*basis.APV(z),eye(2),AbsTol=2e-12)
            selected = basis.projectionOnGrid(z,weights,variable="APV",columns=[2 1]);
            testCase.verifyEqual(selected.physicalDualMatrix*basis.APV(z),[0 1;1 0],AbsTol=2e-12)
        end

        function endpointResponsesKeepCanonicalScaling(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
            k = 1e-4;
            problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=background.N2,zDomain=background.zDomain,k=k,rotationRate=background.rotationRate,latitude=background.latitude,surfaceBoundary=background.surfaceCondition,endpoints=["surface","bottom"]);
            source = IMSolverSpectral(nEVP=64,coordinateKind="z").solveGeostrophicZeroAPVModes(problem);
            basis = source;
            z = linspace(-1000,0,31).';
            testCase.verifyEqual(basis.F(z),source.F(z))
            testCase.verifyEqual(basis.G(z),source.G(z))
            endpointF = basis.F(background.zDomain.');
            endpointG = basis.G(background.zDomain.');
            testCase.verifyEqual([endpointG(2,:)-endpointF(2,:);endpointG(1,:)],eye(2),AbsTol=2e-11)
            testCase.verifyEqual(basis.Gz(z),-background.g*k^2/background.f0^2*basis.F(z),AbsTol=1e-14)
        end

        function intrinsicAndEmptyBasesHaveTheSameShapeContract(testCase)
            rotationRate = 7.2921e-5;
            latitude = 33;
            f0 = 2*rotationRate*sind(latitude);
            background = IMBackground(N2=@(z) f0^2*ones(size(z)),zDomain=[-1000 0],rotationRate=rotationRate,latitude=latitude,surfaceCondition="freeSurface");
            basis = IMZeroDepthWaveModesBasis(background,1e-3,criticalIntervals=[-1000 0],nModesPerInterval=1);
            z = linspace(-1000,0,31).';
            testCase.verifyEqual(basis.F(z),1000*basis.Gz(z),AbsTol=1e-14)
            testCase.verifyEqual(basis.h,0)
            emptyBasis = IMZeroDepthWaveModesBasis(background,1e-3,criticalIntervals=[-1000 0],nModesPerInterval=0);
            testCase.verifySize(emptyBasis.G(z),[31 0])
            testCase.verifySize(emptyBasis.F(z),[31 0])
            testCase.verifyEqual(emptyBasis.nColumns,0)
            nullBasis = basis.addNormalization("null",@(~,~) 0);
            testCase.verifyError(@() nullBasis.F(z,normalization="null"),"IMZeroDepthWaveModesBasis:NullNormalization")
        end

        function maximumNormalizationResolvesInteriorExtrema(testCase)
            background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="rigidLid");
            evp = IMInternalModes.waveModesAtWavenumber(N2=background.N2,zDomain=background.zDomain,k=1e-3,surfaceBoundary=IMBoundaryCondition.dirichlet());
            numerical = IMSolverSpectral(nEVP=64,coordinateKind="z").solveEVP(evp,nModes=3);
            solution = background.analyticalSolution();
            analytical = solution.internalModes(evp,nModes=3);
            for source = {numerical,analytical}
                maxF = testCase.wrap(background,source{1},"maxAbsF",true,1e-3);
                maxG = testCase.wrap(background,source{1},"maxAbsG",true,1e-3);
                testCase.verifyEqual(abs(maxF.F(0)),ones(1,3),AbsTol=2e-12)
                for j = 1:3
                    values = maxG.G(-1000+1000/(2*j));
                    testCase.verifyEqual(abs(values(j)),1,AbsTol=2e-12)
                end
            end
        end

        function intrinsicNormalizationSplitsExactSupportIntervals(testCase)
            rotationRate = 7.2921e-5;
            latitude = 33;
            f0 = 2*rotationRate*sind(latitude);
            background = IMBackground(N2=@(z) f0^2*ones(size(z)),zDomain=[-1000 0],rotationRate=rotationRate,latitude=latitude,surfaceCondition="freeSurface");
            intervals = [-812 -631;-111 -34];
            basis = IMZeroDepthWaveModesBasis(background,1e-3,criticalIntervals=intervals,nModesPerInterval=[1 1],normalization="rmsG");
            peaks = mean(intervals,2);
            expected = diag(sqrt(2000./diff(intervals,1,2)));
            testCase.verifyEqual(basis.G(peaks),expected,AbsTol=2e-13)
            testCase.verifyEqual(basis.F(peaks),1000*basis.Gz(peaks),AbsTol=1e-13)
        end
    end

    methods (Static, Access = private)
        function basis = wrap(~,source,normalization,~,~)
            basis = source;
            basis.normalization = normalization;
        end
    end
end

function values = waveEnergyDensity(basis,z,j,k)
F = basis.F(z(:));
G = basis.G(z(:));
values = reshape(F(:,j).^2+(k*diff(basis.zDomain))^2*G(:,j).^2,size(z));
end
