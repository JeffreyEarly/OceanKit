classdef IMNeutralGeostrophicModesTests < matlab.unittest.TestCase
    % Check regular balanced pairs independently of diagnostic ODE recovery.

    properties (Access = private)
        originalPath
    end

    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.originalPath = path;
            addpath(repoRoot);
            addpath(fullfile(repoRoot,"UnitTestsV2","TestSupport"));
        end
    end

    methods (TestClassTeardown)
        function restoreRepositoryPath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function positiveProfilesRetainAnalyticalAPVModes(testCase)
            N0 = 5e-3;
            zDomain = [-1000 0];
            solutions = {IMConstantStratificationSolution(N0=N0,zDomain=zDomain), ...
                IMExponentialStratificationSolution(N0=N0,b=700,zDomain=zDomain)};
            z = linspace(zDomain(1),zDomain(2),257).';
            for iProfile = 1:numel(solutions)
                solution = solutions{iProfile};
                evp = IMInternalModes.geostrophicAPVModes(N2=@(z) solution.N2(z),zDomain=zDomain,g0=-0.02,gd=Inf,surfaceBoundary="freeSurface");
                numerical = IMSolverSpectral(nEVP=128,coordinateKind="z").solveEVP(evp,nModes=6);
                analytical = solution.internalModes(evp,nModes=6);
                testCase.verifyClass(numerical,"IMGeostrophicModesBasis")
                testCase.verifySize(numerical.nativeModes,[128 6])
                testCase.verifyEqual(numerical.modeNumber,analytical.modeNumber)
                testCase.verifyEqual(numerical.h,analytical.h,RelTol=1e-7)
                numericalScale = testCase.independentDepthScales(numerical);
                analyticalScale = testCase.independentDepthScales(analytical);
                testCase.verifyEqual(numerical.F(z)./numericalScale,analytical.F(z)./analyticalScale,RelTol=1e-7,AbsTol=1e-7)
                testCase.verifyEqual(numerical.G(z)./numericalScale,analytical.G(z)./analyticalScale,RelTol=1e-7,AbsTol=1e-7)
                testCase.verifyEqual(numerical.h,1./numerical.inverseDepth,AbsTol=0)
                testCase.verifyPairEquations(numerical,1e-7)
            end
        end

        function weakConstantStratificationRetainsDimensionallyBalancedPairs(testCase)
            N0 = 1e-8;
            zDomain = [-1000 0];
            solution = IMConstantStratificationSolution(N0=N0,zDomain=zDomain);
            z = linspace(-1000,0,257).';
            endpoints = [0 0;0.02 0.03];
            for iEndpoint = 1:size(endpoints,1)
                evp = IMInternalModes.geostrophicAPVModes(N2=@(z) N0*N0*ones(size(z)), ...
                    zDomain=zDomain,g0=endpoints(iEndpoint,1),gd=endpoints(iEndpoint,2),surfaceBoundary="freeSurface");
                basis = IMSolverSpectral(nEVP=64,coordinateKind="z").solveEVP(evp,nModes=4);
                analytical = solution.internalModes(evp,nModes=4);
                numericalScale = testCase.independentDepthScales(basis);
                analyticalScale = testCase.independentDepthScales(analytical);
                numericalF = basis.F(z)./numericalScale;
                numericalG = basis.G(z)./numericalScale;
                exactF = analytical.F(z)./analyticalScale;
                exactG = analytical.G(z)./analyticalScale;
                testCase.verifyNumElements(basis.h,4)
                testCase.verifyFalse(any(basis.inverseDepth==0))
                testCase.verifyEqual(basis.modeNumber,analytical.modeNumber)
                testCase.verifyEqual(basis.h,analytical.h,RelTol=1e-7)
                testCase.verifyLessThan(max(vecnorm(numericalF-exactF)./vecnorm(exactF)),1e-7)
                testCase.verifyLessThan(max(vecnorm(numericalG-exactG)./vecnorm(exactG)),1e-7)
                testCase.verifyPairEquations(basis,1e-7)
            end
        end

        function positiveEnstrophyModesMatchIndependentScalarPencil(testCase)
            N2 = @(z) (5e-3)^2*exp(2*z/700);
            zDomain = [-1000 0];
            k = 2*pi/10e3;
            pairedEVP = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes( ...
                N2=N2,zDomain=zDomain,k=k,latitude=45,alpha0=Inf,alphaD=Inf);
            f0 = pairedEVP.f0;
            g = pairedEVP.g;
            scalarEVP = IMInternalModes(name=pairedEVP.name,N2=N2,zDomain=zDomain,formulation="F",modeFamily="none", ...
                p=@(z,~) f0*f0./N2(z),q=k*k,r=1, ...
                surfaceBoundary=pairedEVP.surfaceBoundary,bottomBoundary=pairedEVP.bottomBoundary, ...
                rotationRate=pairedEVP.rotationRate,latitude=pairedEVP.latitude,g=g, ...
                hFromEigenvalue=@(lambda) f0*f0./(g*(lambda-k*k)),parameters=pairedEVP.parameters);
            paired = IMSolverSpectral(nEVP=128,coordinateKind="z").solveEVP(pairedEVP,nModes=6);
            scalar = IMSolverSpectral(nEVP=128,coordinateKind="z").solveEVP(scalarEVP,nModes=6);
            z = linspace(-1000,0,257).';
            testCase.verifyEqual(paired.eigenvalues,scalar.eigenvalues,RelTol=1e-7)
            testCase.verifyEqual(paired.h,scalar.h,RelTol=1e-7)
            testCase.verifyEqual(paired.F(z),scalar.F(z),RelTol=1e-7,AbsTol=1e-7)
            testCase.verifyEqual(paired.G(z),scalar.G(z),RelTol=1e-7,AbsTol=1e-7)
        end

        function neutralCapRetainsSignedModesAndAffineLayer(testCase)
            evp = IMInternalModes.geostrophicAPVModes(N2=@IMNeutralGeostrophicModesTests.neutralCap,zDomain=[-1000 0],g0=0,gd=-0.01,surfaceBoundary="freeSurface");
            basis = IMSolverSpectral(nEVP=512,coordinateKind="z").solveEVP(evp,nModes=6);
            testCase.verifyTrue(any(basis.inverseDepth < 0))
            testCase.verifyPairEquations(basis,1e-7)
            testCase.verifyAPVBoundaryConditions(basis,1e-7)
            testCase.verifyEqual(basis.gramMatrix(variable="F")/1000,eye(6),AbsTol=2e-7)
            testCase.verifyEqual(basis.gramMatrix(variable="G")/1000,diag(basis.inverseDepth),AbsTol=2e-7,RelTol=2e-7)

            z = [-150;-100;-50;0];
            F = basis.F(z);
            G = basis.G(z);
            expectedG = G(end,:)+(z-z(end)).*(basis.inverseDepth.*F(end,:));
            testCase.verifyLessThan(max(abs(F-F(end,:)),[],"all"),1e-7)
            testCase.verifyLessThan(max(abs(G-expectedG),[],"all")/max(abs(G),[],"all"),1e-7)
        end

        function interiorNeutralIntervalAndIsolatedZeroAreRegular(testCase)
            profiles = {@(z) (5e-3)^2*(max(abs(z+500)-100,0)/400).^4, ...
                @(z) (5e-3)^2*((z+500)/500).^2};
            for iProfile = 1:numel(profiles)
                evp = IMInternalModes.geostrophicAPVModes(N2=profiles{iProfile},zDomain=[-1000 0],g0=Inf,gd=Inf,surfaceBoundary="rigidLid");
                basis = IMSolverSpectral(nEVP=256,coordinateKind="z").solveEVP(evp,nModes=5);
                testCase.verifyPairEquations(basis,1e-7)
                testCase.verifyAPVBoundaryConditions(basis,1e-7)
                testCase.verifyTrue(all(isfinite(basis.G([-600;-500;-400])),"all"))
                zero = basis.modeNumber == 0;
                testCase.verifyEqual(basis.inverseDepth(zero),0,AbsTol=0)
                testCase.verifyEqual(basis.h(zero),Inf,AbsTol=0)
            end
        end

        function generalizedEnstrophyPairUsesCompanionEndpointGram(testCase)
            evp = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes( ...
                N2=@IMNeutralGeostrophicModesTests.neutralCap,zDomain=[-1000 0],k=2*pi/10e3,latitude=45,alpha0=2e-11,alphaD=3e-11);
            basis = IMSolverSpectral(nEVP=512,coordinateKind="z").solveEVP(evp,nModes=6);
            testCase.verifyPairEquations(basis,1e-7)
            testCase.verifyEnstrophyBoundaryConditions(basis,1e-7)
            testCase.verifyEqual(basis.gramMatrix(variable="F")/1000,eye(6),AbsTol=2e-7)
            testCase.verifyEqual(basis.gramMatrix(variable="G")/1000,diag(basis.inverseDepth),AbsTol=2e-7,RelTol=2e-7)
            testCase.verifyGreaterThan(min(basis.eigenvalues),0)
        end

        function enstrophyZeroInverseDepthHasNonzeroPublicEigenvalue(testCase)
            k = 2*pi/500e3;
            g = 9.81;
            evp = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes( ...
                N2=@IMNeutralGeostrophicModesTests.neutralCap,zDomain=[-1000 0],k=k,g=g,alpha0=g*k*k,alphaD=Inf);
            basis = IMSolverSpectral(nEVP=128,coordinateKind="z").solveEVP(evp,nModes=4);
            zero = basis.inverseDepth == 0;
            testCase.verifyNumElements(find(zero),1)
            testCase.verifyEqual(basis.eigenvalues(zero),k*k,AbsTol=0)
            testCase.verifyEqual(basis.h(zero),Inf,AbsTol=0)
            testCase.verifyGreaterThan(min(basis.eigenvalues),0)
            G = basis.G([-1000;-500;0]);
            F = basis.F([-1000;-500;0]);
            testCase.verifyLessThan(max(abs(G(:,zero)),[],"all"),1e-10)
            testCase.verifyGreaterThan(min(F(:,zero)),0)
            [A,B] = evp.assemble(basis.solver);
            testCase.verifySize(A,[256 256])
            testCase.verifySize(B,[256 256])
            [~,~,map] = evp.prepareSolvePencil(basis.solver,A,B);
            testCase.verifyEqual(map.nativeFieldSizes,[128 128])
            testCase.verifyTrue(map.hasInverseDepth)
            T = basis.solver.physicalDerivativeMatrix(0);
            nativeG = T\basis.rawVariable("G",basis.solver.zNative);
            V = [basis.nativeModes;nativeG];
            residual = vecnorm(A*V-(B*V).*basis.eigenvalues,2,1);
            scale = (norm(A,1)+abs(basis.eigenvalues)*norm(B,1)).*vecnorm(V,2,1);
            testCase.verifyLessThan(max(residual./scale),1e-7)
            testCase.verifyFalse(isfield(basis.modeSelectionDiagnostics,"geostrophicInverseDepth"))
        end

        function bottomTrappedSignedModeHasAResolvedPhaseFallback(testCase)
            rotationRate = 7.2921e-5;
            f0 = 2*rotationRate*sind(45);
            evp = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes( ...
                N2=@IMNeutralGeostrophicModesTests.neutralCap,zDomain=[-1000 0],k=2*pi/1e3, ...
                rotationRate=rotationRate,latitude=45,alpha0=2*f0^2/1000,alphaD=3*f0^2/1000);
            basis = IMSolverSpectral(nEVP=1024,coordinateKind="z").solveEVP(evp,nModes=4);
            testCase.verifyTrue(any(basis.metadata.dominantGOrientationFallback))
            testCase.verifyTrue(all(basis.inverseDepth(basis.metadata.dominantGOrientationFallback)<0))
            testCase.verifyGreaterThan(min(basis.eigenvalues),0)
            z = basis.solver.zNative;
            G = basis.G(z);
            for j = find(basis.metadata.dominantGOrientationFallback)
                [~,index] = max(abs(G(:,j)));
                testCase.verifyGreaterThan(G(index,j),0)
            end
            repeated = basis.orientModeSigns();
            testCase.verifyEqual(repeated.F(z),basis.F(z),AbsTol=0)
            testCase.verifyEqual(repeated.G(z),G,AbsTol=0)
            testCase.verifyEnstrophyBoundaryConditions(basis,1e-7)
        end

        function scalarWavePencilKeepsDefaultMapping(testCase)
            N2 = @(z) (5e-3)^2*ones(size(z));
            evp = IMInternalModes.waveModesAtWavenumber(N2=N2,zDomain=[-1000 0],k=2*pi/10e3,latitude=45,surfaceBoundary=IMBoundaryCondition.dirichlet());
            numerical = IMSolverSpectral(nEVP=64,coordinateKind="z").solveEVP(evp,nModes=4);
            solution = IMConstantStratificationSolution(N0=5e-3,zDomain=[-1000 0],latitude=45);
            analytical = solution.internalModes(evp,nModes=4);
            [A,B] = evp.assemble(numerical.solver);
            [solveA,solveB,map] = evp.prepareSolvePencil(numerical.solver,A,B);
            testCase.verifyEqual(solveA,A,AbsTol=0)
            testCase.verifyEqual(solveB,B,AbsTol=0)
            testCase.verifyEqual(map.nativeFieldSizes,64)
            testCase.verifyFalse(map.hasInverseDepth)
            testCase.verifyEqual(numerical.h,analytical.h,RelTol=1e-7)
        end

        function nearbyEnstrophyWeightsRetainSolvedInverseDepth(testCase)
            k = 2*pi/500e3;
            g = 9.81;
            for perturbation = [-1e-5 1e-5]
                evp = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes( ...
                    N2=@(z) (5e-3)^2*ones(size(z)),zDomain=[-1000 0],k=k,g=g,alpha0=g*k*k*(1+perturbation),alphaD=Inf);
                basis = IMSolverSpectral(nEVP=96,coordinateKind="z").solveEVP(evp,nModes=4);
                [~,iNear] = min(abs(basis.inverseDepth));
                testCase.verifyNotEqual(basis.inverseDepth(iNear),0)
                testCase.verifyEqual(sign(basis.inverseDepth(iNear)),sign(perturbation))
                testCase.verifyEqual(basis.h,1./basis.inverseDepth,AbsTol=0)
                testCase.verifyPairEquations(basis,1e-7)
                testCase.verifyEnstrophyBoundaryConditions(basis,1e-7)
            end
        end

        function fullyNeutralAPVUsesOnlyItsFiniteFamily(testCase)
            N2 = @(z) zeros(size(z));
            evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=[-1000 0],g0=0.02,gd=0.03,surfaceBoundary="freeSurface");
            solver = IMSolverSpectral(nEVP=48,coordinateKind="z");
            testCase.verifyError(@() solver.solveEVP(evp,nModes=6),"IMSolver:InsufficientModes")
            basis = solver.solveEVP(evp,nModes=1);
            expected = evp.g*(1/evp.g+1/0.02+1/0.03)/1000;
            testCase.verifySize(basis.nativeModes,[48 1])
            testCase.verifyEqual(basis.inverseDepth,expected,RelTol=1e-12)
            testCase.verifyPairEquations(basis,1e-10)
            testCase.verifyAPVBoundaryConditions(basis,1e-10)

            zeroEVP = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=[-1000 0],g0=Inf,gd=Inf,surfaceBoundary="rigidLid");
            testCase.verifyError(@() solver.solveEVP(zeroEVP,nModes=6),"IMSolver:InsufficientModes")
            zeroBasis = solver.solveEVP(zeroEVP,nModes=1);
            testCase.verifyEqual(zeroBasis.inverseDepth,0,AbsTol=0)
            testCase.verifyEqual(zeroBasis.h,Inf,AbsTol=0)
            testCase.verifyLessThan(max(abs(zeroBasis.G([-1000;0]))),1e-10)
        end

        function fullyNeutralEnstrophyRequiresExplicitWeights(testCase)
            N2 = @(z) zeros(size(z));
            k = 2*pi/10e3;
            alpha0 = 2e-11;
            alphaD = 3e-11;
            evp = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes( ...
                N2=N2,zDomain=[-1000 0],k=k,alpha0=alpha0,alphaD=alphaD);
            solver = IMSolverSpectral(nEVP=48,coordinateKind="z");
            testCase.verifyError(@() solver.solveEVP(evp,nModes=6),"IMSolver:InsufficientModes")
            basis = solver.solveEVP(evp,nModes=1);
            weightSum = 1/alpha0+1/alphaD;
            expectedR = (1-evp.g*k*k*weightSum)/(1000+evp.f0^2*weightSum);
            expectedLambda = (k*k*1000+evp.f0^2/evp.g)/(1000+evp.f0^2*weightSum);
            testCase.verifySize(basis.nativeModes,[48 1])
            testCase.verifyEqual(basis.inverseDepth,expectedR,RelTol=1e-12)
            testCase.verifyEqual(basis.eigenvalues,expectedLambda,RelTol=1e-12)
            testCase.verifyPairEquations(basis,1e-10)
            testCase.verifyEnstrophyBoundaryConditions(basis,1e-10)
            testCase.verifyError(@() IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes( ...
                N2=N2,zDomain=[-1000 0],k=k),"IMInternalModes:UndefinedStratificationScale")
        end

        function fullyNeutralEnstrophyZeroDepthUsesExactFiniteFamily(testCase)
            k = 1e-4;
            g = 9.81;
            evp = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes( ...
                N2=@(z) zeros(size(z)),zDomain=[-1000 0],k=k,g=g,alpha0=g*k*k,alphaD=Inf);
            solver = IMSolverSpectral(nEVP=64,coordinateKind="z");
            testCase.verifyError(@() solver.solveEVP(evp,nModes=4),"IMSolver:InsufficientModes")
            basis = solver.solveEVP(evp,nModes=1);
            testCase.verifyNumElements(basis.h,1)
            testCase.verifyEqual(basis.inverseDepth,0,AbsTol=0)
            testCase.verifyEqual(basis.h,Inf,AbsTol=0)
            testCase.verifyEqual(basis.eigenvalues,k*k,AbsTol=0)
            testCase.verifyLessThan(max(abs(basis.G([-1000;-500;0]))),1e-10)
            testCase.verifyPairEquations(basis,1e-10)
            testCase.verifyEnstrophyBoundaryConditions(basis,1e-10)
        end

        function missingAPVPhysicalSectorsCannotBeReplacedByPositiveModes(testCase)
            N2 = @(z) (5e-3)^2*ones(size(z));
            evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=[-1000 0],g0=-0.02,gd=Inf,surfaceBoundary="freeSurface");
            solver = IMSolverSpectral(nEVP=64,coordinateKind="z");
            [A,~] = evp.assemble(solver);
            testCase.verifyError(@() evp.selectModes([1;2;3],2,solver,A),"IMGeostrophicModes:UnresolvedPhysicalSector")

            evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=[-1000 0],g0=Inf,gd=Inf,surfaceBoundary="rigidLid");
            [A,~] = evp.assemble(solver);
            testCase.verifyError(@() evp.selectModes([1;2;3],2,solver,A),"IMGeostrophicModes:UnresolvedPhysicalSector")
            selected = evp.selectModes([0;1;2],3,solver,A);
            testCase.verifyEqual(selected.modeNumber,[0 1 2])
        end

        function unqualifiedRequestedModeCannotBeSkipped(testCase)
            evp = IMUnqualifiedSecondModeTestProblem();
            solver = IMSolverSpectral(nEVP=64,coordinateKind="z");
            first = solver.solveEVP(evp,nModes=1);
            testCase.verifyEqual(first.modeNumber,1)
            testCase.verifyEqual(first.eigenvalues,pi*pi,RelTol=1e-10)
            testCase.verifyError(@() solver.solveEVP(evp,nModes=3),"IMGeostrophicModes:UnresolvedMode")
        end

        function invalidStratificationAndFNullFamiliesAreRejected(testCase)
            profiles = {@(z) -ones(size(z)),@(z) 1i*ones(size(z)),@(z) Inf(size(z)),@(z) NaN(size(z))};
            for iProfile = 1:numel(profiles)
                testCase.verifyError(@() IMInternalModes.geostrophicAPVModes(N2=profiles{iProfile},zDomain=[-1000 0],g0=Inf,gd=Inf,surfaceBoundary="freeSurface"), ...
                    "IMGeostrophicModes:InvalidStratification")
            end
            evp = IMInternalModes.geostrophicAPVModes(N2=@(z) zeros(size(z)),zDomain=[-1000 0],g0=0,gd=Inf,surfaceBoundary="freeSurface");
            testCase.verifyError(@() IMSolverSpectral(nEVP=48,coordinateKind="z").solveEVP(evp,nModes=4),"IMSolver:NoValidEigenvalues")
        end

        function storedPairsShareNormalizationsAndOrientation(testCase)
            evp = IMInternalModes.geostrophicAPVModes(N2=@(z) (5e-3)^2*((z+500)/500).^2,zDomain=[-1000 0],g0=Inf,gd=Inf,surfaceBoundary="rigidLid");
            basis = IMSolverSpectral(nEVP=128,coordinateKind="z").solveEVP(evp,nModes=5);
            z = linspace(-1000,0,65).';
            raw = basis.addNormalization("raw",@(~,~) 1);
            scale = basis.normalizationFactors();
            testCase.verifyEqual(basis.F(z).*scale,raw.F(z,normalization="raw"),RelTol=1e-12,AbsTol=1e-12)
            testCase.verifyEqual(basis.G(z).*scale,raw.G(z,normalization="raw"),RelTol=1e-12,AbsTol=1e-12)
            testCase.verifyEqual(basis.u(z),basis.F(z),AbsTol=0)
            oriented = basis.orientModeSigns();
            testCase.verifyEqual(oriented.F(z),basis.F(z),AbsTol=0)
            testCase.verifyEqual(oriented.G(z),basis.G(z),AbsTol=0)
            zNative = basis.solver.zNative;
            expectedGz = basis.solver.differentiateGridValues(basis.G(zNative),1);
            testCase.verifyEqual(basis.Gz(zNative),expectedGz,AbsTol=1e-10,RelTol=1e-10)
        end

        function neutralPairSupportsFProjectionAndRootGrid(testCase)
            evp = IMInternalModes.geostrophicAPVModes(N2=@(z) (5e-3)^2*((z+500)/500).^2,zDomain=[-1000 0],g0=Inf,gd=Inf,surfaceBoundary="rigidLid");
            basis = IMSolverSpectral(nEVP=128,coordinateKind="z").solveEVP(evp,nModes=8);
            z = linspace(-1000,0,257).';
            weights = [0.5;ones(255,1);0.5]*(1000/256);
            transform = basis.discreteTransform(z=z,weights=weights,nModes=5,variables="F",gramTolerance=1);
            coefficients = reshape(1:20,5,4)+1i*reshape(20:-1:1,5,4);
            values = transform.transformBack(coefficients,variable="F");
            testCase.verifyEqual(transform.transformForward(values,variable="F"),coefficients,RelTol=1e-9,AbsTol=1e-9)
            [roots,design] = basis.modeRootGrid(nPoints=6);
            testCase.verifyNumElements(roots,6)
            testCase.verifyGreaterThanOrEqual(min(roots),-1000)
            testCase.verifyLessThanOrEqual(max(roots),0)
            testCase.verifyGreaterThan(design.representedModeCount,0)
        end
    end

    methods (Access = private)
        function scales = independentDepthScales(testCase,basis)
            scales = zeros(1,numel(basis.h));
            for iMode = 1:numel(basis.h)
                scales(iMode) = sqrt(integral(@(z) testCase.modeSquare(basis,z,iMode), ...
                    basis.zDomain(1),basis.zDomain(2),AbsTol=1e-12,RelTol=1e-12)/diff(basis.zDomain));
            end
        end

        function values = modeSquare(~,basis,z,iMode)
            F = basis.F(z(:));
            values = reshape(F(:,iMode).^2,size(z));
        end

        function verifyPairEquations(testCase,basis,tolerance)
            z = basis.solver.zNative;
            F = basis.F(z);
            G = basis.G(z);
            Fz = basis.uz(z);
            Gz = basis.Gz(z);
            N2 = basis.N2(z);
            first = Fz+(N2(:)/basis.evp.g).*G;
            second = Gz-F.*basis.inverseDepth;
            interior = true(size(z));
            interior([basis.solver.boundaryIndex("surface") basis.solver.boundaryIndex("bottom")]) = false;
            depth = diff(basis.zDomain);
            firstScale = max(abs(Fz),[],1)+max(abs((N2(:)/basis.evp.g).*G),[],1)+max(abs(F),[],1)/depth;
            secondScale = max(abs(Gz),[],1)+max(abs(F.*basis.inverseDepth),[],1)+max(abs(G),[],1)/depth;
            testCase.verifyLessThan(max(abs(first(interior,:))./firstScale,[],"all"),tolerance)
            testCase.verifyLessThan(max(abs(second(interior,:))./max(secondScale,realmin),[],"all"),tolerance)
        end

        function verifyAPVBoundaryConditions(testCase,basis,tolerance)
            p = basis.evp.parameters;
            F = basis.F(basis.zDomain.');
            G = basis.G(basis.zDomain.');
            sigma = double(string(p.surfaceBoundary) == "freeSurface");
            if isinf(p.g0)
                surface = G(2,:)-sigma*F(2,:);
            else
                surface = p.g0*G(2,:)-(basis.evp.g+sigma*p.g0)*F(2,:);
            end
            if isinf(p.gd)
                bottom = G(1,:);
            else
                bottom = basis.evp.g*F(1,:)+p.gd*G(1,:);
            end
            scale = max(abs(F),[],1)+max(abs(G),[],1);
            testCase.verifyLessThan(max(abs(surface)./scale),tolerance)
            testCase.verifyLessThan(max(abs(bottom)./scale),tolerance)
        end

        function verifyEnstrophyBoundaryConditions(testCase,basis,tolerance)
            p = basis.evp.parameters;
            F = basis.F(basis.zDomain.');
            G = basis.G(basis.zDomain.');
            if isinf(p.alpha0)
                surface = G(2,:)-F(2,:);
            else
                surface = G(2,:)-F(2,:)+(basis.evp.g/p.alpha0)*basis.eigenvalues.*F(2,:);
            end
            if isinf(p.alphaD)
                bottom = G(1,:);
            else
                bottom = G(1,:)-(basis.evp.g/p.alphaD)*basis.eigenvalues.*F(1,:);
            end
            scale = max(abs(F),[],1)+max(abs(G),[],1);
            testCase.verifyLessThan(max(abs(surface)./scale),tolerance)
            testCase.verifyLessThan(max(abs(bottom)./scale),tolerance)
        end
    end

    methods (Static,Access = private)
        function N2 = neutralCap(z)
            t = min(1,max(0,(-z-200)/200));
            N2 = (5e-3)^2*(10*t.^3-15*t.^4+6*t.^5);
        end
    end
end
