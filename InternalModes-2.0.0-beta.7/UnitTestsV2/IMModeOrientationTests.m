classdef IMModeOrientationTests < matlab.unittest.TestCase

    properties (Access = private)
        originalPath
    end

    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.originalPath = path;
            addpath(repoRoot,fullfile(repoRoot,"UnitTestsV2","TestSupport"));
        end
    end

    methods (TestClassTeardown)
        function restoreRepositoryPath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function numericalOrientationIsInvariantToRawEigenvectorSigns(testCase)
            z = linspace(-1,0,65).';
            N2 = @(z) ones(size(z));
            evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-1 0]);
            solver = IMSolverFiniteDifference(z=z).configuredForEVP(evp);
            GValues = [-ones(size(solver.zNative)), solver.zNative];

            first = IMInternalModesBasis(solver=solver,evp=evp,nativeModes=GValues, ...
                eigenvalues=[1 2],modeNumber=[1 2]);
            second = IMInternalModesBasis(solver=solver,evp=evp,nativeModes=-GValues, ...
                eigenvalues=[1 2],modeNumber=[1 2]);
            first = first.orientModeSigns();
            second = second.orientModeSigns();

            zCheck = linspace(-1,0,101).';
            testCase.verifyEqual(first.rawVariable("G",zCheck),second.rawVariable("G",zCheck),AbsTol=1e-12)
            testCase.verifyGreaterThan(first.rawVariable("G",-1e-6),zeros(1,2))
            testCase.verifyEqual(first.metadata.modeOrientation,"shallowInteriorGPositive-v1")
        end

        function FFormNullModeUsesDeterministicFallback(testCase)
            z = linspace(-1,0,65).';
            N2 = @(z) ones(size(z));
            evp = IMInternalModes.hydrostaticFModes(N2=N2,zDomain=[-1 0]);
            solver = IMSolverFiniteDifference(z=z).configuredForEVP(evp);
            basisSet = IMInternalModesBasis(solver=solver,evp=evp,nativeModes=-ones(size(z)), ...
                eigenvalues=0,modeNumber=0);

            basisSet = basisSet.orientModeSigns();

            testCase.verifyGreaterThan(basisSet.rawVariable("F",0),0)
            testCase.verifyEqual(basisSet.rawVariable("G",z),zeros(size(z)),AbsTol=1e-12)
        end

        function APVAndMDAFamiliesUseTheSameShallowInteriorConvention(testCase)
            D = 4000;
            N0 = 5.2e-3;
            N2 = @(z) N0*N0*ones(size(z));
            solver = IMSolverSpectral(nEVP=128);
            apvEVP = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=[-D 0],g0=0.02,gd=0.03);
            mdaEVP = IMInternalModes.meanDensityAnomalyModes(N2=N2,zDomain=[-D 0],g0=0.02,gd=0.03);

            apvBasis = solver.solveEVP(apvEVP,nModes=6);
            mdaBasis = solver.solveEVP(mdaEVP,nModes=6);
            zShallow = -1e-6*D;

            testCase.verifyGreaterThan(apvBasis.G(zShallow,normalization="wMax"),zeros(1,6))
            testCase.verifyGreaterThan(mdaBasis.G(zShallow,normalization="wMax"),zeros(1,6))
            testCase.verifyEqual(apvBasis.metadata.modeOrientation,IMModeOrientationTools.convention)
            testCase.verifyEqual(mdaBasis.metadata.modeOrientation,IMModeOrientationTools.convention)
        end

        function analyticalAndNumericalRigidSurfaceModesHaveMatchingPhase(testCase)
            D = 3000;
            N0 = 5.2e-3;
            N2 = @(z) N0*N0*ones(size(z));
            evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-D 0]);
            analytical = IMConstantStratificationSolution(N0=N0,zDomain=[-D 0]).internalModes(evp,nModes=4);
            numerical = IMSolverSpectral(nEVP=128).solveEVP(evp,nModes=4);
            z = linspace(-D,0,2001).';
            analyticalG = analytical.G(z,normalization="wMax");
            numericalG = numerical.G(z,normalization="wMax");

            correlation = sum(analyticalG.*numericalG,1);
            testCase.verifyGreaterThan(correlation,zeros(1,4))
            testCase.verifyGreaterThan(analytical.G(-1e-6*D,normalization="wMax"),zeros(1,4))
            testCase.verifyEqual(analytical.metadata.modeOrientation,IMModeOrientationTools.convention)
        end

        function optimizedSpectralOrientationMatchesLegacyAcrossProfilesAndBoundaries(testCase)
            D = 1000;
            constantN2 = @(z) 1e-4*ones(size(z));
            exponentialN2 = @(z) 1e-4*exp(2*z/700);
            transition = @(z) 0.5*(1+tanh((z+40)/10));
            sharpN2 = @(z) transition(z)*(8e-4)^2+(1-transition(z)).*((5e-4)^2+((3e-3)^2-(5e-4)^2)*exp((z+40)/80));
            profiles = {constantN2,exponentialN2,sharpN2};
            coordinateKinds = ["z","wkb","density"];
            surfaceBoundaries = {IMBoundaryCondition.dirichlet(),IMBoundaryCondition(a=0,b=1,c=1,d=0)};

            for iProfile = 1:numel(profiles)
                for iSurface = 1:numel(surfaceBoundaries)
                    evp = IMInternalModes.hydrostaticGModes(N2=profiles{iProfile},zDomain=[-D 0],surfaceBoundary=surfaceBoundaries{iSurface});
                    source = IMSolverSpectral(nEVP=96,coordinateKind=coordinateKinds(iProfile)).solveEVP(evp,nModes=6);
                    scale = (-1).^(1:numel(source.eigenvalues)).*10.^linspace(-8,8,numel(source.eigenvalues));
                    testCase.verifyLegacyEquivalentOrientation(testCase.scaledCopy(source,scale));
                    testCase.verifyLegacyEquivalentOrientation(testCase.scaledCopy(source,-scale));
                end
            end
        end

        function longExternalAndNearInertialWavesMatchLegacyOrientation(testCase)
            D = 4000;
            f0 = 1e-4;
            N2 = @(z) 1e-4*exp(2*z/1400);
            freeSurface = IMBoundaryCondition(a=0,b=1,c=1,d=0);
            longExternal = IMInternalModes.waveModesAtWavenumber(N2=N2,zDomain=[-D 0],k=2*pi/1e7,f0=f0,surfaceBoundary=freeSurface);
            zeroWavenumber = IMInternalModes.waveModesAtWavenumber(N2=N2,zDomain=[-D 0],k=0,f0=f0,surfaceBoundary=freeSurface);
            nearInertial = IMInternalModes.waveModesAtFrequency(N2=N2,zDomain=[-D 0],omega=1.01*f0,f0=f0);
            problems = {longExternal,zeroWavenumber,nearInertial};

            for iProblem = 1:numel(problems)
                source = IMSolverSpectral(nEVP=96,coordinateKind="wkb").solveEVP(problems{iProblem},nModes=6);
                scale = [1e-9 -1e-3 1 -7 2e4 -3e8];
                testCase.verifyLegacyEquivalentOrientation(testCase.scaledCopy(source,scale));
                testCase.verifyLegacyEquivalentOrientation(testCase.scaledCopy(source,-scale));
            end
        end

        function spectralSurfaceAndDerivativeThresholdsMatchLegacyDecisions(testCase)
            relativeTolerance = 1e-10;
            evp = IMInternalModes.hydrostaticGModes(N2=@(z) ones(size(z)),zDomain=[-1 0]);
            solver = IMSolverSpectral(nEVP=33).configuredForEVP(evp);
            coefficients = zeros(solver.nEVP,4);
            resolvedSurface = 2*relativeTolerance;
            unresolvedSurface = 0.5*relativeTolerance;
            resolvedDerivative = 2*relativeTolerance;
            coefficients(1:2,1) = [resolvedSurface-0.5;0.5];
            coefficients(1:2,2) = [unresolvedSurface-0.5;0.5];
            coefficients(1:3,3) = [3/8-resolvedDerivative/2;-0.5+resolvedDerivative/2;1/8];
            coefficients(1:2,4) = 1e-100*[-0.5;0.5];
            basisSet = IMInternalModesBasis(solver=solver,evp=evp,nativeModes=coefficients,eigenvalues=ones(1,4),modeNumber=1:4);

            legacySigns = testCase.legacyOrientationSigns(basisSet);
            [actual,actualSigns] = basisSet.orientModeSigns();

            testCase.verifyEqual(legacySigns,[1 -1 -1 -1])
            testCase.verifyEqual(actualSigns,legacySigns)
            testCase.verifyEqual(actual.nativeModes,basisSet.nativeModes.*legacySigns,AbsTol=0)
        end

        function unresolvedSurfaceAndDerivativeRemainIndeterminate(testCase)
            relativeTolerance = 1e-10;
            derivative = 0.5*relativeTolerance;
            evp = IMInternalModes.hydrostaticGModes(N2=@(z) ones(size(z)),zDomain=[-1 0]);
            solver = IMSolverSpectral(nEVP=33).configuredForEVP(evp);
            coefficients = zeros(solver.nEVP,1);
            coefficients(1:3) = [3/8-derivative/2;-0.5+derivative/2;1/8];
            basisSet = IMInternalModesBasis(solver=solver,evp=evp,nativeModes=coefficients,eigenvalues=1,modeNumber=1);

            testCase.verifyError(@() testCase.legacyOrientationSigns(basisSet),"IMModeOrientationTools:IndeterminateOrientation")
            testCase.verifyError(@() basisSet.orientModeSigns(),"IMModeOrientationTools:IndeterminateOrientation")
        end

        function nearDerivativeCutoffMatchesSampledLegacyOutcome(testCase)
            relativeTolerance = 1e-10;
            derivativeFactors = 1+[-1e-3 -1e-4 0 1e-4 1e-3];
            evp = IMInternalModes.hydrostaticGModes(N2=@(z) ones(size(z)),zDomain=[-1 0]);
            solver = IMSolverSpectral(nEVP=33).configuredForEVP(evp);

            for derivativeFactor = derivativeFactors
                derivative = relativeTolerance*derivativeFactor;
                coefficients = zeros(solver.nEVP,4);
                coefficients(1,1) = 1;
                coefficients(1:3,2) = [3/8-derivative/2;-0.5+derivative/2;1/8];
                coefficients(1,3) = -1;
                coefficients(1:2,4) = 1e-100*[0.75;0.25];
                for columns = {2,1:4}
                    selected = columns{1};
                    basisSet = IMInternalModesBasis(solver=solver,evp=evp,nativeModes=coefficients(:,selected),eigenvalues=ones(1,numel(selected)),modeNumber=selected);
                    legacyOutcome = testCase.orientationOutcome(basisSet,"legacy");
                    actualOutcome = testCase.orientationOutcome(basisSet,"optimized");
                    testCase.verifyEqual(actualOutcome,legacyOutcome)
                end
            end
        end

        function spectralFFormNearZeroGUsesLegacyBarotropicFallback(testCase)
            relativeTolerance = 1e-10;
            g = 9.81;
            evp = IMInternalModes.hydrostaticFModes(N2=@(z) ones(size(z)),zDomain=[-1 0],g=g);
            solver = IMSolverSpectral(nEVP=33).configuredForEVP(evp);
            lowSlope = 0.5*relativeTolerance/g;
            highSlope = 2*relativeTolerance/g;
            coefficients = zeros(solver.nEVP,3);
            coefficients(1,1) = -1;
            coefficients(1:2,2) = [-1+lowSlope/2;-lowSlope/2];
            coefficients(1:2,3) = [-1+highSlope/2;-highSlope/2];
            basisSet = IMInternalModesBasis(solver=solver,evp=evp,nativeModes=coefficients,eigenvalues=zeros(1,3),modeNumber=zeros(1,3));

            legacySigns = testCase.legacyOrientationSigns(basisSet);
            [actual,actualSigns] = basisSet.orientModeSigns();

            testCase.verifyEqual(legacySigns,[-1 -1 1])
            testCase.verifyEqual(actualSigns,legacySigns)
            testCase.verifyEqual(actual.nativeModes,basisSet.nativeModes.*legacySigns,AbsTol=0)
        end

        function customSolverAndBasisOverridesRetainLegacyDispatch(testCase)
            evp = IMInternalModes.hydrostaticGModes(N2=@(z) ones(size(z)),zDomain=[-1 0]);
            customSolver = IMOrientationDerivativeOverrideSolver(nEVP=33).configuredForEVP(evp);
            coefficients = zeros(customSolver.nEVP,1);
            coefficients(1:2) = [-0.5;0.5];
            solverBasis = IMInternalModesBasis(solver=customSolver,evp=evp,nativeModes=coefficients,eigenvalues=1,modeNumber=1);
            standardSolver = IMSolverSpectral(nEVP=33).configuredForEVP(evp);
            standardBasis = IMInternalModesBasis(solver=standardSolver,evp=evp,nativeModes=coefficients,eigenvalues=1,modeNumber=1);
            customBasis = IMOrientationRawGOverrideBasis(standardBasis);

            testCase.verifyEqual(testCase.legacyOrientationSigns(solverBasis),1)
            testCase.verifyEqual(testCase.legacyOrientationSigns(customBasis),1)
            testCase.verifyLegacyEquivalentOrientation(solverBasis);
            testCase.verifyLegacyEquivalentOrientation(customBasis);
        end

        function overResolvedNativeCoefficientsRetainLegacyAliasedDerivative(testCase)
            evp = IMInternalModes.hydrostaticGModes(N2=@(z) ones(size(z)),zDomain=[-1 0]);
            solver = IMSolverSpectral(nEVP=4).configuredForEVP(evp);
            coefficients = [1;0;-2;0;1];
            basisSet = IMInternalModesBasis(solver=solver,evp=evp,nativeModes=coefficients,eigenvalues=1,modeNumber=1);
            legacySigns = testCase.legacyOrientationSigns(basisSet);
            directDerivative = solver.evaluatePhysicalDerivative(coefficients,solver.zNative,1);

            testCase.verifyEqual(legacySigns,1)
            testCase.verifyGreaterThan(directDerivative(solver.boundaryIndex("surface")),0)
            testCase.verifyLegacyEquivalentOrientation(basisSet);
        end
    end

    methods (Access = private)
        function basisSet = scaledCopy(~,source,scale)
            basisSet = IMInternalModesBasis(solver=source.solver,evp=source.evp,nativeModes=source.nativeModes.*scale, ...
                eigenvalues=source.eigenvalues,modeNumber=source.modeNumber,modeSelectionDiagnostics=source.modeSelectionDiagnostics);
        end

        function signs = legacyOrientationSigns(~,basisSet)
            zNative = basisSet.solver.zNative;
            surfaceIndex = basisSet.solver.boundaryIndex("surface");
            GValues = basisSet.rawVariable("G",zNative);
            GzValues = basisSet.solver.differentiateGridValues(GValues,1);
            FValues = basisSet.rawVariable("F",zNative);
            allowFFallback = basisSet.evp.formulation == "F" & basisSet.eigenvalues == 0;
            signs = IMModeOrientationTools.shallowInteriorGPositive(GValues=GValues,GzSurface=GzValues(surfaceIndex,:), ...
                FValues=FValues,depth=diff(basisSet.zDomain),surfaceIndex=surfaceIndex,allowFFallback=allowFFallback);
        end

        function verifyLegacyEquivalentOrientation(testCase,basisSet)
            legacySigns = testCase.legacyOrientationSigns(basisSet);
            [actual,actualSigns] = basisSet.orientModeSigns();
            testCase.verifyEqual(actualSigns,legacySigns)
            testCase.verifyEqual(actual.nativeModes,basisSet.nativeModes.*legacySigns,AbsTol=0)
            testCase.verifyEqual(actual.metadata.modeOrientation,IMModeOrientationTools.convention)
        end

        function outcome = orientationOutcome(testCase,basisSet,implementation)
            try
                if implementation == "legacy"
                    signs = testCase.legacyOrientationSigns(basisSet);
                else
                    [~,signs] = basisSet.orientModeSigns();
                end
                outcome = struct("status","resolved","signs",signs);
            catch exception
                outcome = struct("status",string(exception.identifier),"signs",zeros(1,0));
            end
        end
    end
end
