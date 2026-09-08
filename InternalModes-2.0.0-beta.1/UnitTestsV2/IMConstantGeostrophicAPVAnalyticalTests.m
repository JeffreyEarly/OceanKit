classdef IMConstantGeostrophicAPVAnalyticalTests < matlab.unittest.TestCase

    properties (Access = private)
        originalPath
    end

    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.originalPath = path;
            addpath(repoRoot);
            addpath(fullfile(repoRoot, "UnitTestsV2"));
        end
    end

    methods (TestClassTeardown)
        function restoreRepositoryPath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function positiveAndNegativeRootsMatchIndependentEquations(testCase)
            [N0, zDomain, g, N2] = testCase.profile();
            positiveParameters = testCase.parametersFromStrengths(N0, zDomain, g, 1, 1, "rigidLid");
            negativeParameters = testCase.parametersFromStrengths(N0, zDomain, g, -4, 0, "rigidLid");
            positiveEVP = IMInternalModes.geostrophicAPVModes(N2=N2, zDomain=zDomain, g=g, g0=positiveParameters.g0, gd=positiveParameters.gd, surfaceBoundary="rigidLid");
            negativeEVP = IMInternalModes.geostrophicAPVModes(N2=N2, zDomain=zDomain, g=g, g0=negativeParameters.g0, gd=negativeParameters.gd, surfaceBoundary="rigidLid");
            solution = IMConstantStratificationSolution(N0=N0, zDomain=zDomain, g=g);
            positiveBasis = solution.internalModes(positiveEVP, nModes=4);
            negativeBasis = solution.internalModes(negativeEVP, nModes=4);

            positiveRoots = testCase.independentRoots(N0, zDomain, positiveParameters.betaSurface, positiveParameters.betaBottom, "active", "active", "positive", 4);
            negativeRoots = testCase.independentRoots(N0, zDomain, negativeParameters.betaSurface, negativeParameters.betaBottom, "active", "active", "negative", 1);
            positiveScale = N0*N0*diff(zDomain)^2/g;

            testCase.verifyEqual(positiveBasis.metadata.roots, positiveRoots, RelTol=2e-11)
            testCase.verifyEqual(positiveBasis.h, positiveScale./(positiveRoots.^2), RelTol=2e-11)
            testCase.verifyEqual(negativeBasis.metadata.roots(1), negativeRoots(1), RelTol=2e-11)
            testCase.verifyEqual(negativeBasis.h(1), -positiveScale/(negativeRoots(1)^2), RelTol=2e-11)
            testCase.verifyEqual(negativeBasis.metadata.modeKinds, ["negative" "positive" "positive" "positive"])
            testCase.verifyLessThan(max(abs([positiveBasis.metadata.rootResiduals negativeBasis.metadata.rootResiduals])), 2e-10)
        end

        function affineZeroBranchMatchesIndependentThresholds(testCase)
            [N0, zDomain, g, N2] = testCase.profile();
            IN = N0*N0*diff(zDomain);
            cases = {
                "freeSurface", 0, -IN
                "rigidLid", 0, -IN
                "freeSurface", testCase.g0FromStrength(N0, zDomain, g, -1, "freeSurface"), 0
                "rigidLid", -IN, 0};
            solution = IMConstantStratificationSolution(N0=N0, zDomain=zDomain, g=g);
            z = linspace(zDomain(1), zDomain(2), 101).';
            for iCase = 1:size(cases,1)
                surfaceConvention = string(cases{iCase,1});
                evp = IMInternalModes.geostrophicAPVModes(N2=N2, zDomain=zDomain, g=g, g0=cases{iCase,2}, gd=cases{iCase,3}, surfaceBoundary=surfaceConvention);
                basisSet = solution.internalModes(evp, nModes=3);
                zeroIndex = find(basisSet.modeNumber == 0);
                F = basisSet.F(z, normalization="uMax");
                G = basisSet.G(z, normalization="uMax");

                testCase.verifyNumElements(zeroIndex, 1)
                testCase.verifyEqual(basisSet.eigenvalues(zeroIndex), 0, AbsTol=0)
                testCase.verifyEqual(basisSet.h(zeroIndex), Inf, AbsTol=0)
                testCase.verifyLessThan(max(abs(diff(F(:,zeroIndex),2))), 2e-13)
                testCase.verifyLessThan(max(abs(diff(G(:,zeroIndex)))), 2e-13)
            end
        end

        function endpointInertiaCoversBothConventionsAndAllSigns(testCase)
            [N0, zDomain, g, N2] = testCase.profile();
            strengths = [1 1; -4 0; -8 -4; -0.5 1];
            expectedNegative = [0 1 2 0];
            expectedZero = [0 0 0 1];
            solution = IMConstantStratificationSolution(N0=N0, zDomain=zDomain, g=g);
            for surfaceConvention = ["freeSurface", "rigidLid"]
                for iCase = 1:size(strengths,1)
                    parameters = testCase.parametersFromStrengths(N0, zDomain, g, strengths(iCase,1), strengths(iCase,2), surfaceConvention);
                    evp = IMInternalModes.geostrophicAPVModes(N2=N2, zDomain=zDomain, g=g, g0=parameters.g0, gd=parameters.gd, surfaceBoundary=surfaceConvention);
                    basisSet = solution.internalModes(evp, nModes=4);

                    testCase.verifyEqual(nnz(basisSet.h < 0), expectedNegative(iCase))
                    testCase.verifyEqual(nnz(isinf(basisSet.h)), expectedZero(iCase))
                    testCase.verifyEqual(basisSet.metadata.negativeModeCount, expectedNegative(iCase))
                    testCase.verifyTrue(all(diff(basisSet.eigenvalues) >= 0))
                end
            end
        end

        function zeroAndInfiniteEndpointLimitsAreReducedExactly(testCase)
            [N0, zDomain, g, N2] = testCase.profile();
            solution = IMConstantStratificationSolution(N0=N0, zDomain=zDomain, g=g);
            cases = {
                "freeSurface", -g, Inf, true
                "rigidLid", Inf, Inf, true
                "freeSurface", 0, 0, false
                "freeSurface", Inf, Inf, false
                "rigidLid", 0, Inf, false
                "rigidLid", Inf, 0, false};
            for iCase = 1:size(cases,1)
                evp = IMInternalModes.geostrophicAPVModes(N2=N2, zDomain=zDomain, g=g, g0=cases{iCase,2}, gd=cases{iCase,3}, surfaceBoundary=string(cases{iCase,1}));
                basisSet = solution.internalModes(evp, nModes=3);
                testCase.verifyEqual(any(isinf(basisSet.h)), cases{iCase,4})
            end
        end

        function exactModesSatisfyODEDiagnosticRelationAndDepthNormalization(testCase)
            [N0, zDomain, g, N2] = testCase.profile();
            parameters = testCase.parametersFromStrengths(N0, zDomain, g, -4, 1, "freeSurface");
            evp = IMInternalModes.geostrophicAPVModes(N2=N2, zDomain=zDomain, g=g, g0=parameters.g0, gd=parameters.gd, surfaceBoundary="freeSurface");
            basisSet = IMConstantStratificationSolution(N0=N0, zDomain=zDomain, g=g).internalModes(evp, nModes=4);
            solver = IMSolverSpectral(nEVP=160, coordinateKind="z").configuredForEVP(evp);
            z = solver.zNative;
            F = basisSet.F(z);
            G = basisSet.G(z);
            Fz = basisSet.uz(z);
            Fzz = solver.differentiateGridValues(Fz,1);
            odeResidual = -Fzz/(N0*N0) - F.*reshape(1./(g*basisSet.h),1,[]);
            odeScale = max(1e-12,max(abs(F/(N0*N0*diff(zDomain)^2)),[],1));
            zDepth = linspace(zDomain(1), zDomain(2), 4001).';
            FDepth = basisSet.F(zDepth);

            testCase.verifyEqual(string(basisSet.normalization), "depth")
            testCase.verifyEqual(trapz(zDepth,FDepth.*FDepth,1)/diff(zDomain), ones(1,4), RelTol=2e-5, AbsTol=2e-5)
            testCase.verifyLessThan(max(abs(G + (g/(N0*N0))*Fz),[],"all"), 2e-12)
            testCase.verifyLessThan(max(abs(odeResidual(2:end-1,:))./odeScale,[],"all"), 2e-6)
            testCase.verifyTrue(all(basisSet.normalizationFactors("depth") > 0))
        end

        function refinedSpectralSolveMatchesRootsAndClusteredSubspace(testCase)
            [N0, zDomain, g, N2] = testCase.profile();
            parameters = testCase.parametersFromStrengths(N0, zDomain, g, -8, -8, "rigidLid");
            evp = IMInternalModes.geostrophicAPVModes(N2=N2, zDomain=zDomain, g=g, g0=parameters.g0, gd=parameters.gd, surfaceBoundary="rigidLid");
            solution = IMConstantStratificationSolution(N0=N0, zDomain=zDomain, g=g);
            analytical = solution.internalModes(evp, nModes=4);
            lowOrder = IMSolverSpectral(nEVP=48, coordinateKind="z").solveEVP(evp, nModes=4);
            highOrder = IMSolverSpectral(nEVP=96, coordinateKind="z").solveEVP(evp, nModes=4);
            lowError = norm(lowOrder.eigenvalues - analytical.eigenvalues)/norm(analytical.eigenvalues);
            highError = norm(highOrder.eigenvalues - analytical.eigenvalues)/norm(analytical.eigenvalues);
            independentNegativeRoots = [fzero(@(root) testCase.independentDeterminant(N0, zDomain, parameters.betaSurface, parameters.betaBottom, "active", "active", root, "negative"), [7.98 8]), fzero(@(root) testCase.independentDeterminant(N0, zDomain, parameters.betaSurface, parameters.betaBottom, "active", "active", root, "negative"), [8 8.02])];

            z = linspace(zDomain(1), zDomain(2), 801).';
            analyticalF = analytical.F(z, normalization="uMax");
            numericalF = highOrder.F(z, normalization="uMax");
            subspaceError = IMGeostrophicAPVTestSupport.subspaceError(z, analyticalF(:,1:2), numericalF(:,1:2));

            testCase.verifyEqual(sort(analytical.metadata.roots(1:2)), sort(independentNegativeRoots), RelTol=2e-10)
            testCase.verifyLessThan(lowError, 2e-9)
            testCase.verifyLessThan(highError, 2e-9)
            testCase.verifyLessThan(subspaceError, 2e-6)
        end

        function metadataAdvertisesPublicAPVContract(testCase)
            [N0, zDomain, g, N2] = testCase.profile();
            evp = IMInternalModes.geostrophicAPVModes(N2=N2, zDomain=zDomain, g=g, g0=-0.02, gd=Inf, surfaceBoundary="freeSurface");
            solution = IMConstantStratificationSolution(N0=N0, zDomain=zDomain, g=g);
            basisSet = solution.internalModes(evp, nModes=3);

            testCase.verifyTrue(ismember("depth", basisSet.normalizationNames()))
            testCase.verifyEqual(basisSet.metadata.analyticalFamily, "generalizedEnergyAPV")
            testCase.verifyEqual(basisSet.metadata.g0, evp.parameters.g0, AbsTol=0)
            testCase.verifyEqual(basisSet.metadata.gd, evp.parameters.gd, AbsTol=0)
            testCase.verifyEqual(basisSet.metadata.surfaceBoundary, "freeSurface")
        end
    end

    methods (Static, Access = private)
        function [N0, zDomain, g, N2] = profile()
            N0 = 5.2e-3;
            zDomain = [-4000 0];
            g = 9.81;
            N2 = @(z) N0*N0*ones(size(z));
        end

        function parameters = parametersFromStrengths(N0, zDomain, g, surfaceStrength, bottomStrength, surfaceConvention)
            IN = N0*N0*diff(zDomain);
            betaSurface = surfaceStrength/IN;
            betaBottom = bottomStrength/IN;
            if surfaceConvention == "freeSurface"
                denominator = betaSurface - 1/g;
            else
                denominator = betaSurface;
            end
            if denominator == 0
                g0 = Inf;
            else
                g0 = 1/denominator;
            end
            if betaBottom == 0
                gd = Inf;
            else
                gd = 1/betaBottom;
            end
            parameters = struct("g0", g0, "gd", gd, "betaSurface", betaSurface, "betaBottom", betaBottom);
        end

        function g0 = g0FromStrength(N0, zDomain, g, surfaceStrength, surfaceConvention)
            parameters = IMConstantGeostrophicAPVAnalyticalTests.parametersFromStrengths(N0, zDomain, g, surfaceStrength, 0, surfaceConvention);
            g0 = parameters.g0;
        end

        function roots = independentRoots(N0, zDomain, betaSurface, betaBottom, surfaceKind, bottomKind, branch, count)
            residual = @(root) IMConstantGeostrophicAPVAnalyticalTests.independentDeterminant(N0, zDomain, betaSurface, betaBottom, surfaceKind, bottomKind, root, branch);
            if branch == "positive"
                samples = linspace(1e-7,max(20,(count + 4)*pi),30000).';
            else
                samples = linspace(1e-7,20,30000).';
            end
            values = arrayfun(residual,samples);
            roots = zeros(1,0);
            for iSample = 1:(length(samples)-1)
                if sign(values(iSample)) ~= sign(values(iSample+1))
                    roots(end+1) = fzero(residual,[samples(iSample) samples(iSample+1)]); %#ok<AGROW>
                    if length(roots) == count
                        break;
                    end
                end
            end
            roots = unique(round(roots,12),"stable");
        end

        function value = independentDeterminant(N0, zDomain, betaSurface, betaBottom, surfaceKind, bottomKind, root, branch)
            depth = diff(zDomain);
            N2 = N0*N0;
            if branch == "positive"
                bottomValue = [1 0];
                bottomFlux = root/(depth*N2)*[0 1];
                surfaceValue = [cos(root) sin(root)];
                surfaceFlux = root/(depth*N2)*[-sin(root) cos(root)];
            else
                bottomValue = [1 0];
                bottomFlux = root/(depth*N2)*[0 1];
                surfaceValue = [cosh(root) sinh(root)];
                surfaceFlux = root/(depth*N2)*[sinh(root) cosh(root)];
            end
            if bottomKind == "dirichlet"
                bottomRow = bottomValue;
            else
                bottomRow = bottomFlux - betaBottom*bottomValue;
            end
            if surfaceKind == "dirichlet"
                surfaceRow = surfaceValue;
            else
                surfaceRow = surfaceFlux + betaSurface*surfaceValue;
            end
            value = det([bottomRow;surfaceRow]);
        end
    end
end
