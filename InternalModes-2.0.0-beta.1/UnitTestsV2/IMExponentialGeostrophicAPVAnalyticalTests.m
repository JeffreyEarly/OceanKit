classdef IMExponentialGeostrophicAPVAnalyticalTests < matlab.unittest.TestCase

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
        function localizedSurfaceRootsMatchSpectralModes(testCase)
            N0=5.2e-3; b=1300; g=9.81; zDomain=[-4000 0];
            N2=@(z)N0^2*exp(2*z/b); gd=integral(N2,zDomain(1),0);
            solution=IMExponentialStratificationSolution(N0=N0,b=b,zDomain=zDomain,g=g);
            for surfaceLength=[5 10 20 50]
                g0=-N0^2*surfaceLength;
                evp=IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=zDomain,g=g,g0=g0,gd=gd,surfaceBoundary="freeSurface");
                exact=solution.internalModes(evp,nModes=3);
                numerical=IMSolverSpectral(nEVP=257,coordinateKind="z").solveEVP(evp,nModes=3);
                testCase.verifyEqual(nnz(exact.h<0),1)
                testCase.verifyEqual(exact.h,numerical.h,RelTol=2e-6)
                F=exact.F(0,normalization="uMax"); Fz=exact.uz(0,normalization="uMax");
                residual=Fz/N0^2+(1/g+1/g0)*F;
                scale=abs(Fz/N0^2)+abs((1/g+1/g0)*F);
                testCase.verifyLessThan(max(abs(residual)./scale),1e-10)
                testCase.verifyLessThan(max(abs(exact.metadata.rootResiduals)),1e-10)
            end
        end
        function canonicalBenchmarkMatchesIndependentRootEquations(testCase)
            [N0, b, zDomain, g, N2] = testCase.canonicalProfile();
            g0 = -N0*N0*b;
            evp = IMInternalModes.geostrophicAPVModes( ...
                N2=N2, zDomain=zDomain, g=g, g0=g0, gd=Inf, ...
                surfaceBoundary="freeSurface");
            solution = IMExponentialStratificationSolution( ...
                N0=N0, b=b, zDomain=zDomain, g=g);

            basisSet = solution.internalModes(evp, nModes=4);
            expectedH = testCase.canonicalEigendepths(N0, b, zDomain, g, g0);

            testCase.verifyTrue(ismember("depth", basisSet.normalizationNames()))
            testCase.verifyEqual(basisSet.h, expectedH, RelTol=5e-11)
            testCase.verifyEqual(basisSet.modeNumber, [-1 1 2 3])
            testCase.verifyEqual(string(basisSet.normalization), "depth")
            testCase.verifyEqual(basisSet.metadata.g0, g0, AbsTol=0)
            testCase.verifyEqual(basisSet.metadata.gd, Inf, AbsTol=0)
            testCase.verifyEqual(basisSet.metadata.surfaceBoundary, "freeSurface")
            testCase.verifyEqual(basisSet.metadata.modeKinds, ...
                ["negative" "positive" "positive" "positive"])
            testCase.verifyLessThan(max(abs(basisSet.metadata.rootResiduals)), 1e-10)
        end

        function endpointInertiaControlsBothSurfaceConventions(testCase)
            [N0, b, zDomain, g, N2] = testCase.canonicalProfile();
            solution = IMExponentialStratificationSolution( ...
                N0=N0, b=b, zDomain=zDomain, g=g);
            strengths = [1 1; -4 0; -8 -4; -0.5 1];
            expectedNegative = [0 1 2 0];
            expectedZero = [0 0 0 1];

            for surfaceConvention = ["freeSurface", "rigidLid"]
                for iCase = 1:size(strengths,1)
                    parameters = IMGeostrophicAPVTestSupport.parametersFromStrengths( ...
                        N0=N0, b=b, zDomain=zDomain, g=g, ...
                        surfaceStrength=strengths(iCase,1), ...
                        bottomStrength=strengths(iCase,2), ...
                        surfaceBoundary=surfaceConvention);
                    evp = IMInternalModes.geostrophicAPVModes( ...
                        N2=N2, zDomain=zDomain, g=g, ...
                        g0=parameters.g0, gd=parameters.gd, ...
                        surfaceBoundary=surfaceConvention);
                    basisSet = solution.internalModes(evp, nModes=4);

                    testCase.verifyEqual(nnz(basisSet.h < 0), expectedNegative(iCase))
                    testCase.verifyEqual(nnz(isinf(basisSet.h)), expectedZero(iCase))
                    testCase.verifyTrue(all(diff(basisSet.eigenvalues) >= 0))
                    finiteResiduals = basisSet.metadata.rootResiduals(isfinite(basisSet.metadata.rootResiduals));
                    testCase.verifyLessThan(max(abs(finiteResiduals)), 1e-9)
                    if expectedZero(iCase) == 1
                        testCase.verifyEqual(basisSet.h(expectedNegative(iCase) + 1), Inf)
                        testCase.verifyEqual(basisSet.modeNumber(expectedNegative(iCase) + 1), 0)
                    end
                end
            end
        end

        function zeroAndInfiniteEndpointLimitsAreExplicit(testCase)
            [N0, b, zDomain, g, N2] = testCase.canonicalProfile();
            solution = IMExponentialStratificationSolution( ...
                N0=N0, b=b, zDomain=zDomain, g=g);
            depth = diff(zDomain);
            IN = N0*N0*b*(1 - exp(-2*depth/b))/2;

            freeNeumann = IMInternalModes.geostrophicAPVModes( ...
                N2=N2, zDomain=zDomain, g=g, g0=-g, gd=Inf, ...
                surfaceBoundary="freeSurface");
            rigidNeumann = IMInternalModes.geostrophicAPVModes( ...
                N2=N2, zDomain=zDomain, g=g, g0=Inf, gd=Inf, ...
                surfaceBoundary="rigidLid");
            doubleDirichlet = IMInternalModes.geostrophicAPVModes( ...
                N2=N2, zDomain=zDomain, g=g, g0=0, gd=0, ...
                surfaceBoundary="freeSurface");
            surfaceDirichletThreshold = IMInternalModes.geostrophicAPVModes( ...
                N2=N2, zDomain=zDomain, g=g, g0=0, gd=-IN, ...
                surfaceBoundary="freeSurface");
            rigidBottomDirichletThreshold = IMInternalModes.geostrophicAPVModes( ...
                N2=N2, zDomain=zDomain, g=g, g0=-IN, gd=0, ...
                surfaceBoundary="rigidLid");
            freeInfinite = IMInternalModes.geostrophicAPVModes( ...
                N2=N2, zDomain=zDomain, g=g, g0=Inf, gd=Inf, ...
                surfaceBoundary="freeSurface");

            freeNeumannBasis = solution.internalModes(freeNeumann, nModes=3);
            rigidNeumannBasis = solution.internalModes(rigidNeumann, nModes=3);
            doubleDirichletBasis = solution.internalModes(doubleDirichlet, nModes=3);
            surfaceThresholdBasis = solution.internalModes(surfaceDirichletThreshold, nModes=3);
            bottomThresholdBasis = solution.internalModes(rigidBottomDirichletThreshold, nModes=3);
            freeInfiniteBasis = solution.internalModes(freeInfinite, nModes=3);

            testCase.verifyEqual(freeNeumannBasis.h(1), Inf)
            testCase.verifyEqual(rigidNeumannBasis.h(1), Inf)
            testCase.verifyFalse(any(isinf(doubleDirichletBasis.h)))
            testCase.verifyEqual(surfaceThresholdBasis.h(1), Inf)
            testCase.verifyEqual(bottomThresholdBasis.h(1), Inf)
            testCase.verifyFalse(any(isinf(freeInfiniteBasis.h)))

            z = linspace(zDomain(1), zDomain(2), 17).';
            constantF = freeNeumannBasis.F(z, normalization="uMax");
            constantG = freeNeumannBasis.G(z, normalization="uMax");
            testCase.verifyEqual(constantF(:,1), ones(size(z)), AbsTol=1e-12)
            testCase.verifyEqual(constantG(:,1), zeros(size(z)), AbsTol=1e-12)
        end

        function exactModesSatisfyODEDiagnosticRelationAndEndpoints(testCase)
            [N0, b, zDomain, g, N2] = testCase.canonicalProfile();
            parameters = IMGeostrophicAPVTestSupport.parametersFromStrengths( ...
                N0=N0, b=b, zDomain=zDomain, g=g, ...
                surfaceStrength=-0.5, bottomStrength=1, ...
                surfaceBoundary="freeSurface");
            evp = IMInternalModes.geostrophicAPVModes( ...
                N2=N2, zDomain=zDomain, g=g, ...
                g0=parameters.g0, gd=parameters.gd, ...
                surfaceBoundary="freeSurface");
            solution = IMExponentialStratificationSolution( ...
                N0=N0, b=b, zDomain=zDomain, g=g);
            basisSet = solution.internalModes(evp, nModes=4);
            solver = IMSolverSpectral(nEVP=160, coordinateKind="z");
            solver = solver.configuredForEVP(evp);
            z = solver.zNative;
            F = basisSet.F(z, normalization="uMax");
            G = basisSet.G(z, normalization="uMax");
            Fz = basisSet.uz(z, normalization="uMax");
            flux = Fz ./ N2(z);
            fluxDerivative = solver.differentiateGridValues(flux, 1);
            rhs = F .* reshape(1./(g*basisSet.h), 1, []);
            odeResidual = -fluxDerivative - rhs;
            referenceScale = max(abs(F),[],1)/(min(N2(z))*diff(zDomain)^2);
            odeScale = max([referenceScale; max(abs(fluxDerivative),[],1); max(abs(rhs),[],1)], [], 1);

            testCase.verifyLessThan(max(abs(odeResidual)./odeScale, [], "all"), 2e-8)
            diagnosticScale = max(1e-12, max(abs(G),[],1));
            diagnosticResidual = G + (g./N2(z)).*Fz;
            testCase.verifyLessThan(max(abs(diagnosticResidual)./diagnosticScale, [], "all"), 5e-13)

            [~, surfaceIndex] = max(z);
            [~, bottomIndex] = min(z);
            surfaceResidual = -evp.surfaceBoundary.a*F(surfaceIndex,:) ...
                + evp.surfaceBoundary.b*flux(surfaceIndex,:);
            bottomResidual = -evp.bottomBoundary.a*F(bottomIndex,:) ...
                + evp.bottomBoundary.b*flux(bottomIndex,:);
            endpointScale = max(1e-12, max(abs([F(surfaceIndex,:); F(bottomIndex,:); ...
                flux(surfaceIndex,:); flux(bottomIndex,:)]),[],1));
            testCase.verifyLessThan(max(abs(surfaceResidual)./endpointScale), 2e-11)
            testCase.verifyLessThan(max(abs(bottomResidual)./endpointScale), 2e-11)
        end

        function analyticalModesConvergeToNumericalEVP(testCase)
            [N0, b, zDomain, g, N2] = testCase.canonicalProfile();
            parameters = IMGeostrophicAPVTestSupport.parametersFromStrengths( ...
                N0=N0, b=b, zDomain=zDomain, g=g, ...
                surfaceStrength=-8, bottomStrength=-4, ...
                surfaceBoundary="rigidLid");
            evp = IMInternalModes.geostrophicAPVModes( ...
                N2=N2, zDomain=zDomain, g=g, ...
                g0=parameters.g0, gd=parameters.gd, ...
                surfaceBoundary="rigidLid");
            solution = IMExponentialStratificationSolution( ...
                N0=N0, b=b, zDomain=zDomain, g=g);
            analytical = solution.internalModes(evp, nModes=4);
            lowOrder = IMSolverSpectral(nEVP=48, coordinateKind="z").solveEVP(evp, nModes=4);
            highOrder = IMSolverSpectral(nEVP=96, coordinateKind="z").solveEVP(evp, nModes=4);
            lowError = norm(lowOrder.eigenvalues - analytical.eigenvalues) ...
                /norm(analytical.eigenvalues);
            highError = norm(highOrder.eigenvalues - analytical.eigenvalues) ...
                /norm(analytical.eigenvalues);

            testCase.verifyLessThan(lowError, 2e-8)
            testCase.verifyLessThan(highError, 2e-8)

            z = linspace(zDomain(1), zDomain(2), 801).';
            analyticalF = analytical.F(z, normalization="uMax");
            numericalF = highOrder.F(z, normalization="uMax");
            [analyticalF, ~] = IMGeostrophicAPVTestSupport.orientColumns( ...
                analyticalF, zeros(size(analyticalF)));
            [numericalF, ~] = IMGeostrophicAPVTestSupport.orientColumns( ...
                numericalF, zeros(size(numericalF)));
            subspaceError = IMGeostrophicAPVTestSupport.subspaceError( ...
                z, analyticalF(:,3:4), numericalF(:,3:4));
            testCase.verifyLessThan(subspaceError, 2e-5)
        end

        function negativeModesUsePositiveRealNormalization(testCase)
            [N0, b, zDomain, g, N2] = testCase.canonicalProfile();
            parameters = IMGeostrophicAPVTestSupport.parametersFromStrengths( ...
                N0=N0, b=b, zDomain=zDomain, g=g, ...
                surfaceStrength=-8, bottomStrength=-4, ...
                surfaceBoundary="freeSurface");
            evp = IMInternalModes.geostrophicAPVModes( ...
                N2=N2, zDomain=zDomain, g=g, ...
                g0=parameters.g0, gd=parameters.gd, ...
                surfaceBoundary="freeSurface");
            solution = IMExponentialStratificationSolution( ...
                N0=N0, b=b, zDomain=zDomain, g=g);
            basisSet = solution.internalModes(evp, nModes=4);

            factors = basisSet.normalizationFactors("geostrophic");
            signedGram = basisSet.gramMatrix(variable="G");
            majorantGram = basisSet.majorantGramMatrix(variable="G");
            testCase.verifyTrue(all(isreal(factors)))
            testCase.verifyTrue(all(isfinite(factors)))
            testCase.verifyTrue(all(factors > 0))
            testCase.verifyTrue(all(basisSet.h(1:2) < 0))
            testCase.verifyEqual(nnz(eig(signedGram) < 0),2)
            testCase.verifyGreaterThan(min(eig(majorantGram)),0)
            testCase.verifyEqual(basisSet.majorantNorm(ones(4,1),variable="G"), ...
                sqrt(real(ones(1,4)*majorantGram*ones(4,1))),RelTol=2e-13)
        end

        function mismatchedStratificationIsRejected(testCase)
            [N0, b, zDomain, g, ~] = testCase.canonicalProfile();
            wrongN2 = @(z) N0*N0*ones(size(z));
            evp = IMInternalModes.geostrophicAPVModes( ...
                N2=wrongN2, zDomain=zDomain, g=g, g0=-g, gd=Inf, ...
                surfaceBoundary="freeSurface");
            solution = IMExponentialStratificationSolution( ...
                N0=N0, b=b, zDomain=zDomain, g=g);

            testCase.verifyError(@() solution.internalModes(evp,nModes=4),"IMExponentialStratificationSolution:UnsupportedStratification")
        end
    end

    methods (Static, Access = private)
        function [N0, b, zDomain, g, N2] = canonicalProfile()
            N0 = 5.2e-3;
            b = 1300;
            zDomain = [-4000 0];
            g = 9.81;
            N2 = @(z) N0*N0*exp(2*z/b);
        end

        function h = canonicalEigendepths(N0, b, zDomain, g, g0)
            depth = diff(zDomain);
            betaSurface = 1/g + 1/g0;
            bottomFactor = exp(-depth/b);
            positiveResidual = @(s) (s/b).*( ...
                bessely(0,bottomFactor*s).*besselj(0,s) ...
                - besselj(0,bottomFactor*s).*bessely(0,s)) ...
                + N0*N0*betaSurface.*( ...
                bessely(0,bottomFactor*s).*besselj(1,s) ...
                - besselj(0,bottomFactor*s).*bessely(1,s));
            negativeResidual = @(s) (s/b).*( ...
                besselk(0,bottomFactor*s).*besseli(0,s) ...
                - besseli(0,bottomFactor*s).*besselk(0,s)) ...
                + N0*N0*betaSurface.*( ...
                besselk(0,bottomFactor*s).*besseli(1,s) ...
                + besseli(0,bottomFactor*s).*besselk(1,s));
            negativeRoot = fzero(negativeResidual, [0.6 0.8]);
            positiveRoots = [ ...
                fzero(positiveResidual, [2 3]) ...
                fzero(positiveResidual, [6 7]) ...
                fzero(positiveResidual, [9 10])];
            h = [-(b*N0/negativeRoot)^2/g (b*N0./positiveRoots).^2/g];
        end
    end
end
