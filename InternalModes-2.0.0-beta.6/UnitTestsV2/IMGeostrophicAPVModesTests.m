classdef IMGeostrophicAPVModesTests < matlab.unittest.TestCase

    properties (Access = private)
        originalPath
    end

    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            repoRoot = fileparts(fileparts(mfilename("fullpath")));
            testCase.originalPath = path;
            addpath(repoRoot);
        end
    end

    methods (TestClassTeardown)
        function restoreRepositoryPath(testCase)
            path(testCase.originalPath);
        end
    end

    methods (Test)
        function factoryMapsEveryEndpointLimitAndStoresMetadata(testCase)
            [N2, zDomain, g] = testCase.profile();
            g0Values = [2 -3 -g 0 Inf];
            gdValues = [4 -5 0 Inf];
            for surfaceConvention = ["freeSurface", "rigidLid"]
                for g0 = g0Values
                    for gd = gdValues
                        evp = IMInternalModes.geostrophicAPVModes(N2=N2, zDomain=zDomain, g=g, g0=g0, gd=gd, surfaceBoundary=surfaceConvention);
                        expectedSurface = testCase.expectedSurfaceCoefficients(g, g0, surfaceConvention);
                        expectedBottom = testCase.expectedBottomCoefficients(gd);

                        testCase.verifyEqual(string(evp.name), "geostrophicAPVModes")
                        testCase.verifyEqual(evp.formulation, "F")
                        testCase.verifyEqual(evp.modeFamily, "hydrostatic")
                        testCase.verifyEqual([evp.surfaceBoundary.a evp.surfaceBoundary.b evp.surfaceBoundary.c evp.surfaceBoundary.d], expectedSurface, AbsTol=10*eps)
                        testCase.verifyEqual([evp.bottomBoundary.a evp.bottomBoundary.b evp.bottomBoundary.c evp.bottomBoundary.d], expectedBottom, AbsTol=10*eps)
                        testCase.verifyEqual(evp.parameters.g0, g0, AbsTol=0)
                        testCase.verifyEqual(evp.parameters.gd, gd, AbsTol=0)
                        testCase.verifyEqual(string(evp.parameters.surfaceBoundary), surfaceConvention)
                    end
                end
            end
        end

        function factoryCoefficientsAndDiagnosticRelationAreCanonical(testCase)
            [N2, zDomain, g] = testCase.profile();
            evp = IMInternalModes.geostrophicAPVModes(N2=N2, zDomain=zDomain, g=g, g0=2, gd=3);
            z = linspace(zDomain(1), zDomain(2), 7).';
            context = struct("N2", N2, "g", g);
            p = IMEigenvalueProblem.evaluateCoefficient(evp.p, z, context);
            q = IMEigenvalueProblem.evaluateCoefficient(evp.q, z, context);
            r = IMEigenvalueProblem.evaluateCoefficient(evp.r, z, context);
            Fz = reshape(1:numel(z),[],1);

            testCase.verifyEqual(p, 1./N2(z), RelTol=1e-14)
            testCase.verifyEqual(q, zeros(size(z)), AbsTol=0)
            testCase.verifyEqual(r, ones(size(z))/g, RelTol=1e-14)
            testCase.verifyEqual(evp.GfromFz(z, Fz, 1, context), -(g./N2(z)).*Fz, RelTol=1e-14)
        end

        function endpointAccelerationsAreValidatedAndF0IsNotAnOption(testCase)
            [N2, zDomain] = testCase.profile();
            testCase.verifyError(@() IMInternalModes.geostrophicAPVModes(N2=N2, zDomain=zDomain, g0=NaN, gd=1), "IMInternalModes:InvalidSurfaceAcceleration")
            testCase.verifyError(@() IMInternalModes.geostrophicAPVModes(N2=N2, zDomain=zDomain, g0=-Inf, gd=1), "IMInternalModes:InvalidSurfaceAcceleration")
            testCase.verifyError(@() IMInternalModes.geostrophicAPVModes(N2=N2, zDomain=zDomain, g0=1, gd=NaN), "IMInternalModes:InvalidBottomAcceleration")
            testCase.verifyError(@() IMInternalModes.geostrophicAPVModes(N2=N2, zDomain=zDomain, g0=1, gd=-Inf), "IMInternalModes:InvalidBottomAcceleration")
            testCase.verifyError(@() IMInternalModes.geostrophicAPVModes(N2=N2, zDomain=zDomain, g0=1, gd=1, f0=1e-4), "MATLAB:TooManyInputs")
        end

        function innerProductCatalogMatchesPhysicalEndpointEnergy(testCase)
            [N2, zDomain, g] = testCase.profile();
            cases = {
                "freeSurface", 2, 3
                "freeSurface", -3, -4
                "freeSurface", -g, Inf
                "freeSurface", 0, 0
                "freeSurface", Inf, Inf
                "rigidLid", 2, 3
                "rigidLid", -3, -4
                "rigidLid", 0, 0
                "rigidLid", Inf, Inf};
            z = linspace(zDomain(1), zDomain(2), 5).';
            context = struct("N2", N2, "g", g);
            for iCase = 1:size(cases,1)
                surfaceConvention = string(cases{iCase,1});
                g0 = cases{iCase,2};
                gd = cases{iCase,3};
                evp = IMInternalModes.geostrophicAPVModes(N2=N2, zDomain=zDomain, g=g, g0=g0, gd=gd, surfaceBoundary=surfaceConvention);
                FSpec = evp.innerProduct("F");
                GSpec = evp.innerProduct("G");

                testCase.verifyTrue(FSpec.hasInnerProduct)
                testCase.verifyTrue(GSpec.hasInnerProduct)
                testCase.verifyEmpty(FSpec.endpointInnerProductTerms)
                testCase.verifyEqual(FSpec.interiorWeight(z, context), ones(size(z)), AbsTol=0)
                testCase.verifyEqual(GSpec.interiorWeight(z, context), N2(z)/g, RelTol=1e-14)

                expectedSurface = testCase.expectedGSurfaceCoefficient(g, g0, surfaceConvention);
                expectedBottom = testCase.expectedGBottomCoefficient(g, gd);
                testCase.verifyEndpointCoefficient(GSpec.endpointInnerProductTerms, "surface", expectedSurface)
                testCase.verifyEndpointCoefficient(GSpec.endpointInnerProductTerms, "bottom", expectedBottom)
            end
        end

        function exactModesSatisfyIndependentPhysicalEndpointEquations(testCase)
            [N2, zDomain, g, N0] = testCase.profile();
            cases = {
                "freeSurface", 2, 3
                "freeSurface", -g, Inf
                "freeSurface", 0, -4
                "freeSurface", Inf, 0
                "rigidLid", -3, 4
                "rigidLid", Inf, Inf
                "rigidLid", 0, 0};
            solution = IMConstantStratificationSolution(N0=N0, zDomain=zDomain, g=g);
            for iCase = 1:size(cases,1)
                surfaceConvention = string(cases{iCase,1});
                g0 = cases{iCase,2};
                gd = cases{iCase,3};
                evp = IMInternalModes.geostrophicAPVModes(N2=N2, zDomain=zDomain, g=g, g0=g0, gd=gd, surfaceBoundary=surfaceConvention);
                basisSet = solution.internalModes(evp, nModes=3);
                FSurface = basisSet.F(zDomain(2));
                GSurface = basisSet.G(zDomain(2));
                FBottom = basisSet.F(zDomain(1));
                GBottom = basisSet.G(zDomain(1));
                surfaceResidual = testCase.physicalSurfaceResidual(FSurface, GSurface, g, g0, surfaceConvention);
                bottomResidual = testCase.physicalBottomResidual(FBottom, GBottom, g, gd);
                scale = max(1,max(abs([FSurface GSurface FBottom GBottom]),[],2));

                testCase.verifyLessThan(max(abs(surfaceResidual)./scale), 2e-11)
                testCase.verifyLessThan(max(abs(bottomResidual)./scale), 2e-11)
                testCase.verifyEqual(basisSet.metadata.g0, g0, AbsTol=0)
                testCase.verifyEqual(basisSet.metadata.gd, gd, AbsTol=0)
                testCase.verifyEqual(basisSet.metadata.surfaceBoundary, surfaceConvention)
            end
        end

        function depthNormalizationIsDefaultPositiveAndCoupled(testCase)
            [N2, zDomain, g] = testCase.profile();
            evp = IMInternalModes.geostrophicAPVModes(N2=N2, zDomain=zDomain, g=g, g0=-0.02, gd=Inf);
            solver = IMSolverSpectral(nEVP=64, coordinateKind="z");
            basisSet = solver.solveEVP(evp, nModes=4);
            factors = basisSet.normalizationFactors(Normalization.depth);
            z = basisSet.solver.innerProductGrid(zDomain);
            F = basisSet.F(z);
            gram = zeros(4,4);
            for iMode = 1:4
                for jMode = iMode:4
                    value = basisSet.solver.integrateInnerProduct(z, F(:,iMode).*F(:,jMode), zDomain)/diff(zDomain);
                    gram(iMode,jMode) = value;
                    gram(jMode,iMode) = value;
                end
            end

            rawBasis = basisSet.addNormalization("raw", @(~,~) 1);
            rawF = rawBasis.F(z, normalization="raw");
            rawG = rawBasis.G(z, normalization="raw");
            testCase.verifyEqual(string(basisSet.normalization), "depth")
            testCase.verifyTrue(ismember("depth", basisSet.normalizationNames()))
            testCase.verifyTrue(all(isfinite(factors) & factors > 0))
            testCase.verifyTrue(any(basisSet.h < 0))
            testCase.verifyEqual(gram, eye(4), RelTol=2e-7, AbsTol=2e-7)
            testCase.verifyEqual(F.*factors, rawF, RelTol=2e-12, AbsTol=2e-12)
            testCase.verifyEqual(basisSet.G(z).*factors, rawG, RelTol=2e-12, AbsTol=2e-12)
            testCase.verifyEqual(basisSet.metadata.g0, evp.parameters.g0, AbsTol=0)
            testCase.verifyEqual(basisSet.metadata.gd, evp.parameters.gd, AbsTol=0)
            testCase.verifyEqual(basisSet.metadata.surfaceBoundary, "freeSurface")
        end

        function numericalAndAnalyticalMajorantsGivePositiveStateNorms(testCase)
            [N2,zDomain,g] = testCase.profile();
            evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=zDomain,g=g,g0=-0.02,gd=Inf);
            numerical = IMSolverSpectral(nEVP=128,coordinateKind="z").solveEVP(evp,nModes=4);
            analytical = IMConstantStratificationSolution(N0=sqrt(N2(0)),zDomain=zDomain,g=g).internalModes(evp,nModes=4);
            signedGram = numerical.gramMatrix(variable="G");
            majorantGram = numerical.majorantGramMatrix(variable="G");
            analyticalMajorant = analytical.majorantGramMatrix(variable="G");
            numericalRecipe = numerical.majorantInnerProduct(variable="G");
            analyticalRecipe = analytical.majorantInnerProduct(variable="G");
            [vectors,values] = eig(0.5*(signedGram+signedGram.'));
            signedEigenvalues = diag(values);
            negativeIndex = find(signedEigenvalues < 0,1);
            positiveIndex = find(signedEigenvalues > 0,1);
            cancellation = vectors(:,negativeIndex)/sqrt(-signedEigenvalues(negativeIndex)) ...
                + vectors(:,positiveIndex)/sqrt(signedEigenvalues(positiveIndex));
            expectedNorm = sqrt(real(cancellation'*(majorantGram*cancellation)));

            testCase.assertNotEmpty(negativeIndex)
            testCase.assertNotEmpty(positiveIndex)
            testCase.verifyEqual(numericalRecipe.kind,"inducedHilbertMajorant")
            testCase.verifyEqual(analyticalRecipe.kind,"inducedHilbertMajorant")
            testCase.verifyGreaterThan(min(eig(majorantGram)),0)
            testCase.verifyEqual(majorantGram,analyticalMajorant,RelTol=1e-5,AbsTol=1e-5)
            testCase.verifyLessThan(abs(real(cancellation'*(signedGram*cancellation))),1e-10)
            testCase.verifyEqual(numerical.majorantNorm(cancellation,variable="G"),expectedNorm,RelTol=2e-13)
            testCase.verifyGreaterThan(expectedNorm,0)

            positiveEVP = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=zDomain,g=g,g0=0.02,gd=0.03,surfaceBoundary="rigidLid");
            positiveBasis = IMSolverSpectral(nEVP=96).solveEVP(positiveEVP,nModes=4);
            testCase.verifyEqual(positiveBasis.majorantGramMatrix(variable="G"), ...
                positiveBasis.gramMatrix(variable="G"),RelTol=2e-10,AbsTol=2e-10)
        end

        function otherHydrostaticBasesGainDepthWithoutDefaultChange(testCase)
            [N2, zDomain] = testCase.profile();
            evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain);
            numerical = IMSolverSpectral(nEVP=48).solveEVP(evp, nModes=2);
            analytical = IMConstantStratificationSolution(N0=sqrt(N2(0)), zDomain=zDomain).internalModes(evp, nModes=2);

            testCase.verifyTrue(ismember("depth", numerical.normalizationNames()))
            testCase.verifyTrue(ismember("depth", analytical.normalizationNames()))
            testCase.verifyEqual(string(numerical.normalization), "geostrophic")
            testCase.verifyEqual(string(analytical.normalization), "geostrophic")
        end

        function numericalZeroModeIsRetainedExactly(testCase)
            [N2, zDomain, g] = testCase.profile();
            evp = IMInternalModes.geostrophicAPVModes(N2=N2, zDomain=zDomain, g=g, g0=Inf, gd=Inf, surfaceBoundary="rigidLid");
            basisSet = IMSolverSpectral(nEVP=48, coordinateKind="z").solveEVP(evp, nModes=4);
            zeroIndex = find(basisSet.modeNumber == 0);

            testCase.verifyNumElements(zeroIndex, 1)
            testCase.verifyEqual(basisSet.eigenvalues(zeroIndex), 0, AbsTol=0)
            testCase.verifyEqual(basisSet.h(zeroIndex), Inf, AbsTol=0)
            testCase.verifyTrue(all(diff(basisSet.eigenvalues) >= 0))
        end

        function defaultExponentialNegativeModeConvergesThroughHighResolution(testCase)
            D = 4000;
            N0 = 5.2e-3;
            b = 1300;
            N2 = @(z) N0*N0*exp(2*z/b);
            g0 = -integral(N2,-D,0);
            evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=[-D 0],g0=g0,gd=Inf);
            expected = -0.393119878698696;

            for nEVP = [128 256 512]
                basisSet = IMSolverSpectral(nEVP=nEVP).solveEVP(evp,nModes=4);
                testCase.verifyEqual(basisSet.modeNumber(1),-1,AbsTol=0)
                testCase.verifyEqual(basisSet.eigenvalues(1),expected,RelTol=1e-5)
                testCase.verifyGreaterThan(basisSet.eigenvalues(2),0)
                testCase.verifyLessThan(max(abs(basisSet.eigenvalues)),1e3)
            end
        end
    end

    methods
        function verifyEndpointCoefficient(testCase, terms, location, expected)
            locations = string({terms.location});
            index = find(locations == string(location));
            if isnan(expected)
                testCase.verifyEmpty(index)
            else
                testCase.verifyNumElements(index, 1)
                testCase.verifyEqual(terms(index).coefficient, expected, RelTol=1e-13, AbsTol=1e-13)
                testCase.verifyEqual(string(terms(index).variable), "G")
            end
        end
    end

    methods (Static, Access = private)
        function [N2, zDomain, g, N0] = profile()
            N0 = 5.2e-3;
            zDomain = [-4000 0];
            g = 9.81;
            N2 = @(z) N0*N0*ones(size(z));
        end

        function coefficients = expectedSurfaceCoefficients(g, g0, surfaceConvention)
            if g0 == 0
                coefficients = [1 0 0 0];
            elseif surfaceConvention == "freeSurface"
                if isinf(g0)
                    coefficients = [-1/g 1 0 0];
                else
                    coefficients = [-(1/g + 1/g0) 1 0 0];
                end
            elseif isinf(g0)
                coefficients = [0 1 0 0];
            else
                coefficients = [-1/g0 1 0 0];
            end
        end

        function coefficients = expectedBottomCoefficients(gd)
            if gd == 0
                coefficients = [1 0 0 0];
            elseif isinf(gd)
                coefficients = [0 1 0 0];
            else
                coefficients = [1/gd 1 0 0];
            end
        end

        function coefficient = expectedGSurfaceCoefficient(g, g0, surfaceConvention)
            coefficient = NaN;
            if g0 == 0 || (surfaceConvention == "freeSurface" && g0 == -g) || (surfaceConvention == "rigidLid" && isinf(g0))
                return;
            end
            if surfaceConvention == "freeSurface"
                if isinf(g0)
                    coefficient = 1;
                else
                    coefficient = g0/(g + g0);
                end
            elseif ~isinf(g0)
                coefficient = g0/g;
            end
        end

        function coefficient = expectedGBottomCoefficient(g, gd)
            coefficient = NaN;
            if gd ~= 0 && ~isinf(gd)
                coefficient = gd/g;
            end
        end

        function residual = physicalSurfaceResidual(F, G, g, g0, surfaceConvention)
            if g0 == 0
                residual = F;
            elseif surfaceConvention == "freeSurface"
                if isinf(g0)
                    residual = G - F;
                else
                    residual = G - (1 + g/g0)*F;
                end
            elseif isinf(g0)
                residual = G;
            else
                residual = G - (g/g0)*F;
            end
        end

        function residual = physicalBottomResidual(F, G, g, gd)
            if gd == 0
                residual = F;
            elseif isinf(gd)
                residual = G;
            else
                residual = G + (g/gd)*F;
            end
        end
    end
end
