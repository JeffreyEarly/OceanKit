classdef IMInternalModesInnerProductCatalogTests < matlab.unittest.TestCase

    properties (Access = private)
        originalPath
    end

    methods (TestClassSetup)
        function addRepositoryPath(testCase)
            repoRoot = fileparts(fileparts(mfilename('fullpath')));
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
        function simpleHydrostaticDiagnosticRowsAreKnown(testCase)
            [N2, zDomain] = testCase.profile();
            gEVP = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain);
            fEVP = IMInternalModes.hydrostaticFModes(N2=N2, zDomain=zDomain);

            fSpec = gEVP.innerProduct("F");
            gSpec = fEVP.innerProduct("G");

            testCase.verifyTrue(fSpec.hasInnerProduct)
            testCase.verifyTrue(gSpec.hasInnerProduct)
            testCase.verifyEmpty(fSpec.endpointInnerProductTerms)
            testCase.verifyEmpty(gSpec.endpointInnerProductTerms)
            testCase.verifyFalse(isfield(fSpec, "status"))
            testCase.verifyFalse(isfield(fSpec, "endpointFunctionals"))
            testCase.verifyTrue(contains(fSpec.reason, "interior-only"))
        end

        function diagnosticEndpointValueRowsAreKnown(testCase)
            [N2, zDomain, ~, g] = testCase.profile();
            gEVP = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=2, b=3));
            fRobinEVP = IMInternalModes.hydrostaticFModes(N2=N2, zDomain=zDomain, g=g, ...
                surfaceBoundary=IMBoundaryCondition(a=2, b=3));
            fLinkedEVP = IMInternalModes.hydrostaticFModes(N2=N2, zDomain=zDomain, g=g, ...
                surfaceBoundary=IMBoundaryCondition(a=2, b=4, c=1, d=0));

            fSpec = gEVP.innerProduct("F");
            gRobinSpec = fRobinEVP.innerProduct("G");
            gLinkedSpec = fLinkedEVP.innerProduct("G");

            testCase.verifyTrue(fSpec.hasInnerProduct)
            testCase.verifyEqual(fSpec.endpointInnerProductTerms(1).catalogCase, "G-P5")
            testCase.verifyEqual(fSpec.endpointInnerProductTerms(1).variable, "F")
            testCase.verifyEqual(fSpec.endpointInnerProductTerms(1).coefficient, -3/2, AbsTol=0)
            testCase.verifyTrue(gRobinSpec.hasInnerProduct)
            testCase.verifyEqual(gRobinSpec.endpointInnerProductTerms(1).catalogCase, "F-P2")
            testCase.verifyEqual(gRobinSpec.endpointInnerProductTerms(1).variable, "G")
            testCase.verifyEqual(gRobinSpec.endpointInnerProductTerms(1).coefficient, -3/(2*g), RelTol=1e-12)
            testCase.verifyTrue(gLinkedSpec.hasInnerProduct)
            testCase.verifyEqual(gLinkedSpec.endpointInnerProductTerms(1).catalogCase, "F-T4")
            testCase.verifyEqual(gLinkedSpec.endpointInnerProductTerms(1).variable, "F")
            testCase.verifyEqual(gLinkedSpec.endpointInnerProductTerms(1).coefficient, -g/2, RelTol=1e-12)
        end

        function majorantRecipeUsesAbsoluteSolvedAndCatalogEndpointCoefficients(testCase)
            [N2,zDomain,~,g] = testCase.profile();
            catalogEVP = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=2,b=3));
            signedCatalog = catalogEVP.innerProduct("F");
            majorantCatalog = catalogEVP.majorantInnerProduct("F");
            apvEVP = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=zDomain,g=g, ...
                g0=-0.02,gd=-0.03,surfaceBoundary="rigidLid");
            signedAPV = apvEVP.innerProduct("F");
            majorantAPV = apvEVP.majorantInnerProduct("F");

            testCase.verifyEqual(signedCatalog.kind,"signedPontryagin")
            testCase.verifyEqual(majorantCatalog.kind,"inducedHilbertMajorant")
            testCase.verifyEqual([majorantCatalog.endpointInnerProductTerms.coefficient], ...
                abs([signedCatalog.endpointInnerProductTerms.coefficient]),AbsTol=0)
            testCase.verifyEqual([majorantAPV.surfaceWeights.coefficient],abs([signedAPV.surfaceWeights.coefficient]),AbsTol=0)
            testCase.verifyEqual([majorantAPV.bottomWeights.coefficient],abs([signedAPV.bottomWeights.coefficient]),AbsTol=0)
            z = linspace(zDomain(1),zDomain(2),5).';
            context.N2 = N2;
            context.g = g;
            testCase.verifyEqual(majorantAPV.interiorWeight(z,context),signedAPV.interiorWeight(z,context),AbsTol=0)
        end

        function diagnosticInteriorRowsMatchCatalog(testCase)
            [N2, zDomain] = testCase.profile();
            fP6 = IMInternalModes.hydrostaticFModes(N2=N2, zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=0, b=2, c=3, d=0));
            gP2 = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=0, b=2, c=3, d=0));
            gP4 = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=0, b=2, c=0, d=3));

            testCase.verifyInteriorOnly(fP6.innerProduct("G"))
            testCase.verifyInteriorOnly(gP2.innerProduct("F"))
            testCase.verifyInteriorOnly(gP4.innerProduct("F"))
        end

        function gFormGEqualsAFBoundaryHasDiagnosticFInnerProduct(testCase)
            [N2, zDomain] = testCase.profile();
            A = 7;
            evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=0, b=A, c=1, d=0));

            testCase.verifyInteriorOnly(evp.innerProduct("F"))
        end

        function endpointInnerProductTermsEnterGramMatrix(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=2, b=3));
            basisSet = solver.solveEVP(evp, nModes=2);
            spec = evp.innerProduct("F");

            z = solver.configuredForEVP(evp).innerProductGrid(zDomain);
            F = basisSet.F(z);
            gramInterior = zeros(2,2);
            for iMode = 1:2
                for jMode = iMode:2
                    value = basisSet.solver.integrateInnerProduct(z, F(:,iMode).*F(:,jMode), zDomain);
                    gramInterior(iMode,jMode) = value;
                    gramInterior(jMode,iMode) = value;
                end
            end
            FSurface = basisSet.F(zDomain(2));
            endpointTerms = basisSet.endpointGramTerms(variable="F");
            rawEndpointTerms = basisSet.endpointGramTerms(variable="F", useNormalized=false);
            factors = basisSet.normalizationFactors(basisSet.normalization);
            expected = gramInterior + endpointTerms(1).coefficient*(endpointTerms(1).values(:)*endpointTerms(1).values(:).');
            partialBounds = [zDomain(1) mean(zDomain)];
            partialEndpointTerms = basisSet.endpointGramTerms(variable="F", zBounds=partialBounds);
            zPartial = solver.configuredForEVP(evp).innerProductGrid(partialBounds);
            FPartial = basisSet.F(zPartial);
            expectedPartial = zeros(2,2);
            for iMode = 1:2
                for jMode = iMode:2
                    value = basisSet.solver.integrateInnerProduct(zPartial, FPartial(:,iMode).*FPartial(:,jMode), partialBounds);
                    expectedPartial(iMode,jMode) = value;
                    expectedPartial(jMode,iMode) = value;
                end
            end

            testCase.verifyNumElements(endpointTerms, 1)
            testCase.verifyEqual(endpointTerms(1).kind, "endpointInnerProductTerm")
            testCase.verifyEqual(endpointTerms(1).location, "surface")
            testCase.verifyEqual(endpointTerms(1).coefficient, spec.endpointInnerProductTerms(1).coefficient, AbsTol=0)
            testCase.verifyEqual(endpointTerms(1).values, FSurface, RelTol=1e-12, AbsTol=1e-12)
            testCase.verifyEqual(endpointTerms(1).values, rawEndpointTerms(1).values./factors, RelTol=1e-12, AbsTol=1e-12)
            testCase.verifyEmpty(partialEndpointTerms)
            testCase.verifyEqual(basisSet.gramMatrix(variable="F"), expected, RelTol=1e-10, AbsTol=1e-10)
            testCase.verifyEqual(basisSet.gramMatrix(zBounds=partialBounds, variable="F"), ...
                expectedPartial, RelTol=1e-10, AbsTol=1e-10)
        end

        function linkedEndpointInnerProductTermsEnterGramMatrix(testCase)
            [N2, zDomain, nEVP, g] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMInternalModes.hydrostaticFModes(N2=N2, zDomain=zDomain, g=g, ...
                surfaceBoundary=IMBoundaryCondition(a=2, b=4, c=1, d=0));
            basisSet = solver.solveEVP(evp, nModes=2);
            spec = evp.innerProduct("G");

            z = solver.configuredForEVP(evp).innerProductGrid(zDomain);
            G = basisSet.G(z);
            weight = N2(z)/g;
            gramInterior = zeros(2,2);
            for iMode = 1:2
                for jMode = iMode:2
                    value = basisSet.solver.integrateInnerProduct(z, weight.*G(:,iMode).*G(:,jMode), zDomain);
                    gramInterior(iMode,jMode) = value;
                    gramInterior(jMode,iMode) = value;
                end
            end
            FSurface = basisSet.F(zDomain(2));
            expected = gramInterior + spec.endpointInnerProductTerms(1).coefficient*(FSurface.'*FSurface);

            testCase.verifyEqual(basisSet.gramMatrix(variable="G"), expected, RelTol=1e-10, AbsTol=1e-10)
        end

        function unsupportedHydrostaticRowsHaveNoInnerProduct(testCase)
            [N2, zDomain] = testCase.profile();
            fP5 = IMInternalModes.hydrostaticFModes(N2=N2, zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=1, b=0, c=0, d=1));
            fT3 = IMInternalModes.hydrostaticFModes(N2=N2, zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=1, b=1, c=0, d=1));
            fT2 = IMInternalModes.hydrostaticFModes(N2=N2, zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=1, b=0, c=1, d=1));
            gP1 = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=1, b=0, c=0, d=1));
            gT1 = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=1, b=0, c=1, d=1));
            gT2 = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=1, b=1, c=0, d=1));
            gT3 = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=1, b=1, c=1, d=0));
            allFour = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=1, b=1, c=1, d=1));

            testCase.verifyUnavailable(fP5.innerProduct("G"), "mixed identity")
            testCase.verifyUnavailable(fT3.innerProduct("G"), "mixed identity")
            testCase.verifyUnavailable(gT3.innerProduct("F"), "derivative endpoint terms")
            testCase.verifyUnavailable(gP1.innerProduct("F"), "eigenvalue-dependent")
            testCase.verifyUnavailable(gT1.innerProduct("F"), "eigenvalue-dependent")
            testCase.verifyUnavailable(gT2.innerProduct("F"), "eigenvalue-dependent")
            testCase.verifyUnavailable(fT2.innerProduct("G"), "eigenvalue-dependent")
            testCase.verifyUnavailable(allFour.innerProduct("F"), "outside the value-only")
        end

        function unavailableDiagnosticGramMatrixThrows(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain, ...
                surfaceBoundary=IMBoundaryCondition(a=1, b=0, c=0, d=1));
            basisSet = solver.solveEVP(evp, nModes=2);

            testCase.verifyError(@() basisSet.gramMatrix(variable="F"), "IMInternalModesBasis:UnavailableInnerProduct")
        end

        function waveDiagnosticInnerProductsRemainUnavailable(testCase)
            [N2, zDomain, nEVP, g] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMInternalModes.waveModesAtFrequency(N2=N2, zDomain=zDomain, omega=1e-3, g=g);
            basisSet = solver.solveEVP(evp, nModes=2);
            spec = evp.innerProduct("F");

            testCase.verifyFalse(spec.hasInnerProduct)
            testCase.verifyTrue(contains(spec.reason, "hydrostatic"))
            testCase.verifyError(@() basisSet.gramMatrix(variable="F"), "IMInternalModesBasis:UnavailableInnerProduct")
        end

        function customHydrostaticFamilyUsesCatalogWithoutHydrostaticName(testCase)
            [N2, zDomain, ~, g] = testCase.profile();
            surface = IMBoundaryCondition(a=2, b=3);
            bottom = IMBoundaryCondition(a=-4, b=5);
            evp = IMInternalModes(name="unforced-APV-modes", formulation="F", ...
                modeFamily="hydrostatic", N2=N2, zDomain=zDomain, ...
                p=@(z,ctx) 1./ctx.N2(z), q=@(z,~) zeros(size(z)), ...
                r=@(z,ctx) ones(size(z))/ctx.g, g=g, ...
                surfaceBoundary=surface, bottomBoundary=bottom);

            spec = evp.innerProduct("G");

            testCase.verifyEqual(evp.name, "unforced-APV-modes")
            testCase.verifyEqual(evp.modeFamily, "hydrostatic")
            testCase.verifyTrue(spec.hasInnerProduct)
            testCase.verifyEqual([spec.endpointInnerProductTerms.catalogCase], ["F-P2" "F-P2"])
            testCase.verifyEqual(spec.endpointInnerProductTerms(1).coefficient, -3/(2*g), RelTol=1e-12)
            testCase.verifyEqual(spec.endpointInnerProductTerms(2).coefficient, -5/(4*g), RelTol=1e-12)
        end

        function geostrophicNormalizationCouplesFAndGForGFormulatedModes(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain);
            basisSet = solver.solveEVP(evp, nModes=2);
            basisSet.normalization = "geostrophic";

            testCase.verifyCoupledGeostrophicNormalization(basisSet)
        end

        function geostrophicNormalizationCouplesFAndGForFFormulatedModes(testCase)
            [N2, zDomain, nEVP] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMInternalModes.hydrostaticFModes(N2=N2, zDomain=zDomain);
            basisSet = solver.solveEVP(evp, nModes=3);
            basisSet.normalization = "geostrophic";

            testCase.verifyCoupledGeostrophicNormalization(basisSet)
        end

        function geostrophicNormalizationUsesEndpointValueDiagnosticTerms(testCase)
            [N2, zDomain, nEVP, g] = testCase.profile();
            solver = IMSolverSpectral(nEVP=nEVP);
            evp = IMInternalModes.hydrostaticFModes(N2=N2, zDomain=zDomain, g=g, ...
                surfaceBoundary=IMBoundaryCondition(a=2, b=3));
            basisSet = solver.solveEVP(evp, nModes=2);
            basisSet = basisSet.addNormalization("raw", @(~,~) 1);
            basisSet.normalization = "raw";
            rawGGram = basisSet.gramMatrix(variable="G");
            rawFactor = sqrt(abs(rawGGram(1,1)));
            factors = basisSet.normalizationFactors("geostrophic");

            testCase.verifyEqual( ...
                factors(1), ...
                rawFactor, RelTol=2e-8, AbsTol=2e-8)
        end
    end

    methods
        function verifyUnavailable(testCase, spec, reasonText)
            testCase.verifyFalse(spec.hasInnerProduct)
            testCase.verifyTrue(contains(spec.reason, reasonText))
            testCase.verifyFalse(isfield(spec, "status"))
        end

        function verifyInteriorOnly(testCase, spec)
            testCase.verifyTrue(spec.hasInnerProduct)
            testCase.verifyEmpty(spec.endpointInnerProductTerms)
            testCase.verifyFalse(isfield(spec, "status"))
            testCase.verifyTrue(contains(spec.reason, "interior-only"))
        end

        function verifyCoupledGeostrophicNormalization(testCase, basisSet)
            baroclinic = basisSet.modeNumber ~= 0;
            gramG = basisSet.gramMatrix(variable="G");
            gramF = basisSet.gramMatrix(variable="F");

            testCase.verifyEqual(diag(gramG(baroclinic,baroclinic)).', ...
                ones(1,nnz(baroclinic)), RelTol=2e-8, AbsTol=2e-8)
            testCase.verifyEqual(diag(gramF(baroclinic,baroclinic)).', ...
                basisSet.h(baroclinic), RelTol=2e-7, AbsTol=2e-7)
        end
    end

    methods (Static, Access = private)
        function [N2, zDomain, nEVP, g] = profile()
            N0 = 5.2e-3;
            zDomain = [-1000 0];
            nEVP = 24;
            g = 9.81;
            N2 = @(z) N0*N0*(1 + 0.1*z/abs(zDomain(1)));
        end
    end
end
