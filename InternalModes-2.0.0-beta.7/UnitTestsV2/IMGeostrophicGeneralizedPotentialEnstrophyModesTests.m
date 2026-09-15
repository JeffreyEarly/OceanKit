classdef IMGeostrophicGeneralizedPotentialEnstrophyModesTests < matlab.unittest.TestCase

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
        function factoryMapsCanonicalProblemAndEndpointLimits(testCase)
            [N2,zDomain,f0,g,k] = testCase.profile();
            alpha0 = 2e-11;
            alphaD = 3e-11;
            evp = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes( ...
                N2=N2,zDomain=zDomain,k=k,f0=f0,g=g,alpha0=alpha0,alphaD=alphaD);
            z = linspace(zDomain(1),zDomain(2),7).';
            context = struct("N2",N2,"f0",f0,"k",k);

            testCase.verifyEqual(string(evp.name),"geostrophicGeneralizedPotentialEnstrophyModes")
            testCase.verifyEqual(evp.formulation,"F")
            testCase.verifyEqual(evp.modeFamily,"none")
            testCase.verifyEqual(IMEigenvalueProblem.evaluateCoefficient(evp.p,z,context),f0*f0./N2(z),RelTol=1e-14)
            testCase.verifyEqual(IMEigenvalueProblem.evaluateCoefficient(evp.q,z,context),k*k*ones(size(z)),RelTol=1e-14)
            testCase.verifyEqual(IMEigenvalueProblem.evaluateCoefficient(evp.r,z,context),ones(size(z)),AbsTol=0)
            testCase.verifyEqual([evp.surfaceBoundary.a evp.surfaceBoundary.b evp.surfaceBoundary.c evp.surfaceBoundary.d], ...
                [-f0*f0/g 1 f0*f0/alpha0 0],RelTol=1e-14,AbsTol=1e-14)
            testCase.verifyEqual([evp.bottomBoundary.a evp.bottomBoundary.b evp.bottomBoundary.c evp.bottomBoundary.d], ...
                [0 1 -f0*f0/alphaD 0],RelTol=1e-14,AbsTol=1e-14)
            testCase.verifyEqual(evp.parameters.alpha0Source,"explicit")

            inactiveEVP = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes( ...
                N2=N2,zDomain=zDomain,k=k,f0=f0,g=g,alpha0=Inf,alphaD=Inf);
            testCase.verifyEqual([inactiveEVP.surfaceBoundary.a inactiveEVP.surfaceBoundary.b inactiveEVP.surfaceBoundary.c inactiveEVP.surfaceBoundary.d], ...
                [-f0*f0/g 1 0 0],RelTol=1e-14,AbsTol=1e-14)
            testCase.verifyEqual([inactiveEVP.bottomBoundary.a inactiveEVP.bottomBoundary.b inactiveEVP.bottomBoundary.c inactiveEVP.bottomBoundary.d], ...
                [0 1 0 0],AbsTol=0)
        end

        function defaultSurfaceWeightUsesStratificationScale(testCase)
            [N2,zDomain,f0,~,k] = testCase.profile();
            depth = diff(zDomain);
            evp = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes( ...
                N2=N2,zDomain=zDomain,k=k,f0=f0);

            testCase.verifyEqual(evp.parameters.bEffective,depth/4,RelTol=1e-12)
            testCase.verifyEqual(evp.parameters.alpha0,4*f0*f0/depth,RelTol=1e-12)
            testCase.verifyEqual(evp.parameters.alphaD,Inf)
            testCase.verifyEqual(evp.parameters.alpha0Source,"stratification")
        end

        function invalidPhysicalParametersAreRejected(testCase)
            [N2,zDomain,f0,~,k] = testCase.profile();
            testCase.verifyError(@() IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes( ...
                N2=N2,zDomain=zDomain,k=0,f0=f0),"MATLAB:validators:mustBePositive")
            testCase.verifyError(@() IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes( ...
                N2=N2,zDomain=zDomain,k=k,f0=0),"IMInternalModes:InvalidCoriolisParameter")
            testCase.verifyError(@() IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes( ...
                N2=N2,zDomain=zDomain,k=k,f0=f0,alpha0=0),"IMInternalModes:InvalidSurfaceEnstrophyWeight")
            testCase.verifyError(@() IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes( ...
                N2=N2,zDomain=zDomain,k=k,f0=f0,alphaD=-1),"IMInternalModes:InvalidBottomEnstrophyWeight")
        end

        function numericalModesSatisfyBoundaryNormAndEigendepthIdentities(testCase)
            [N2,zDomain,f0,g,k] = testCase.profile();
            alpha0 = 2e-11;
            alphaD = 3e-11;
            evp = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes( ...
                N2=N2,zDomain=zDomain,k=k,f0=f0,g=g,alpha0=alpha0,alphaD=alphaD);
            basisSet = IMSolverSpectral(nEVP=96,coordinateKind="z").solveEVP(evp,nModes=6);
            z = basisSet.solver.zNative;
            F = basisSet.F(z);
            Fz = basisSet.solver.differentiateGridValues(F,1);
            surfaceIndex = basisSet.solver.boundaryIndex("surface");
            bottomIndex = basisSet.solver.boundaryIndex("bottom");
            Lambda = basisSet.eigenvalues;
            surfaceResidual = (f0*f0/N2(z(surfaceIndex)))*Fz(surfaceIndex,:) ...
                +(f0*f0/g)*F(surfaceIndex,:)-Lambda.*(f0*f0/alpha0).*F(surfaceIndex,:);
            bottomResidual = (f0*f0/N2(z(bottomIndex)))*Fz(bottomIndex,:) ...
                +Lambda.*(f0*f0/alphaD).*F(bottomIndex,:);
            residualScale = max(1,max(abs(F),[],1));

            testCase.verifyEqual(string(basisSet.normalization),"generalizedPotentialEnstrophy")
            testCase.verifyTrue(ismember("generalizedPotentialEnstrophy",basisSet.normalizationNames()))
            testCase.verifyFalse(ismember("depth",basisSet.normalizationNames()))
            testCase.verifyGreaterThan(min(Lambda),0)
            testCase.verifyLessThan(max(abs(surfaceResidual)./residualScale),2e-7)
            testCase.verifyLessThan(max(abs(bottomResidual)./residualScale),2e-7)
            testCase.verifyEqual(basisSet.h,f0*f0./(g*(Lambda-k*k)),RelTol=1e-13)
            testCase.verifyEqual(basisSet.gramMatrix(variable="F")/diff(zDomain),eye(6),RelTol=2e-7,AbsTol=2e-7)
        end

        function existingDepthNormalizationRemainsVolumeOnly(testCase)
            [N2,zDomain] = testCase.profile();
            evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=zDomain,g0=-0.02,gd=Inf);
            basisSet = IMSolverSpectral(nEVP=64).solveEVP(evp,nModes=4);
            z = basisSet.solver.innerProductGrid(zDomain);
            F = basisSet.F(z,normalization=Normalization.depth);
            weights = basisSet.solver.innerProductWeights(z,zDomain);
            volumeGram = F.'*(weights.*F)/diff(zDomain);

            testCase.verifyEqual(volumeGram,eye(4),RelTol=2e-7,AbsTol=2e-7)
        end

        function modeRootQuadratureVariesWithWavenumber(testCase)
            D = 1200;
            N2 = @(z) (5.2e-3)^2*exp(z/700);
            zDomain = [-D 0];
            f0 = 1e-4;
            kValues = 2*pi./[500e3 20e3];
            grids = cell(size(kValues));
            for iK = 1:numel(kValues)
                evp = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes( ...
                    N2=N2,zDomain=zDomain,k=kValues(iK),f0=f0);
                basisSet = IMSolverSpectral(nEVP=96,coordinateKind="wkb").solveEVP(evp,nModes=16);
                [grids{iK},gridDesign] = basisSet.modeRootGrid(nPoints=12);
                [weights,weightFit] = basisSet.quadratureWeightsForPoints( ...
                    z=grids{iK},nModes=gridDesign.representedModeCount,variables="F");

                testCase.verifyNumElements(grids{iK},12)
                testCase.verifyEqual(sum(weights),D,AbsTol=1e-8)
                testCase.verifyGreaterThanOrEqual(min(weights),-1e-10)
                testCase.verifyLessThan(weightFit.transform.relativeGramOperatorError(variable="F"),0.2)
            end
            testCase.verifyGreaterThan(norm(grids{1}-grids{2}),1e-3)
        end
    end

    methods (Static,Access = private)
        function [N2,zDomain,f0,g,k] = profile()
            N0 = 5.2e-3;
            zDomain = [-4000 0];
            N2 = @(z) N0*N0*ones(size(z));
            f0 = 1e-4;
            g = 9.81;
            k = 2*pi/100e3;
        end
    end
end
