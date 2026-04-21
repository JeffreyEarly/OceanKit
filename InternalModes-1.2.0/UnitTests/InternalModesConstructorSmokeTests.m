classdef InternalModesConstructorSmokeTests < matlab.unittest.TestCase

    properties (Access = private)
        originalPath
        pathWarningState
    end

    methods (TestClassSetup)
        function addRepositoryPaths(testCase)
            repoRoot = fileparts(fileparts(mfilename('fullpath')));
            parentRoot = fileparts(repoRoot);

            testCase.originalPath = path;
            testCase.pathWarningState = warning('off', 'MATLAB:path:reorderPackageFolders');

            addpath(repoRoot);
            addpath(fullfile(parentRoot, 'spline-core'));
            addpath(fullfile(parentRoot, 'distributions'));
            addpath(fullfile(parentRoot, 'chebfun'));
        end
    end

    methods (TestClassTeardown)
        function restoreRepositoryPaths(testCase)
            path(testCase.originalPath);
            warning(testCase.pathWarningState);
        end
    end

    methods (Test)
        function spectralConstructorSupportsRhoInitialization(testCase)
            [rhoFunction, ~, zIn, zOut, N0] = testCase.exponentialProfile();

            im = InternalModesSpectral(rho=rhoFunction, zIn=zIn, zOut=zOut, latitude=33, nEVP=33, nModes=4);
            [F, G, h] = im.modesAtFrequency(0.8*N0);

            testCase.verifyModes(F, G, h, zOut, 4)
        end

        function spectralConstructorSupportsN2Initialization(testCase)
            [~, N2Function, zIn, zOut, N0] = testCase.exponentialProfile();

            im = InternalModesSpectral(N2=N2Function, zIn=zIn, zOut=zOut, latitude=33, nEVP=33, nModes=4);
            [F, G, h] = im.modesAtFrequency(0.8*N0);

            testCase.verifyModes(F, G, h, zOut, 4)
        end

        function wkbSpectralConstructorSupportsRhoInitialization(testCase)
            [rhoFunction, ~, zIn, zOut, N0] = testCase.exponentialProfile();

            im = InternalModesWKBSpectral(rho=rhoFunction, zIn=zIn, zOut=zOut, latitude=33, nEVP=33, nModes=4);
            [F, G, h] = im.modesAtFrequency(0.8*N0);

            testCase.verifyModes(F, G, h, zOut, 4)
        end

        function wkbSpectralConstructorSupportsN2Initialization(testCase)
            [~, N2Function, zIn, zOut, N0] = testCase.exponentialProfile();

            im = InternalModesWKBSpectral(N2=N2Function, zIn=zIn, zOut=zOut, latitude=33, nEVP=33, nModes=4);
            [F, G, h] = im.modesAtFrequency(0.8*N0);

            testCase.verifyModes(F, G, h, zOut, 4)
        end

        function densitySpectralConstructorComputesModes(testCase)
            [rhoFunction, ~, zIn, zOut, N0] = testCase.exponentialProfile();

            im = InternalModesDensitySpectral(rho=rhoFunction, zIn=zIn, zOut=zOut, latitude=33, nEVP=33, nModes=4);
            [F, G, h] = im.modesAtFrequency(0.8*N0);

            testCase.verifyModes(F, G, h, zOut, 4)
        end

        function adaptiveSpectralConstructorComputesModes(testCase)
            [rhoFunction, ~, zIn, zOut, N0] = testCase.exponentialProfile();

            im = InternalModesAdaptiveSpectral(rho=rhoFunction, zIn=zIn, zOut=zOut, latitude=33, nEVP=33, nModes=4);
            [F, G, h] = im.modesAtFrequency(0.8*N0);

            testCase.verifyModes(F, G, h, zOut, 4)
        end

        function finiteDifferenceConstructorComputesModes(testCase)
            [rhoFunction, ~, ~, zOut, ~] = testCase.exponentialProfile();

            im = InternalModesFiniteDifference(rho=rhoFunction, zIn=zOut, zOut=zOut, latitude=33, nModes=4, orderOfAccuracy=4);
            [F, G, h] = im.modesAtWavenumber(1e-4);

            testCase.verifyModes(F, G, h, zOut, 4)
        end

        function wkbConstructorComputesModes(testCase)
            [rhoFunction, ~, zIn, zOut, N0] = testCase.exponentialProfile();

            im = InternalModesWKB(rho=rhoFunction, zIn=zIn, zOut=zOut, latitude=33, nEVP=33, nModes=4);
            [F, G, h] = im.modesAtFrequency(0.8*N0);

            testCase.verifyModes(F, G, h, zOut, 4)
        end

        function wkbHydrostaticConstructorComputesModes(testCase)
            [rhoFunction, ~, zIn, zOut, N0] = testCase.exponentialProfile();

            im = InternalModesWKBHydrostatic(rho=rhoFunction, zIn=zIn, zOut=zOut, latitude=33, nEVP=33, nModes=4);
            [F, G, h] = im.modesAtFrequency(0.8*N0);

            testCase.verifyModes(F, G, h, zOut, 4)
        end

        function exponentialStratificationConstructorComputesModes(testCase)
            [~, ~, zIn, zOut, N0, b] = testCase.exponentialProfile();

            im = InternalModesExponentialStratification(N0=N0, b=b, zIn=zIn, zOut=zOut, latitude=33, nModes=4);
            [F, G, h] = im.modesAtFrequency(0.8*N0);

            testCase.verifyModes(F, G, h, zOut, 4)
        end
        
        function spectralSurfaceModeComputesWithLowerCamelApi(testCase)
            [rhoFunction, ~, zIn, zOut, ~] = testCase.exponentialProfile();
            
            im = InternalModesSpectral(rho=rhoFunction, zIn=zIn, zOut=zOut, latitude=33, nEVP=33, nModes=4);
            psi = im.surfaceModesAtWavenumber(1e-4);
            
            testCase.verifyEqual(size(psi, ndims(psi)), numel(zOut))
        end
        
        function spectralEigenmatricesForFrequencyComputesWithLowerCamelApi(testCase)
            [rhoFunction, ~, zIn, zOut, N0] = testCase.exponentialProfile();
            
            im = InternalModesSpectral(rho=rhoFunction, zIn=zIn, zOut=zOut, latitude=33, nEVP=33, nModes=4);
            [A, B] = im.eigenmatricesForFrequency(0.8*N0);
            
            testCase.verifyEVPMatrices(A, B, im.nEVP)
        end
        
        function spectralEigenmatricesForWavenumberComputesWithLowerCamelApi(testCase)
            [rhoFunction, ~, zIn, zOut] = testCase.exponentialProfile();
            
            im = InternalModesSpectral(rho=rhoFunction, zIn=zIn, zOut=zOut, latitude=33, nEVP=33, nModes=4);
            [A, B] = im.eigenmatricesForWavenumber(1e-4);
            
            testCase.verifyEVPMatrices(A, B, im.nEVP)
        end
        
        function adaptiveEigenmatricesForFrequencyComputesWithLowerCamelApi(testCase)
            [rhoFunction, ~, zIn, zOut, N0] = testCase.exponentialProfile();
            omega = 0.8*N0;
            
            im = InternalModesAdaptiveSpectral(rho=rhoFunction, zIn=zIn, zOut=zOut, latitude=33, nEVP=33, nModes=4);
            im.modesAtFrequency(omega);
            [A, B] = im.eigenmatricesForFrequency(omega);
            
            testCase.verifyEVPMatrices(A, B, im.nEVP)
        end

        function wrapperDetectsConstantStratification(testCase)
            [rhoFunction, ~, zIn] = InternalModes.StratificationProfileWithName('constant');
            zOut = linspace(min(zIn), max(zIn), 33)';

            im = InternalModes(rhoFunction, zIn, zOut, 33, 'nModes', 4);
            [F, G, h] = im.modesAtFrequency(0.8*5.2e-3);

            testCase.verifyClass(im.internalModes, 'InternalModesConstantStratification')
            testCase.verifyModes(F, G, h, zOut, 4)
        end

        function wrapperDetectsExponentialStratification(testCase)
            [rhoFunction, ~, zIn] = InternalModes.StratificationProfileWithName('exponential');
            zOut = linspace(min(zIn), max(zIn), 33)';

            im = InternalModes(rhoFunction, zIn, zOut, 33, 'nModes', 4);
            [F, G, h] = im.modesAtFrequency(0.8*5.2e-3);

            testCase.verifyClass(im.internalModes, 'InternalModesExponentialStratification')
            testCase.verifyModes(F, G, h, zOut, 4)
        end

        function wrapperBuildsGenericAdaptivePath(testCase)
            [rhoFunction, ~, zIn] = InternalModes.StratificationProfileWithName('pycnocline-constant');
            zOut = linspace(min(zIn), max(zIn), 33)';

            im = InternalModes(rhoFunction, zIn, zOut, 33, 'nEVP', 33, 'nModes', 4);
            omega = im.f0 + 0.5*(max(sqrt(im.N2)) - im.f0);
            [F, G, h] = im.modesAtFrequency(omega);

            testCase.verifyClass(im.internalModes, 'InternalModesAdaptiveSpectral')
            testCase.verifyModes(F, G, h, zOut, 4)
        end

        function wrapperAppliesWrapperOnlyConstructorOptions(testCase)
            [rhoFunction, ~, zIn] = InternalModes.StratificationProfileWithName('pycnocline-constant');
            zOut = linspace(min(zIn), max(zIn), 33)';

            im = InternalModes(rhoFunction, zIn, zOut, 33, ...
                'method', 'spectral', 'nEVP', 33, 'nModes', 4, ...
                'shouldShowDiagnostics', 1, ...
                'upperBoundary', UpperBoundary.freeSurface, ...
                'lowerBoundary', LowerBoundary.noSlip, ...
                'normalization', Normalization.omegaConstant);
            [F, G, h] = im.modesAtWavenumber(1e-4);

            testCase.verifyEqual(im.shouldShowDiagnostics, 1)
            testCase.verifyEqual(im.upperBoundary, UpperBoundary.freeSurface)
            testCase.verifyEqual(im.lowerBoundary, LowerBoundary.noSlip)
            testCase.verifyEqual(im.normalization, Normalization.omegaConstant)
            testCase.verifyModes(F, G, h, zOut, 4)
        end
        
        function spectralModesAtFrequencyCompatibilityAliasMatchesLowerCamel(testCase)
            [rhoFunction, ~, zIn, zOut, N0] = testCase.exponentialProfile();
            omega = 0.8*N0;
            
            im = InternalModesSpectral(rho=rhoFunction, zIn=zIn, zOut=zOut, latitude=33, nEVP=33, nModes=4);
            [FNew, GNew, hNew, kNew] = im.modesAtFrequency(omega);
            [FOld, GOld, hOld, kOld] = im.ModesAtFrequency(omega);
            
            testCase.verifyEqual(FNew, FOld)
            testCase.verifyEqual(GNew, GOld)
            testCase.verifyEqual(hNew, hOld)
            testCase.verifyEqual(kNew, kOld)
        end
        
        function spectralEigenmatricesForFrequencyCompatibilityAliasMatchesLowerCamel(testCase)
            [rhoFunction, ~, zIn, zOut, N0] = testCase.exponentialProfile();
            omega = 0.8*N0;
            
            im = InternalModesSpectral(rho=rhoFunction, zIn=zIn, zOut=zOut, latitude=33, nEVP=33, nModes=4);
            [ANew, BNew] = im.eigenmatricesForFrequency(omega);
            [AOld, BOld] = im.EigenmatricesForFrequency(omega);
            
            testCase.verifyEqual(ANew, AOld)
            testCase.verifyEqual(BNew, BOld)
        end
        
        function spectralModesAtWavenumberCompatibilityAliasMatchesLowerCamel(testCase)
            [rhoFunction, ~, zIn, zOut, ~] = testCase.exponentialProfile();
            k = 1e-4;
            
            im = InternalModesSpectral(rho=rhoFunction, zIn=zIn, zOut=zOut, latitude=33, nEVP=33, nModes=4);
            [FNew, GNew, hNew, omegaNew] = im.modesAtWavenumber(k);
            [FOld, GOld, hOld, omegaOld] = im.ModesAtWavenumber(k);
            
            testCase.verifyEqual(FNew, FOld)
            testCase.verifyEqual(GNew, GOld)
            testCase.verifyEqual(hNew, hOld)
            testCase.verifyEqual(omegaNew, omegaOld)
        end
        
        function spectralEigenmatricesForWavenumberCompatibilityAliasMatchesLowerCamel(testCase)
            [rhoFunction, ~, zIn, zOut] = testCase.exponentialProfile();
            k = 1e-4;
            
            im = InternalModesSpectral(rho=rhoFunction, zIn=zIn, zOut=zOut, latitude=33, nEVP=33, nModes=4);
            [ANew, BNew] = im.eigenmatricesForWavenumber(k);
            [AOld, BOld] = im.EigenmatricesForWavenumber(k);
            
            testCase.verifyEqual(ANew, AOld)
            testCase.verifyEqual(BNew, BOld)
        end
        
        function wrapperSurfaceModesCompatibilityAliasMatchesLowerCamel(testCase)
            [rhoFunction, ~, zIn] = InternalModes.StratificationProfileWithName('exponential');
            zOut = linspace(min(zIn), max(zIn), 33)';
            k = 1e-4;
            
            im = InternalModes(rhoFunction, zIn, zOut, 33, 'method', 'spectral', 'nEVP', 33, 'nModes', 4);
            psiNew = im.surfaceModesAtWavenumber(k);
            psiOld = im.SurfaceModesAtWavenumber(k);
            
            testCase.verifyEqual(psiNew, psiOld)
        end
    end

    methods (Access = private)
        function verifyModes(testCase, F, G, h, zOut, minModeCount)
            testCase.verifyEqual(size(F, 1), numel(zOut))
            testCase.verifyEqual(size(G, 1), numel(zOut))
            testCase.verifyGreaterThanOrEqual(size(F, 2), minModeCount)
            testCase.verifyGreaterThanOrEqual(size(G, 2), minModeCount)
            testCase.verifyGreaterThanOrEqual(numel(h), minModeCount)
            testCase.verifyGreaterThan(h(1), 0)
        end
        
        function verifyEVPMatrices(testCase, A, B, nEVP)
            testCase.verifySize(A, [nEVP nEVP])
            testCase.verifySize(B, [nEVP nEVP])
            testCase.verifyFalse(any(isnan(A(:))))
            testCase.verifyFalse(any(isnan(B(:))))
        end
    end

    methods (Static, Access = private)
        function [rhoFunction, N2Function, zIn, zOut, N0, b] = exponentialProfile()
            N0 = 5.2e-3;
            b = 1300;
            rho0 = 1025;
            g = 9.81;
            zIn = [-5000 0];
            zOut = linspace(zIn(1), zIn(2), 33)';

            rhoFunction = @(z) rho0*(1 + b*N0*N0/(2*g)*(1 - exp(2*z/b)));
            N2Function = @(z) N0*N0*exp(2*z/b);
        end
    end

end
