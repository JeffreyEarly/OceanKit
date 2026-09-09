classdef InternalModesBasisTransformTests < matlab.unittest.TestCase

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
            addpath(fullfile(parentRoot, 'class-annotations'));
            addpath(fullfile(parentRoot, 'netcdf'));
        end
    end

    methods (TestClassTeardown)
        function restoreRepositoryPaths(testCase)
            path(testCase.originalPath);
            warning(testCase.pathWarningState);
        end
    end

    methods (Test)
        function geostrophicBasisHasCanonicalFAndG(testCase)
            basis = testCase.geostrophicBasis(12);
            transform = basis.nativeTransform(component="both",projectionMethod="weightedPseudoinverse");

            testCase.verifyTrue(transform.forwardProjectionAvailableF)
            testCase.verifyTrue(transform.forwardProjectionAvailableG)
            testCase.verifyEqual(transform.componentRoleF, "eigenfunction")
            testCase.verifyEqual(transform.componentRoleG, "eigenfunction")
            testCase.verifyLessThan(transform.gramErrorF, 1e-8)
            testCase.verifyLessThan(transform.gramErrorG, 1e-3)
        end

        function modeAdaptedRigidLidGUsesInteriorRows(testCase)
            basis = testCase.geostrophicBasis(16);
            transform = basis.nativeTransform(component="G",projectionMethod="auto");

            inactiveRows = testCase.inactiveRowsForG(basis);
            expectedRetained = length(basis.z) - length(inactiveRows);

            testCase.verifyEqual(length(transform.retainedModesG), expectedRetained)
            testCase.verifyEqual(transform.rejectedModesG, size(basis.G,2))
            testCase.verifyEqual(size(transform.forwardG,2), length(basis.z))
            testCase.verifyEqual(size(transform.inverseG,1), length(basis.z))
            testCase.verifyEqual(transform.forwardG(:,inactiveRows), zeros(size(transform.forwardG,1),length(inactiveRows)), AbsTol=1e-14)
            testCase.verifyEqual(transform.inverseG(inactiveRows,:), zeros(length(inactiveRows),size(transform.inverseG,2)), AbsTol=1e-12)
            testCase.verifyLessThan(transform.gramErrorG, 1e-12)
        end

        function highRequestedModelGRetainsInteriorResolvableModes(testCase)
            basis = testCase.geostrophicBasis(16);
            transform = basis.modelTransform(component="G",nModes=size(basis.G,2),nonlinearAliasingPolicy="none",projectionTolerance=1e-2);

            expectedRetained = length(basis.z) - length(testCase.inactiveRowsForG(basis));
            testCase.verifyEqual(length(transform.retainedModesG), expectedRetained)
            testCase.verifyEqual(transform.rejectedModesG, size(basis.G,2))
            testCase.verifyEqual(transform.selectionReason, "projectionQualityLimit")
            testCase.verifyLessThan(transform.gramErrorG, 1e-12)
        end

        function quadraticAliasingStillCapsModeAdaptedG(testCase)
            nModes = 30;
            basis = testCase.geostrophicBasis(nModes);
            transform = basis.modelTransform(component="G",nModes=nModes,nonlinearAliasingPolicy="quadratic",projectionTolerance=1e-2);

            expectedLimit = floor(2*nModes/3);
            testCase.verifyEqual(length(transform.retainedModesG), expectedLimit)
            testCase.verifyEqual(transform.nonlinearAliasLimit, expectedLimit)
            testCase.verifyEqual(transform.selectionReason, "nonlinearAliasingLimit")
            testCase.verifyLessThan(transform.gramErrorG, 1e-12)
        end

        function canonicalGeostrophicGUsesProjectionNormalization(testCase)
            basis = testCase.geostrophicBasis(16);
            transform = basis.nativeTransform(component="G",projectionMethod="canonical");

            testCase.verifyEqual(transform.transformStatusG, "canonical")
            testCase.verifyLessThan(transform.gramErrorG, 1e-2)
        end

        function canonicalWaveGDoesNotDivideBySpectrumWeight(testCase)
            transform = testCase.waveBasis(12, 2*pi/150).nativeTransform(component="G",projectionMethod="canonical");

            activeRows = find(max(abs(transform.inverseG),[],2) > 0);
            expectedForward = zeros(size(transform.forwardG));
            expectedForward(:,activeRows) = transform.inverseG(activeRows,:).' .* transform.weightsG(activeRows).';

            testCase.verifyEqual(transform.transformStatusG, "canonical")
            testCase.verifyEqual(transform.forwardG, expectedForward, AbsTol=1e-12)
        end

        function waveBasisRejectsCanonicalFProjection(testCase)
            basis = testCase.waveBasis(12, 2*pi/150);
            transform = basis.nativeTransform(component="both",projectionMethod="weightedPseudoinverse",allowNoncanonical=true);

            testCase.verifyFalse(transform.forwardProjectionAvailableF)
            testCase.verifyTrue(transform.forwardProjectionAvailableG)
            testCase.verifyEqual(transform.componentRoleF, "diagnostic")
            testCase.verifyError(@() transform.spectrum(ones(length(transform.spectralWeightsF),1),component="F"), ...
                "InternalModesTransform:NoncanonicalSpectrum")
            testCase.verifyWarningFree(@() transform.forward(component="F",allowNoncanonical=true));
        end

        function spectraUseParsevalWeights(testCase)
            basis = testCase.geostrophicBasis(10);
            transform = basis.nativeTransform(component="both",projectionMethod="weightedPseudoinverse");

            aF = (1:length(transform.spectralWeightsF)).';
            bF = 0.25*aF;
            SF = transform.crossSpectrum(aF,bF,component="F",horizontalMultiplicity=2);
            expectedF = 2*transform.spectralWeightsF .* real(aF .* conj(bF));

            aG = (1:length(transform.spectralWeightsG)).';
            bG = -0.5*aG;
            SG = transform.crossSpectrum(aG,bG,component="G");
            expectedG = transform.g * real(aG .* conj(bG));

            testCase.verifyEqual(SF, expectedF, RelTol=1e-14)
            testCase.verifyEqual(SG, expectedG, RelTol=1e-14)
            testCase.verifyEqual(transform.spectrum(aG,component="G"), transform.crossSpectrum(aG,aG,component="G"))
        end

        function quadraticAliasingCapsConstantStratification(testCase)
            nModes = 30;
            z = linspace(-1000,0,nModes+1).';
            im = InternalModesConstantStratification(N0=5.2e-3,zIn=[-1000 0],zOut=z,latitude=33,nModes=nModes);
            im.normalization = Normalization.kConstant;

            basis = InternalModesBasis.fromSolverAtFrequency(im,0,nModes=nModes);
            transform = basis.modelTransform(component="F",nModes=nModes+1,nonlinearAliasingPolicy="quadratic",projectionTolerance=1);

            expectedLimit = floor(2*nModes/3);
            testCase.verifyEqual(length(transform.retainedModesF), expectedLimit)
            testCase.verifyEqual(transform.nonlinearAliasLimit, expectedLimit)
            testCase.verifyEqual(transform.selectionReason, "nonlinearAliasingLimit")
        end

        function fixedGridWaveTransformZerosRejectedRows(testCase)
            hydroBasis = testCase.geostrophicBasis(24);
            hydroTransform = hydroBasis.modelTransform(component="G",nModes=16,projectionTolerance=1e-2);
            waveBasis = testCase.waveBasis(24, 2*pi/110);

            waveTransform = waveBasis.fixedGridTransform(hydroTransform,component="G",projectionTolerance=1e-2,preserveSize=true);

            rejected = waveTransform.rejectedModesG;
            testCase.verifyNotEmpty(rejected)
            testCase.verifyEqual(max(vecnorm(waveTransform.forwardG(rejected,:),2,2)), 0, AbsTol=1e-14)
            testCase.verifyEqual(max(vecnorm(waveTransform.inverseG(:,rejected),2,1)), 0, AbsTol=1e-14)
        end

        function observationProjectionSelectsResolvableModes(testCase)
            basis = testCase.geostrophicBasis(20);
            zObs = linspace(-900,-100,14).';
            projection = basis.observationProjection(zObs,component="G",nModes=20,maxConditionNumber=20);

            testCase.verifyGreaterThan(length(projection.retainedModes), 0)
            testCase.verifyLessThanOrEqual(length(projection.retainedModes), 14)
            coefficients = ones(length(projection.retainedModes),1);
            testCase.verifySize(projection.spectrum(coefficients), size(coefficients))
        end

        function transformPersistenceRoundTripsRequiredState(testCase)
            basis = testCase.geostrophicBasis(8);
            transform = basis.nativeTransform(component="both",projectionMethod="weightedPseudoinverse");
            outputPath = fullfile(tempdir, "InternalModesTransformRoundTrip.nc");
            cleanup = onCleanup(@() deleteIfPresent(outputPath));

            transform.writeToFile(outputPath, shouldOverwriteExisting=true);
            restored = CAAnnotatedClass.annotatedClassFromFile(outputPath);

            testCase.verifyClass(restored, "InternalModesTransform")
            testCase.verifyEqual(restored.forwardG, transform.forwardG, AbsTol=0)
            testCase.verifyEqual(restored.inverseG, transform.inverseG, AbsTol=0)
            testCase.verifyEqual(restored.spectralWeightsG, transform.spectralWeightsG, AbsTol=0)
            testCase.verifyEqual(restored.componentRoleG, transform.componentRoleG)
        end
    end

    methods (Access = private)
        function basis = geostrophicBasis(~, nModes)
            [N2,zDomain,z] = localProfile(nModes);
            im = InternalModesWKBSpectral(N2=N2,zIn=zDomain,zOut=z,latitude=31,nEVP=max(128,ceil(2.1*(nModes+1))),nModes=nModes,g=9.81);
            im.normalization = Normalization.geostrophic;
            im.upperBoundary = UpperBoundary.rigidLid;
            basis = InternalModesBasis.fromSolverAtFrequency(im,0,nModes=nModes,useModeAdaptedGrid=true,g=9.81);
        end

        function basis = waveBasis(~, nModes, kappa)
            [N2,zDomain,z] = localProfile(nModes);
            im = InternalModesWKBSpectral(N2=N2,zIn=zDomain,zOut=z,latitude=31,nEVP=max(128,ceil(2.1*(nModes+1))),nModes=nModes,g=9.81);
            im.normalization = Normalization.kConstant;
            im.upperBoundary = UpperBoundary.rigidLid;
            basis = InternalModesBasis.fromSolverAtWavenumber(im,kappa,nModes=nModes,g=9.81);
        end

        function rows = inactiveRowsForG(~, basis)
            rowScale = max(abs(basis.G),[],2);
            tolerance = sqrt(eps(max(1,max(rowScale))));
            rows = find(rowScale <= tolerance);
        end
    end
end

function [N2,zDomain,z] = localProfile(nModes)
Lz = 1200;
N0 = 3*2*pi/3600;
L_gm = 900;
N2 = @(z) N0*N0*exp(2*z/L_gm);
zDomain = [-Lz 0];
z = linspace(zDomain(1),zDomain(2),nModes+1).';
end

function deleteIfPresent(path)
if isfile(path)
    delete(path);
end
end
