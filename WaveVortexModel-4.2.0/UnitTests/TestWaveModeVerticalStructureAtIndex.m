classdef TestWaveModeVerticalStructureAtIndex < matlab.unittest.TestCase
    % Verify wave vertical-structure factors across transform geometries.

    properties
        wvt
    end

    properties (ClassSetupParameter)
        transformType = struct( ...
            boussinesq="boussinesq", ...
            constantNonhydrostatic="constantNonhydrostatic", ...
            constantHydrostatic="constantHydrostatic", ...
            hydrostatic="hydrostatic")
        shouldAntialias = struct(full=false,dealiased=true)
    end

    methods (TestClassSetup)
        function createTransform(testCase,transformType,shouldAntialias)
            Lxyz = [60e3 40e3 1e3];
            Nxyz = [8 6 7];
            switch transformType
                case "boussinesq"
                    testCase.wvt = WVTransformBoussinesq(Lxyz,Nxyz,N2=@(z) 2e-5*exp(z/4000),latitude=45,shouldAntialias=shouldAntialias);
                case "constantNonhydrostatic"
                    testCase.wvt = WVTransformConstantStratification(Lxyz,Nxyz,N0=sqrt(2e-5),latitude=45,shouldAntialias=shouldAntialias,isHydrostatic=false);
                case "constantHydrostatic"
                    testCase.wvt = WVTransformConstantStratification(Lxyz,Nxyz,N0=sqrt(2e-5),latitude=45,shouldAntialias=shouldAntialias,isHydrostatic=true);
                case "hydrostatic"
                    testCase.wvt = WVTransformHydrostatic(Lxyz,Nxyz,N2=@(z) 2e-5*exp(z/4000),latitude=45,shouldAntialias=shouldAntialias);
            end
        end
    end

    methods (Test)
        function factorsMatchNativeTransforms(testCase)
            previousRandomState = rng;
            randomStateCleanup = onCleanup(@()rng(previousRandomState));
            rng(19374,"twister")

            wvt = testCase.wvt;
            waveMask = logical(wvt.waveComponent.maskAp) | logical(wvt.waveComponent.maskAm);
            Apm = (randn(wvt.spectralMatrixSize)+1i*randn(wvt.spectralMatrixSize)).*waveMask;
            FField = wvt.transformToSpatialDomainWithF(Apm=Apm);
            FFourier = wvt.transformFromSpatialDomainWithFourier(FField);
            GField = wvt.transformToSpatialDomainWithG(Apm=Apm);
            dGdzFourier = wvt.transformFromSpatialDomainWithFourier(wvt.diffZG(GField));

            for iZ = unique([1 ceil(wvt.Nz/2) wvt.Nz])
                [F,dGdz] = wvt.waveModeVerticalStructureAtIndex(iZ);
                testCase.verifySize(F,wvt.spectralMatrixSize)
                testCase.verifySize(dGdz,wvt.spectralMatrixSize)
                testCase.verifyTrue(all(isfinite(F),"all"))
                testCase.verifyTrue(all(isfinite(dGdz),"all"))
                testCase.verifyLessThanOrEqual(testCase.relativeError(sum(F.*Apm,1),FFourier(iZ,:)),5e-12)
                testCase.verifyLessThanOrEqual(testCase.relativeError(sum(dGdz.*Apm,1),dGdzFourier(iZ,:)),5e-12)

                FOnly = wvt.waveModeVerticalStructureAtIndex(iZ);
                testCase.verifyEqual(FOnly,F)
            end
            clear randomStateCleanup
        end

        function invalidVerticalIndexUsesStableError(testCase)
            for iZ = [0 testCase.wvt.Nz+1 1.5 Inf NaN]
                testCase.verifyError(@()testCase.wvt.waveModeVerticalStructureAtIndex(iZ),"WVStratification:InvalidVerticalIndex")
            end
        end

        function constantResolutionConversionsPreserveConfiguration(testCase)
            source = testCase.wvt;
            if ~isa(source,"WVTransformConstantStratification")
                return
            end

            source.t0 = 17;
            source.t = 43;
            source.A0 = complex(reshape(1:numel(source.A0),size(source.A0)));
            source.Ap = (2-1i)*source.A0;
            source.Am = (-1+3i)*source.A0;
            target = source.waveVortexTransformWithResolution([10 8 9]);
            testCase.verifyEqual(target.N0,source.N0)
            testCase.verifyEqual(target.isHydrostatic,source.isHydrostatic)
            testCase.verifyEqual(target.shouldAntialias,source.shouldAntialias)
            testCase.verifyEqual([target.Lx target.Ly target.Lz],[source.Lx source.Ly source.Lz])
            testCase.verifyEqual(target.t0,source.t0)
            testCase.verifyEqual(target.t,source.t)
            [expectedA0,expectedAp,expectedAm] = source.spectralVariableWithResolution(target,source.A0,source.Ap,source.Am);
            testCase.verifyEqual(target.A0,expectedA0)
            testCase.verifyEqual(target.Ap,expectedAp)
            testCase.verifyEqual(target.Am,expectedAm)

            if ~source.shouldAntialias
                return
            end
            explicit = source.waveVortexTransformWithExplicitAntialiasing();
            testCase.verifyEqual(explicit.N0,source.N0)
            testCase.verifyEqual(explicit.isHydrostatic,source.isHydrostatic)
            testCase.verifyFalse(explicit.shouldAntialias)
            testCase.verifyEqual([explicit.Lx explicit.Ly explicit.Lz],[source.Lx source.Ly source.Lz])
            testCase.verifyEqual(explicit.t0,source.t0)
            testCase.verifyEqual(explicit.t,source.t)
            [expectedA0,expectedAp,expectedAm] = source.spectralVariableWithResolution(explicit,source.A0,source.Ap,source.Am);
            testCase.verifyEqual(explicit.A0,expectedA0)
            testCase.verifyEqual(explicit.Ap,expectedAp)
            testCase.verifyEqual(explicit.Am,expectedAm)
        end
    end

    methods (Static, Access = private)
        function errorValue = relativeError(actual,expected)
            errorValue = norm(actual(:)-expected(:))/max(norm(expected(:)),eps);
        end
    end
end
