classdef TestWVPseudoTopographicWaveGeneration < matlab.unittest.TestCase
    % Verify prescribed bottom wave generation and model integration.

    methods (Test, TestTags = "full")
        function darwinSymbolsSelectFrequency(testCase)
            wvt = TestWVPseudoTopographicWaveGeneration.createTransform(false);
            terrain = TestWVPseudoTopographicWaveGeneration.sinusoidalTopography(wvt,50);
            symbols = ["M2" "S2" "N2" "K1" "O1"];
            periodHours = [12.420602 12.000000 12.65834751 23.93447213 25.81933871];
            for iConstituent = 1:numel(symbols)
                forcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0],darwinSymbol=symbols(iConstituent));
                testCase.verifyEqual(forcing.darwinSymbol,symbols(iConstituent))
                testCase.verifyEqual(forcing.frequency,2*pi/(periodHours(iConstituent)*3600),RelTol=10*eps)
            end

            defaultForcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0]);
            testCase.verifyEqual(defaultForcing.darwinSymbol,"M2")
            testCase.verifyEqual(defaultForcing.name,"pseudo-topographic wave generation")

            customFrequency = 1.37e-4;
            customForcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0],frequency=customFrequency);
            testCase.verifyEqual(customForcing.frequency,customFrequency)
            testCase.verifyEqual(customForcing.darwinSymbol,"")

            testCase.verifyError(@()WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0],frequency=customFrequency,darwinSymbol="M2"),"WVPseudoTopographicWaveGeneration:ConflictingFrequencyOptions")
            testCase.verifyError(@()WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0],darwinSymbol="m2"),"WVPseudoTopographicWaveGeneration:UnknownConstituent")
            testCase.verifyError(@()WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0],darwinSymbol=""),"WVPseudoTopographicWaveGeneration:UnknownConstituent")
        end

        function supportedWaveTransformsProduceFiniteForcing(testCase)
            transforms = {
                WVTransformConstantStratification([4e3 4e3 2e3],[8 4 5],N0=sqrt(2e-5),latitude=45,isHydrostatic=false,shouldAntialias=true)
                WVTransformHydrostatic([4e3 4e3 2e3],[8 4 5],N2=@(z)2e-5*ones(size(z)),latitude=45,shouldAntialias=true)
                WVTransformBoussinesq([4e3 4e3 2e3],[8 4 5],N2=@(z)2e-5*ones(size(z)),latitude=45,shouldAntialias=true)
                };
            for iTransform = 1:numel(transforms)
                wvt = transforms{iTransform};
                terrain = TestWVPseudoTopographicWaveGeneration.sinusoidalTopography(wvt,50);
                forcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0],darwinSymbol="M2");
                [Fp,Fm,F0] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
                testCase.verifyTrue(all(isfinite([Fp(:); Fm(:)])))
                testCase.verifyGreaterThan(norm([Fp(:); Fm(:)]),0)
                testCase.verifyEqual(F0,zeros(size(wvt.A0)))
            end
        end

        function constructorValidatesScientificContract(testCase)
            for shouldAntialias = [false true]
                wvt = TestWVPseudoTopographicWaveGeneration.createTransform(shouldAntialias);
                terrains = {
                    zeros(wvt.Nx,wvt.Ny)
                    100*ones(wvt.Nx,wvt.Ny)
                    TestWVPseudoTopographicWaveGeneration.sinusoidalTopography(wvt,50)
                    TestWVPseudoTopographicWaveGeneration.bandLimitedTopography(wvt)
                    };
                for iTerrain = 1:numel(terrains)
                    forcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrains{iTerrain},barotropicVelocityAmplitude=[0.05; 0]);
                    testCase.verifyEqual(forcing.topographicHeight,terrains{iTerrain})
                end
            end

            wvt = TestWVPseudoTopographicWaveGeneration.createTransform(false);
            terrain = zeros(wvt.Nx,wvt.Ny);
            forcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05+0.01i; -0.02i],name="terrain waves");
            testCase.verifyEqual(forcing.name,"terrain waves")
            testCase.verifyEqual(forcing.frequency,2*pi/(12.420602*3600),"RelTol",10*eps)
            testCase.verifyEqual(forcing.darwinSymbol,"M2")
            testCase.verifyEqual(forcing.rampDuration,0)
            testCase.verifyEqual(forcing.startTime,wvt.t)
            testCase.verifyTrue(forcing.shouldAvoidAdaptiveDamping)
            testCase.verifyEqual(forcing.maximumForcedHorizontalWavenumber,Inf)
            testCase.verifyEqual(forcing.maximumForcedVerticalMode,Inf)

            testCase.verifyError(@()WVPseudoTopographicWaveGeneration(wvt,topographicHeight=zeros(wvt.Nx-1,wvt.Ny),barotropicVelocityAmplitude=[0.05; 0]),"WVPseudoTopographicWaveGeneration:InvalidTopographicHeightSize")
            testCase.verifyError(@()WVPseudoTopographicWaveGeneration(wvt,topographicHeight=1i*ones(wvt.Nx,wvt.Ny),barotropicVelocityAmplitude=[0.05; 0]),"WVPseudoTopographicWaveGeneration:InvalidTopographicHeight")
            invalidTerrain = terrain;
            invalidTerrain(1) = NaN;
            testCase.verifyError(@()WVPseudoTopographicWaveGeneration(wvt,topographicHeight=invalidTerrain,barotropicVelocityAmplitude=[0.05; 0]),"WVPseudoTopographicWaveGeneration:InvalidTopographicHeight")
            testCase.verifyError(@()WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05 0]),"WVPseudoTopographicWaveGeneration:InvalidBarotropicVelocityAmplitude")
            testCase.verifyError(@()WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[Inf; 0]),"WVPseudoTopographicWaveGeneration:InvalidBarotropicVelocityAmplitude")
            testCase.verifyError(@()WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0],frequency=0),"WVPseudoTopographicWaveGeneration:InvalidFrequency")
            testCase.verifyError(@()WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0],rampDuration=-1),"WVPseudoTopographicWaveGeneration:InvalidRampDuration")
            testCase.verifyError(@()WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0],startTime=NaN),"WVPseudoTopographicWaveGeneration:InvalidStartTime")
            testCase.verifyError(@()WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0],maximumForcedHorizontalWavenumber=-1),"WVPseudoTopographicWaveGeneration:InvalidMaximumForcedHorizontalWavenumber")
            testCase.verifyError(@()WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0],maximumForcedVerticalMode=NaN),"WVPseudoTopographicWaveGeneration:InvalidMaximumForcedVerticalMode")
            testCase.verifyError(@()WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0],name=""),"WVPseudoTopographicWaveGeneration:InvalidName")

            barotropic = WVTransformBarotropicQG([wvt.Lx wvt.Ly],[wvt.Nx wvt.Ny],latitude=45,shouldAntialias=false);
            testCase.verifyError(@()WVPseudoTopographicWaveGeneration(barotropic,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0]),"WVPseudoTopographicWaveGeneration:UnsupportedTransform")
            variableN2 = WVTransformBoussinesq([wvt.Lx wvt.Ly wvt.Lz],[wvt.Nx wvt.Ny wvt.Nz],N2=@(z)2e-5*exp(z/4000),latitude=45,shouldAntialias=false);
            variableForcing = WVPseudoTopographicWaveGeneration(variableN2,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0]);
            testCase.verifyClass(variableForcing,"WVPseudoTopographicWaveGeneration")
        end

        function manualSpectralBoundsRestrictGeneration(testCase)
            wvt = WVTransformBoussinesq([40e3 30e3 2e3],[8 8 9],N2=@(z)2e-5*exp(z/4000),latitude=45,shouldAntialias=false);
            terrain = TestWVPseudoTopographicWaveGeneration.bandLimitedTopography(wvt);
            horizontalBound = 1.5*min(wvt.dk,wvt.dl);
            verticalBound = 2;
            forcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05+0.01i; -0.02],maximumForcedHorizontalWavenumber=horizontalBound,maximumForcedVerticalMode=verticalBound);
            unrestricted = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=forcing.barotropicVelocityAmplitude,shouldAvoidAdaptiveDamping=false);

            [mask,components] = forcing.spectralGenerationMask();
            expectedHorizontal = wvt.Kh <= horizontalBound;
            expectedVertical = wvt.J <= verticalBound;
            expectedWaveValidity = logical(wvt.waveComponent.maskAp) | logical(wvt.waveComponent.maskAm);
            testCase.verifyEqual(components.horizontalBound,expectedHorizontal)
            testCase.verifyEqual(components.verticalBound,expectedVertical)
            testCase.verifyEqual(components.waveValidity,expectedWaveValidity)
            testCase.verifyEqual(components.adaptiveDamping,true(size(mask)))
            testCase.verifyEqual(mask,expectedWaveValidity & expectedHorizontal & expectedVertical)
            testCase.verifyEqual(components.effectivePositive,mask & logical(wvt.waveComponent.maskAp))
            testCase.verifyEqual(components.effectiveNegative,mask & logical(wvt.waveComponent.maskAm))

            wvt.t = 317;
            [Fp,Fm] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            [unrestrictedFp,unrestrictedFm] = unrestricted.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            testCase.verifyEqual(Fp,mask.*unrestrictedFp,AbsTol=1e-15)
            testCase.verifyEqual(Fm,mask.*unrestrictedFm,AbsTol=1e-15)
        end

        function adaptiveDampingMaskTracksForcingChanges(testCase)
            wvt = WVTransformBoussinesq([40e3 30e3 2e3],[8 8 9],N2=@(z)2e-5*exp(z/4000),latitude=45,shouldAntialias=false);
            terrain = TestWVPseudoTopographicWaveGeneration.bandLimitedTopography(wvt);
            forcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; -0.02]);
            legacyForcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=forcing.barotropicVelocityAmplitude,shouldAvoidAdaptiveDamping=false);
            wvt.removeAllForcing();
            wvt.addForcing(forcing);
            wvt.t = 421;

            [legacyFp,legacyFm] = legacyForcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            [beforeFp,beforeFm] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            testCase.verifyEqual(beforeFp,legacyFp)
            testCase.verifyEqual(beforeFm,legacyFm)

            damping = WVAdaptiveDamping(wvt);
            wvt.addForcing(damping);
            [mask,components] = forcing.spectralGenerationMask();
            testCase.verifyEqual(components.adaptiveDamping,damping.damp == 0)
            testCase.verifyFalse(any(mask(damping.damp ~= 0),"all"))
            testCase.verifyGreaterThan(nnz(mask),0)
            [maskedFp,maskedFm] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            testCase.verifyEqual(maskedFp,mask.*legacyFp,AbsTol=1e-15)
            testCase.verifyEqual(maskedFm,mask.*legacyFm,AbsTol=1e-15)
            testCase.verifyEqual(maskedFp(damping.damp ~= 0),zeros(nnz(damping.damp ~= 0),1))
            testCase.verifyEqual(maskedFm(damping.damp ~= 0),zeros(nnz(damping.damp ~= 0),1))

            originalDamping = damping.damp;
            antialiasing = WVAntialiasing(wvt,Nj=4);
            wvt.addForcing(antialiasing);
            [rebuiltMask,rebuiltComponents] = forcing.spectralGenerationMask();
            testCase.verifyNotEqual(damping.damp,originalDamping)
            testCase.verifyEqual(rebuiltComponents.adaptiveDamping,damping.damp == 0)
            testCase.verifyFalse(any(rebuiltMask(damping.damp ~= 0),"all"))
            wvt.removeForcing(antialiasing);

            wvt.removeForcing(damping);
            [afterFp,afterFm] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            testCase.verifyEqual(afterFp,legacyFp)
            testCase.verifyEqual(afterFm,legacyFm)
        end

        function maskedBottomWorkMatchesModalInjection(testCase)
            seedRandomNumberGenerator(testCase,92841);
            wvt = WVTransformBoussinesq([40e3 30e3 2e3],[8 8 9],N2=@(z)2e-5*exp(z/4000),latitude=45,shouldAntialias=true);
            terrain = TestWVPseudoTopographicWaveGeneration.bandLimitedTopography(wvt);
            forcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; -0.02],maximumForcedVerticalMode=2);
            wvt.addForcing(forcing);
            wvt.addForcing(WVAdaptiveDamping(wvt));
            wvt.t = 813;
            wvt.Ap = (randn(size(wvt.Ap))+1i*randn(size(wvt.Ap))).*wvt.waveComponent.maskAp;
            wvt.Am = (randn(size(wvt.Am))+1i*randn(size(wvt.Am))).*wvt.waveComponent.maskAm;
            wvt.A0(:) = 0;

            diagnostics = TestWVPseudoTopographicWaveGeneration.sourceDiagnostics(wvt,forcing);
            testCase.verifyLessThanOrEqual(diagnostics.powerError,5e-12)
        end

        function prescribedVelocityAndRampAreCorrect(testCase)
            wvt = TestWVPseudoTopographicWaveGeneration.createTransform(true);
            terrain = TestWVPseudoTopographicWaveGeneration.bandLimitedTopography(wvt);
            amplitude = [0.08+0.02i; -0.03+0.04i];
            frequency = 1.4e-4;
            forcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=amplitude,frequency=frequency,rampDuration=200,startTime=100);

            testCase.verifyEqual(forcing.barotropicVelocityAtTime(99),zeros(2,1))
            testCase.verifyEqual(forcing.barotropicVelocityAtTime(100),zeros(2,1))
            expectedMidRamp = 0.5*real(amplitude*exp(-1i*frequency*100));
            testCase.verifyEqual(forcing.barotropicVelocityAtTime(200),expectedMidRamp,"AbsTol",10*eps)
            expectedAfterRamp = real(amplitude*exp(-1i*frequency*350));
            testCase.verifyEqual(forcing.barotropicVelocityAtTime(450),expectedAfterRamp,"AbsTol",10*eps)

            velocity = forcing.barotropicVelocityAtTime(200);
            expectedBottom = velocity(1)*wvt.diffX(terrain)+velocity(2)*wvt.diffY(terrain);
            testCase.verifyEqual(forcing.bottomVelocityAtTime(200),expectedBottom,"AbsTol",1e-14)

            immediate = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=amplitude,startTime=100);
            testCase.verifyEqual(immediate.barotropicVelocityAtTime(100),real(amplitude))
        end

        function directOracleMatchesFourierSelectionAndKernel(testCase)
            for shouldAntialias = [false true]
                wvt = TestWVPseudoTopographicWaveGeneration.createTransform(shouldAntialias);
                terrains = {
                    TestWVPseudoTopographicWaveGeneration.sinusoidalTopography(wvt,50)
                    TestWVPseudoTopographicWaveGeneration.bandLimitedTopography(wvt)
                    };
                for iTerrain = 1:numel(terrains)
                    amplitude = [0.05+0.01i; -0.02+0.015i];
                    forcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrains{iTerrain},barotropicVelocityAmplitude=amplitude,frequency=1.405e-4);
                    for t = [0 813]
                        wvt.t = t;
                        velocity = forcing.barotropicVelocityAtTime(t);
                        [directFp,directFm,fourierFp,fourierFm,gBottom,gBottomFourier] = referenceBoundaryProjection(wvt,terrains{iTerrain},velocity);
                        [Fp,Fm,F0] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
                        testCase.verifyLessThanOrEqual(TestWVPseudoTopographicWaveGeneration.relativeError(directFp,fourierFp),5e-12)
                        testCase.verifyLessThanOrEqual(TestWVPseudoTopographicWaveGeneration.relativeError(directFm,fourierFm),5e-12)
                        testCase.verifyLessThanOrEqual(TestWVPseudoTopographicWaveGeneration.relativeError(Fp,directFp),5e-12)
                        testCase.verifyLessThanOrEqual(TestWVPseudoTopographicWaveGeneration.relativeError(Fm,directFm),5e-12)
                        testCase.verifyEqual(F0,zeros(size(wvt.A0)))

                        reconstructed = wvt.transformToSpatialDomainWithFourier(repmat(gBottomFourier,wvt.Nz,1));
                        testCase.verifyLessThanOrEqual(norm(reconstructed(:,:,1)-gBottom,"fro")/max(norm(gBottom,"fro"),eps),5e-12)
                        testCase.verifyLessThanOrEqual(max(abs(imag(reconstructed)),[],"all"),1e-12*max(1,max(abs(gBottom),[],"all")))
                    end
                end
            end
        end

        function kernelSelectionLimitsAndScaling(testCase)
            wvt = TestWVPseudoTopographicWaveGeneration.createTransform(true);
            terrain = TestWVPseudoTopographicWaveGeneration.sinusoidalTopography(wvt,50);
            xForcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0]);
            [Fp,Fm,F0] = xForcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            expectedHorizontal = abs(abs(wvt.K)-2*pi/wvt.Lx) < 100*eps(2*pi/wvt.Lx) & abs(wvt.L) < eps;
            testCase.verifyGreaterThan(norm([Fp(:); Fm(:)]),0)
            testCase.verifyLessThanOrEqual(norm(Fp(~expectedHorizontal)),1e-18)
            testCase.verifyLessThanOrEqual(norm(Fm(~expectedHorizontal)),1e-18)
            testCase.verifyTrue(all(Fp(~wvt.waveComponent.maskAp) == 0))
            testCase.verifyTrue(all(Fm(~wvt.waveComponent.maskAm) == 0))
            testCase.verifyEqual(F0,zeros(size(wvt.A0)))

            transverse = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0; 0.05]);
            [FpTransverse,FmTransverse] = transverse.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            testCase.verifyLessThanOrEqual(norm([FpTransverse(:); FmTransverse(:)]),1e-15)

            flat = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=zeros(wvt.Nx,wvt.Ny),barotropicVelocityAmplitude=[0.05; 0]);
            uniform = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=50*ones(wvt.Nx,wvt.Ny),barotropicVelocityAmplitude=[0.05; 0]);
            zeroCurrent = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=zeros(2,1));
            TestWVPseudoTopographicWaveGeneration.verifyZeroForcing(testCase,wvt,flat)
            TestWVPseudoTopographicWaveGeneration.verifyZeroForcing(testCase,wvt,uniform)
            TestWVPseudoTopographicWaveGeneration.verifyZeroForcing(testCase,wvt,zeroCurrent)

            doubledTerrain = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=2*terrain,barotropicVelocityAmplitude=[0.05; 0]);
            doubledCurrent = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.10; 0]);
            [FpTerrain,FmTerrain] = doubledTerrain.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            [FpCurrent,FmCurrent] = doubledCurrent.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            testCase.verifyEqual(FpTerrain,2*Fp,"RelTol",1e-12,"AbsTol",1e-15)
            testCase.verifyEqual(FmTerrain,2*Fm,"RelTol",1e-12,"AbsTol",1e-15)
            testCase.verifyEqual(FpCurrent,2*Fp,"RelTol",1e-12,"AbsTol",1e-15)
            testCase.verifyEqual(FmCurrent,2*Fm,"RelTol",1e-12,"AbsTol",1e-15)
        end

        function antialiasLayoutsAgreeOnCommonModes(testCase)
            wvtFull = TestWVPseudoTopographicWaveGeneration.createTransform(false);
            wvtDealiased = TestWVPseudoTopographicWaveGeneration.createTransform(true);
            terrainFull = TestWVPseudoTopographicWaveGeneration.bandLimitedTopography(wvtFull);
            terrainDealiased = TestWVPseudoTopographicWaveGeneration.bandLimitedTopography(wvtDealiased);
            amplitude = [0.05+0.01i; -0.02];
            forcingFull = WVPseudoTopographicWaveGeneration(wvtFull,topographicHeight=terrainFull,barotropicVelocityAmplitude=amplitude);
            forcingDealiased = WVPseudoTopographicWaveGeneration(wvtDealiased,topographicHeight=terrainDealiased,barotropicVelocityAmplitude=amplitude);
            [FpFull,FmFull] = forcingFull.addSpectralForcing(wvtFull,zeros(size(wvtFull.Ap)),zeros(size(wvtFull.Am)),zeros(size(wvtFull.A0)));
            [FpDealiased,FmDealiased] = forcingDealiased.addSpectralForcing(wvtDealiased,zeros(size(wvtDealiased.Ap)),zeros(size(wvtDealiased.Am)),zeros(size(wvtDealiased.A0)));

            for index = reshape(find(wvtDealiased.waveComponent.maskAp),1,[])
                match = find(wvtFull.waveComponent.maskAp & wvtFull.K == wvtDealiased.K(index) & wvtFull.L == wvtDealiased.L(index) & wvtFull.J == wvtDealiased.J(index));
                testCase.assertNumElements(match,1)
                testCase.verifyEqual(FpDealiased(index),FpFull(match),"RelTol",1e-12,"AbsTol",1e-15)
                testCase.verifyEqual(FmDealiased(index),FmFull(match),"RelTol",1e-12,"AbsTol",1e-15)
            end
        end

        function callbackAccumulatesWithoutMutatingState(testCase)
            wvt = TestWVPseudoTopographicWaveGeneration.createTransform(false);
            terrain = TestWVPseudoTopographicWaveGeneration.bandLimitedTopography(wvt);
            forcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05+0.01i; -0.02i]);
            wvt.t = 647;
            wvt.Ap(2,2) = 0.3+0.2i;
            wvt.Am(2,2) = -0.1+0.4i;
            wvt.A0(2,2) = 0.2-0.3i;
            tBefore = wvt.t;
            ApBefore = wvt.Ap;
            AmBefore = wvt.Am;
            A0Before = wvt.A0;
            [contributionFp,contributionFm] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));

            baseFp = complex(reshape(1:numel(wvt.Ap),size(wvt.Ap)),reshape(numel(wvt.Ap):-1:1,size(wvt.Ap)));
            baseFm = -2*baseFp;
            baseF0 = 3*baseFp;
            [Fp,Fm,F0] = forcing.addSpectralForcing(wvt,baseFp,baseFm,baseF0);
            testCase.verifyEqual(Fp,baseFp+contributionFp)
            testCase.verifyEqual(Fm,baseFm+contributionFm)
            testCase.verifyEqual(F0,baseF0)
            testCase.verifyEqual(wvt.t,tBefore)
            testCase.verifyEqual(wvt.Ap,ApBefore)
            testCase.verifyEqual(wvt.Am,AmBefore)
            testCase.verifyEqual(wvt.A0,A0Before)

            otherTransform = TestWVPseudoTopographicWaveGeneration.createTransform(false);
            testCase.verifyError(@()forcing.addSpectralForcing(otherTransform,baseFp,baseFm,baseF0),"WVPseudoTopographicWaveGeneration:TransformMismatch")
            convertedForcing = forcing.forcingWithResolutionOfTransform(otherTransform);
            testCase.verifyClass(convertedForcing,"WVPseudoTopographicWaveGeneration")
        end

        function noPVAndEnergyWorkIdentities(testCase)
            seedRandomNumberGenerator(testCase,67291);
            for shouldAntialias = [false true]
                wvt = TestWVPseudoTopographicWaveGeneration.createTransform(shouldAntialias);
                terrain = TestWVPseudoTopographicWaveGeneration.sinusoidalTopography(wvt,50);
                amplitude = [0.05+0.01i; -0.02+0.015i];
                forcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=amplitude,frequency=1.405e-4);
                wvt.t = 813;
                [forcingFp,forcingFm] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
                randomAp = (randn(size(wvt.Ap))+1i*randn(size(wvt.Ap))).*wvt.waveComponent.maskAp;
                randomAm = (randn(size(wvt.Am))+1i*randn(size(wvt.Am))).*wvt.waveComponent.maskAm;
                states = {
                    1e4*forcingFp,zeros(size(wvt.Am))
                    zeros(size(wvt.Ap)),1e4*forcingFm
                    1e-3*randomAp+1e4*forcingFp,1e-3*randomAm+1e4*forcingFm
                    };
                for iState = 1:size(states,1)
                    wvt.Ap = states{iState,1};
                    wvt.Am = states{iState,2};
                    wvt.A0(:) = 0;
                    diagnostics = TestWVPseudoTopographicWaveGeneration.sourceDiagnostics(wvt,forcing);
                    testCase.verifyEqual(diagnostics.F0,zeros(size(wvt.A0)))
                    testCase.verifyEqual(diagnostics.modalQGPVSource,zeros(size(wvt.A0)))
                    testCase.verifyLessThanOrEqual(diagnostics.normalizedQGPV,1e-10)
                    testCase.verifyLessThanOrEqual(diagnostics.powerError,5e-12)
                end

                wvt.Ap = states{3,1};
                wvt.Am = states{3,2};
                baseDiagnostics = TestWVPseudoTopographicWaveGeneration.sourceDiagnostics(wvt,forcing);
                doubledTerrain = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=2*terrain,barotropicVelocityAmplitude=amplitude,frequency=forcing.frequency,startTime=forcing.startTime);
                doubledCurrent = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=2*amplitude,frequency=forcing.frequency,startTime=forcing.startTime);
                terrainDiagnostics = TestWVPseudoTopographicWaveGeneration.sourceDiagnostics(wvt,doubledTerrain);
                currentDiagnostics = TestWVPseudoTopographicWaveGeneration.sourceDiagnostics(wvt,doubledCurrent);
                testCase.verifyEqual(terrainDiagnostics.modalPower,2*baseDiagnostics.modalPower,RelTol=1e-12,AbsTol=1e-15)
                testCase.verifyEqual(terrainDiagnostics.bottomPower,2*baseDiagnostics.bottomPower,RelTol=1e-12,AbsTol=1e-15)
                testCase.verifyEqual(currentDiagnostics.modalPower,2*baseDiagnostics.modalPower,RelTol=1e-12,AbsTol=1e-15)
                testCase.verifyEqual(currentDiagnostics.bottomPower,2*baseDiagnostics.bottomPower,RelTol=1e-12,AbsTol=1e-15)

                flat = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=zeros(wvt.Nx,wvt.Ny),barotropicVelocityAmplitude=amplitude,frequency=forcing.frequency,startTime=forcing.startTime);
                zeroCurrent = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=zeros(2,1),frequency=forcing.frequency,startTime=forcing.startTime);
                TestWVPseudoTopographicWaveGeneration.verifyZeroDiagnostics(testCase,wvt,flat)
                TestWVPseudoTopographicWaveGeneration.verifyZeroDiagnostics(testCase,wvt,zeroCurrent)
            end

            wvt = TestWVPseudoTopographicWaveGeneration.createTransform(true);
            terrain = TestWVPseudoTopographicWaveGeneration.sinusoidalTopography(wvt,50);
            forcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0]);
            wvt.removeAllForcing();
            wvt.addForcing(forcing);
            warningState = warning;
            warningCleanup = onCleanup(@()warning(warningState));
            warning("off","all")
            model = WVModel(wvt);
            clear warningCleanup
            model.setupIntegrator(integratorType="adaptive",absTolerance=1e-10,relTolerance=1e-8);
            model.integrateToTime(600,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            evolvedDiagnostics = TestWVPseudoTopographicWaveGeneration.sourceDiagnostics(wvt,forcing);
            testCase.verifyEqual(evolvedDiagnostics.F0,zeros(size(wvt.A0)))
            testCase.verifyEqual(evolvedDiagnostics.modalQGPVSource,zeros(size(wvt.A0)))
            testCase.verifyLessThanOrEqual(evolvedDiagnostics.normalizedQGPV,1e-10)
            testCase.verifyLessThanOrEqual(evolvedDiagnostics.powerError,5e-12)
        end

        function adaptiveModelSmoke(testCase)
            wvt = TestWVPseudoTopographicWaveGeneration.createTransform(true);
            terrain = TestWVPseudoTopographicWaveGeneration.sinusoidalTopography(wvt,50);
            forcing = WVPseudoTopographicWaveGeneration(wvt,topographicHeight=terrain,barotropicVelocityAmplitude=[0.05; 0]);
            wvt.removeAllForcing();
            wvt.addForcing(forcing);
            [Fp,Fm,F0] = wvt.nonlinearFlux();
            testCase.verifyGreaterThan(norm([Fp(:); Fm(:)]),0)
            testCase.verifyEqual(F0,zeros(size(wvt.A0)))

            warningState = warning;
            warningCleanup = onCleanup(@()warning(warningState));
            warning("off","all")
            model = WVModel(wvt);
            clear warningCleanup
            model.setupIntegrator(integratorType="adaptive",absTolerance=1e-10,relTolerance=1e-8);
            model.integrateToTime(300,shouldShowIntegrationDiagnostics=false,callback=@(~)[]);
            testCase.verifyEqual(wvt.t,300,"AbsTol",1e-10)
            testCase.verifyTrue(all(isfinite([wvt.Ap(:); wvt.Am(:); wvt.A0(:)])))
            testCase.verifyGreaterThan(norm([wvt.Ap(:); wvt.Am(:)]),0)
            testCase.verifyEqual(wvt.A0,zeros(size(wvt.A0)))
        end

        function goffTopographyIsDeterministicAndLocal(testCase)
            wvt = TestWVPseudoTopographicWaveGeneration.createTransform(false);
            previousRandomState = rng;
            [virtualDepth,terrain,diagnostics] = WVPseudoTopographicWaveGeneration.goffAbyssalHillTopography(wvt,rmsHeight=30,cornerWavenumber=1e-3,minimumWavelength=4e3,randomSeed=2023);
            randomStateAfter = rng;
            [secondVirtualDepth,secondTerrain,secondDiagnostics] = WVPseudoTopographicWaveGeneration.goffAbyssalHillTopography(wvt,rmsHeight=30,cornerWavenumber=1e-3,minimumWavelength=4e3,randomSeed=2023);

            testCase.verifyEqual(randomStateAfter,previousRandomState)
            testCase.verifyEqual(secondTerrain,terrain)
            testCase.verifyEqual(secondVirtualDepth,virtualDepth)
            testCase.verifyEqual(secondDiagnostics.fourierCoefficients,diagnostics.fourierCoefficients)
            testCase.verifyEqual(mean(terrain,"all"),0,AbsTol=1e-12)
            testCase.verifyEqual(sqrt(mean(terrain.^2,"all")),30,RelTol=1e-14)
            testCase.verifyEqual(virtualDepth,wvt.Lz-terrain)
            testCase.verifyGreaterThan(min(virtualDepth,[],"all"),0)
            testCase.verifyError(@()WVPseudoTopographicWaveGeneration.goffAbyssalHillTopography(wvt,minimumWavelength=100),"WVPseudoTopographicWaveGeneration:UnresolvedCutoff")
            testCase.verifyError(@()WVPseudoTopographicWaveGeneration.goffAbyssalHillTopography(wvt,randomSeed=double(intmax("uint32"))+1),"WVPseudoTopographicWaveGeneration:InvalidRandomSeed")
        end
    end

    methods (Static, Access = private)
        function wvt = createTransform(shouldAntialias)
            wvt = WVTransformBoussinesq([4e3 4e3 2e3],[8 4 5],N2=@(z)2e-5*ones(size(z)),latitude=45,shouldAntialias=shouldAntialias);
        end

        function terrain = sinusoidalTopography(wvt,amplitude)
            x = reshape(wvt.x,[],1);
            terrain = amplitude*cos(2*pi*x/wvt.Lx).*ones(1,wvt.Ny);
        end

        function terrain = bandLimitedTopography(wvt)
            x = reshape(wvt.x,[],1);
            y = reshape(wvt.y,1,[]);
            previousRandomState = rng;
            randomStateCleanup = onCleanup(@()rng(previousRandomState));
            rng(48271,"twister")
            modes = [1 0; 0 1; 1 1; 2 1];
            amplitudes = 5+20*rand(size(modes,1),1);
            phases = 2*pi*rand(size(modes,1),1);
            terrain = zeros(wvt.Nx,wvt.Ny);
            for iMode = 1:size(modes,1)
                terrain = terrain+amplitudes(iMode)*cos(2*pi*(modes(iMode,1)*x/wvt.Lx+modes(iMode,2)*y/wvt.Ly)+phases(iMode));
            end
        end

        function errorValue = relativeError(actual,expected)
            errorValue = norm(actual(:)-expected(:))/max(norm(expected(:)),eps);
        end

        function verifyZeroForcing(testCase,wvt,forcing)
            [Fp,Fm,F0] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            testCase.verifyEqual(Fp,zeros(size(wvt.Ap)))
            testCase.verifyEqual(Fm,zeros(size(wvt.Am)))
            testCase.verifyEqual(F0,zeros(size(wvt.A0)))
        end

        function diagnostics = sourceDiagnostics(wvt,forcing)
            [Fp,Fm,F0] = forcing.addSpectralForcing(wvt,zeros(size(wvt.Ap)),zeros(size(wvt.Am)),zeros(size(wvt.A0)));
            [Fu,Fv,~,Feta] = wvt.transformWaveVortexToUVWEta(Fp,Fm,F0,wvt.t);
            vorticityX = wvt.diffX(Fv);
            vorticityY = wvt.diffY(Fu);
            stretching = wvt.f*wvt.diffZG(Feta);
            qgpvSource = vorticityX-vorticityY-stretching;
            qgpvScale = norm(vorticityX(:))+norm(vorticityY(:))+norm(stretching(:));
            normalizedQGPV = norm(qgpvSource(:))/max(qgpvScale,eps);
            modalPower = 2*sum(wvt.Apm_TE_factor(:).*real(Fp(:).*conj(wvt.Ap(:))+Fm(:).*conj(wvt.Am(:))));
            [~,maskComponents] = forcing.spectralGenerationMask();
            kinematicPressure = wvt.g*wvt.transformToSpatialDomainWithF(Apm=wvt.NAp.*wvt.Apt.*maskComponents.effectivePositive+wvt.NAm.*wvt.Amt.*maskComponents.effectiveNegative);
            [~,iBottom] = min(wvt.z);
            bottomPower = mean(kinematicPressure(:,:,iBottom).*forcing.bottomVelocityAtTime(wvt.t),"all");
            powerError = abs(modalPower-bottomPower)/max([abs(modalPower) abs(bottomPower) eps]);
            modalQGPVSource = wvt.A0_QGPV_factor.*F0;
            diagnostics = struct(Fp=Fp,Fm=Fm,F0=F0,modalQGPVSource=modalQGPVSource,qgpvSource=qgpvSource,normalizedQGPV=normalizedQGPV,modalPower=modalPower,bottomPower=bottomPower,powerError=powerError);
        end

        function verifyZeroDiagnostics(testCase,wvt,forcing)
            diagnostics = TestWVPseudoTopographicWaveGeneration.sourceDiagnostics(wvt,forcing);
            testCase.verifyEqual(diagnostics.Fp,zeros(size(wvt.Ap)))
            testCase.verifyEqual(diagnostics.Fm,zeros(size(wvt.Am)))
            testCase.verifyEqual(diagnostics.F0,zeros(size(wvt.A0)))
            testCase.verifyEqual(diagnostics.modalQGPVSource,zeros(size(wvt.A0)))
            testCase.verifyEqual(diagnostics.qgpvSource,zeros(size(diagnostics.qgpvSource)))
            testCase.verifyEqual(diagnostics.modalPower,0)
            testCase.verifyEqual(diagnostics.bottomPower,0)
        end
    end
end
